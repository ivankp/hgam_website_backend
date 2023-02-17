#include <cstdlib>
#include <cstring>
#include <cmath>

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <chrono>
#include <vector>

#include <unistd.h>
#include <fcntl.h>
// #include <sys/stat.h>
// #include <sys/types.h>

#include "strings.hh"
#include "debug.hh"

using std::cout, std::endl;
using ivan::cat;

template <typename... T>
[[noreturn]] void error(const T&... args) {
  throw std::runtime_error(cat(args...));
}

struct uniform_axis {
  const double min, max, width;
  const unsigned nbins;

  uniform_axis(unsigned nbins, double min, double max)
  : min(min), max(max), width((max-min)/nbins), nbins(nbins) { }

  unsigned index(double x) const noexcept {
    // no overflow or underflow
    if (x < min || max <= x) return -1;
    return (x-min)/width;
  }

  double edge(unsigned i) const noexcept {
    return min + i*width;
  }
};

struct nonuniform_axis {
  const double* const edges;
  const unsigned nbins;

  nonuniform_axis(unsigned nbins, const double* edges)
  : edges(edges), nbins(nbins) { }

  unsigned index(double x) const noexcept {
    // edge at index must be less than or equal to x
    // only works for nbins >= 1
    // no overflow or underflow
    const double *a = edges, *p = a;
    unsigned n = nbins+1, d;
    while (n > 0) {
      d = n/2;
      p = a + d;
      if (x < *p) {
        n = d;
      } else {
        a = ++p;
        n -= d+1;
      }
    }
    const unsigned i = (p - edges) - 1;
    return i == nbins ? -1 : i;
  }

  double edge(unsigned i) const noexcept {
    return edges[i];
  }
};

template <auto f, typename T>
auto guard(T x) {
  class wrapper {
    T x;
  public:
    wrapper(T x): x(x) { }
    ~wrapper() { f(x); }
    operator T() const noexcept { return x; }
    T operator+() const noexcept { return x; }
    auto&& operator=(T o) noexcept { return x = o; }
  };
  return wrapper(x);
}

template <
  typename event_t,
  typename H,
  typename E
>
void read_via_buffer(
  const char* filename,
  char* buffer,
  size_t buffer_size,
  H&& header_f,
  E&& event_f
) {
  constexpr size_t event_size = sizeof(event_t);

  auto fd = guard<close>( ::open(filename,O_RDONLY) );
  if (fd < 0) error("failed to open file ",filename);

  bool header = true;
  size_t nbuf = 0; // data left in the buffer
  for (;;) { // read data file in chunks of buffer size
    const ssize_t nread = ::read(+fd,buffer+nbuf,buffer_size-nbuf);
    if (nread < 0) error("failed to read file ",filename);
    if (nread == 0) {
      if (nbuf) error("unusable data at the end of file ",filename);
      break;
    }
    nbuf += nread;
    event_t* e = reinterpret_cast<event_t*>(buffer);
    if (header) {
      const auto n = header_f(buffer,buffer+nbuf) - buffer;
      nbuf -= n;
      // assume that header size is a multiple of event_size
      e = reinterpret_cast<event_t*>( reinterpret_cast<char*>(e) - n );
      header = false;
    }
    for (;;) { // loop over events
      if (nbuf < event_size) {
        if (nbuf) memmove(buffer,e,nbuf);
        break;
      }
      event_f(*e); // call event function
      nbuf -= event_size;
      ++e;
    }
  }
}

int main(int argc, char** argv) {
  using clock = std::chrono::system_clock;
  using time  = std::chrono::time_point<clock>;
  const time start = clock::now();
try {

  float lumi = 0;
  std::string_view var;
  std::vector<double> edges;
  uint64_t n_events_data = 0, n_events_mc = 0;
  double signal_region[2] { 121, 129 };

  // parse query string =============================================
  for (char *a=argv[1], *b=a, *c{};;) {
    b = strchr(a,'&');
    if (b) *b = '\0';

    c = strchr(a,'=');
    if (c) { *c = '\0'; ++c; }
    if (!strcmp(a,"lumi")) {
      if (c && c!=b) {
        lumi = atof(c);
      }
    } else if (!strcmp(a,"var")) {
      if (c && c!=b) {
        var = c;
      }
    } else if (!strcmp(a,"edges")) {
      if (c && c!=b) {
        unsigned max_nedges = 0;
        for (char *b=c;; ++b) { // find number of edges
          char k = *b;
          if (!k) break;
          if (k=='+') ++max_nedges;
        }
        edges.reserve(++max_nedges);
        for (char *a=c, *b=a;; ++b) { // parse edges
          char k = *b;
          if (!k || k=='+') {
            *b = '\0';
            if (b!=a) edges.push_back(atof(a));
            a = b+1;
            if (!k) break;
          }
        }
      }
    }

    if (!b) break;
    a = b+1;
  }

  // validate input =================================================
  if (var.empty()) error("missing var");
  if (edges.empty()) error("missing edges");
  if (edges.size() < 2) error("fewer than 2 edges");

  std::sort(edges.begin(),edges.end());
  edges.erase( std::unique(edges.begin(),edges.end()), edges.end() );

  // binning ========================================================
  nonuniform_axis var_axis(edges.size()-1,edges.data());
  uniform_axis    myy_axis(55,105,160);

  struct data_bin { unsigned long n=0; };
  struct   mc_bin { double w=0, w2=0; };

  const unsigned nbins = var_axis.nbins * myy_axis.nbins;
  std::vector<data_bin> data_hist(nbins);
  std::vector<  mc_bin>   mc_hist(var_axis.nbins);

  std::vector<char> buffer(1 << 15);

#define READ(X) \
  memcpy(&X,m,sizeof(X)); m += sizeof(X);

  struct data_event_t { float myy, var; };
  read_via_buffer<data_event_t>(
    cat("data/",var,"_data.dat").c_str(),
    buffer.data(), buffer.size(),
    [&](char* m, char* end){ // read header
      if ((end-m) < (ssize_t)var.size()) goto too_short;
      if (var != m) error("var name ",var," doesn't match in file");
      m += var.size();

      while (!*++m) { if (m==end) goto too_short; }
      if (*m != 'd') error("wrong file type (not d)");
      ++m;

      if ((end-m) < (ssize_t)sizeof(lumi)) goto too_short;
      READ(lumi)
      if ((end-m) < (ssize_t)sizeof(n_events_data)) goto too_short;
      READ(n_events_data)

      return m;
too_short:
      error("file header is too short");
    },
    [&](data_event_t& e){ // read event
      if (signal_region[0] <= e.myy && e.myy <= signal_region[1]) return;
      // TODO: don't use empty bins

      const unsigned myy_i = myy_axis.index(e.myy);
      if (myy_i == unsigned(-1)) return;
      const unsigned var_i = var_axis.index(e.var);
      if (var_i == unsigned(-1)) return;
      const unsigned i = var_i * myy_axis.nbins + myy_i;

      ++data_hist.at(i).n;
    }
  );

  struct mc_event_t { float myy, var, truth, weight; };
  read_via_buffer<mc_event_t>(
    cat("data/",var,"_mc.dat").c_str(),
    // "data/pT_yy_mc.dat",
    buffer.data(), buffer.size(),
    [&](char* m, char* end){ // read header
      if ((end-m) < (ssize_t)var.size()) goto too_short;
      if (var != m) error("var name ",var," doesn't match in file");
      m += var.size();

      while (!*++m) { if (m==end) goto too_short; }
      if (*m != 'm') error("wrong file type (not m)");
      ++m;

      if ((end-m) < (ssize_t)sizeof(n_events_mc)) goto too_short;
      READ(n_events_mc)

      return m;
too_short:
      error("file header is too short");
    },
    [&](mc_event_t& e){ // read event
      if (!(signal_region[0] <= e.myy && e.myy <= signal_region[1])) return;

      // const unsigned myy_i = myy_axis.index(e.myy);
      // if (myy_i == unsigned(-1)) return;
      const unsigned var_i = var_axis.index(e.var);
      if (var_i == unsigned(-1)) return;
      const unsigned i = var_i; // * myy_axis.nbins + myy_i;

      auto& bin = mc_hist.at(i);
      bin.w += e.weight;
      bin.w2 += e.weight * e.weight;

      // TODO: migration (compare truth & reco)
    }
  );

  // JSON response ==================================================
  cout << "{"
    "\"n_events_data\":" << n_events_data << ","
    "\"n_events_mc\":" << n_events_mc << ","
    "\"lumi\":" << lumi << ","
    "\"var\":\"" << var << "\","
    "\"signal_region\":[" << signal_region[0] <<','<< signal_region[1] << "],"
    "\"var_bins\":[";

  for (unsigned i=0, n=edges.size(); i<n; ++i) {
    if (i) cout << ',';
    const auto x = edges[i];
    if (std::isinf(x)) cout << (x < 0 ? "\"-inf\"" : "\"inf\"");
    else cout << x;
  }

  cout << "],"
    "\"myy_bins\":["
      << myy_axis.nbins << ','
      << myy_axis.min << ','
      << myy_axis.max << "],"
    "\"sig\":[";

  for (unsigned i=0, n=var_axis.nbins; i<n; ++i) {
    if (i) cout << ',';
    cout << mc_hist[i].w * lumi;
  }

  cout << "],"
    "\"sig_sys\":[";

  for (unsigned i=0, n=var_axis.nbins; i<n; ++i) {
    if (i) cout << ',';
    cout << std::sqrt(mc_hist[i].w2) * lumi;
  }

  cout << "],"
    "\"time\":"
      << std::setprecision(3)
      << std::chrono::duration<double>(clock::now() - start).count()*1e3
      << "}";

} catch (const std::exception& e) {
  cout << "{\"error\":" << std::quoted(e.what()) << "}";
}}

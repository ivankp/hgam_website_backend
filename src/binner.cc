#include <iostream>
#include <iomanip>
#include <algorithm>
#include <chrono>
#include <vector>

#include <cstdlib>
#include <cstring>
#include <cmath>

#include <unistd.h>
#include <fcntl.h>

#include "error.hh"
// #include "debug.hh"

using std::cout;
using ivan::cat, ivan::error;

struct uniform_axis {
  double min, max, width;
  unsigned nbins;

  uniform_axis(unsigned nbins, double min, double max) noexcept
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
  const double* edges;
  unsigned nbins;

  nonuniform_axis(unsigned nbins, const double* edges) noexcept
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
  // static_assert(event_size <= buffer_size);

  auto fd = guard<::close>( ::open(filename,O_RDONLY) );
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

unsigned is_var(const char* s) {
  if (!s) return 0;
  if (*s != 'x') return 0;
  const char c = *++s;
  if (c < '1' || '9' < c) return 0;
  if (*++s != '\0') return 0;
  return c - '0';
}

int main(int argc, char** argv) { try {
  using clock = std::chrono::system_clock;
  using time  = std::chrono::time_point<clock>;
  const time start = clock::now();

  std::ios_base::sync_with_stdio(false);

  if (argc!=2) error("missing argument");

  float lumi = 0;
  struct var_t: nonuniform_axis {
    std::string_view name;
    var_t(): nonuniform_axis(0,nullptr) { }
    ~var_t() { if (edges) free(const_cast<double*>(edges)); }
  };
  std::vector<var_t> vars; vars.reserve(9);
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
    } else if (const unsigned xi = is_var(a)) {
      if (c && c!=b) {
        if (vars.size() < xi) vars.resize(xi);
        auto& var = vars[xi-1];

        int i = 0;
        for (char *b=c;; ++b) { // find number of edges
          char k = *b;
          if (!k) break;
          if (k=='+') ++i;
        }
        double* edges = static_cast<double*>(malloc(i*sizeof(double)));
        i = -1;
        for (char *a=c, *b=a;; ++b) { // parse edges
          char k = *b;
          if (!k || k=='+') {
            *b = '\0';
            if (b!=a) {
              if (i < 0) var.name = a;
              else edges[i] = atof(a);
              ++i;
            }
            a = b+1;
            if (!k) break;
          }
        }

        if (i < 0) error("missing variable name");
        if (i < 1) error("missing edges for ",var.name);

        std::sort( edges, edges+i );
        i = std::unique( edges, edges+i ) - edges;

        if (i < 2) error("fewer than 2 edges for ",var.name);

        var.edges = edges;
        var.nbins = i-1;
      }
    }

    if (!b) break;
    a = b+1;
  }

  // binning ========================================================
  const uniform_axis myy_axis(55,105,160);

  struct data_bin { unsigned long n=0; };
  struct   mc_bin { double w=0, w2=0; };

  unsigned n_data_bins = myy_axis.nbins;
  unsigned n_mc_bins = 1;
  for (const auto& var : vars) {
    n_data_bins *= var.nbins;
    n_mc_bins *= var.nbins;
  }
  std::vector<data_bin> data_hist(n_data_bins);
  std::vector<  mc_bin>   mc_hist(n_mc_bins);

  std::vector<char> buffer(1 << 15);

#define READ(X) \
  memcpy(&X,m,sizeof(X)); m += sizeof(X);

  // TODO: read multiple files in parallel

  read_via_buffer(
    [&]{
      std::vector<std::string> files;
      files.reserve(vars.size());
      for (const auto& var : vars)
        files.emplace_back(cat("data/",var.name,"_data.dat"));
      return files;
    }(),
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
    [&](double m_yy, double* x){ // read event
      if (signal_region[0] <= m_yy && m_yy <= signal_region[1]) return;
      // TODO: don't use empty bins

      unsigned B = 0;
      for (unsigned i=vars.size(); i; ) {
        const auto& var = vars[--i];
        const unsigned b = var.index(x[i]);
        if (b == unsigned(-1)) return;
        B = B * var.nbins + b;
      }
      { const unsigned b = myy_axis.index(m_yy);
        if (b == unsigned(-1)) return;
        B = B * myy_axis.nbins + b;
      }

      ++data_hist[B].n;
    }
  );

  // struct mc_event_t { float m_yy, var, truth, weight; };
  read_via_buffer/*<mc_event_t>*/(
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
    [&](double m_yy, double* x){ // read event
      if (!(signal_region[0] <= m_yy && m_yy <= signal_region[1])) return;

      unsigned B = 0;
      for (unsigned i=vars.size(); i; ) {
        const auto& var = vars[--i];
        const unsigned b = var.index(x[i]);
        if (b == unsigned(-1)) return;
        B = B * var.nbins + b;
      }

      auto& bin = mc_hist[B];
      bin.w += e.weight;
      bin.w2 += e.weight * e.weight;

      // TODO: migration (compare truth & reco)
    }
  );

  // JSON response ==================================================
  cout << "{"
    "\"lumi\":" << lumi << ","
    "\"data_events\":" << n_events_data << ","
    "\"mc_events\":" << n_events_mc << ","
    "\"signal_region\":[" << signal_region[0] <<','<< signal_region[1] << "],"
    "\"m_yy\":["
      << myy_axis.nbins << ','
      << myy_axis.min << ','
      << myy_axis.max << "],"
    "\"vars\":[";

  bool first_var = true;
  for (const auto& var : vars) {
    if (first_var) first_var = false;
    else cout << ',';
    cout << "[\"" << var.name << "\",[";
    for (unsigned i=0, n=var.nbins+1; i<n; ++i) {
      if (i) cout << ',';
      const auto x = var.edges[i];
      if (std::isinf(x)) cout << (x < 0 ? "\"-inf\"" : "\"inf\"");
      else cout << x;
    }
    cout << "]]";
  }

  cout << "],"
    "\"sig\":[";

  for (unsigned i=0; i<n_mc_bins; ++i) {
    if (i) cout << ',';
    cout << mc_hist[i].w * lumi;
  }

  cout << "],"
    "\"sig_sys\":[";

  for (unsigned i=0; i<n_mc_bins; ++i) {
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

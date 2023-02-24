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
#include "debug.hh"

using std::cout;
using ivan::cat, ivan::error;

size_t round_down(size_t x, size_t n) noexcept { return x-(x%n); }

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

// template <auto f, typename T>
// auto guard(T x) {
//   class wrapper {
//     T x;
//   public:
//     wrapper(T x): x(x) { }
//     ~wrapper() { f(x); }
//     operator T() const noexcept { return x; }
//     T operator+() const noexcept { return x; }
//     auto&& operator=(T o) noexcept { return x = o; }
//   };
//   return wrapper(x);
// }

unsigned is_var(const char* s) {
  if (!s) return 0;
  if (*s != 'x') return 0;
  const char c = *++s;
  if (c < '1' || '9' < c) return 0;
  if (*++s != '\0') return 0;
  return c - '0';
}
constexpr unsigned maxnvars = 9;

struct Var: nonuniform_axis {
  std::string_view name;
  Var(): nonuniform_axis(0,nullptr) { }
  ~Var() { if (edges) free((void*)edges); }
} vars[maxnvars];
unsigned nvars = 0;

double lumi = 0;
uint64_t n_events_data = 0, n_events_mc = 0;
double signal_region[2] { 121, 129 };

template <bool mc, typename F>
void read_events(F&& event_f) {
  struct file {
    char *buf = nullptr;
    unsigned buflen, nbuf = 0, event_size;
    int fd = -1;
    char type;
    std::string path;

    ~file() {
      if (fd >= 0) ::close(fd);
      if (buf) ::free(buf);
    }
  } files[maxnvars];

  // open files and read headers
  for (unsigned i=0; i<nvars; ++i) {
    auto& var = vars[i];
    auto& f = files[i];

    // open file for reading
    f.fd = ::open((f.path = cat(
      "data/", var.name, mc ? "_mc.dat" : "_data.dat"
    )).c_str(), O_RDONLY);
    if (f.fd < 0) error("failed to open ",f.path);

    TEST(f.path)

    // read headers -----------------------------------------------
    char buf[128]; // assume this is enough for any header
    const char* m = buf;

    float _lumi;
    uint64_t _nevents;

    const size_t header_tail_len = 2 + (mc
      ? sizeof(_nevents)
      : sizeof(_lumi) + sizeof(_nevents)
    );
    const size_t header_len = round_down(
      var.name.size()+8 + header_tail_len, 8
    );
    if (sizeof(buf) < header_len)
      error("insufficient header buffer size for ",f.path);

    const ssize_t nread = ::read(f.fd,buf,header_len);
    if (nread < 0) error("failed to read ",f.path);
    if (nread != ssize_t(header_len))
      error("failed to read whole header in ",f.path);

    // check that variable name is as expected
    if (var.name != m) error(f.path," does not start with ",var.name);

    // next byte should be d or m
    m += (header_len-header_tail_len);
    if (*m != (mc ? 'm' : 'd')) error("wrong data marker in ",f.path);
    // next byte specifies the variable type
    f.type = *++m;
    ++m;

    // read constants
    if (!mc) { memcpy(&_lumi,m,sizeof(_lumi)); m += sizeof(_lumi); }
    memcpy(&_nevents,m,sizeof(_nevents));

    if (i) {
      if (!mc) {
        if (_lumi != lumi) error("different lumi in ",f.path);
        if (_nevents != n_events_data)
          error("different number of events in ",f.path);
      } else {
        if (_nevents != n_events_mc)
          error("different number of events in ",f.path);
      }
    } else {
      if (!mc) {
        lumi = _lumi;
        n_events_data = _nevents;
      } else {
        n_events_mc = _nevents;
      }
    }

    // allocate buffer
    switch (f.type) {
      case 'f': f.event_size = 4; break;
      case 'i': f.event_size = 4; break;
      case 'B': f.event_size = 1; break;
      case 'c': f.event_size = 1; break;
      default: error("unexpected type marker \'",f.type,"\' in ",f.path);
    }
    if (mc) f.event_size *= 2;
    f.buflen = f.event_size * (1 << 16);
    f.buf = static_cast<char*>(malloc(f.buflen));
  }

  double event[sizeof(double)*2*maxnvars];
  for (;;) { // read data file in chunks of buffer size
    unsigned nevents = -1;
    for (unsigned i=0; i<nvars; ++i) {
      auto& f = files[i];
      const ssize_t nread = ::read(f.fd,f.buf+f.nbuf,f.buflen-f.nbuf);
      if (nread < 0) error("failed to read ",f.path);
      f.nbuf += nread;
      const unsigned read_events = f.nbuf / f.event_size;
      if (read_events < nevents) nevents = read_events;
    } // for files
    if (nevents==0) break;
    for (unsigned e=0; e<nevents; ++e) { // loop over events
      double* x = event;
      // combine events from multiple files
      for (unsigned i=0; i<nvars; ++i) {
        auto& f = files[i];
        const char* m = f.buf + f.event_size*e;

        for (unsigned n=(mc?2:1); n--; ++x) {
          switch (f.type) {
#define CASE(t,type) \
            case t: \
              *x = *reinterpret_cast<const type*>(m); \
              m += sizeof(type); \
              break;

            CASE('f',float)
            CASE('i',int32_t)
            CASE('B',uint8_t)
            CASE('c',char)
            default: ;

#undef CASE
          }
        }
      } // for files
      // event_f(event); // call event processing function
    } // for events
    for (unsigned i=0; i<nvars; ++i) {
      auto& f = files[i];
      unsigned processed = f.event_size*nevents;
      if (f.nbuf -= processed)
        memmove(f.buf,f.buf+processed,f.nbuf);
    }
  } // for ;;
}

int main(int argc, char** argv) {
  std::ios_base::sync_with_stdio(false);
  if (argc!=2) {
    std::cerr << "usage: " << argv[0] << " query_string\n";
    return 1;
  }

try {
  using clock = std::chrono::system_clock;
  using time  = std::chrono::time_point<clock>;
  const time start = clock::now();

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

  nvars = std::remove_if(vars,vars+maxnvars,[](const auto& var){
    return var.name.empty();
  }) - vars;

  if (nvars==0) error("no variables specified");
  TEST(nvars)

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

  read_events<false>( // read data
    [&](double* x){ // read event
      const double m_yy = x[0];
      if (signal_region[0] <= m_yy && m_yy <= signal_region[1]) return;
      // TODO: don't use empty bins

      unsigned B = 0;
      for (unsigned i=nvars; i; ) {
        const auto& var = vars[--i];
        const unsigned b = i ? var.index(x[i]) : myy_axis.index(m_yy);
        if (b == unsigned(-1)) return;
        B = B * var.nbins + b;
      }

      ++data_hist[B].n;
    }
  );

  read_events<true>( // read mc
    [&](double* x){ // read event
      const double m_yy = x[0];
      if (!(signal_region[0] <= m_yy && m_yy <= signal_region[1])) return;
      const double weight = x[1];

      unsigned B = 0;
      for (unsigned i=nvars; i; ) {
        const auto& var = vars[--i];
        const unsigned b = var.index(x[i*2]);
        if (b == unsigned(-1)) return;
        B = B * var.nbins + b;
      }

      auto& bin = mc_hist[B];
      bin.w += weight;
      bin.w2 += weight * weight;

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

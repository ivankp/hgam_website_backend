#include <iostream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <numeric>
#include <chrono>
#include <limits>
#include <vector>

#include <cstdlib>
#include <cstring>
#include <cmath>

#include <unistd.h>
#include <fcntl.h>

#include "least_squares.h"
#include "pool.hh"
#include "numconv.hh"
#include "error.hh"
#include "debug.hh"

using std::cout;
using ivan::cat, ivan::error;

// Global variables =================================================
double data_lumi = 0, lumi = 0;
uint64_t nevents_data = 0, nevents_mc = 0;
double fiducial_myy[] { 105, 160 };
double   signal_myy[] { 121, 129 };
double myy_binwidth = 1;
constexpr unsigned maxnvars = 9;
unsigned nvars = 0;
int ExpPolyN = 0;
// ==================================================================

size_t round_down(size_t x, size_t n) noexcept { return x-(x%n); }

template <typename T>
constexpr T sq(T x) noexcept { return x*x; }

constexpr unsigned overflow = std::numeric_limits<unsigned>::max();

struct uniform_axis {
  double min, max, width;
  unsigned nbins;

  uniform_axis(unsigned nbins, double min, double max) noexcept
  : min(min), max(max), width((max-min)/nbins), nbins(nbins) { }

  unsigned index(double x) const noexcept {
    // no overflow or underflow
    if (x < min || max <= x) return overflow;
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
    return i == nbins ? overflow : i;
  }

  double edge(unsigned i) const noexcept {
    return edges[i];
  }
};

unsigned is_var(const char* s) {
  if (!s) return 0;
  if (*s != 'x') return 0;
  const char c = *++s;
  if (c < '1' || '9' < c) return 0;
  if (*++s != '\0') return 0;
  return c - '0';
}

struct Var: nonuniform_axis {
  std::string_view name;
  Var(): nonuniform_axis(0,nullptr) { }
  Var& operator=(Var& o) = delete;
  Var& operator=(Var&& o) noexcept {
    nbins = o.nbins;
    edges = o.edges; o.edges = nullptr;
    name = o.name;
    return *this;
  }
  ~Var() { if (edges) { free((void*)edges); } }
} vars[maxnvars];

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
  } files[maxnvars+1];

  // open files and read headers
  for (unsigned v=0; v<=nvars; ++v) {
    const auto var_name = v
      ? std::string_view(vars[v-1].name)
      : std::string_view("m_yy");
    auto& f = files[v];

    // open file for reading
    f.fd = ::open((f.path = cat(
      "data/", var_name, mc ? "_mc.dat" : "_data.dat"
    )).c_str(), O_RDONLY);
    if (f.fd < 0) error("failed to open ",f.path);

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
      var_name.size()+8 + header_tail_len, 8
    );
    if (sizeof(buf) < header_len)
      error("insufficient header buffer size for ",f.path);

    const ssize_t nread = ::read(f.fd,buf,header_len);
    if (nread < 0) error("failed to read ",f.path);
    if (nread != ssize_t(header_len))
      error("failed to read whole header in ",f.path);

    // check that variable name is as expected
    if (var_name != m) error(f.path," does not start with ",var_name);

    // next byte should be d or m
    m += (header_len-header_tail_len);
    if (*m != (mc ? 'm' : 'd')) error("wrong data marker in ",f.path);
    // next byte specifies the variable type
    f.type = *++m;
    ++m;

    // read constants
    if (!mc) { memcpy(&_lumi,m,sizeof(_lumi)); m += sizeof(_lumi); }
    memcpy(&_nevents,m,sizeof(_nevents));

    if (v) {
      if (!mc) {
        if (_lumi != data_lumi) error("different lumi in ",f.path);
        if (_nevents != nevents_data)
          error("different number of events in ",f.path);
      } else {
        if (_nevents != nevents_mc)
          error("different number of events in ",f.path);
      }
    } else {
      if (!mc) {
        data_lumi = _lumi;
        nevents_data = _nevents;
      } else {
        nevents_mc = _nevents;
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

  uint64_t nevents_total = 0;
  double event[sizeof(double)*2*maxnvars];
  for (;;) { // read data file in chunks of buffer size
    unsigned nevents = -1;
    for (unsigned v=0; v<=nvars; ++v) {
      auto& f = files[v];
      const ssize_t nread = ::read(f.fd,f.buf+f.nbuf,f.buflen-f.nbuf);
      if (nread < 0) error("failed to read ",f.path);
      f.nbuf += nread;
      const unsigned nevents_batch = f.nbuf / f.event_size;
      if (nevents_batch < nevents) nevents = nevents_batch;
    } // for files
    if (nevents==0) break;
    nevents_total += nevents;
    for (unsigned e=0; e<nevents; ++e) { // loop over events
      double* x = event;
      // combine events from multiple files
      for (unsigned v=0; v<=nvars; ++v) {
        auto& f = files[v];
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
      event_f(event); // call event processing function
    } // for events
    for (unsigned v=0; v<nvars; ++v) {
      auto& f = files[v];
      unsigned processed = f.event_size*nevents;
      if (f.nbuf -= processed)
        memmove(f.buf,f.buf+processed,f.nbuf);
    }
  } // for ;;
  if (nevents_total != (mc ? nevents_mc : nevents_data)) error(
    (mc?"mc":"data")," contains ",(mc ? nevents_mc : nevents_data)," events"
    ", but ",nevents_total," were read"
  );
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
  bool fit_passed = false;
  for (char *a=argv[1], *b=a, *c{};;) {
    b = strchr(a,'&');
    if (b) *b = '\0';

    c = strchr(a,'=');
    if (c) { *c = '\0'; ++c; }
    if (!strcmp(a,"lumi")) {
      if (c && c!=b) {
        lumi = atof(c);
        if (lumi <= 0) lumi = 0;
      }
    } else if (!strcmp(a,"fit")) {
      if (c && c!=b) {
        if (ivan::starts_with(c,"ExpPoly")) {
          c += 7;
          const char d = *c;
          if (d < '0' || '5' < d || c[1]!='\0') error(
            "invalid ExpPoly fit function degree ",c
          );
          ExpPolyN = d - '0';
          fit_passed = true;
        } else error("invalid fit function ",c);
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

  // Binning ========================================================
  const uniform_axis myy_axis(
    round((fiducial_myy[1]-fiducial_myy[0])/myy_binwidth),
    fiducial_myy[0],
    fiducial_myy[1]
  );
  if (myy_binwidth != myy_axis.width) error(
    "m_yy bin width does not fit into the fiducial region "
    "a whole number of times"
  );

  // make sure signal region matches bit edges
  const unsigned myy_nbins_left =
    floor( (signal_myy[0]-fiducial_myy[0])/myy_binwidth );
  const unsigned myy_nbins_right =
    ceil ( (fiducial_myy[1]-signal_myy[1])/myy_binwidth );

  signal_myy[0] = fiducial_myy[0] + myy_nbins_left *myy_binwidth;
  signal_myy[1] = fiducial_myy[1] - myy_nbins_right*myy_binwidth;

  // exclude signal region bins
  const unsigned myy_nbins_sides  = myy_nbins_left + myy_nbins_right;
  const unsigned myy_nbins_signal = myy_axis.nbins - myy_nbins_sides;

  struct data_bin { unsigned long n=0; };
  struct   mc_bin { double w=0, w2=0; };

  unsigned nbins_vars = 1;
  for (unsigned v=nvars; v; ) {
    nbins_vars *= vars[--v].nbins;
  }
  const unsigned nbins_data = myy_nbins_sides * nbins_vars;

  std::vector<data_bin> data_hist(nbins_data);
  std::vector<  mc_bin>   mc_hist(nbins_vars);
  std::vector<double> migration_hist(sq(nbins_vars));

  read_events<false>( // read data
    [&](double* x){ // read event
      const double myy = x[0];
      if (signal_myy[0] <= myy && myy <= signal_myy[1]) return;
      ++x;

      unsigned B = 0;
      { for (unsigned v=nvars; v; ) {
          const auto& var = vars[--v];
          const unsigned b = var.index(x[v]);
          if (b == overflow) return;
          B = B * var.nbins + b;
        }
        unsigned b = myy_axis.index(myy);
        if (b == overflow) return;
        if (b > myy_nbins_left) // safe because of myy cut above
          b -= myy_nbins_signal;
        B = B * myy_nbins_sides + b;
      }

      ++data_hist[B].n;
    }
  );

  read_events<true>( // read mc
    [&](double* x){ // read event
      const double myy = x[0];
      if (!(signal_myy[0] <= myy && myy <= signal_myy[1])) return;
      const double weight = x[1];
      x += 2;

      unsigned B = 0, Bt/*truth*/ = 0;
      for (unsigned v=nvars; v; ) {
        const auto& var = vars[--v];
        { const unsigned b = var.index(x[v*2]);
          if (b == overflow) return;
          B = B * var.nbins + b;
        }
        if (Bt != overflow) {
          const unsigned b = var.index(x[v*2+1]);
          Bt = ( b==overflow ? overflow : (Bt * var.nbins + b) );
        }
      }

      auto& bin = mc_hist[B];
      bin.w += weight;
      bin.w2 += weight * weight;

      if (Bt != overflow) {
        migration_hist[B*nbins_vars+Bt] += weight;
      }
    }
  );

  // Fit ============================================================
  std::stringstream fits;

  if (!fit_passed) ExpPolyN = 2;
  { const unsigned np = ExpPolyN + 1;
    const unsigned nb = myy_nbins_sides;
    const unsigned N = np*(np+1)/2;

    const auto [
      y, u, A, W, P, cov
    ] = ivan::pool<double>(
      nb, nb, nb*np, N+nb, np, N
    );

    { double* a = A;
      for (unsigned i=0; i<nb; ++i, ++a) {
        *a = 1;
      }
      for (unsigned p=1; p<np; ++p) {
        for (unsigned i=0; i<nb; ++i, ++a) {
          const double x =
            ( ( i<myy_nbins_left ? i : i+myy_nbins_signal ) + 0.5 ) // center
            * myy_binwidth - (125. - fiducial_myy[0]);
          *a = *(a-nb) * x;
        }
      }
    }

    const auto* h = data_hist.data();
    for (unsigned b=0; b<nbins_vars; ++b) {
      for (unsigned i=0; i<nb; ++i, ++h) {
        const double n = h->n;
        if (n > 0) {
          y[i] = std::log(n);
          u[i] = 1./std::sqrt(n);
        } else {
          y[i] = -1;
          u[i] = 10;
        }
      }

      linear_least_squares(A, y, u, nb, np, W, P, cov);

      if (b) fits << ',';
      fits << '[';
      for (unsigned p=0; p<np; ++p) {
        if (p) fits << ',';
        fits << P[p];
      }
      fits << ']';
    }
  }

  // JSON response ==================================================
  if (lumi == 0) lumi = data_lumi;
  cout << "{"
    "\"lumi\":[" << lumi;
  if (lumi != data_lumi) cout << ',' << data_lumi;
  cout << "],"
    "\"nevents_data\":" << nevents_data << ","
    "\"nevents_mc\":" << nevents_mc << ","
    "\"fiducial_myy\":[" << fiducial_myy[0] <<','<< fiducial_myy[1] << "],"
    "\"signal_myy\":[" << signal_myy[0] <<','<< signal_myy[1] << "],"
    "\"m_yy\":["
      << myy_axis.nbins << ','
      << myy_axis.min << ','
      << myy_axis.max << "],"
    "\"vars\":[";

  for (unsigned v=0; v<nvars; ++v) {
    const auto& var = vars[v];
    if (v) cout << ',';
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

  for (unsigned i=0; i<nbins_vars; ++i) {
    if (i) cout << ',';
    cout << mc_hist[i].w * lumi;
  }

  cout << "],"
    "\"sig_sys\":[";

  for (unsigned i=0; i<nbins_vars; ++i) {
    if (i) cout << ',';
    cout << std::sqrt(mc_hist[i].w2) * lumi;
  }

  cout << "],"
    "\"bkg\":[";
  {  const auto* h = data_hist.data();
    for (unsigned i=0; i<nbins_vars; ++i) {
      if (i) cout << ',';
      const auto* h2 = h + myy_nbins_left;
      cout << '['
        << std::accumulate(h,h2,0ul,[](auto a, auto x){ return a+x.n; })
        << ',';
      h += myy_nbins_sides;
      cout << 0 << ',';
      cout << std::accumulate(h2,h,0ul,[](auto a, auto x){ return a+x.n; })
        << ']';
    }
  }

  cout << "],"
    "\"fit\":[" << fits.rdbuf();

  cout << "],"
    "\"migration\":[";

  for (unsigned i=0, n=sq(nbins_vars); i<n; ++i) {
    if (i) cout << ',';
    cout << migration_hist[i] * lumi;
  }

  cout << "],"
    "\"time\":"
      << std::setprecision(3)
      << std::chrono::duration<double>(clock::now() - start).count()*1e3
      << "}";

} catch (const std::exception& e) {
  cout << "{\"error\":" << std::quoted(e.what()) << "}";
}}

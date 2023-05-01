#include <iostream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <numeric>
#include <chrono>
#include <limits>

#include <cstdlib>
#include <cmath>

#include <unistd.h>
#include <fcntl.h>

#include "serial_parsing.hh"
#include "pool.hh"
#include "numconv.hh"
#include "error.hh"
#include "debug.hh"

using std::cout;
using ivan::cat, ivan::error, ivan::cnt_of;

size_t round_down(size_t x, size_t n) noexcept { return x-(x%n); }

template <typename It, typename T>
auto find(It begin, It end, T& x) {
  const auto it = std::lower_bound(begin, end, x);
  return (it == end || *it != x) ? end : it;
}

uint64_t nevents;

constexpr unsigned maxnvars = 16;
unsigned nvars = 0;

constexpr unsigned maxncuts = 16;
unsigned ncuts = 0;

constexpr unsigned maxnallvars = maxnvars + maxncuts;
unsigned nallvars = 0;

const char* dataset { };
struct {
  std::string_view name;
  unsigned i;
} vars[maxnvars];
struct {
  std::string_view name;
  unsigned i;
  bool gt;
  double val;
} cuts[maxncuts];
std::string_view allvars[maxnallvars];

template <typename F>
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
  for (unsigned v=0; v<nallvars; ++v) {
    const auto var_name = allvars[v];
    auto& f = files[v];

    // open file for reading
    f.fd = ::open((f.path = cat(
      "data/", dataset, '/', var_name, "_data.dat"
    )).c_str(), O_RDONLY);
    if (f.fd < 0) error("failed to open ",f.path);

    // read headers -----------------------------------------------
    char buf[128]; // assume this is enough for any header
    const char* m = buf;

    uint64_t _nevents;

    const size_t header_tail_len = 2 + 4/*lumi*/ + sizeof(_nevents);
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
    if (*m != 'd') error("wrong data marker in ",f.path);
    // next byte specifies the variable type
    f.type = *++m;
    ++m;

    // read constants
    m += 4; // skip lumi
    memcpy(&_nevents,m,sizeof(_nevents));

    if (v) {
      if (_nevents != nevents)
        error("different number of events in ",f.path);
    } else {
      nevents = _nevents;
    }

    // allocate buffer
    switch (f.type) {
      case 'f': f.event_size = 4; break;
      case 'i': f.event_size = 4; break;
      case 'B': f.event_size = 1; break;
      case 'c': f.event_size = 1; break;
      default: error("unexpected type marker \'",f.type,"\' in ",f.path);
    }
    f.buflen = f.event_size * (1 << 16);
    f.buf = static_cast<char*>(malloc(f.buflen));
  }

  uint64_t nevents_total = 0;
  double event[maxnallvars];
  for (;;) { // read data file in chunks of buffer size
    unsigned nevents = -1;
    for (unsigned v=0; v<nallvars; ++v) {
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
      for (unsigned v=0; v<nallvars; ++v) {
        auto& f = files[v];
        const char* m = f.buf + f.event_size*e;

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
        ++x;
      } // for files
      if (event_f(event)) return; // call event processing function
    } // for events
    for (unsigned v=0; v<nallvars; ++v) {
      auto& f = files[v];
      unsigned processed = f.event_size*nevents;
      if (f.nbuf -= processed)
        memmove(f.buf,f.buf+processed,f.nbuf);
    }
  } // for ;;
  if (nevents_total != nevents) error(
    "data was expected to contain ",nevents," events, "
    "but ",nevents_total," were read"
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

  ivan::serial_split(argv[1],'&',[&](char* key, char* end){
    char* val = strchr(key,'=');
    if (val) {
      if (val == end) return; // assumes all keys must have values
      *val = '\0';
      ++val;
    } else return;
    if (!strcmp(key,"ds")) {
      dataset = val;
    } else {
      const bool is_vars = !strcmp(key,"vars");
      const bool is_cuts = !strcmp(key,"cuts");
      if (is_vars || is_cuts) {
        unsigned i = 0;
        auto* cut = cuts;
        ivan::serial_split(val,'+',[&](char* val, char* end){
          if (is_vars) {
            vars[nvars].name = { val, end ? end-val : strlen(val) };
            if (++nvars > maxncuts) error("too many cuts");
          } else { // is_cuts
            switch (i%3) {
              case 0:
                if (++ncuts > maxncuts) error("too many cuts");
                cut->name = { val, end ? end-val : strlen(val) };
                break;
              case 1:
                cut->gt = *val != '0';
                break;
              case 2:
                cut->val = atof(val);
                ++cut;
                break;
            }
            ++i;
          }
        });
      }
    }
  });

  if (!dataset) error("no dataset specified");
  if (nvars==0) error("no variables specified");

  // collect distinct variables =====================================
  for (unsigned i=0; i<nvars; ++i)
    allvars[i] = vars[i].name;
  for (unsigned i=0; i<ncuts; ++i)
    allvars[nvars+i] = cuts[i].name;
  nallvars = nvars + ncuts;

  std::sort(allvars,allvars+nallvars);
  nallvars = std::unique(allvars,allvars+nallvars) - allvars;

  for (unsigned i=0; i<nvars; ++i)
    vars[i].i = find(allvars,allvars+nallvars,vars[i].name) - allvars;
  for (unsigned i=0; i<ncuts; ++i)
    cuts[i].i = find(allvars,allvars+nallvars,cuts[i].name) - allvars;

  std::stringstream tab;
  tab.precision(8);

  { bool first = true;
    unsigned noutevents = 1000;
    unsigned event_i = 0;
    read_events( // read data
      [&](double* xs) -> bool { // read event
        ++event_i;
        for (unsigned i=0; i<ncuts; ++i) {
          const auto& cut = cuts[i];
          const double x = xs[cut.i];
          const double v = cut.val;
          if (!(cut.gt ? x > v : x < v)) return false;
        }
        if (first) first = false;
        else tab << ',';
        tab << '[' << event_i;
        for (unsigned i=0; i<nvars; ++i) {
          tab << ',';
          const double x = xs[vars[i].i];
          if (std::isnan(x)) tab << "null";
          else tab << x;
        }
        tab << ']';
        return !--noutevents;
      }
    );
  }

  // JSON response ==================================================
  cout << "{"
    "\"selection\":{"
      "\"ds\":\"" << dataset << "\","
      "\"vars\":[";
  for (unsigned i=0; i<nvars; ++i) {
    if (i) cout << ',';
    cout << std::quoted(vars[i].name);
  }
  cout <<
      "],"
      "\"cuts\":[";
  for (unsigned i=0; i<ncuts; ++i) {
    if (i) cout << ',';
    const auto& cut = cuts[i];
    cout << '[' << std::quoted(cut.name) << ','
         << char('0'+cut.gt) << ','
         << cut.val << ']';
  }
  cout <<
      "]"
    "},"
    "\"events\":["
      // << tab.rdbuf() <<
      << tab.view() <<
    "],"
    "\"time\":"
      << std::setprecision(3)
      << std::chrono::duration<double>(clock::now() - start).count()*1e3
      << "}";

} catch (const std::exception& e) {
  cout << "{\"error\":" << std::quoted(e.what()) << "}";
}}

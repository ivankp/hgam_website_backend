#include <iostream>

#include <unistd.h>
#include <fcntl.h>

#include "error.hh"
#include "debug.hh"

using std::cout, std::cerr;
using ivan::error;

size_t round_down(size_t x, size_t n) noexcept { return x-(x%n); }

int main(int argc, char** argv) {
  std::ios_base::sync_with_stdio(false);
  if (argc != 2) {
    cerr << "usage: " << argv[0] << " file.dat\n";
    return 1;
  }

  char buf[128];
  const int fd = ::open(argv[1],O_RDONLY);
  if (fd < 0) error("failed to open ",argv[1]);

  const ssize_t nread = ::read(fd,buf,sizeof(buf));
  if (nread < 0) error("failed to read ",argv[1]);

  const char *m = static_cast<const char*>(memchr(buf,'\0',nread));

  if (!m) error("no null character in header");

  cout.write(buf,m-buf) << '\n';

  while (!*m) ++m;

  switch (*m) {
    case 'd':
      cout << "d: Data events\n";
      break;
    case 'm':
      cout << "m: Monte Carlo simulated events\n";
      break;
    default:
      error("unexpected source byte \'",*m,'\'');
  }
  const bool mc = *m == 'm';
  ++m;

  switch (*m) {
    case 'f':
      cout << "f: float, 4 bytes\n";
      break;
    case 'i':
      cout << "i: int, 4 bytes\n";
      break;
    case 'B':
      cout << "B: unsigned int, 1 byte\n";
      break;
    case 'c':
      cout << "c: char, 1 byte\n";
      break;
    default:
      error("unexpected type byte \'",*m,'\'');
  }
  ++m;

  if (!mc) {
    cout << "Luminosity: "
      << (*reinterpret_cast<const float*>(m)) << " ifb\n";
    m += sizeof(float);
  }

  cout << "N events: "
    << (*reinterpret_cast<const uint64_t*>(m)) << '\n';
  // m += sizeof(uint64_t);

  ::close(fd);
}

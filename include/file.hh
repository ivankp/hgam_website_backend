#ifndef IVAN_FILE_HH
#define IVAN_FILE_HH

#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <cstdlib>
#include <cerrno>
#include <stdexcept>
#include "strings.hh"

namespace ivan {

struct file {
  int fd;
  file(const char* name): fd(::open(name,O_RDONLY)) {
    if (fd < 0) throw std::runtime_error(cat(
      '\"',name,"\": open(): ",strerror(errno)
    ));
  }
  file(int i) noexcept: fd(i) { }
  ~file() { ::close(fd); }
  operator int() const noexcept { return fd; }
};

struct regular_file: file {
  size_t size;
  regular_file(const char* name): file(name) {
    struct stat sb;
    if (fstat(fd,&sb) < 0) throw std::runtime_error(cat(
      '\"',name,"\": fstat(): ",strerror(errno)
    ));
    if (!S_ISREG(sb.st_mode)) throw std::runtime_error(cat(
      '\"',name,"\": not a regular file"
    ));
    size = (size_t)sb.st_size;
  }
};

template <typename T = char>
std::pair<T*,size_t> read_whole_file(const char* name) {
  regular_file f(name);
  T* m = static_cast<T*>(malloc(f.size));
  if (read(f,m,f.size) < 0) {
    free(m);
    throw std::runtime_error(cat(
      '\"',name,"\": read(): ",strerror(errno)
    ));
  }
  return { m, f.size };
}
inline
std::pair<char*,size_t> read_whole_file(const char* name, bool null) {
  regular_file f(name);
  char* m = static_cast<char*>(malloc(f.size+null));
  if (read(f,m,f.size) < 0) {
    free(m);
    throw std::runtime_error(cat(
      '\"',name,"\": read(): ",strerror(errno)
    ));
  }
  if (null) m[f.size] = '\0';
  return { m, f.size };
}

inline
bool file_exists(const char* fname) {
  struct stat sb;
  return stat(fname,&sb) == 0;
}

inline
void write_loop(int fd, const char* r, const char* end) {
  while (r<end) {
    const ssize_t w = write(fd,r,end-r);
    if (w < 0) throw std::runtime_error(cat(
      "write(): ",strerror(errno)
    ));
    r += w;
  }
}

}

#endif

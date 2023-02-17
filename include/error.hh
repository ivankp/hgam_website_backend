#ifndef IVAN_ERROR_HH
#define IVAN_ERROR_HH

#include "strings.hh"
#include <stdexcept>

namespace ivan {

[[noreturn]] void error(const char* str) {
  throw std::runtime_error(str);
}
[[noreturn]] void error(const std::string& str) {
  throw std::runtime_error(str);
}
template <typename... T>
[[noreturn]] void error(const T&... args) {
  throw std::runtime_error(cat(args...));
}

}

#endif

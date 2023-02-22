#ifndef IVAN_ERROR_HH
#define IVAN_ERROR_HH

#include "strings.hh"
#include <stdexcept>

namespace ivan {

[[noreturn]]
[[gnu::always_inline]]
inline void error(const char* str) {
  throw std::runtime_error(str);
}
[[noreturn]]
[[gnu::always_inline]]
inline void error(const std::string& str) {
  throw std::runtime_error(str);
}
template <typename... T>
[[gnu::always_inline]]
[[noreturn]]
inline void error(const T&... args) {
  throw std::runtime_error(cat(args...));
}

}

#endif

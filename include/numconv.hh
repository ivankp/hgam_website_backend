#ifndef IVAN_NUMCONV_HH
#define IVAN_NUMCONV_HH

#include <string_view>
#include <type_traits>

#if __has_include(<charconv>)
#  include <charconv>
#  define IVAN_CHARCONV_EXISTS
#else
#  include <string>
#endif

namespace ivan {

#if __cplusplus < 202002L
template <typename T, typename=void>
#else
template <typename T>
#endif
struct xtos { };

template <typename T>
#if __cplusplus < 202002L
struct xtos<T,std::enable_if_t<std::is_arithmetic_v<T>>>
#else
requires std::is_arithmetic_v<T>
struct xtos<T>
#endif
{
#ifdef IVAN_CHARCONV_EXISTS
  unsigned char n;
  char s[
    std::is_integral_v<T>
    ? ((616*sizeof(T)) >> 8) + 1 + std::is_signed_v<T>
    : std::is_same_v<T,double> ? 24 : 16
  ];
  xtos(T x) noexcept
  : n(std::to_chars(s,s+sizeof(s),x).ptr-s) { }
  explicit operator std::string_view() const noexcept { return { s, n }; }
#else
  std::string s;
  xtos(T x): s(std::to_string(x)) { }
  explicit operator std::string_view() const noexcept { return { s }; }
#endif
};

template <>
struct xtos<bool> {
  bool x;
  constexpr xtos(bool x) noexcept: x(x) { }
  explicit constexpr operator std::string_view() const noexcept
  { return x ? "true" : "false"; }
};

template <typename T>
xtos(T x) -> xtos<T>;

}

// ------------------------------------------------------------------

#include "error.hh"

namespace ivan {

template <typename T>
#if __cplusplus < 202002L
std::enable_if_t<std::is_arithmetic_v<T>,T>
#else
requires std::is_arithmetic_v<T> T
#endif

stox(std::string_view s) {

#ifdef IVAN_CHARCONV_EXISTS
  T x;
  const char* end = s.data() + s.size();
  const auto [p,e] = std::from_chars(s.data(),end,x);
  switch (e) {
    case std::errc::invalid_argument:
      error("invalid value: \"",s,'\"');
    case std::errc::result_out_of_range:
      error("value out-of-range: \"",s,'\"');
    default: ;
  }
  if (const auto r = end-p)
    error(r," bytes unconverted: \"",s,'\"');
  return x;
#else
  try {
    if        constexpr (std::is_same_v<T,int>) {
      return std::stoi  (std::string(s));
    } else if constexpr (std::is_same_v<T,long>) {
      return std::stol  (std::string(s));
    } else if constexpr (std::is_same_v<T,long long>) {
      return std::stoll (std::string(s));
    } else if constexpr (std::is_same_v<T,unsigned long>) {
      return std::stoul (std::string(s));
    } else if constexpr (std::is_same_v<T,unsigned long long>) {
      return std::stoull(std::string(s));
    } else if constexpr (std::is_same_v<T,float>) {
      return std::stof  (std::string(s));
    } else if constexpr (std::is_same_v<T,double>) {
      return std::stod  (std::string(s));
    } else if constexpr (std::is_same_v<T,long double>) {
      return std::stold (std::string(s));
    }
  } catch (const std::invalid_argument&) {
    error("invalid value: \"",s,'\"');
  } catch (const std::out_of_range&) {
    error("value out-of-range: \"",s,'\"');
  }
#endif

}

}

#endif

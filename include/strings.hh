#ifndef IVAN_STRINGS_HH
#define IVAN_STRINGS_HH

#include <string>
#include <cstring>
#include <string_view>
#include <type_traits>

namespace ivan {

[[nodiscard]]
inline const char* cstr(const char* x) noexcept { return x; }
[[nodiscard]]
inline const char* cstr(const std::string& x) noexcept { return x.c_str(); }

namespace impl {

template <typename>
struct failed_string_view_conversion;

template <typename T>
[[nodiscard]]
[[gnu::always_inline]]
inline auto to_string_view(const T& x) noexcept {
  if constexpr (std::is_same_v<T,char>)
    return std::string_view(&x,1);
#ifdef IVAN_NUMCONV_HH
  else if constexpr (std::is_constructible_v<xtos<T>,const T&>)
    return xtos<T>(x);
#endif
  else if constexpr (std::is_convertible_v<T,std::string_view>)
    return std::string_view(x);
  else
    return failed_string_view_conversion<T>{};
}

template <typename T>
inline constexpr bool is_sv = std::is_same_v<T,std::string_view>;

}

[[nodiscard]]
inline std::string cat() noexcept { return { }; }
[[nodiscard]]
inline std::string cat(std::string x) noexcept { return x; }
[[nodiscard]]
inline std::string cat(const char* x) noexcept { return x; }
[[nodiscard]]
inline std::string cat(char x) noexcept { return std::string(1,x); }
[[nodiscard]]
inline std::string cat(std::string_view x) noexcept { return std::string(x); }

template <typename... T>
[[nodiscard]]
[[gnu::always_inline]]
inline std::enable_if_t<
  ( (sizeof...(T) > 1) && (impl::is_sv<T> && ...) ),
  std::string
> cat(T... x) {
  std::string s((x.size() + ...),{});
  char* p = s.data();
  ( ( memcpy(p, x.data(), x.size()), p += x.size() ), ...);
  return s;
}

template <typename... T>
[[nodiscard]]
[[gnu::always_inline]]
inline std::enable_if_t<
  ( (sizeof...(T) > 1) && !(impl::is_sv<T> && ...) ),
  std::string
> cat(const T&... x) {
  return cat(static_cast<std::string_view>(impl::to_string_view(x))...);
}

inline bool starts_with(const char* a, const char* b) noexcept {
  for (;;) {
    const char c = *b;
    if (!c) return true;
    if (*a != c) return false;
    ++a; ++b;
  }
}
inline bool ends_with(const char* a, const char* b) noexcept {
  ssize_t i = ssize_t(strlen(a)) - ssize_t(strlen(b));
  if (i < 0) return false;
  return !strcmp(a+i,b);
}

inline bool starts_with(const char* a, const char* b, size_t na) noexcept {
  for (;;) {
    if (!na) return false;
    const char c = *b;
    if (!c) return true;
    if (*a != c) return false;
    ++a; ++b; --na;
  }
}

} // end namespace ivan

#endif

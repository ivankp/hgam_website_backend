#ifndef IVAN_PULL_HH
#define IVAN_PULL_HH

#include <cstdlib>
#include <type_traits>
#include <tuple>

// pool<0> uses calloc, allocated memory is zeroed out
// pool<1> uses malloc, allocated memory is NOT initialized

namespace ivan {

template <typename T>
struct type_constant { using type = T; };

template <typename T = double>
struct cnt_of {
  using type = T;
  const size_t n;
  cnt_of(size_t n) noexcept: n(n) { }
  operator size_t() const noexcept { return n; }
};

template <typename T>
struct cnt_of_trait {
  using type = double;
  static constexpr bool is = false;
  constexpr operator bool() const noexcept { return is; }
};
template <typename T>
struct cnt_of_trait<cnt_of<T>> {
  using type = T;
  static constexpr bool is = true;
  constexpr operator bool() const noexcept { return is; }
};
template <typename T>
using cnt_of_t = typename cnt_of_trait<T>::type;

template <typename... T>
struct pool_pointers {
  using tuple = std::tuple<T*...>;
  const tuple m;
  ~pool_pointers() {
    if constexpr (sizeof...(T) > 0)
      free(std::get<0>(m));
  }
  template <size_t I>
  friend std::tuple_element_t<I,tuple> get(const pool_pointers& x) {
    return std::get<I>(x.m);
  }
};

template <bool use_malloc, typename... T>
auto pool(T... n) -> std::enable_if_t<
  !(cnt_of_trait<T>{} && ...),
  pool_pointers<cnt_of_t<T>...>
> {
  return pool<use_malloc,cnt_of_t<T>...>( cnt_of<cnt_of_t<T>>{n}... );
}

template <bool use_malloc, typename... T>
auto pool(cnt_of<T>... n) -> pool_pointers<T...> {
  static constexpr size_t N = sizeof...(T);
  size_t len[N+1] { 0, (sizeof(T)*n)... };
  size_t i = 0;
  ([&](auto t){
    if (i) {
      using type = typename decltype(t)::type;
      size_t r = len[i] % alignof(type);
      if (r) r = alignof(type) - r;
      len[i+1] += (len[i] += r);
    }
    ++i;
  }(type_constant<T>{}),...);
  char* m;
  if constexpr (use_malloc)
    m = static_cast<char*>(malloc(len[N]));
  else
    m = static_cast<char*>(calloc(len[N],1));
  i = 0;
  return {{ reinterpret_cast<T*>(m+len[i++])... }};
}

}

namespace std {
  template <typename... T>
  struct tuple_size<ivan::pool_pointers<T...>>
  : std::integral_constant<size_t,sizeof...(T)> { };

  template <size_t I, typename... T>
  struct tuple_element<I,ivan::pool_pointers<T...>>
  : std::tuple_element<I,std::tuple<T*...>> { };
}

#endif

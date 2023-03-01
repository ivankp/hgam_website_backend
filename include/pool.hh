#ifndef IVAN_PULL_HH
#define IVAN_PULL_HH

#include <cstdlib>
#include <type_traits>
#include <tuple>

namespace ivan {

template <typename T, size_t N>
struct pool_array {
  T* const m[N];
  ~pool_array() { free(*m); }
  template <size_t I>
  friend T* get(const pool_array& x) { return x.m[I]; }
};

template <typename T, bool nzero=true, typename... N>
auto pool(N... n) -> pool_array<T,sizeof...(n)> {
  T* m;
  if constexpr (nzero)
    m = static_cast<T*>(malloc( (size_t(n) + ...)*sizeof(T) ));
  else
    m = static_cast<T*>(calloc( (size_t(n) + ...),sizeof(T) ));
  return {{ ( m+=size_t(n), m-size_t(n) ) ... }};
}

}

namespace std {
  template <typename T, size_t N>
  struct tuple_size<ivan::pool_array<T,N>>: std::integral_constant<size_t,N> { };

  template <size_t I, typename T, size_t N>
  struct tuple_element<I,ivan::pool_array<T,N>> { using type = T*; };
}

#endif

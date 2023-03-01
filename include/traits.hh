#ifndef IVAN_TRAITS_HH
#define IVAN_TRAITS_HH

namespace ivan {

template <typename T, typename...> struct pack_head { using type = T; };
template <typename... T> using head_t = typename pack_head<T...>::type;

template <typename T>
struct type_constant { using type = T; };

template <typename...>
struct type_sequence { };

}

#endif

#ifndef IVAN_PROGRAM_OPTIONS_HH
#define IVAN_PROGRAM_OPTIONS_HH

#include <utility> // declval
#include <tuple>
#include "numconv.hh"

namespace ivan {
namespace po {

template <typename>
class ref;

// ref specializations ----------------------------------------------
template <>
class ref<bool> {
  bool& x;
public:
  explicit ref(bool& x): x(x) { }
  void operator()() { x = !x; }
};

template <typename T>
requires std::is_assignable_v<T&,const char*>
class ref<T> {
  T& x;
public:
  explicit ref(T& x): x(x) { }
  void operator()(const char* val) { x = val; }
};

template <typename T>
requires std::is_arithmetic_v<T>
class ref<T> {
  T& x;
public:
  explicit ref(T& x): x(x) { }
  void operator()(const char* val) { x = stox<T>(val); }
};

// ref concepts -----------------------------------------------------
template <typename T>
concept Ref = requires { ref<T>(std::declval<T&>()); };
template <typename T>
concept CRef = Ref<const T>;

template <typename F>
concept OptFcnNoArg = std::is_invocable_v<F>;
template <typename F>
concept OptFcnArg = std::is_invocable_v<F,char*>;
template <typename F>
concept OptFcn = OptFcnNoArg<F> || OptFcnArg<F>;

// factory ----------------------------------------------------------
template <typename T, typename... Args>
class factory {
  T& x;
  std::tuple<Args&&...> args;
public:
  explicit factory(T& x, Args&&... args)
  : x(x), args(std::forward<Args>(args)...) { }
  void operator()() { x = std::make_from_tuple<T>(args); }
};

// opt --------------------------------------------------------------
template <OptFcn F>
struct opt {
  const char* s;
  F f;

  static constexpr bool no_arg = OptFcnNoArg<F>;
  static constexpr bool with_arg = OptFcnArg<F>;

  template <typename T>
  requires std::is_same_v<
    std::remove_cvref_t<F>,
    std::remove_cvref_t<T> >
  opt(const char* s, T&& f): s(s), f(std::forward<T>(f)) { }

  template <CRef T>
  requires std::is_same_v<F,ref<const T>>
  opt(const char* s, const T& x): s(s), f(ref<const T>(x)) { }

  template <Ref T>
  requires std::is_same_v<F,ref<T>>
  opt(const char* s, T& x): s(s), f(ref<T>(x)) { }

  template <typename T, typename... Args>
  requires std::is_same_v<F,factory<T,Args...>>
  opt(const char* s, T& x, Args&&... args)
  : s(s), f(factory<T,Args...>(x,std::forward<Args>(args)...)) { }
};

template <OptFcn F>
opt(const char*, F&&) -> opt<F>;

template <CRef T>
opt(const char*, const T&) -> opt<ref<const T>>;

template <Ref T>
opt(const char*, T&) -> opt<ref<T>>;

template <typename T, typename... Args>
opt(const char*, T&, Args&&...) -> opt<factory<T,Args...>>;

// consumer ---------------------------------------------------------
struct consumer_struct {
  template <typename F>
  static bool wrap(void* f, char* val) {
    return (*reinterpret_cast<F*>(f))(val);
  }
  bool(*w)(void*,char*) = nullptr;
  void* f;
  const char* s;
  operator bool() const noexcept { return w; }
  consumer_struct() = default;
  template <typename F>
  consumer_struct(opt<F>& o): w(&wrap<F>), f(&o.f), s(o.s) { }
  bool operator()(char* val) const { return w(f,val); }
};

} // end namespace po

// NOTE: short option must not be '-', long option must not contain =

// program_options function -----------------------------------------
template <typename... F>
void program_options(int argc, char** argv, po::opt<F>&&... opts) {
  po::consumer_struct consumer;
  for (int i=0; i<argc; ++i) {
    char* arg = argv[i];
    if (*arg == '-') { // option
      if (consumer) error(
        "option -",(consumer.s[1] ? "-" : ""),consumer.s,
        " expected more values"
      );
      ++arg;
      if (*arg == '-') { // long option (--)
        ++arg;
        if (*arg=='\0') { // --
          // TODO: find a way to elegantly handle --
          error("unexpected option --");
        } else {
          if (!([&](auto& o){
            const size_t n = strcspn(arg,"=");
            const bool eq = arg[n];
            if (memcmp(arg,o.s,n)) // skip other opts
              return false;
            if constexpr(!o.with_arg) {
              if (eq) error("option --",o.s," does not accept a value");
              o.f();
            } else if (eq) { // provided value
              arg += n+1;
              goto use_arg;
            } else if (i+1 < argc && argv[i+1][0]!='-') {
              arg = argv[++i];
use_arg:
              using ret_t = decltype(o.f(arg));
              if constexpr (std::is_same_v<ret_t,bool>) {
                if (o.f(arg)) consumer = o;
              } else o.f(arg);
            } else {
              if constexpr (o.no_arg) o.f();
              else error("option --",o.s," requires a value");
            }
            return true;
          }(opts) || ... )) error("unexpected option --",arg);
        }
      } else { // short option (-)
        if (*arg=='\0') { // -
          // TODO: find a way to elegantly handle -
          error("unexpected option -");
        } else {
          for (;;) {
            const char c = *arg++;
            if (!([&](auto& o){
              if (!(o.s[0]==c && o.s[1]=='\0')) // skip other opts
                return false;
              if constexpr(!o.with_arg) o.f();
              else if (*arg!='\0') { // provided value
                goto use_arg;
              } else if (i+1 < argc && argv[i+1][0]!='-') {
                arg = argv[++i];
use_arg:
                using ret_t = decltype(o.f(arg));
                if constexpr (std::is_same_v<ret_t,bool>) {
                  if (o.f(arg)) consumer = o;
                } else o.f(arg);
                arg = nullptr;
              } else {
                if constexpr (o.no_arg) o.f();
                else error("option -",c," requires a value");
              }
              return true;
            }(opts) || ... )) error("unexpected option -",c);
            if (!arg || !*arg) break;
          }
        }
      }
    } else { // not an option
      if (consumer) {
        if (!consumer(arg)) consumer = { };
      } else {
        std::cout << "free argument: \"" << arg << '"' << std::endl;
      }
    }
  }
}

} // end namespace ivan

#endif

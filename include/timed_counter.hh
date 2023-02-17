#ifndef IVAN_TIMED_COUNTER_HH
#define IVAN_TIMED_COUNTER_HH

#include <chrono>
#include <sstream>
#include <locale>
#include <cstdio>

namespace ivan {

template <typename Counter = unsigned>
class timed_counter {
public:
  using value_type   = Counter;
  using clock_type   = std::chrono::steady_clock;
  using time_type    = std::chrono::time_point<clock_type>;
  using sec_type     = std::chrono::duration<double>;

private:
  value_type cnt{}, cnt_start{}, cnt_end{};
  time_type t_start = clock_type::now(), t_last = t_start;
  std::stringstream ss;
  unsigned width[3]{};

public:
  void print(bool check=true) {
    const auto now = clock_type::now();
    if (check && sec_type(now-t_last).count() < 1) return;
    t_last = now;
    const auto dt = sec_type(now-t_start).count();
    const int hours   = dt/3600;
    const int minutes = (dt-hours*3600)/60;
    const int seconds = (dt-hours*3600-minutes*60);

    ss.str({});
    ss.width(width[0]);
    ss << cnt << " | ";
    ss.width(width[1]);
    ss << ( (100.*(cnt-cnt_start))/(cnt_end-cnt_start) ) << "% | ";
    bool only_sec = false;
    if (hours)   goto prt_hours;
    if (minutes) goto prt_min;
    if (seconds) {
      only_sec = true;
      goto prt_sec;
    }
    ss << int(dt*1000) << "ms";
    goto prt_ms;
prt_hours:
    ss.width(5);
    ss << hours << ':';
prt_min:
    ss.width(2);
    ss << minutes << ':';
    ss.fill('0');
prt_sec:
    ss.width(2);
    ss << seconds;
    if (only_sec) ss << 's';
    else ss.fill(' ');
prt_ms: ;

    const auto& s = ss.str();
    printf(ss.str().c_str());
    for (int i=width[2]-s.size(); i>0; --i) putchar(' ');
    putchar('\r');
    fflush(stdout);
    width[2] = s.size();
  }

private:
  void init() {
    ss.imbue(std::locale(""));
    ss.precision(2);
    ss << std::right << std::fixed;
    ss << cnt_start;
    unsigned len = ss.str().size();
    if (len > width[0]) width[0] = len;
    ss.str({});
    ss << cnt_end;
    len = ss.str().size();
    if (len > width[0]) width[0] = len;
    width[0] += 2;
    ss.str({});
    ss << 100.f;
    len = ss.str().size();
    if (len > width[1]) width[1] = len;

    print(false);
  }

public:
  timed_counter() noexcept = default;
  timed_counter(value_type i, value_type n) noexcept
  : cnt(i), cnt_start(i), cnt_end(n) { init(); }
  timed_counter(value_type n) noexcept
  : cnt{}, cnt_start{}, cnt_end(n) { init(); }
  ~timed_counter() {
    print(false);
    puts("\n");
  }

  void reset(value_type i, value_type n) noexcept {
    cnt = i;
    cnt_start = i;
    cnt_end = n;
    t_start = clock_type::now();
    t_last = t_start;
  }
  void reset(value_type n) noexcept { reset({},n); }

  bool done() const noexcept { return !(cnt < cnt_end); }
  bool operator!() const noexcept { return done(); }
  operator bool() const noexcept { return !done(); }

  operator value_type() const noexcept { return cnt; }
  value_type operator*() const noexcept { return cnt; }

  // prefix
  value_type operator++() { print(); return ++cnt; }
  value_type operator--() { print(); return --cnt; }

  // postfix
  value_type operator++(int) { print(); return cnt++; }
  value_type operator--(int) { print(); return cnt--; }

  template <typename T>
  value_type operator+=(T i) { print(); return cnt += i; }
  template <typename T>
  value_type operator-=(T i) { print(); return cnt -= i; }

  template <typename T>
  bool operator==(T i) const noexcept { return cnt == i; }
  template <typename T>
  bool operator!=(T i) const noexcept { return cnt != i; }
  template <typename T>
  bool operator< (T i) const noexcept { return cnt <  i; }
  template <typename T>
  bool operator<=(T i) const noexcept { return cnt <= i; }
  template <typename T>
  bool operator> (T i) const noexcept { return cnt >  i; }
  template <typename T>
  bool operator>=(T i) const noexcept { return cnt >= i; }
};

template <typename T> timed_counter(T end) -> timed_counter<T>;
template <typename T> timed_counter(T start, T end) -> timed_counter<T>;

} // end namespace

#endif

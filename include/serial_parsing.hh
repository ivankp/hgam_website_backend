#ifndef IVAN_SERIAL_PARSING_HH
#define IVAN_SERIAL_PARSING_HH

#include <tuple>
#include <cstring>

namespace ivan {

template <typename F>
void serial_split(char* a, char d, F&& f) {
  for (;;) {
    char* b = strchr(a,d);
    if (b) *b = '\0';

    f(a,b);

    if (!b) break;
    a = b+1;
  }
}

}

#endif

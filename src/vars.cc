#include <iostream>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <cmath>
#include "strings.hh"

#include "debug.hh"

using std::cout;
using ivan::cat; //, ivan::error, ivan::cnt_of;

std::string_view lstrip(std::string_view s) noexcept {
  s.remove_prefix(s.find_first_not_of(" \t\r"));
  return s;
}
std::string_view rstrip(std::string_view s) noexcept {
  s.remove_suffix(s.size()-(s.find_last_not_of(" \t\r")+1));
  return s;
}
std::string_view strip(std::string_view s) noexcept {
  return lstrip(rstrip(s));
}

int main(int argc, char** argv) {
  std::ios_base::sync_with_stdio(false);
  if (argc!=2) {
    std::cerr << "usage: " << argv[0] << " vars.txt\n";
    return 1;
  }

  std::ifstream f(argv[1]);
  std::string_view name;
  bool first_var = true;

  cout << '[';
  for (std::string line; std::getline(f,line); ) {
    if (line.empty()) continue;

    auto line_tail = lstrip(line);
    if (line_tail.empty() || line_tail[0]=='#') continue;

    if (line_tail.size() == line.size()) { // var name
      if (!name.empty()) cout << ']';
      name = rstrip(line);
      if (!(
        std::filesystem::exists(cat("data/",name,"_data.dat")) &&
        std::filesystem::exists(cat("data/",name,"_mc.dat"))
      )) {
        name = { };
        continue;
      }
      if (first_var) first_var = false;
      else cout << ',';
      cout << "[\"" << name << "\"";
    } else {
      if (name.empty()) continue;
      auto col = line_tail.find_first_of(':');
      cout << ",[\"";
      if (col != std::string_view::npos)
        cout << rstrip(line_tail.substr(0,col));
      cout << "\"";

      std::stringstream ss(std::string(line_tail.substr(col+1)));
      for (std::string s; ss >> s; ) {
        double x = std::stod(s);
        cout << ',';
        if (std::isinf(x)) cout << (x < 0 ? "\"-inf\"" : "\"inf\"");
        else cout << x;
      }

      cout << ']';
    }

  }
  cout << ']';
}

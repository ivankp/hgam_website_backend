#include <iostream>
#include <iomanip>
#include <vector>
#include <array>
#include <algorithm>
#include <cstring>

#include <TFile.h>
#include <TTree.h>
#include <TLeaf.h>

using std::cout, std::cerr, std::endl;

#define ERROR(...) ivan::error(__FILE__ ":" STR(__LINE__) ": ",__VA_ARGS__)

const char* const var_pref = "HGamEventInfoAuxDyn.";
const size_t var_pref_len = strlen(var_pref);

int main(int argc, char** argv) {
  if (argc != 2) {
    cout << "usage: " << argv[0] << " mxaod.root\n";
    return 1;
  }

  TFile file(argv[1]);
  TTree* tree = dynamic_cast<TTree*>(file.Get("CollectionTree"));

  std::vector<std::array<std::string_view,2>> branches;

  for (auto* _branch : *tree->GetListOfBranches()) {
    TBranch* branch = dynamic_cast<TBranch*>(_branch);
    const char* name = branch->GetName();

    if (name && !strncmp(var_pref, name, var_pref_len)) {
      TLeaf* leaf = static_cast<TLeaf*>(branch->GetListOfLeaves()->First());
      const char* type = leaf->GetTypeName();

      if (type)
        branches.push_back({ name+var_pref_len, type });
    }
  }

  { using namespace std;
    sort(
      begin(branches), end(branches),
      [](const auto& a, const auto& b) {
        return lexicographical_compare(
          begin(a[0]), end(a[0]),
          begin(b[0]), end(b[0]),
          [](char a, char b) { return tolower(a) < tolower(b); }
        );
      }
    );
  }

  size_t ntype = 0;
  for (const auto& [name,type] : branches) {
    if (ntype < type.size()) ntype = type.size();
  }

  // cout << '\n';
  // for (const auto& [name,type] : branches) {
  //   cout << "output_file<" << std::setw(ntype) << type << ">"
  //     " f_" << name << "(\"" << name << "\");\n";
  // }

  // cout << '\n';
  // for (const auto& [name,type] : branches) {
  //   cout << "branch_reader<" << std::setw(ntype) << type << ">"
  //     " b_" << name << "(reader, VAR_PREF \"" << name << "\");\n";
  // }

  // cout << '\n';
  // for (const auto& [name,type] : branches) {
  //   cout << "const " << type << ' ' << name << " = *b_" << name;
  //   cout << ";\n"
  //     "f_" << name << ".write(" << name << ");\n\n";
  // }

  cout << '\n';
  for (const auto& [name,type] : branches) {
    cout << "VAR(" << std::setw(ntype) << type << ","
      << name << ") \\\n";
  }
}

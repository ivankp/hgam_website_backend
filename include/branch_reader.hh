// ------------------------------------------------------------------
// Class template for reading TTree branches of any of specified types
// Works similarly to std::variant
// Writted by Ivan Pogrebnyak
// ------------------------------------------------------------------

#ifndef BRANCH_READER_HH
#define BRANCH_READER_HH

#include <tuple>
#include <algorithm>

#include <TLeaf.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

#include "traits.hh"
#include "error.hh"

namespace ivan {

template <typename>
constexpr const char* root_type_str();
template <>
constexpr const char* root_type_str<Char_t>() { return "Char_t"; }
template <>
constexpr const char* root_type_str<UChar_t>() { return "UChar_t"; }
template <>
constexpr const char* root_type_str<Short_t>() { return "Short_t"; }
template <>
constexpr const char* root_type_str<UShort_t>() { return "UShort_t"; }
template <>
constexpr const char* root_type_str<Int_t>() { return "Int_t"; }
template <>
constexpr const char* root_type_str<UInt_t>() { return "UInt_t"; }
template <>
constexpr const char* root_type_str<Float_t>() { return "Float_t"; }
template <>
constexpr const char* root_type_str<Double_t>() { return "Double_t"; }
template <>
constexpr const char* root_type_str<Long64_t>() { return "Long64_t"; }
template <>
constexpr const char* root_type_str<ULong64_t>() { return "ULong64_t"; }
template <>
constexpr const char* root_type_str<Bool_t>() { return "Bool_t"; }

template <typename... Ts>
class branch_reader {
public:
  using value_type = std::common_type_t<std::remove_extent_t<Ts>...>;

  static constexpr bool is_array = std::is_array<head_t<Ts...>>::value;

private:
  size_t index;

  template <size_t I>
  using type = typename std::tuple_element<I,std::tuple<Ts...>>::type;

  template <typename T>
  using reader_type = std::conditional_t< is_array,
    TTreeReaderArray<std::remove_extent_t<T>>,
    TTreeReaderValue<std::remove_extent_t<T>>
  >;

  char data[ std::max({sizeof(reader_type<Ts>)...}) ];

  template <size_t I=0>
  std::enable_if_t<(I<sizeof...(Ts)),size_t>
  get_index(const char* type_name) {
    if (!strcmp(root_type_str<std::remove_extent_t<type<I>>>(),type_name))
      return I;
    else return get_index<I+1>(type_name);
  }
  template <size_t I=0>
  std::enable_if_t<(I==sizeof...(Ts)),size_t>
  get_index [[noreturn]] (const char* type_name) {
    error("this branch_reader cannot read ",type_name);
  }

  template <size_t I>
  auto cast() { return reinterpret_cast<reader_type<type<I>>*>(data); }

  template <typename F, size_t I=0>
  std::enable_if_t<(sizeof...(Ts)-I>1)> call(F&& f) {
    if (index==I) f(cast<I>());
    else call<F,I+1>(std::forward<F>(f));
  }
  template <typename F, size_t I=0>
  std::enable_if_t<(sizeof...(Ts)-I==1)> call(F&& f) {
    f(cast<I>());
  }

  static auto* get_leaf(TTree* tree, const char* name) {
    auto* leaf = tree->GetLeaf(name);
    if (leaf) return leaf;
    error("no leaf \"",name,"\" in tree \"",tree->GetName(),'\"');
  }

public:
  branch_reader(TTreeReader& reader, const char* branch_name)
  : index(get_index(
      get_leaf(reader.GetTree(),branch_name)->GetTypeName()
    ))
  {
    call([&](auto* p) {
      using T = std::decay_t<decltype(*p)>;
      new(p) T(reader,branch_name);
    });
  }
  branch_reader(TTreeReader& reader, const std::string& branch_name)
  : branch_reader(reader,branch_name.c_str()) { }
  ~branch_reader() {
    call([](auto* p) {
      using T = std::decay_t<decltype(*p)>;
      p->~T();
    });
  }

  value_type operator*() {
    value_type x;
    call([&](auto* p){ x = **p; });
    return x;
  }

  value_type operator[](size_t i) {
    value_type x;
    call([&,i](auto* p){ x = (*p)[i]; });
    return x;
  }

  const char* GetBranchName() {
    const char* x;
    call([&](auto* p){ x = p->GetBranchName(); });
    return x;
  }
};

template <typename T>
class branch_reader<T> {
public:
  using value_type = std::remove_extent_t<T>;

  static constexpr bool is_array = std::is_array<T>::value;

private:
  using impl_type = std::conditional_t< is_array,
    TTreeReaderArray<value_type>,
    TTreeReaderValue<value_type>
  >;
  impl_type impl;

public:
  branch_reader(TTreeReader& reader, const char* branch_name)
  : impl(reader,branch_name) { }
  branch_reader(TTreeReader& reader, const std::string& branch_name)
  : branch_reader(reader,branch_name.c_str()) { }

  const char* GetBranchName() const noexcept { return impl.GetBranchName(); }

  value_type& operator*() { return *impl; }

  template <bool A = is_array, typename V = value_type>
  std::enable_if_t<A,V&>
  operator[](size_t i) { return impl[i]; }

  template <bool A = is_array, typename V = value_type>
  std::enable_if_t<!A, decltype(std::declval<V&>()[0])>
  operator[](size_t i) { return (*impl)[i]; }
};

}

#endif

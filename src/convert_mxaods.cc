#include <iostream>
#include <filesystem>
#include <optional>
#include <regex>
#include <cstdlib>
#include <cerrno>

#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <TFile.h>
#include <TH1.h>
#include <TKey.h>
#include <TLorentzVector.h>

#include "numconv.hh"
#include "strings.hh"
#include "timed_counter.hh"
#include "branch_reader.hh"
#include "debug.hh"

using std::cout, std::cerr, std::endl;
using ivan::cat, ivan::branch_reader;

static_assert( std::is_same_v<float,Float_t> );

#define ERROR(...) ivan::error(__FILE__ ":" STR(__LINE__) ": ",__VA_ARGS__)

constexpr float fnan = std::numeric_limits<float>::quiet_NaN();

std::filesystem::path dir;
bool mc = false;
uint64_t nevents = 0;
float lumi = 0;
double mc_factor = 0;

template <typename T> constexpr char type_byte = '\0';
// https://docs.python.org/3/library/struct.html#format-characters
template <> constexpr char type_byte<float   > = 'f';
template <> constexpr char type_byte<char    > = 'c';
template <> constexpr char type_byte< uint8_t> = 'B';
template <> constexpr char type_byte< int32_t> = 'i';
template <> constexpr char type_byte<uint32_t> = 'I';
template <> constexpr char type_byte< int64_t> = 'q';
template <> constexpr char type_byte<uint64_t> = 'Q';

template <typename To, typename From>
To ceiling_cast(From x) noexcept {
  static constexpr To max = std::numeric_limits<To>::max();
  return x < max ? To(x) : max;
}

template <typename T>
char* bufcpy(char* m, T x) noexcept {
  memcpy(m,&x,sizeof(x));
  return m + sizeof(x);
}

template <typename T>
struct output_file {
  std::filesystem::path path;
  int fd = -1;
  char* buf = nullptr;
  char* m = nullptr;
  off_t skipped = 0;
  bool absent_var = false;

  static constexpr size_t buflen = 1 << 20;

  output_file() { }
  void open(std::string_view name) {
    path = dir / cat(name, mc ? "_mc" : "_data", ".dat");
    fd = ::open(path.c_str(), O_CREAT|O_WRONLY|O_TRUNC, 0644);
    if (fd < 0) ERROR('\"',path.native(),"\": open(): ",strerror(errno));
    buf = static_cast<char*>( ::malloc(buflen) );
    m = buf;

    const size_t skip = 2ul + (!mc ? sizeof(lumi) : 0ul) + sizeof(nevents);
    const size_t nulls = 8-((name.size()+skip)%8);
    const size_t bufneed = name.size()+nulls+skip;
    if (buflen < bufneed)
      ERROR('\"',path.native(),"\": buffer length ",buflen," , need ",bufneed);

    cout << path.native() << '\n';

    // write header
    memcpy(m,name.data(),name.size()); m += name.size();
    char tmp[8] { };
    memcpy(m,tmp,nulls); m += nulls;
    m[0] = mc ? 'm' : 'd';
    static_assert(type_byte<T> != '\0');
    m[1] = type_byte<T>;

    skipped = (m+2) - buf;
    m += skip;
  }
  void remove() {
    if (fd >= 0) {
      cerr << "\033[31m" "removing output file " << path << "\033[0m" << endl;
      ::close(fd);
      fd = -1;
      ::free(buf);
      buf = nullptr;
      m = nullptr;
      std::filesystem::remove(path);
      path.clear();
    }
  }
  operator bool() const noexcept {
    return fd >= 0;
  }
  ~output_file() {
    if (fd >= 0) {
      flush();
      if (!mc) m = bufcpy(m,lumi);
      m = bufcpy(m,nevents);
      if ((
        ::lseek(fd,skipped,SEEK_SET)
      ) < 0) ERROR('\"',path.native(),"\": lseek(): ",strerror(errno));
      if ((
        ::write(fd,buf,m-buf)
      ) < 0) ERROR('\"',path.native(),"\": write(): ",strerror(errno));
      ::close(fd);
      ::free(buf);
    }
  }

  void write(const char* s, size_t n) {
    if (n < buflen) {
      if (buf + buflen < m + n) flush();
      memcpy(m,s,n); m += n;
    } else {
      flush();
      if ((
        ::write(fd,s,n)
      ) < 0) ERROR('\"',path.native(),"\": write(): ",strerror(errno));
    }
  }
  template <typename U>
  void write(U x) requires(std::is_same_v<U,T>) {
    static_assert( sizeof(x) < buflen );
    if (buf + buflen < m + sizeof(x)) flush();
    m = bufcpy(m,x);
  }
  void flush() {
    if (m != buf) {
      if ((
        ::write(fd,buf,m-buf)
      ) < 0) ERROR('\"',path.native(),"\": write(): ",strerror(errno));
      m = buf;
    }
  }
};

int main(int argc, char** argv) {
  std::ios_base::sync_with_stdio(false);
  if (argc < 3) {
    cout << "usage: " << argv[0] << " output_dir input.root ...\n";
    return 1;
  }

try {
  cout << "Input files:\n";
  for (int i=2; i<argc; ++i) {
    cout << argv[i] << '\n';

    // Luminosity
    static const std::regex re(R"((\d+(?:\.\d*)?)i([fp])b\b)");
    std::cmatch match;
    if (std::regex_search(argv[i],match,re)) {
      if (i>2 &&  mc) ERROR("data file after mc files");
      mc = false;
      double this_lumi = atof(match.str(1).c_str());
      const char unit = match.str(2)[0];
      switch (unit) {
        case 'p': this_lumi *= 1e-3; break;
        case 'f': break;
        default: ERROR("unexpected lumi units \"i",unit,"b\"");
      }
      lumi += this_lumi;
    } else {
      if (i>2 && !mc) ERROR("mc file after data files");
      mc = true;
    }
  }
  cout << endl;

  if (mc) {
    cout << "Monte Carlo\n\n";
  } else {
    cout << "Data\nLumi = " << lumi << " ifb\n\n";
  }

  // Open output files ==============================================
  std::filesystem::create_directories(dir = argv[1]);
  cout << "Output files:\n";

  output_file<uint8_t> f_N_j_30;
  output_file<Float_t> f_m_yy;

  output_file<uint32_t> f_runNumber;
  output_file<uint64_t> f_eventNumber;

// https://en.wikipedia.org/wiki/X_Macro
#define VARS \
  VAR(Float_t,abs_Zepp) \
  VAR( Char_t,catXS_lepton) \
  VAR( Char_t,catXS_MET) \
  VAR( Char_t,catXS_nbjet) \
  VAR( Char_t,catXS_ttH) \
  VAR( Char_t,catXS_VBF) \
  VAR(Float_t,Dphi_j_j_30) \
  VAR(Float_t,Dphi_j_j_30_signed) \
  VAR(Float_t,Dphi_yy_jj_30) \
  VAR(Float_t,HT_30) \
  VAR(Float_t,m_jj_30) \
  VAR(Float_t,m_yyj_30) \
  VAR(Float_t,maxTau_yyj_30) \
  VAR(Float_t,pT_j1_30) \
  VAR(Float_t,pT_yy) \
  VAR(Float_t,pT_yy_JV_30) \
  VAR(Float_t,pT_yy_JV_40) \
  VAR(Float_t,pT_yy_JV_50) \
  VAR(Float_t,pT_yy_JV_60) \
  VAR(Float_t,pT_yyj_30) \
  VAR(Float_t,pT_yyjj_30) \
  VAR(Float_t,rel_DpT_y_y) \
  VAR(Float_t,rel_pT_y1) \
  VAR(Float_t,rel_pT_y2) \
  VAR(Float_t,rel_sumpT_y_y) \
  VAR(Float_t,sumTau_yyj_30) \
  VAR(Float_t,yAbs_yy)

#define WEIGHTS \
    VAR(weightJvt_30) \
    VAR(weightCatXS_lepton) \
    VAR(weightCatXS_nbjet) \
    VAR(weightCatXS_ttH) \
    VAR(weightCatXS_VBF)

#define VAR(TYPE,NAME) \
  output_file<TYPE> f_##NAME;

    VARS

#undef VAR

#define VAR(NAME) \
  output_file<Float_t> f_##NAME;

    WEIGHTS

#undef VAR

  cout << endl;

  for (int i=2; i<argc; ++i) { // loop over input files
    cout << argv[i] << endl;
    TFile file(argv[i]);

    if (mc) {
      for (auto* key : *file.GetListOfKeys()) {
        const char* name = key->GetName();
        if (!ivan::starts_with(name,"CutFlow_") ||
            !ivan::ends_with(name,"_noDalitz_weighted")) continue;
        TH1 *h = dynamic_cast<TH1*>(static_cast<TKey*>(key)->ReadObj());
        cout << name << endl;
        const double n_all = h->GetBinContent(3);
        cout << h->GetXaxis()->GetBinLabel(3) << " = " << n_all << endl;
        mc_factor = 1e3/n_all;
        break;
      }
    }

    TTreeReader reader("CollectionTree",&file);
    auto* const tree = reader.GetTree();

#define VAR_PREF "HGamEventInfoAuxDyn."
#define VAR_PREF_TRUTH "HGamTruthEventInfoAuxDyn."

    branch_reader<Char_t> b_isPassed(reader, VAR_PREF "isPassed");

    std::optional<branch_reader< Char_t>> b_isFiducial;
    std::optional<branch_reader<Float_t>> b_cs_br_fe;
    std::optional<branch_reader<Float_t>> b_weight;

    if (mc) {
      b_isFiducial.emplace(reader, VAR_PREF_TRUTH "isFiducial");
      b_cs_br_fe.emplace(reader, VAR_PREF "crossSectionBRfilterEff");
      b_weight.emplace(reader, VAR_PREF "weight");
    }

    branch_reader<Float_t> b_m_yy(reader, VAR_PREF "m_yy");
    f_m_yy.open("m_yy");

#define VAR(TYPE,NAME) \
  std::optional<branch_reader<TYPE>> b_##NAME, b_##NAME##_truth; \
  if (!f_##NAME.absent_var) { \
    auto& f = f_##NAME; \
    bool absent_branches = false; \
    const char* name[] = { VAR_PREF #NAME, VAR_PREF_TRUTH #NAME }; \
    for (int i=0; i<(mc?2:1); ++i) \
      if (!tree->FindBranch(name[i])) { \
        cerr << "\033[31m" "missing branch " << name[i] << "\033[0m" << endl; \
        absent_branches = true; \
      } \
    if (!absent_branches) { \
      b_##NAME.emplace(reader, name[0]); \
      if (mc) b_##NAME##_truth.emplace(reader, name[1]); \
      if (!f) f.open(#NAME); \
    } else { \
      if (f) f.remove(); \
      f.absent_var = true; \
    } \
  }

    VAR(Int_t,N_j_30)

    VARS

#undef VAR
#undef VARS

#define VAR(TYPE,NAME) \
  std::optional<branch_reader<TYPE>> b_##NAME; \
  if (!mc && !f_##NAME.absent_var) { \
    auto& f = f_##NAME; \
    const char* name = "EventInfoAux." #NAME; \
    if (tree->FindBranch(name)) { \
      b_##NAME.emplace(reader, name); \
      if (!f) f.open(#NAME); \
    } else { \
      cerr << "\033[31m" "missing branch " << name << "\033[0m" << endl; \
      if (f) f.remove(); \
      f.absent_var = true; \
    } \
  }

    static_assert( sizeof(UInt_t) == sizeof(uint32_t) );
    static_assert( sizeof(ULong64_t) == sizeof(uint64_t) );

    VAR(UInt_t,runNumber)
    VAR(ULong64_t,eventNumber)

#undef VAR

#define VAR(NAME) \
  std::optional<branch_reader<Float_t>> b_##NAME; \
  if (mc && !f_##NAME.absent_var) { \
    auto& f = f_##NAME; \
    const char* name = VAR_PREF #NAME; \
    if (tree->FindBranch(name)) { \
      b_##NAME.emplace(reader, name); \
      if (!f) f.open(#NAME); \
    } else { \
      cerr << "\033[31m" "missing branch " << name << "\033[0m" << endl; \
      if (f) f.remove(); \
      f.absent_var = true; \
    } \
  }

    WEIGHTS

#undef VAR

    // LOOP over events =============================================
    double weight;
    for (
      ivan::timed_counter<Long64_t> ent(reader.GetEntries(true));
      reader.Next();
      ++ent
    ) {
      // selection cut
      if (!*b_isPassed) continue;

      if (mc) {
        weight = **b_weight;
        if (weight==0) continue;
      }

      // diphoton mass cut
      const float m_yy = *b_m_yy*1e-3;
      if (m_yy < 105 || 160 < m_yy) continue;

      ++nevents;

      f_m_yy.write(m_yy);
      if (mc) {
        f_m_yy.write<float>( // weight
          weight * double(**b_cs_br_fe) * mc_factor
        );

#define VAR(NAME) \
  if (auto& f = f_##NAME) f.write( **b_##NAME );

        WEIGHTS

#undef VAR
      } else {
        if (f_runNumber  ) f_runNumber  .write<uint32_t>( **b_runNumber   );
        if (f_eventNumber) f_eventNumber.write<uint64_t>( **b_eventNumber );
      }
      // TODO: weightCatXS

      // ============================================================
      const uint8_t N_j_30 = ceiling_cast<uint8_t>(**b_N_j_30);
      f_N_j_30.write(N_j_30);
      if (mc) f_N_j_30.write(ceiling_cast<uint8_t>(**b_N_j_30_truth));

#define _1e3(X) (X)*1e-3

#define VAR(TYPE,NAME,F) \
  if (auto& f = f_##NAME) { \
    f.write<TYPE>(F(**b_##NAME)); \
    if (mc) f.write<TYPE>(F(**b_##NAME##_truth)); \
  }

      VAR(Float_t,abs_Zepp,)
      VAR( Char_t,catXS_lepton,)
      VAR( Char_t,catXS_MET,)
      VAR( Char_t,catXS_nbjet,)
      VAR( Char_t,catXS_ttH,)
      VAR( Char_t,catXS_VBF,)
      VAR(Float_t,yAbs_yy,)
      VAR(Float_t,HT_30,_1e3)
      VAR(Float_t,pT_yy,_1e3)
      VAR(Float_t,pT_yy_JV_30,_1e3)
      VAR(Float_t,pT_yy_JV_40,_1e3)
      VAR(Float_t,pT_yy_JV_50,_1e3)
      VAR(Float_t,pT_yy_JV_60,_1e3)
      VAR(Float_t,rel_DpT_y_y,_1e3)
      VAR(Float_t,rel_pT_y1,_1e3)
      VAR(Float_t,rel_pT_y2,_1e3)
      VAR(Float_t,rel_sumpT_y_y,_1e3)

#undef VAR

      // 1 jet ======================================================
      bool enough_jets = N_j_30 >= 1;

#define VAR(NAME,F) \
  if (auto& f = f_##NAME) { \
    f.write(enough_jets ? float(F(**b_##NAME)) : fnan); \
    if (mc) f.write(enough_jets ? float(F(**b_##NAME##_truth)) : fnan); \
  }

      VAR(pT_j1_30,_1e3)
      VAR(pT_yyj_30,_1e3)
      VAR(m_yyj_30,_1e3)
      VAR(maxTau_yyj_30,_1e3)
      VAR(sumTau_yyj_30,_1e3)

      // 2 jets =====================================================
      enough_jets = N_j_30 >= 2;

      VAR(Dphi_j_j_30,)
      VAR(Dphi_j_j_30_signed,)
      VAR(Dphi_yy_jj_30,)

      VAR(pT_yyjj_30,_1e3)
      VAR(m_jj_30,_1e3)

    } // event loop
  } // file loop
} catch (const std::exception& e) {
  cerr << "\033[0;31m" << e.what() << "\033[0m\n";
  return 1;
}} // main

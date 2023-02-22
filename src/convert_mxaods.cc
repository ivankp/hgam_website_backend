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

#define ERROR(...) ivan::error(__FILE__ ":" STR(__LINE__) ": ",__VA_ARGS__)

constexpr float fnan = std::numeric_limits<float>::quiet_NaN();

std::filesystem::path dir;
bool mc = false;
uint64_t nevents = 0;
float lumi = 0, weight = 0;
double mc_factor = 0;

template <typename T> constexpr char type_byte = '\0';
// https://docs.python.org/3/library/struct.html#format-characters
template <> constexpr char type_byte<float  > = 'f';
template <> constexpr char type_byte<char   > = 'c';
template <> constexpr char type_byte<int32_t> = 'i';
template <> constexpr char type_byte<uint8_t> = 'B';

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
class output_file {
  std::filesystem::path path;
  const int fd;
  char* const buf;
  char* m;
  off_t skipped;

  static constexpr size_t buflen = 1 << 20;

public:
  output_file(std::string_view name)
  : path(dir/cat(name, mc ? "_mc" : "_data", ".dat")),
    fd([this]{
      int fd = ::open(path.c_str(), O_CREAT|O_WRONLY|O_TRUNC, 0644);
      if (fd < 0) ERROR('\"',path.native(),"\": open(): ",strerror(errno));
      return fd;
    }()),
    buf( static_cast<char*>( ::malloc(buflen) ) ),
    m(buf)
  {
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
  ~output_file() {
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

int main(int argc, char** argv) { try {
  if (argc < 3) {
    cout << "usage: " << argv[0] << " output_dir input.root ...\n";
    return 1;
  }

  cout << "Input files:\n";
  for (int i=2; i<argc; ++i) {
    cout << argv[i] << '\n';

    // Luminosity
    static const std::regex re(R"((\d+(?:\.\d*)?)ipb\b)");
    std::cmatch match;
    if (std::regex_search(argv[i],match,re)) {
      if (i>2 &&  mc) ERROR("data file after mc files");
      mc = false;
      lumi += ivan::stox<double>(match[1].str());
    } else {
      if (i>2 && !mc) ERROR("mc file after data files");
      mc = true;
    }
  }
  cout << endl;

  if (mc) {
    cout << "Monte Carlo\n\n";
  } else {
    cout << "Data\nLumi = " << lumi << " ipb\n\n";
  }

  // Open output files ==============================================
  std::filesystem::create_directories(dir = argv[1]);
  cout << "Output files:\n";

  output_file<uint8_t> f_N_j_30("N_j_30");

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
  VAR(Float_t,m_yy) \
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

#define VAR(TYPE,NAME) \
  output_file<TYPE> f_##NAME(#NAME);

    VARS

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

#define VAR(TYPE,NAME) \
  branch_reader<TYPE> b_##NAME(reader, VAR_PREF #NAME); \
  std::optional<branch_reader<TYPE>> b_##NAME##_truth; \
  if (mc) b_##NAME##_truth.emplace(reader, VAR_PREF_TRUTH #NAME);

    VAR(Int_t,N_j_30)

    VARS

#undef VAR
#undef VARS

    // LOOP over events =============================================
    for (
      ivan::timed_counter<Long64_t> ent(reader.GetEntries(true));
      reader.Next();
      ++ent
    ) {
      // selection cut
      if (!*b_isPassed) continue;

      // diphoton mass cut
      const float m_yy = *b_m_yy*1e-3;
      if (m_yy < 105 || 160 < m_yy) continue;

      ++nevents;

      f_m_yy.write(m_yy);

      if (mc) {
        weight = double(**b_weight) * double(**b_cs_br_fe) * mc_factor;
        f_m_yy.write(weight);
      }

      // ============================================================
      const uint8_t N_j_30 = ceiling_cast<uint8_t>(*b_N_j_30);
      f_N_j_30.write(N_j_30);
      if (mc) f_N_j_30.write(ceiling_cast<uint8_t>(**b_N_j_30_truth));

#define VAR(TYPE,NAME) \
  const TYPE NAME = *b_##NAME; \
  f_##NAME.write(NAME); \
  if (mc) f_##NAME.write<TYPE>(**b_##NAME##_truth);

#define VAR_1e3(NAME) \
  const float NAME = *b_##NAME * 1e-3; \
  f_##NAME.write(NAME); \
  if (mc) f_##NAME.write<float>(**b_##NAME##_truth * 1e-3);

      VAR(Float_t,abs_Zepp)
      VAR( Char_t,catXS_lepton)
      VAR( Char_t,catXS_MET)
      VAR( Char_t,catXS_nbjet)
      VAR( Char_t,catXS_ttH)
      VAR( Char_t,catXS_VBF)
      VAR(Float_t,HT_30)
      VAR(Float_t,yAbs_yy)

      VAR_1e3(pT_yy)
      VAR_1e3(pT_yy_JV_30)
      VAR_1e3(pT_yy_JV_40)
      VAR_1e3(pT_yy_JV_50)
      VAR_1e3(pT_yy_JV_60)
      VAR_1e3(rel_DpT_y_y)
      VAR_1e3(rel_pT_y1)
      VAR_1e3(rel_pT_y2)
      VAR_1e3(rel_sumpT_y_y)

#undef VAR
#undef VAR_1e3

      // 1 jet ======================================================
      bool enough_jets = N_j_30 >= 1;

#define VAR(NAME) \
  const float NAME = enough_jets ? float(*b_##NAME) : fnan; \
  f_##NAME.write(NAME); \
  if (mc) f_##NAME.write(enough_jets ? float(**b_##NAME##_truth) : fnan);

#define VAR_1e3(NAME) \
  const float NAME = enough_jets ? float(*b_##NAME * 1e-3) : fnan; \
  f_##NAME.write(NAME); \
  if (mc) f_##NAME.write(enough_jets ? float(**b_##NAME##_truth * 1e-3) : fnan);

      VAR(pT_j1_30)
      VAR(pT_yyj_30)
      VAR(m_yyj_30)
      VAR(maxTau_yyj_30)
      VAR(sumTau_yyj_30)

      // 2 jets =====================================================
      enough_jets = N_j_30 >= 2;

      VAR(Dphi_j_j_30)
      VAR(Dphi_j_j_30_signed)
      VAR(Dphi_yy_jj_30)

      VAR_1e3(pT_yyjj_30)
      VAR_1e3(m_jj_30)

    } // event loop
  } // file loop
} catch (const std::exception& e) {
  cerr << "\033[0;31m" << e.what() << "\033[0m\n";
  return 1;
}} // main

#include <iostream>
#include <filesystem>
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

#include "program_options.hh"
#include "timed_counter.hh"
#include "branch_reader.hh"
#include "debug.hh"

using std::cout, std::cerr, std::endl;
using ivan::cat, ivan::branch_reader;

#define VAR_PREF "HGamEventInfoAuxDyn."
#define VAR_PREF_TRUTH "HGamTruthEventInfoAuxDyn."

#define ERROR(...) ivan::error(__FILE__ ":" STR(__LINE__) ": ",__VA_ARGS__)

constexpr float fnan = std::numeric_limits<float>::quiet_NaN();

std::filesystem::path dir;
bool is_mc = false;
uint64_t nevents = 0;
float lumi = 0;

template <typename T> constexpr char type_byte = '\0';
// https://docs.python.org/3/library/struct.html#format-characters
template <> constexpr char type_byte<float> = 'f';
template <> constexpr char type_byte<char > = 'c';
template <> constexpr char type_byte<int  > = 'i';

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
  : path(dir/cat(name, is_mc ? "_mc" : "_data", ".dat")),
    fd([this]{
      int fd = ::open(path.c_str(), O_CREAT|O_WRONLY|O_TRUNC, 0644);
      if (fd < 0) ERROR('\"',path.native(),"\": open(): ",strerror(errno));
      return fd;
    }()),
    buf( static_cast<char*>( ::malloc(buflen) ) ),
    m(buf)
  {
    static constexpr size_t skip = sizeof(lumi) + sizeof(nevents);
    TEST(name)
    TEST(name.size())
    TEST(skip+2)
    const size_t nulls = 8-((name.size()+2+skip)%8);
    TEST(nulls)
    const size_t bufneed = name.size()+nulls+2+skip;
    if (buflen < bufneed)
      ERROR('\"',path.native(),"\": buffer length ",buflen," , need ",bufneed);

    cout << path.native() << '\n';

    // write header
    memcpy(m,name.data(),name.size()); m += name.size();
    char tmp[8] { };
    memcpy(m,tmp,nulls); m += nulls;
    *m++ = is_mc ? 'm' : 'd';
    static_assert(type_byte<T> != '\0');
    *m++ = type_byte<T>;

    skipped = m - buf;
    m += skip;
  }
  ~output_file() {
    flush();
    m = bufcpy(buf,lumi);
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

int main(int argc, char** argv) {
  if (argc < 3) {
    cout << "usage: " << argv[0] << " output_dir input.root ...\n";
    return 1;
  }

  cout << "Input files:\n";
  for (int i=2; i<argc; ++i) {
    cout << argv[i] << '\n';
  }
  cout << endl;

  // Open output files ==============================================
  std::filesystem::create_directories(dir = argv[1]);
  cout << "Output files:\n";

  output_file<Float_t> f_abs_Zepp("abs_Zepp");
  output_file< Char_t> f_catXS_lepton("catXS_lepton");
  output_file< Char_t> f_catXS_MET("catXS_MET");
  output_file< Char_t> f_catXS_nbjet("catXS_nbjet");
  output_file< Char_t> f_catXS_ttH("catXS_ttH");
  output_file< Char_t> f_catXS_VBF("catXS_VBF");
  // output_file<  Int_t> f_cutFlow("cutFlow");
  output_file<Float_t> f_Dphi_j_j_30("Dphi_j_j_30");
  output_file<Float_t> f_Dphi_j_j_30_signed("Dphi_j_j_30_signed");
  output_file<Float_t> f_Dphi_yy_jj_30("Dphi_yy_jj_30");
  output_file<Float_t> f_HT_30("HT_30");
  // output_file< Char_t> f_isPassed("isPassed");
  output_file<Float_t> f_m_jj_30("m_jj_30");
  output_file<Float_t> f_m_yy("m_yy");
  output_file<Float_t> f_m_yyj_30("m_yyj_30");
  output_file<Float_t> f_maxTau_yyj_30("maxTau_yyj_30");
  // output_file<Float_t> f_mu("mu");
  output_file<  Int_t> f_N_j_30("N_j_30");
  output_file<Float_t> f_pT_j1_30("pT_j1_30");
  output_file<Float_t> f_pT_yy("pT_yy");
  output_file<Float_t> f_pT_yy_JV_30("pT_yy_JV_30");
  output_file<Float_t> f_pT_yy_JV_40("pT_yy_JV_40");
  output_file<Float_t> f_pT_yy_JV_50("pT_yy_JV_50");
  output_file<Float_t> f_pT_yy_JV_60("pT_yy_JV_60");
  output_file<Float_t> f_pT_yyj_30("pT_yyj_30");
  output_file<Float_t> f_pT_yyjj_30("pT_yyjj_30");
  output_file<Float_t> f_rel_DpT_y_y("rel_DpT_y_y");
  output_file<Float_t> f_rel_pT_y1("rel_pT_y1");
  output_file<Float_t> f_rel_pT_y2("rel_pT_y2");
  output_file<Float_t> f_rel_sumpT_y_y("rel_sumpT_y_y");
  output_file<Float_t> f_sumTau_yyj_30("sumTau_yyj_30");
  output_file<Float_t> f_yAbs_yy("yAbs_yy");

  cout << endl;

  for (int i=2; i<argc; ++i) { // loop over input files
    cout << argv[i] << endl;
    TFile file(argv[i]);
    TTreeReader reader("CollectionTree",&file);

    branch_reader< Char_t> b_isPassed(reader, VAR_PREF "isPassed");

    branch_reader<Float_t> b_abs_Zepp(reader, VAR_PREF "abs_Zepp");
    branch_reader< Char_t> b_catXS_lepton(reader, VAR_PREF "catXS_lepton");
    branch_reader< Char_t> b_catXS_MET(reader, VAR_PREF "catXS_MET");
    branch_reader< Char_t> b_catXS_nbjet(reader, VAR_PREF "catXS_nbjet");
    branch_reader< Char_t> b_catXS_ttH(reader, VAR_PREF "catXS_ttH");
    branch_reader< Char_t> b_catXS_VBF(reader, VAR_PREF "catXS_VBF");
    // branch_reader<  Int_t> b_cutFlow(reader, VAR_PREF "cutFlow");
    branch_reader<Float_t> b_Dphi_j_j_30(reader, VAR_PREF "Dphi_j_j_30");
    branch_reader<Float_t> b_Dphi_j_j_30_signed(reader, VAR_PREF "Dphi_j_j_30_signed");
    branch_reader<Float_t> b_Dphi_yy_jj_30(reader, VAR_PREF "Dphi_yy_jj_30");
    branch_reader<Float_t> b_HT_30(reader, VAR_PREF "HT_30");
    branch_reader<Float_t> b_m_jj_30(reader, VAR_PREF "m_jj_30");
    branch_reader<Float_t> b_m_yy(reader, VAR_PREF "m_yy");
    branch_reader<Float_t> b_m_yyj_30(reader, VAR_PREF "m_yyj_30");
    branch_reader<Float_t> b_maxTau_yyj_30(reader, VAR_PREF "maxTau_yyj_30");
    // branch_reader<Float_t> b_mu(reader, VAR_PREF "mu");
    branch_reader<  Int_t> b_N_j_30(reader, VAR_PREF "N_j_30");
    branch_reader<Float_t> b_pT_j1_30(reader, VAR_PREF "pT_j1_30");
    branch_reader<Float_t> b_pT_yy(reader, VAR_PREF "pT_yy");
    branch_reader<Float_t> b_pT_yy_JV_30(reader, VAR_PREF "pT_yy_JV_30");
    branch_reader<Float_t> b_pT_yy_JV_40(reader, VAR_PREF "pT_yy_JV_40");
    branch_reader<Float_t> b_pT_yy_JV_50(reader, VAR_PREF "pT_yy_JV_50");
    branch_reader<Float_t> b_pT_yy_JV_60(reader, VAR_PREF "pT_yy_JV_60");
    branch_reader<Float_t> b_pT_yyj_30(reader, VAR_PREF "pT_yyj_30");
    branch_reader<Float_t> b_pT_yyjj_30(reader, VAR_PREF "pT_yyjj_30");
    branch_reader<Float_t> b_rel_DpT_y_y(reader, VAR_PREF "rel_DpT_y_y");
    branch_reader<Float_t> b_rel_pT_y1(reader, VAR_PREF "rel_pT_y1");
    branch_reader<Float_t> b_rel_pT_y2(reader, VAR_PREF "rel_pT_y2");
    branch_reader<Float_t> b_rel_sumpT_y_y(reader, VAR_PREF "rel_sumpT_y_y");
    branch_reader<Float_t> b_sumTau_yyj_30(reader, VAR_PREF "sumTau_yyj_30");
    branch_reader<Float_t> b_yAbs_yy(reader, VAR_PREF "yAbs_yy");

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

      // ============================================================

      const Float_t abs_Zepp = *b_abs_Zepp;
      f_abs_Zepp.write(abs_Zepp);

      const Char_t catXS_lepton = *b_catXS_lepton;
      f_catXS_lepton.write(catXS_lepton);

      const Char_t catXS_MET = *b_catXS_MET;
      f_catXS_MET.write(catXS_MET);

      const Char_t catXS_nbjet = *b_catXS_nbjet;
      f_catXS_nbjet.write(catXS_nbjet);

      const Char_t catXS_ttH = *b_catXS_ttH;
      f_catXS_ttH.write(catXS_ttH);

      const Char_t catXS_VBF = *b_catXS_VBF;
      f_catXS_VBF.write(catXS_VBF);

      const Float_t HT_30 = *b_HT_30*1e-3;
      f_HT_30.write(HT_30);

      const Float_t pT_yy = *b_pT_yy*1e-3;
      f_pT_yy.write(pT_yy);

      const Float_t pT_yy_JV_30 = *b_pT_yy_JV_30*1e-3;
      f_pT_yy_JV_30.write(pT_yy_JV_30);

      const Float_t pT_yy_JV_40 = *b_pT_yy_JV_40*1e-3;
      f_pT_yy_JV_40.write(pT_yy_JV_40);

      const Float_t pT_yy_JV_50 = *b_pT_yy_JV_50*1e-3;
      f_pT_yy_JV_50.write(pT_yy_JV_50);

      const Float_t pT_yy_JV_60 = *b_pT_yy_JV_60*1e-3;
      f_pT_yy_JV_60.write(pT_yy_JV_60);

      const Float_t rel_DpT_y_y = *b_rel_DpT_y_y*1e-3;
      f_rel_DpT_y_y.write(rel_DpT_y_y);

      const Float_t rel_pT_y1 = *b_rel_pT_y1*1e-3;
      f_rel_pT_y1.write(rel_pT_y1);

      const Float_t rel_pT_y2 = *b_rel_pT_y2*1e-3;
      f_rel_pT_y2.write(rel_pT_y2);

      const Float_t rel_sumpT_y_y = *b_rel_sumpT_y_y*1e-3;
      f_rel_sumpT_y_y.write(rel_sumpT_y_y);

      const Float_t yAbs_yy = *b_yAbs_yy;
      f_yAbs_yy.write(yAbs_yy);

      const Int_t N_j_30 = *b_N_j_30;
      f_N_j_30.write(N_j_30);

      // 1 jet ======================================================
      bool enough_jets = N_j_30 >= 1;
#define J(...) \
  enough_jets ? float(__VA_ARGS__) : fnan;

      const Float_t pT_j1_30 = J(*b_pT_j1_30*1e-3);
      f_pT_j1_30.write(pT_j1_30);

      const Float_t pT_yyj_30 = J(*b_pT_yyj_30*1e-3);
      f_pT_yyj_30.write(pT_yyj_30);

      const Float_t m_yyj_30 = J(*b_m_yyj_30*1e-3);
      f_m_yyj_30.write(m_yyj_30);

      const Float_t maxTau_yyj_30 = J(*b_maxTau_yyj_30);
      f_maxTau_yyj_30.write(maxTau_yyj_30);

      const Float_t sumTau_yyj_30 = J(*b_sumTau_yyj_30);
      f_sumTau_yyj_30.write(sumTau_yyj_30);

      // 2 jets =====================================================
      enough_jets = N_j_30 >= 2;

      const Float_t pT_yyjj_30 = J(*b_pT_yyjj_30*1e-3);
      f_pT_yyjj_30.write(pT_yyjj_30);

      const Float_t m_jj_30 = J(*b_m_jj_30*1e-3);
      f_m_jj_30.write(m_jj_30);

      const Float_t Dphi_j_j_30 = J(*b_Dphi_j_j_30);
      f_Dphi_j_j_30.write(Dphi_j_j_30);

      const Float_t Dphi_j_j_30_signed = J(*b_Dphi_j_j_30_signed);
      f_Dphi_j_j_30_signed.write(Dphi_j_j_30_signed);

      const Float_t Dphi_yy_jj_30 = J(*b_Dphi_yy_jj_30);
      f_Dphi_yy_jj_30.write(Dphi_yy_jj_30);
    }
  }
}

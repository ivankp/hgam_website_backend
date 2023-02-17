#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <limits>
#include <regex>
#include <optional>

#include <boost/preprocessor/facilities/expand.hpp>

#include <TFile.h>
#include <TH1.h>
#include <TKey.h>
#include <TLorentzVector.h>

#include "program_options.hh"
#include "timed_counter.hh"
#include "branch_reader.hh"
#include "debug.hh"

using std::cout, std::cerr, std::endl, std::get;
using ivan::cat, ivan::error, ivan::branch_reader;

constexpr float fnan = std::numeric_limits<float>::quiet_NaN();

template <typename T, typename U>
inline bool in(const std::array<T,2>& reg, U x) {
  return (get<0>(reg) < x && x < get<1>(reg));
}

double get_lumi(TFile& f) {
  const auto* lumi = f.Get("Lumi");
  if (lumi) return atof(lumi->GetTitle());
  static const std::regex re("^(?:.*/)?data.*_(\\d+)ipb\\..*\\.root$");
  std::cmatch match;
  return std::regex_match(f.GetName(),match,re)
    ? ivan::stox<double>(match[1].str()) : 0.;
}

bool is_mc = false, is_fiducial = false;
const char* opref;
size_t nevents = 0;
float lumi = 0, weight = 0;
double mc_factor = 0;
const char* catXS_name = nullptr;

class var_file {
  std::ofstream f;
  decltype(f.tellp()) head;
public:
  size_t nevents = 0;
  var_file(const char* var_name)
  : f(cat(opref,var_name,"_",(is_mc?"mc":"data"),".dat"))
  {
    if (catXS_name) f << catXS_name << '_';
    f << var_name;
    head = f.tellp();
    fill_n(std::ostream_iterator<char>(f),
      8-(((unsigned)head+(is_mc?1:5))%8), ' ');
    f << (is_mc ? 'm' : 'd');
    head = f.tellp();
    if (!is_mc) write(lumi);
    write(nevents);
  }
  ~var_file() {
    f.seekp(head);
    if (!is_mc) write((float)(lumi*1e-3));
    write(nevents);
  }
  template <typename T>
  inline var_file& write(const T& x) {
    f.write(reinterpret_cast<const char*>(&x),sizeof(T));
    return *this;
  }
};

template <typename C, typename T>
unsigned find(const C& c, const T& x) {
  const auto _begin = begin(c);
  const auto _end = end(c);
  const auto it = std::find(_begin,_end,x);
  if (it==_end) throw error(x," not found");
  return std::distance(_begin,it);
}

template <typename T>
std::array<T,2> split_range(const char* s) {
  const char* d = strchr(s,':');
  if (!d) error("missing : in range option");
  return {
    ivan::stox<T>({s,d}),
    ivan::stox<T>({d+1,strlen(d+1)})
  };
}

int main(int argc, char* argv[]) {
  std::vector<const char*> ifnames;
  std::array<float,2> dat_reg{105,160}, sig_reg{121,129};

  try {
    using ivan::po::opt;
    ivan::program_options(argc-1,argv+1,
      opt("i",[&](const char* x){ ifnames.push_back(x); }),
      opt("o",opref),
      opt("mc",is_mc),
      opt("d",[&](const char* x){ dat_reg = split_range<float>(x); }),
      opt("s",[&](const char* x){ sig_reg = split_range<float>(x); }),
      opt("cat",catXS_name)
    );
  } catch (const std::exception& e) {
    cerr << e.what() << endl;
    return 1;
  }

  // Output files ===================================================
#define VAR_FILE(X) var_file f_##X(#X);

  VAR_FILE(pT_yy)
  VAR_FILE(HT)
  VAR_FILE(HTall)
  VAR_FILE(yAbs_yy)
  VAR_FILE(cosTS_yy)
  VAR_FILE(pTt_yy)
  VAR_FILE(Dy_y_y)

  VAR_FILE(N_j)

  VAR_FILE(pT_j1)
  VAR_FILE(yAbs_j1)
  VAR_FILE(sumTau_yyj)
  VAR_FILE(maxTau_yyj)

  VAR_FILE(pT_j2)
  VAR_FILE(yAbs_j2)
  VAR_FILE(Dphi_j_j)
  VAR_FILE(Dphi_j_j_signed)
  VAR_FILE(Dy_j_j)
  VAR_FILE(m_jj)
  VAR_FILE(pT_yyjj)
  VAR_FILE(Dphi_yy_jj)

  VAR_FILE(pT_j3)

  VAR_FILE(Hj_mass)
  VAR_FILE(Hj_pT)

  TLorentzVector y1, y2, jet1;

  // Read files =====================================================
  for (const char* ifname : ifnames) {
    cout << (is_mc ? "mc" : "data") << ": \033[32m" << ifname << "\033[0m"
         << endl;

    TFile file(ifname);

    if (is_mc) { // MC
      for (auto* key : *file.GetListOfKeys()) {
        const char* name = key->GetName();
        if (!ivan::starts_with(name,"CutFlow_") ||
            !ivan::ends_with(name,"_noDalitz_weighted")) continue;
        TH1 *h = static_cast<TH1*>(static_cast<TKey*>(key)->ReadObj());
        cout << name << endl;
        const double n_all = h->GetBinContent(3);
        cout << h->GetXaxis()->GetBinLabel(3) << " = " << n_all << endl;
        mc_factor = 1e3/n_all;
        break;
      }
    } else { // Data
      lumi += get_lumi(file);
    }

    TTreeReader reader("CollectionTree",&file);

#define VAR_PREF "HGamEventInfoAuxDyn."
#define VAR_PREF_TRUTH "HGamTruthEventInfoAuxDyn."

#define GET_MACRO(_1,_2,_3,NAME,...) NAME

#define VAR(...) GET_MACRO(__VA_ARGS__, VAR_3, VAR_2, VAR_1)(__VA_ARGS__)
#define VAR_3(NAME,TYPE,ACTUAL) \
    branch_reader<TYPE> _##NAME(reader,VAR_PREF ACTUAL); \
    MCVAR(NAME##_truth,TYPE,VAR_PREF_TRUTH ACTUAL)
#define VAR_2(NAME,TYPE) VAR_3(NAME,TYPE,#NAME)
#define VAR_1(NAME) VAR_2(NAME,Float_t)

#define MCVAR(...) GET_MACRO(__VA_ARGS__, MCVAR_3, MCVAR_2, MCVAR_1)(__VA_ARGS__)
#define MCVAR_3(NAME,TYPE,FULL_ACTUAL) \
    std::optional<branch_reader<TYPE>> _##NAME; \
    if (is_mc) _##NAME.emplace(reader,FULL_ACTUAL);
#define MCVAR_2(NAME,TYPE) MCVAR_3(NAME,TYPE,VAR_PREF #NAME)
#define MCVAR_1(NAME) MCVAR_2(NAME,Float_t)

#define VARJ(...) GET_MACRO(__VA_ARGS__, _3, VARJ_2, VARJ_1)(__VA_ARGS__)
#define VARJ_2(NAME,TYPE) VAR(NAME,TYPE,#NAME"_30")
#define VARJ_1(NAME) VAR(NAME,Float_t,#NAME"_30")

    branch_reader<Char_t> _isPassed(reader,VAR_PREF "isPassed");

    MCVAR(isFiducial,Char_t,VAR_PREF_TRUTH "isFiducial")
    MCVAR(cs_br_fe,Float_t,VAR_PREF "crossSectionBRfilterEff")
    MCVAR(weight)

    VARJ(N_j,Int_t)

    VAR(m_yy)
    VAR(pT_yy) VAR(yAbs_yy) VAR(cosTS_yy) VAR(pTt_yy) VAR(Dy_y_y)

    VARJ(HT)         VARJ(HTall)
    VARJ(pT_j1)      VARJ(pT_j2)      VARJ(pT_j3)
    VARJ(yAbs_j1)    VARJ(yAbs_j2)
    VARJ(Dphi_j_j)   VAR(Dphi_j_j_signed,Float_t,"Dphi_j_j_30_signed")
    VARJ(Dy_j_j)     VARJ(m_jj)
    VARJ(sumTau_yyj) VARJ(maxTau_yyj)
    VARJ(pT_yyjj)    VARJ(Dphi_yy_jj)

    std::optional<branch_reader<Char_t>> _catXS, _catXS_truth;
    if (catXS_name) {
      _catXS.emplace(reader,cat(VAR_PREF,"catXS_",catXS_name).c_str());
      if (is_mc) _catXS_truth.emplace(reader,
          cat(VAR_PREF_TRUTH,"catXS_",catXS_name).c_str());
    }

    std::array<branch_reader<std::vector<float>>,4> _photons {{
      {reader,"HGamPhotonsAuxDyn.pt"},
      {reader,"HGamPhotonsAuxDyn.eta"},
      {reader,"HGamPhotonsAuxDyn.phi"},
      {reader,"HGamPhotonsAuxDyn.m"}
    }};
    std::array<branch_reader<std::vector<float>>,4> _jets {{
      {reader,"HGamAntiKt4EMTopoJetsAuxDyn.pt"},
      {reader,"HGamAntiKt4EMTopoJetsAuxDyn.eta"},
      {reader,"HGamAntiKt4EMTopoJetsAuxDyn.phi"},
      {reader,"HGamAntiKt4EMTopoJetsAuxDyn.m"}
    }};

    float m_yy=0;
    Int_t nj=0, nj_truth=99;

    // LOOP over events =============================================
    for (
      ivan::timed_counter<Long64_t> ent(reader.GetEntries(true));
      reader.Next();
      ++ent
    ) {
      // selection cut
      if (!*_isPassed) continue;
      if (_catXS && !**_catXS) continue;

      // diphoton mass cut
      m_yy = *_m_yy*1e-3;
      if (!in(dat_reg,m_yy)) continue;

      if (is_mc) { // signal from MC
        weight = (**_weight) * (**_cs_br_fe) * mc_factor;
        if (!in(sig_reg,m_yy)) continue;
        const auto m_yy_truth = **_m_yy_truth*1e-3;
        is_fiducial = **_isFiducial && in(sig_reg,m_yy_truth);
        if (_catXS_truth && !**_catXS_truth) is_fiducial = false;
      } else { // background from data
        if (in(sig_reg,m_yy)) continue;
      }
      ++nevents;

      // Write out variables ========================================
#define WRITE(X,EQ) \
      f_##X.write(m_yy).write( (float)(BOOST_PP_EXPAND( EQ(*_##X) )) ); \
      if (is_mc) f_##X.write( is_fiducial \
          ? (float)(BOOST_PP_EXPAND( EQ(**_##X##_truth) )) \
          : fnan \
        ).write(weight); \
      ++f_##X.nevents;

#define F_1e3(X) X*1e-3
#define F_ABS(X) std::abs(X)

      WRITE(pT_yy, F_1e3)
      WRITE(HT, F_1e3)
      WRITE(HTall, F_1e3)
      WRITE(yAbs_yy, )
      WRITE(cosTS_yy, F_ABS)
      WRITE(pTt_yy, F_1e3)
      WRITE(Dy_y_y, F_ABS)

#define F_nj(X) nj
      nj = *_N_j;
      if (is_mc) nj_truth = **_N_j_truth;
      WRITE(N_j,F_nj)

      if (nj < 1) continue; // 1 jet --------------------------------
      if (nj_truth < 1) is_fiducial = false;

      WRITE(pT_j1, F_1e3)
      WRITE(yAbs_j1, )
      WRITE(sumTau_yyj, F_1e3)
      WRITE(maxTau_yyj, F_1e3)

      // ============================================================
      y1.SetPtEtaPhiM(
        _photons[0][0],
        _photons[1][0],
        _photons[2][0],
        _photons[3][0]
      );
      y2.SetPtEtaPhiM(
        _photons[0][1],
        _photons[1][1],
        _photons[2][1],
        _photons[3][1]
      );
      jet1.SetPtEtaPhiM(
        _jets[0][0],
        _jets[1][0],
        _jets[2][0],
        _jets[3][0]
      );

      const TLorentzVector Hj = y1 + y2 + jet1;

      f_Hj_mass.write(m_yy).write( (float)(Hj.M()*1e-3) );
      if (is_mc) f_Hj_mass.write(fnan).write(weight);
      ++f_Hj_mass.nevents;

      f_Hj_pT.write(m_yy).write( (float)(Hj.Pt()*1e-3) );
      if (is_mc) f_Hj_pT.write(fnan).write(weight);
      ++f_Hj_pT.nevents;
      // ============================================================

      if (nj < 2) continue; // 2 jet --------------------------------
      if (nj_truth < 2) is_fiducial = false;

      WRITE(pT_j2, F_1e3)
      WRITE(yAbs_j2, )
      WRITE(Dphi_yy_jj, F_ABS)
      WRITE(Dphi_j_j_signed, )
      WRITE(Dphi_j_j, F_ABS)
      WRITE(Dy_j_j, F_ABS)
      WRITE(m_jj, F_1e3)
      WRITE(pT_yyjj, F_1e3)

      if (nj < 3) continue; // 3 jet --------------------------------
      if (nj_truth < 3) is_fiducial = false;

      WRITE(pT_j3, F_1e3)

    }
  }
  TEST(nevents)
  if (!is_mc) TEST(lumi)
}

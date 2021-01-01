#ifndef _ANA_EFF_HODO__H_
#define _ANA_EFF_HODO__H_
#include <fstream>
#include "TH1D.h"
#include <fun4all/SubsysReco.h>
#include "TGraphAsymmErrors.h"
#include "TEfficiency.h"
class TFile;
class TH1;

class AnaEffHodo: public SubsysReco {
  unsigned int n_evt_all;
  unsigned int n_evt_trig;
  unsigned int n_evt_nhit;
  unsigned int fail;
  unsigned int track_N1;
  unsigned int track_N2;
  unsigned int total[8]={0};
  unsigned int success[8]={0};

  std::ofstream ofs;
  std::ofstream out_File;
  TH1I* hist_all[10];
  TH1I* hist_acc[10];
  TH1D* h2xb_all[16];
  TH1D* h2xb_acc[16];
  TH1D* h2xt_all[16];
  TH1D* h2xt_acc[16];
  TH1D* hist_eff[10];
  TH1D* HodopaddleNamesh2xt[16];
  TH1D* HodopaddleNamesh2xb[16];
  TH1D* h2x_acc[8][16];

  TH1I* h_acc[12][16];
  TH1I* h_all[12][16];
  TH1I* hodo_names[12][16];

  TH1D* h2x_all[8][16];
  TH1I* paddle_diff[10];
  TH1I* hist_allY[8];
  TH1I* hist_accY[10];
  TH1I* paddle_diffY[10];
  TH1D* hist_effY[8];
  TGraphAsymmErrors* graph_eff[8];
  TEfficiency* hodo_eff[10];
  TEfficiency* hodo_effY[8];

   TEfficiency* eff_h2xt[16];
   TEfficiency* eff_h2x[10][16];
   TEfficiency* eff_hx[10][16];
   //vector<TEfficiency*> eff_h2x[16];

  TFile* f_out;
  TH1* h1_eff_all;
  TH1* h1_eff_ok;
  TH1* h1_nhit;
  TH1* h1_ele;
  TH1* h1_time;
  TH1* sqrt_chi2;
  TH1D* n_trial;
  TH1D* n_pass;
  TH1* exp_xH3XT;
  TH1* exp_xH3XB;
  TH1* exp_yDP2TR;
  TH1* TPos_HitPos_diff;

  TH1I* paddle_all[10][16];
  TH1I* paddle_acc[10][16];
  TH1D* paddle_eff[10][16];
 public:
  AnaEffHodo();
  virtual ~AnaEffHodo() {;}
  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);
};

#endif /* _ANA_EFF_HODO__H_ */

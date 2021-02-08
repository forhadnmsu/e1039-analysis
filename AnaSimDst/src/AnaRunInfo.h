#ifndef _ANA_RUN_INFO__H_
#define _ANA_RUN_INFO__H_
#include <fstream>
#include <fun4all/SubsysReco.h>
class PHGenIntegral;

/// An example class to analyze the simulated uDST file.
class AnaRunInfo: public SubsysReco {
  std::string m_fn_lumi_list;
  std::string m_fn_lumi_tot;
  double lumi_inte_tot;
  long   n_evt_gen_tot;
  long   n_evt_pro_tot;
  long   n_run;

  PHGenIntegral* mi_gen_inte;

  std::ofstream ofs_lumi;

 public:
  AnaRunInfo(const char* fn_lumi_list, const char* fn_lumi_tot);
  virtual ~AnaRunInfo() {;}
  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);
};

#endif // _ANA_RUN_INFO__H_

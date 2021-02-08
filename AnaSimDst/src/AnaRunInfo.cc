#include <iomanip>
#include <phhepmc/PHGenIntegral.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/getClass.h>
#include "AnaRunInfo.h"
using namespace std;

AnaRunInfo::AnaRunInfo(const char* fn_lumi_list, const char* fn_lumi_tot)
  : m_fn_lumi_list(fn_lumi_list)
  , m_fn_lumi_tot (fn_lumi_tot )
  , lumi_inte_tot(0)
  , n_evt_gen_tot(0)
  , n_evt_pro_tot(0)
  , n_run        (0)
{
  ;
}

int AnaRunInfo::Init(PHCompositeNode* topNode)
{
  ofs_lumi.open(m_fn_lumi_list.c_str());
  return Fun4AllReturnCodes::EVENT_OK;
}

int AnaRunInfo::InitRun(PHCompositeNode* topNode)
{
  mi_gen_inte = findNode::getClass<PHGenIntegral >(topNode, "PHGenIntegral");
  if (!mi_gen_inte) return Fun4AllReturnCodes::ABORTEVENT;

  double lumi_inte = mi_gen_inte->get_Integrated_Lumi();
  long   n_evt_gen = mi_gen_inte->get_N_Generator_Accepted_Event();
  long   n_evt_pro = mi_gen_inte->get_N_Processed_Event();
  cout << "Integrated luminosity = " << lumi_inte << " /pb\n"
       << "N of gen. acc. events = " << n_evt_gen << "\n"
       << "N of processed events = " << n_evt_pro << endl;
  ofs_lumi << lumi_inte << "\t" << n_evt_gen << "\t" << n_evt_pro << "\n";
  lumi_inte_tot += lumi_inte;
  n_evt_gen_tot += n_evt_gen;
  n_evt_pro_tot += n_evt_pro;
  n_run++;

  return Fun4AllReturnCodes::EVENT_OK;
}

int AnaRunInfo::process_event(PHCompositeNode* topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int AnaRunInfo::End(PHCompositeNode* topNode)
{
  ofs_lumi.close();
  ofstream ofs(m_fn_lumi_tot);
  cout << "Integrated luminosity = " << lumi_inte_tot << " /pb\n"
       << "N of gen. acc. events = " << n_evt_gen_tot << "\n"
       << "N of processed events = " << n_evt_pro_tot << "\n"
       << "N of input runs       = " << n_run         << endl;
  ofs.close();
  return Fun4AllReturnCodes::EVENT_OK;
}

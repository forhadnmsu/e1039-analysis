#ifndef _ANA_Module__H_
#define _ANA_Module__H_

#include <map>
#include <set>
#include <fun4all/SubsysReco.h>
#include <TString.h>
#include <TVector3.h>
#include <ktracker/SRecEvent.h>
#include <ktracker/FastTracklet.h>
#include <geom_svc/GeomSvc.h>
#include <interface_main/SQHit_v1.h>
#include "TH1D.h"
#include "TGraphAsymmErrors.h"
#include "TEfficiency.h"
#include <TCanvas.h>




class TFile;
class TTree;
class SQHitVector;
class SQTrackVector;
class SQDimuonVector;
class TH1;

class AnaModule: public SubsysReco 
{
public:
  AnaModule(const std::string& name = "AnaModule");
  virtual ~AnaModule();

  int Init(PHCompositeNode* topNode);
  int InitRun(PHCompositeNode* topNode);
  int process_event(PHCompositeNode* topNode);
  int End(PHCompositeNode* topNode);

  void set_output_filename(const TString& n) { saveName = n; }
  void registerDetector(TString name);

private:
  int GetNodes(PHCompositeNode* topNode);
  void MakeTree();

  SQHit* findHit(int detectorID, int elementID);
  std::set<int> detectorIDs;

  GeomSvc* p_geomSvc;

  // Input
  SQHitVector*    hitVector;
  TrackletVector* trackletVec;

  // Output
  TString saveName;
  TFile*  saveFile;
  TTree*  saveTree;

  int eventID;
  int detectorID;
  int elementID_exp;
  int elementID_closest;
  double x_exp;
  double y_exp;
  int nHits;
  double chisq;


//eff plot
 int eleHY[8];
 int eleHY_hit[8];
 TH1I* hist_all[8];
 TH1I* hist_acc[8];
 TH1D* hist_eff[8];
 TH1D* paddle_diff[8];
 TEfficiency* hodo_eff[8];
 TEfficiency* hodo_effY[8];
  TH1I* hist_allY[8];
  TH1I* hist_accY[8];
  TH1D* paddle_diffY[8];
  TH1D* hist_effY[8];
  TGraphAsymmErrors* graph_eff[8];
int hitInHX[3]; 
int x_ele[8] ;
int y_ele[8];
int y_eleHit[8];
};

#endif

#ifndef _ANA_SIM_DST__H_
#define _ANA_SIM_DST__H_
#include <map>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <fun4all/SubsysReco.h>
#include "TreeData.h"
#include<vector>
namespace HepMC {
  class GenParticle;
};
class TFile;
class TTree;
class SQEvent;
class SRecEvent;
class SQMCEvent;
class SQTrackVector;
class SQDimuonVector;
class PHHepMCGenEventMap;
class PHGenIntegral;
class SQHitVector;  //this class is for getting the list of spectrometer hits 
class PHG4TruthInfoContainer; //truth container
/// An example class to analyze the simulated uDST file.
class AnaSimDst: public SubsysReco {
  /// Input
  SQEvent       * mi_evt;
  SRecEvent     * mi_srec;
  SQMCEvent     * mi_evt_true;
  SQTrackVector * mi_vec_trk;
  SQDimuonVector* mi_vec_dim;
  SQHitVector*    mi_sqhitvec;
  PHHepMCGenEventMap* genevtmap;
  PHGenIntegral * mi_gen_inte;

  /// Output
  TFile* file;
  TTree* tree;
  EventData  mo_evt;
  TrackList  mo_trk_true;
  TrackList  mo_trk_reco;
  DimuonList mo_dim_true;
  DimuonList mo_dim_reco;
  
  int count=0;
  int mu_ID=0;
  int parent_id_[2];
 public:
  AnaSimDst() {;}
  virtual ~AnaSimDst() {;}
  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

 private:
  int GetNodes(PHCompositeNode *topNode);
  void MakeTree();

  typedef std::map<int, int> IdMap_t; // For now the key is not ID but index.
  void FindTrackRelation (IdMap_t& id_map);
  void FindDimuonRelation(IdMap_t& id_map);
  //void TraceParent(const HepMC::GenParticle* par, const int depth);
  void TraceParent(const HepMC::GenParticle* par, const int depth, int parent_id_[2]);
};

#endif /* _ANA_SIM_DST__H_ */

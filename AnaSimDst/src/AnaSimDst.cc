#include <iomanip>
#include <TFile.h>
#include <TTree.h>
#include <interface_main/SQEvent.h>
#include <interface_main/SQMCEvent.h>
#include <interface_main/SQTrackVector.h>
#include <interface_main/SQTrack.h>
#include <interface_main/SQDimuonVector.h>
#include <ktracker/SRecEvent.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <phool/getClass.h>
#include <UtilAna/UtilDimuon.h>
#include "AnaSimDst.h"
#include <geom_svc/GeomSvc.h>
#include <interface_main/SQHitVector_v1.h>
using namespace std;

int AnaSimDst::Init(PHCompositeNode* topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int AnaSimDst::InitRun(PHCompositeNode* topNode)
{
  int ret = GetNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;
  MakeTree();
  return Fun4AllReturnCodes::EVENT_OK;
}

int AnaSimDst::process_event(PHCompositeNode* topNode)
{
  std::map<int, int> parID_nhits_hodo;
  static unsigned int n_evt = 0;
  if    (++n_evt % 100000 == 0) cout << n_evt << endl;
  else if (n_evt %  10000 == 0) cout << " . " << flush;

  ///
  /// Event info
  ///
  mo_evt.proc_id = mi_evt_true->get_process_id();
  for (int ii = 0; ii < 4; ii++) {
    mo_evt.par_id [ii] = mi_evt_true->get_particle_id      (ii);
    mo_evt.par_mom[ii] = mi_evt_true->get_particle_momentum(ii);
  }
  mo_evt.trig_bits  = mi_evt->get_trigger();
  mo_evt.rec_stat   = mi_srec->getRecStatus();
  mo_evt.n_dim_true = mi_vec_dim->size();
  mo_evt.n_dim_reco = mi_srec->getNDimuons();
  mo_evt.trig_matrix1 = mi_evt->get_trigger(SQEvent::MATRIX1);
  mo_evt.trig_matrix2 = mi_evt->get_trigger(SQEvent::MATRIX2);


 ///
  /// Track info
  ///
  IdMap_t id_trk_t2r;
  FindTrackRelation(id_trk_t2r);
  mo_trk_true.clear();
  mo_trk_reco.clear();
  int tot_trk=0; 
  for (unsigned int ii = 0; ii < mi_vec_trk->size(); ii++) {
	  SQTrack* trk = mi_vec_trk->at(ii);
	  TrackData td;
	  td.charge  = trk->get_charge();
	  td.pos_vtx = trk->get_pos_vtx();
	  //get hodoscope hits
	  GeomSvc::UseDbSvc(true);
	  GeomSvc* geom = GeomSvc::instance();
	  int det_ID1 = geom->getDetectorID("H1B");
	  int det_ID2 = geom->getDetectorID("H4T");
	  int nhit_hodo=0;
	  if(mi_sqhitvec) {
		  for(int ihit=0; ihit<mi_sqhitvec->size(); ++ihit) {
			  SQHit *hit = mi_sqhitvec->at(ihit);
			  //if(_truth) {
			  int track_id = hit->get_track_id();
			  if(hit->get_detector_id() >= det_ID1 and hit->get_detector_id() <=det_ID2) {
				//  if (trk->get_track_id() != hit->get_track_id()) continue;
				  if (trk->get_track_id() == hit->get_track_id()){
				  cout <<"charge: "<<trk->get_charge()<<"     get_track_id: "<<hit->get_track_id() <<"     get_detector_id () : "<< hit->get_detector_id () << "get_tdc_time () :"<<  hit->get_tdc_time ()<< "     getDetectorName: " << geom->getDetectorName( hit->get_detector_id () ) <<endl;

				  nhit_hodo++;
				}
			  }
			 // if (nhit_hodo >=8)tot_trk = 1;
			 // else  continue;
			  td.nhodo = nhit_hodo;
			  }
		  //}
	  }

	  cout << " next track ================================================"<<endl;

	  mo_trk_true.push_back(td);

	  //
	  //    TrackData tdr;
	  //    if (id_trk_t2r[ii] >= 0) {
	  //      SRecTrack* trk_reco = &mi_srec->getTrack(id_trk_t2r[ii]);
	  //      tdr.charge  = trk_reco->getCharge();
	  //      tdr.pos_vtx = trk_reco->getVertex();
	  //      tdr.mom_vtx = trk_reco->getMomentumVertex();
	  //    }
	  //    mo_trk_reco.push_back(tdr);
  }
  ///
  /// Dimuon info
  ///
  IdMap_t id_dim_t2r;
  FindDimuonRelation(id_dim_t2r);
  mo_dim_true.clear();
  mo_dim_reco.clear();
  for (unsigned int ii = 0; ii < mi_vec_dim->size(); ii++) {
	  SQDimuon* dim = mi_vec_dim->at(ii);
	  DimuonData dd;
	  dd.pdg_id  = dim->get_pdg_id();
	  dd.pos     = dim->get_pos();
	  dd.mom     = dim->get_mom();
	  dd.mom_pos = dim->get_mom_pos();
	  dd.mom_neg = dim->get_mom_neg();
	  UtilDimuon::GetX1X2(dim, dd.x1, dd.x2);
	  mo_dim_true.push_back(dd);

	  DimuonData ddr;
	  if (id_dim_t2r[ii] >= 0) {
		  SRecDimuon dim_reco = mi_srec->getDimuon(id_dim_t2r[ii]);
		  ddr.pos     = dim_reco.vtx;
		  ddr.mom     = dim_reco.p_pos + dim_reco.p_neg;
		  ddr.mom_pos = dim_reco.p_pos;
		  ddr.mom_neg = dim_reco.p_neg;
		  ddr.x1      = dim_reco.x1;
		  ddr.x2      = dim_reco.x2;
	  }

	  mo_dim_reco.push_back(ddr);
  }

  tree->Fill();
  return Fun4AllReturnCodes::EVENT_OK;
}

int AnaSimDst::End(PHCompositeNode* topNode)
{
	file->cd();
	file->Write();
	file->Close();
	return Fun4AllReturnCodes::EVENT_OK;
}

int AnaSimDst::GetNodes(PHCompositeNode *topNode)
{
	mi_evt      = findNode::getClass<SQEvent       >(topNode, "SQEvent");
	mi_srec     = findNode::getClass<SRecEvent     >(topNode, "SRecEvent");
	mi_evt_true = findNode::getClass<SQMCEvent     >(topNode, "SQMCEvent");
	mi_vec_trk  = findNode::getClass<SQTrackVector >(topNode, "SQTruthTrackVector");
	mi_vec_dim  = findNode::getClass<SQDimuonVector>(topNode, "SQTruthDimuonVector");
	mi_sqhitvec    = findNode::getClass<SQHitVector     >(topNode, "SQHitVector"); // SQHitVector this class will give the event wise spectrometer hit information.
	if (!mi_evt || !mi_srec || !mi_evt_true || !mi_vec_trk || !mi_vec_dim || !mi_sqhitvec) {
		return Fun4AllReturnCodes::ABORTEVENT;
	}
	return Fun4AllReturnCodes::EVENT_OK;
}

void AnaSimDst::MakeTree()
{
	file = new TFile("sim_tree.root", "RECREATE");
	tree = new TTree("tree", "Created by AnaSimDst");

	//tree->Branch("sqevt"   ,  mi_evt);
	tree->Branch("evt"     , &mo_evt);
	tree->Branch("trk_true", &mo_trk_true);
	tree->Branch("trk_reco", &mo_trk_reco);
	tree->Branch("dim_true", &mo_dim_true);
	tree->Branch("dim_reco", &mo_dim_reco);
}

void AnaSimDst::FindTrackRelation(IdMap_t& id_map)
{

	id_map.clear();
	for (unsigned int i_true = 0; i_true < mi_vec_trk->size(); i_true++) {
		SQTrack* trk_true = mi_vec_trk->at(i_true);
		int     ch_true = trk_true->get_charge();
		double mom_true = trk_true->get_mom_vtx().Mag();

		int i_reco_best = -1;
		double mom_diff_best;
		for (int i_reco = 0; i_reco < mi_srec->getNTracks(); i_reco++) {
			SRecTrack* trk_reco = &mi_srec->getTrack(i_reco);
			if (trk_reco->getCharge() != ch_true) continue;
			double mom_diff = fabs(trk_reco->getMomentumVertex().Mag() - mom_true);
			if (i_reco_best < 0 || mom_diff < mom_diff_best) {
				i_reco_best   = i_reco;
				mom_diff_best = mom_diff;
			}
		}
		id_map[i_true] = i_reco_best;
	}
}

void AnaSimDst::FindDimuonRelation(IdMap_t& id_map)
{
	id_map.clear();
	for (unsigned int i_true = 0; i_true < mi_vec_dim->size(); i_true++) {
		SQDimuon* dim_true = mi_vec_dim->at(i_true);
		double mass_true = dim_true->get_mom().M();

		int i_reco_best = -1;
		double mass_diff_best;
		for (int i_reco = 0; i_reco < mi_srec->getNDimuons(); i_reco++) {
			SRecDimuon dim_reco = mi_srec->getDimuon(i_reco);
			double mass_diff = fabs(dim_reco.mass - mass_true);
			if (i_reco_best < 0 || mass_diff < mass_diff_best) {
				i_reco_best   = i_reco;
				mass_diff_best = mass_diff;
			}
		}
		id_map[i_true] = i_reco_best;
	}
}

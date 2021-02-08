#include "TFile.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TVirtualPad.h"
#include "TTree.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include <UtilAna/UtilDimuon.h>
//if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,0)
R__LOAD_LIBRARY(libana_sim_dst)
//endif
using namespace std;
namespace {
		int nfiles = 4;
		const char* inputs [] = {
/*
 // unoptimized
			"data/main/JPsi.root",    //JPsi data masscut at 1.5 GeV	
			"data/main/JPsi.root",    //Psi' data masscut at 1.5 GeV	
			"data/main/OC.root",	  //Open Charm
			"data/main/DY.root"	  //Drell-Yan masscut at 1.50 GeV
*/

			"opt_present/JPsi_optimized.root",        
                        "opt_present/JPsi_optimized.root",   
                        "opt_present/OC_optimized.root",      //Open Charm 4.0 GeV
                        "opt_present/DY_optimized.root"       //Drell-Yan masscut at 2.00 GeV
		};
}
void ana() {
	gStyle->SetOptFit();

	const int nbin = 14;
	const double xmin = -M_PI;
	const double xmax = M_PI;
	int tot_trk =0;

	TCanvas *c1 = new TCanvas("c1","c1"); c1->SetGrid();
	c1->SetLogy();

	TCanvas *c2 = new TCanvas("c2","c2"); c2->SetGrid();
	c2->SetLogy();

/*
	TCanvas *c3 = new TCanvas("c3","c3"); c3->SetGrid();
	c3->SetLogy();

	TCanvas *c32 = new TCanvas("c32","c32"); c32->SetGrid();
        c32->SetLogy();
	TCanvas *c4 = new TCanvas("c4","c4"); c4->SetGrid();
	c4->SetLogy();
	TCanvas *c5 = new TCanvas("c5","c5"); c5->SetGrid();
        c5->SetLogy();
	TCanvas *c6 = new TCanvas("c6","c6"); c6->SetGrid();
	c6->SetLogy();
	TCanvas *c7 = new TCanvas("c7","c7"); c7->SetGrid();
	c7->SetLogy();
	TCanvas *c8 = new TCanvas("c8","c8"); c8->SetGrid();
        c8->SetLogy();
	TCanvas *c9 = new TCanvas("c9","c9"); c9->SetGrid();
        c9->SetLogy();
	TCanvas *c10 = new TCanvas("c10","c10"); c10->SetGrid();
        c10->SetLogy();
	TCanvas *c11 = new TCanvas("c11","c11"); c11->SetGrid();
        c11->SetLogy();
	TCanvas *c12 = new TCanvas("c12","c12"); c12->SetGrid();
        c12->SetLogy();
        TCanvas *c13 = new TCanvas("c13","c13"); c13->SetGrid();
        c13->SetLogy();

	TCanvas *c14 = new TCanvas("c14","c14"); c14->SetGrid();
        c14->SetLogy();
	TCanvas *c15 = new TCanvas("c15","c15"); c15->SetGrid();
        c15->SetLogy();
	TCanvas *c16 = new TCanvas("c16","c16"); c16->SetGrid();
        c16->SetLogy();
	TCanvas *c17 = new TCanvas("c17","c17"); c17->SetGrid();
        c17->SetLogy();
*/
	int color = kBlack;
	auto legend = new TLegend(0.1,0.7,0.28,0.9);
	auto legend1 = new TLegend(0.1,0.7,0.28,0.9);
	auto legend2 = new TLegend(0.1,0.7,0.28,0.9);
	auto legend3 = new TLegend(0.1,0.7,0.28,0.9);
	auto legend32 = new TLegend(0.1,0.7,0.28,0.9);
	auto legend4 = new TLegend(0.1,0.7,0.28,0.9);
	auto legend5 = new TLegend(0.1,0.7,0.28,0.9);
	auto legend6 = new TLegend(0.1,0.7,0.28,0.9);
	auto legend7 = new TLegend(0.1,0.7,0.28,0.9);
	auto legend8 = new TLegend(0.1,0.7,0.28,0.9);
	auto legend9 = new TLegend(0.1,0.7,0.28,0.9);
	auto legend10 = new TLegend(0.1,0.7,0.28,0.9);
	auto legend11= new TLegend(0.1,0.7,0.28,0.9);
	auto legend12= new TLegend(0.1,0.7,0.28,0.9);
	auto legend13= new TLegend(0.1,0.7,0.28,0.9);
	auto legend14= new TLegend(0.1,0.7,0.28,0.9);
	auto legend15= new TLegend(0.1,0.7,0.28,0.9);
	auto legend16= new TLegend(0.1,0.7,0.28,0.9);
	auto legend17= new TLegend(0.1,0.7,0.28,0.9);

	TH1D *DY_copy;
        TH1D *JPsi_copy;
        TH1D *OC_copy;
	std::ofstream out_File;
	std::ofstream OCout_File;
	out_File.open("yieldP.tsv");
	OCout_File.open("OCout_File.csv");
	vector<int> duplicates[2];
	vector <double> yield_OC;
	vector <double> err_OC;
        vector <double> yield_DY;	
        vector <double> err_DY;	
        vector <double> bcx;	
        vector <double> binErrorNormOC;	
        vector <double> binErrorNormDY;	

	TH1D *trueM_true = new TH1D("trueM_true","trueM_true; trueM_true; percentage", 20, 0.0, 10.0);

	for (int j=0; j<nfiles; ++j) {

		TH1D *M_Px = new TH1D("M_Px","True M_Px; Normalized Yields",  100, -50, 50.0);
		TH1D *mupPz = new TH1D("mupPz","True #mu^{+} Pz; #mu^{+} Pz; Normalized Yields",  10, 0, 120.0);
		TH1D *mumPz = new TH1D("mumPz","True #mu^{-} Pz; #mu^{-}Pz; Normalized Yields",  10, 0, 120.0);
		TH1D *mupPzReco = new TH1D("mupPzReco","Reconstructed #mu^{+} Pz; #mu^{+} Pz; Normalized Yields",  20, 0, 120.0);
                TH1D *mumPzReco = new TH1D("mumPzReco","Reconstructed #mu^{-} Pz; #mu^{-} Pz; Normalized Yields",  20, 0, 120.0);
		TH1D *dimuPz = new TH1D("dimuPz","True Dimuon Pz; Pz; Normalized Yields", 20, 0, 140.0);
		TH1D *dimuPz_reco = new TH1D("dimuPz_reco","Reconstructed Dimuon Pz; Pz; Normalized Yields", 20, 0, 140.0);
		TH1D *dimuqT = new TH1D("dimuqT","True Dimuon q_{T}; q_{T}; Normalized Yields", 50, -1.5, 6.0);
		TH1D *dimuqT_reco = new TH1D("dimuqT_reco","Reconstructed Dimuon q_{T}; q_{T}; Normalized Yields", 50, -1.5, 6.0);
		TH1D *dimuxFtr = new TH1D("dimuxFtr","True Feynman x_{F}; x_{F}; Normalized Yields", 20, -0.5, 1.5);
		TH1D *dimuxFreco = new TH1D("dimuxFreco","Reconstructed Feynman x_{F}; x_{F}; Normalized Yields", 20, -0.5, 1.5);
		TH1D *dimuxB_tr = new TH1D("dimuxB_tr","True Bjorken Beam x_{B}; x_{B}; Normalized Yields", 20, -0.2, 1.5);
		TH1D *dimuxB_reco = new TH1D("dimuxB_reco","Reconstructed Bjorken Beam x_{B}; x_{B}; Normalized Yields", 20, -0.2, 1.5);
		TH1D *dimuxT_tr = new TH1D("dimuxT_tr","True Bjorken Target x_{T}; x_{T}; Normalized Yields", 20, -0.2, 0.8);
		TH1D *dimuxT_reco = new TH1D("dimuxT_reco","Reconstructed Bjorken Target x_{T}; x_{T}; Normalized Yields", 20, -0.2, 0.8);
		TH1D *trueM = new TH1D("trueM","trueM; M_{#mu^{+}#mu^{-}}; Normalized Yields", 20, 0.0, 10.0);
                TH1D *recoM = new TH1D("recoM","recoM; M_{#mu^{+}#mu^{-}}; Normalized Yields", 20, 0.0, 10.0);
		TH1D *costheta_tr = new TH1D("costheta_tr","True Dimuon cos#theta; cos#theta; Normalized Yields", 20, -1.5, 1.5);
		TH1D *costheta_reco = new TH1D("costheta_reco","Reco Dimuon cos#theta; cos#theta; Normalized Yields", 20, -1.5, 1.5);

		TFile *f = TFile::Open(inputs[j],"read");
		TTree *tr = (TTree*) f->Get("tree");

		EventData* ed = new EventData();
		TrackList* trkTrue = new TrackList();
		TrackList* trkReco = new TrackList();
		DimuonList* trueDimu = new DimuonList();
		DimuonList* recoDimu = new DimuonList();  
		tr->SetBranchAddress("evt", &ed);
		tr->SetBranchAddress("trk_true", &trkTrue);
		tr->SetBranchAddress("trk_reco", &trkReco);
		tr->SetBranchAddress("dim_true", &trueDimu);
		tr->SetBranchAddress("dim_reco", &recoDimu);
		TLorentzVector dimu_trk,dimu_trk_reco;
		int n_ent = tr->GetEntries();
		for (int i_ent = 0; i_ent < n_ent; i_ent++) {
			tr->GetEntry(i_ent);
			//making all the unique pairs of parent muons for OC decays
			//
			if(j==2){
				duplicates[0].push_back(ed->parent_id[0]);
				duplicates[1].push_back(ed->parent_id[1]);
			}

			//trigger conditions
			//if (!(ed->fpga1 ==true))continue;

			//particle selections
			if (j==0)if (!(ed->parent_id[0] ==443 && ed->parent_id[1]==443) ) continue;
			if (j==1)if (!(ed->parent_id[0] ==100443 && ed->parent_id[1]==100443) ) continue;

			// reconstruction status 
			if(ed->rec_stat != 0) continue; 
			if(j==0 ||j==1|| j==3)if(!(trueDimu->size() > 0)) continue;   // DY and JPsi dimu size has to be atleast 1

			double costh_tr, mup_pz_tr, mum_pz_tr, mup_E_tr, mum_E_tr, M_tr, Pt_tr,Pz_tr;
			double costh_reco, mup_pz_reco, mum_pz_reco, mup_E_reco, mum_E_reco,M_reco, Pt_reco,Pz_reco;
			double Px_tr;
			TLorentzVector dimu_reco, dimu_true;
			if (j==2){  //For OpenCharm True Dimuon Pair only
				if(!((trkTrue->at(0).charge >0.0 && trkTrue->at(1).charge<0.0) ||(trkTrue->at(0).charge <0.0 && trkTrue->at(1).charge>0.0)))continue;
				dimu_trk = trkTrue->at(0).mom_vtx + trkTrue->at(1).mom_vtx;
				dimu_trk_reco = trkReco->at(0).mom_vtx + trkReco->at(1).mom_vtx;
				dimu_reco = trkReco->at(0).mom_vtx + trkReco->at(1).mom_vtx;
				dimu_true = trkTrue->at(0).mom_vtx + trkTrue->at(1).mom_vtx;
				M_tr = dimu_trk.M();
				Pt_tr = dimu_trk.Pt();
				Pz_tr = dimu_trk.Pz();
				Px_tr = dimu_trk.Px();
				Pz_reco = dimu_trk_reco.Pz();
				M_reco = dimu_trk_reco.M();
				Pt_reco = dimu_trk_reco.Pt();
				if(trkTrue->at(0).charge >0.0){   //true open charm dimuon making 
					mup_pz_tr = trkTrue->at(0).mom_vtx.Pz();
					mum_pz_tr = trkTrue->at(1).mom_vtx.Pz();
					mup_E_tr =  trkTrue->at(0).mom_vtx.E();
					mum_E_tr =  trkTrue->at(1).mom_vtx.E();
				}
				else {
					mup_pz_tr = trkTrue->at(1).mom_vtx.Pz();
					mum_pz_tr = trkTrue->at(0).mom_vtx.Pz();
					mup_E_tr =  trkTrue->at(1).mom_vtx.E();
					mum_E_tr =  trkTrue->at(0).mom_vtx.E();
				}
				if(trkReco->at(0).charge >0.0){    ////reco open charm dimuon making 
					mup_pz_reco = trkReco->at(0).mom_vtx.Pz();
					mum_pz_reco = trkReco->at(1).mom_vtx.Pz();
					mup_E_reco =  trkReco->at(0).mom_vtx.E();
					mum_E_reco =  trkReco->at(1).mom_vtx.E();
				}
				else {
					mup_pz_reco = trkReco->at(1).mom_vtx.Pz();
					mum_pz_reco = trkReco->at(0).mom_vtx.Pz();
					mup_E_reco =  trkReco->at(1).mom_vtx.E();
					mum_E_reco =  trkReco->at(0).mom_vtx.E();
				}
			}
			else{  ////For JPsi/DrellYan True/Reco Dimuon Pair 
				mup_pz_tr = trueDimu->at(0).mom_pos.Pz();             
				mum_pz_tr = trueDimu->at(0).mom_neg.Pz();             
				mup_E_tr =  trueDimu->at(0).mom_pos.E();              
				mum_E_tr =  trueDimu->at(0).mom_neg.E();      
				M_tr =  trueDimu->at(0).mom.M();
				Pt_tr = trueDimu->at(0).mom.Pt();
				dimu_true = trueDimu->at(0).mom;
				Pz_tr = trueDimu->at(0).mom.Pz();
				Px_tr = trueDimu->at(0).mom.Px();

				mup_pz_reco = recoDimu->at(0).mom_pos.Pz();		
				mum_pz_reco = recoDimu->at(0).mom_neg.Pz();		
				mup_E_reco =  recoDimu->at(0).mom_pos.E();		
				mum_E_reco =  recoDimu->at(0).mom_neg.E();		
				M_reco = recoDimu->at(0).mom.M();
				Pt_reco = recoDimu->at(0).mom.Pt();
				dimu_reco = recoDimu->at(0).mom;
				Pz_reco = recoDimu->at(0).mom.Pz();
			}

			//if (M_reco < 4.40) continue; 

			costh_tr = 2.*(mum_E_tr*mup_pz_tr - mup_E_tr*mum_pz_tr)/((M_tr)*sqrt(M_tr*M_tr + Pt_tr*Pt_tr));
			costh_reco = 2.*(mum_E_reco*mup_pz_reco - mup_E_reco*mum_pz_reco)/((M_reco)*sqrt(M_reco*M_reco + Pt_reco*Pt_reco));
			//cout << " true costh "<< costh_tr<< " reco costh"<<costh_reco <<endl;

			//xf, and x_B def are adopted from SRecEvent just to be consistent
			Double_t mp = 0.938;
			Double_t ebeam = 120.;

			TLorentzVector p_beam(0., 0., sqrt(ebeam*ebeam - mp*mp), ebeam);
			TLorentzVector p_target(0., 0., 0., mp);

			TLorentzVector p_cms = p_beam + p_target;

			double x1_tr = (p_target*dimu_true)/(p_target*p_cms);
			double x1_reco = (p_target*dimu_reco)/(p_target*p_cms);
			double x2_tr = (p_beam*dimu_true)/(p_beam*p_cms);
			double x2_reco = (p_beam*dimu_reco)/(p_beam*p_cms);

			Double_t s = p_cms.M2();
			TVector3 bv_cms = p_cms.BoostVector();
			dimu_reco.Boost(-bv_cms);
			dimu_true.Boost(-bv_cms);
			double xF_tr = 2.*dimu_true.Pz()/TMath::Sqrt(s)/(1. - M_reco*M_reco/s);
			double xF_reco = 2.*dimu_reco.Pz()/TMath::Sqrt(s)/(1. - M_reco*M_reco/s);
			//SeaQuest Selections will be applied here: sinlge muon pz,Dimuon Pz  and collin sopper Cos\theta selections		
			//if(!(abs(costh_tr)>=0.5 && abs(costh_reco)>=0.5))continue;
			//if(!(abs(costh_tr)<=0.5 && (mup_pz_tr >=40.0 && mum_pz_tr>=40.0)))continue;
			//if(!(abs(costh_reco)<=0.5 && (mup_pz_reco >=30.0 && mum_pz_reco>=30.0)))continue;
			//if(!(abs(costh_reco)<=0.5 && (mup_pz_reco >=40.0 && mum_pz_reco>=40.0)))continue;
			//if(!(abs(x1)>=0.40))continue;
			//if(!(Pz_tr>=50.0))continue;

			//if(j==3)cout << "true Pz : reco Pz "<< Pz_tr << "\t" << Pz_reco<<endl;



			mupPz->Fill(mup_pz_tr);
			mumPz->Fill(mum_pz_tr);
			mupPzReco->Fill(mup_pz_reco);
			mumPzReco->Fill(mum_pz_reco);
			dimuPz->Fill(Pz_tr);
			dimuPz_reco->Fill(Pz_reco);
			costheta_tr->Fill(costh_tr);
			costheta_reco->Fill(costh_reco);	
			dimuxFtr->Fill(xF_tr);
			dimuxFreco->Fill(xF_reco);
			trueM->Fill(M_tr);
			recoM->Fill(M_reco);
			dimuxB_tr->Fill(x1_tr);
			dimuxB_reco->Fill(x1_reco);
			dimuxT_tr->Fill(x2_tr);
			dimuxT_reco->Fill(x2_reco);
			dimuqT->Fill(Pt_reco);
			dimuqT_reco->Fill(Pt_tr);
		}
		double norm; //scaling all the distributions by the luminosity numbers

		/*
		//begin unoptimized
		if (j==0|| j==1) norm = 1643.09; //JPsi, Psi'
		else if (j==2)norm = 49706.1; //Oopen charm
		else if (j==3) norm= 208474.0; //DY
		//end
		*/ 
		//begin-optimized
		if (j==0|| j==1) norm = 821.166; //JPsi, Psi'
                else if (j==2)norm = 10308.6; //Oopen charm
                else if (j==3) norm= 166416.0; //DY
		//
		dimuPz->Sumw2(kTRUE);
		dimuPz_reco->Sumw2(kTRUE);
		trueM->Sumw2(kTRUE);
		recoM->Sumw2(kTRUE);
		dimuxFtr->Sumw2(kTRUE);
		dimuxFreco->Sumw2(kTRUE);
		dimuxB_tr->Sumw2(kTRUE);
		dimuxB_reco->Sumw2(kTRUE);
		dimuxT_tr->Sumw2(kTRUE);
		dimuxT_reco->Sumw2(kTRUE);
		costheta_tr->Sumw2(kTRUE);
		costheta_reco->Sumw2(kTRUE);
		dimuqT_reco->Sumw2(kTRUE);
		dimuqT->Sumw2(kTRUE);
		mupPz->Sumw2(kTRUE);
		mumPz->Sumw2(kTRUE);
		mupPzReco->Sumw2(kTRUE);
		mumPzReco->Sumw2(kTRUE);

		trueM->Scale(1/norm);
		recoM->Scale(1/norm);
		dimuxFtr->Scale(1/norm);		
		dimuxFreco->Scale(1/norm);		
		dimuxB_tr->Scale(1/norm);
		dimuxB_reco->Scale(1/norm);
		dimuxT_tr->Scale(1/norm);
		dimuxT_reco->Scale(1/norm);
		costheta_tr->Scale(1/norm); 
		costheta_reco->Scale(1/norm); 
		dimuqT_reco->Scale(1/norm);
		dimuqT->Scale(1/norm);
		dimuPz->Scale(1/norm);
		dimuPz_reco->Scale(1/norm);
		mupPz->Scale(1/norm);
		mumPz->Scale(1/norm);
		mupPzReco->Scale(1/norm);
		mumPzReco->Scale(1/norm);


		TPaveStats *st[20];
		int col;
		if(j==0){
			//gStyle->SetOptStat(111111);
			color = kBlack;
			c1->cd();
			trueM->SetLineColor(color);
			trueM->SetMarkerColor(color);
			trueM->SetMarkerStyle(20);
			trueM->Draw("e");
			trueM->SetMinimum(10e-5);
			legend1->AddEntry(trueM, "J/#psi", "l");
			   c2->cd();
			   recoM->SetLineColor(color);
			   recoM->SetMarkerColor(color);
			   recoM->SetMarkerStyle(20);
			   recoM->Draw("e");
			   recoM->SetMinimum(10e-5);
			//trueM->SetMinimum(0.1);
			legend2->AddEntry(recoM, "J/#psi-reco", "l");
/*

			c3->cd();
			costheta_reco->SetLineColor(color);
			costheta_reco->SetMarkerColor(color);
			costheta_reco->SetMarkerStyle(20);
			costheta_reco->Draw("e");
			costheta_reco->SetMinimum(10e-5);
			legend3->AddEntry(costheta_reco, "J/#psi-reco", "l");


			c32->cd();
			costheta_tr->SetLineColor(color);
			costheta_tr->SetMarkerColor(color);
			costheta_tr->SetMarkerStyle(20);
			costheta_tr->Draw("e");
			costheta_tr->SetMinimum(10e-5);
			legend32->AddEntry(costheta_tr, "J/#psi", "l");

			c4->cd();
			dimuxFtr->SetLineColor(color);
			dimuxFtr->SetMarkerColor(color);
			dimuxFtr->SetMarkerStyle(20);
			dimuxFtr->Draw("e");
			dimuxFtr->SetMinimum(10e-5);
			legend4->AddEntry(dimuxFtr, "J/#psi", "l");

			c5->cd();
			dimuxFreco->SetLineColor(color);
			dimuxFreco->SetMarkerColor(color);
			dimuxFreco->SetMarkerStyle(20);
			dimuxFreco->Draw("e");
			dimuxFreco->SetMinimum(10e-5);
			legend5->AddEntry(dimuxFreco, "J/#psi-reco", "l");

			c6->cd();
			dimuxB_tr->SetLineColor(color);
			dimuxB_tr->SetMarkerColor(color);
			dimuxB_tr->SetMarkerStyle(20);
			dimuxB_tr->Draw("e");
			dimuxB_tr->SetMinimum(10e-5);
			legend6->AddEntry(dimuxB_tr, "J/#psi", "l");

			c7->cd();
			dimuxB_reco->SetLineColor(color);
			dimuxB_reco->SetMarkerColor(color);
			dimuxB_reco->SetMarkerStyle(20);
			dimuxB_reco->Draw("e");
			dimuxB_reco->SetMinimum(10e-5);
			legend7->AddEntry(dimuxB_reco, "J/#psi-reco", "l");

			c8->cd();
			dimuxT_reco->SetLineColor(color);
			dimuxT_reco->SetMarkerColor(color);
			dimuxT_reco->SetMarkerStyle(20);
			dimuxT_reco->Draw("e");
			dimuxT_reco->SetMinimum(10e-5);
			legend8->AddEntry(dimuxT_reco, "J/#psi-reco", "l");

			c9->cd();
			dimuxT_tr->SetLineColor(color);
			dimuxT_tr->SetMarkerColor(color);
			dimuxT_tr->SetMarkerStyle(20);
			dimuxT_tr->Draw("e");
			dimuxT_tr->SetMinimum(10e-5);
			legend9->AddEntry(dimuxT_tr, "J/#psi", "l");

			c10->cd();
			dimuqT->SetLineColor(color);
			dimuqT->SetMarkerColor(color);
			dimuqT->SetMarkerStyle(20);
			dimuqT->Draw("e");
			dimuqT->SetMinimum(10e-5);
			legend10->AddEntry(dimuqT, "J/#psi", "l");

			c11->cd();
			dimuqT_reco->SetLineColor(color);
			dimuqT_reco->SetMarkerColor(color);
			dimuqT_reco->SetMarkerStyle(20);
			dimuqT_reco->Draw("e");
			dimuqT_reco->SetMinimum(10e-5);
			legend11->AddEntry(dimuqT_reco, "J/#psi-reco", "l");

			c12->cd();
			dimuPz->SetLineColor(color);
			dimuPz->SetMarkerColor(color);
			dimuPz->SetMarkerStyle(20);
			dimuPz->Draw("e");
			dimuPz->SetMinimum(10e-5);
			legend12->AddEntry(dimuPz, "J/#psi", "l");

			c13->cd();
			dimuPz_reco->SetLineColor(color);
			dimuPz_reco->SetMarkerColor(color);
			dimuPz_reco->SetMarkerStyle(20);
			dimuPz_reco->Draw("e");
			dimuPz_reco->SetMinimum(10e-5);
			legend13->AddEntry(dimuPz_reco, "J/#psi-reco", "l");


			c14->cd();
                        mupPz->SetLineColor(color);
                        mupPz->SetMarkerColor(color);
                        mupPz->SetMarkerStyle(20);
                        mupPz->Draw("e");
                        mupPz->SetMinimum(10e-5);
                        legend14->AddEntry(mupPz, "J/#psi", "l");


			c15->cd();
                        mumPz->SetLineColor(color);
                        mumPz->SetMarkerColor(color);
                        mumPz->SetMarkerStyle(20);
                        mumPz->Draw("e");
                        mumPz->SetMinimum(10e-5);
                        legend15->AddEntry(mumPz, "J/#psi", "l");

			c16->cd();
                        mupPzReco->SetLineColor(color);
                        mupPzReco->SetMarkerColor(color);
                        mupPzReco->SetMarkerStyle(20);
                        mupPzReco->Draw("e");
                        mupPzReco->SetMinimum(10e-5);
                        legend16->AddEntry(mupPzReco, "J/#psi-reco", "l");


                        c17->cd();
                        mumPzReco->SetLineColor(color);
                        mumPzReco->SetMarkerColor(color);
                        mumPzReco->SetMarkerStyle(20);
                        mumPzReco->Draw("e");
                        mumPzReco->SetMinimum(10e-5);
                        legend17->AddEntry(mumPzReco, "J/#psi-reco", "l");
*/

		}
		else {
			//gStyle->SetOptStat(111111);
			color = j+1;
			c1->cd();
			trueM->SetLineColor(color);
			trueM->SetMarkerColor(color);
			trueM->SetMarkerStyle(20);
			trueM->Draw("e sames");
			//trueM->SetMinimum(0.1);
			if(j==1)legend1->AddEntry(trueM, "#psi'", "l");
			if(j==2)legend1->AddEntry(trueM, "Open-Charm", "l");
			if(j==3)legend1->AddEntry(trueM, "Drell-Yan", "l");
			if (st[j] != NULL){
				gPad->Update();
				st[j] = (TPaveStats*)trueM->FindObject("stats");
				st[j]->SetTextColor(color);
				st[j]->SetY1NDC(.695-(j-1)*.200);
				st[j]->SetY2NDC(.495-(j-1)*.200);
			}
			legend1->Draw();
			if(j==3)c1->SaveAs("trueMass.png");


			c2->cd();
			recoM->SetLineColor(color);
			recoM->SetMarkerColor(color);
			recoM->SetMarkerStyle(20);
			recoM->Draw("e sames");
			if(j==1)legend2->AddEntry(recoM, "#psi'-reco", "l");
			if(j==2)legend2->AddEntry(recoM, "Open-Charm-reco", "l");
			if(j==3)legend2->AddEntry(recoM, "Drell-Yan-reco", "l");
			if (st[j] != NULL){
				gPad->Update();
				st[j] = (TPaveStats*)recoM->FindObject("stats");
				st[j]->SetTextColor(color);
				st[j]->SetY1NDC(.695-(j-1)*.200);
				st[j]->SetY2NDC(.495-(j-1)*.200);

			}
			legend2->Draw();	
			if(j==3) c2->SaveAs("recoMass.png");
			/*

			   c3->cd();
			   costheta_reco->SetLineColor(color);
			   costheta_reco->SetMarkerColor(color);
			   costheta_reco->SetMarkerStyle(20);
			   costheta_reco->Draw("esames");
			   if(j==1)legend3->AddEntry(costheta_reco, "#psi'-reco", "l");
			   if(j==2)legend3->AddEntry(costheta_reco, "Open-Charm-reco", "l");
			   if(j==3)legend3->AddEntry(costheta_reco, "Drell-Yan-reco", "l");
			   if (st[j] != NULL){
			   gPad->Update();
			   st[j] = (TPaveStats*)costheta_reco->FindObject("stats");
			   st[j]->SetTextColor(color);
			   st[j]->SetY1NDC(.695-(j-1)*.200);
			   st[j]->SetY2NDC(.495-(j-1)*.200);
			   }
			   legend3->Draw();
			   if(j==3)c3->SaveAs("recoCostheta.png");



			   c32->cd();
			   costheta_tr->SetLineColor(color);
			   costheta_tr->SetMarkerColor(color);
			   costheta_tr->SetMarkerStyle(20);
			   costheta_tr->Draw("esames");
			   if(j==1)legend32->AddEntry(costheta_tr, "#psi'", "l");
			   if(j==2)legend32->AddEntry(costheta_tr, "Open-Charm", "l");
			   if(j==3)legend32->AddEntry(costheta_tr, "Drell-Yan", "l");
			   if (st[j] != NULL){
			   gPad->Update();
			   st[j] = (TPaveStats*)costheta_tr->FindObject("stats");
			   st[j]->SetTextColor(color);
			   st[j]->SetY1NDC(.695-(j-1)*.200);
			   st[j]->SetY2NDC(.495-(j-1)*.200);
			   }
			   legend32->Draw();
			   if(j==3)c32->SaveAs("trueCostheta.png");


			   c4->cd();
			   dimuxFtr->SetLineColor(color);
			   dimuxFtr->SetMarkerColor(color);
			   dimuxFtr->SetMarkerStyle(20);
			   dimuxFtr->Draw("e sames");

			   if(j==1)legend4->AddEntry(dimuxFtr, "#psi'", "l");
			   if(j==2)legend4->AddEntry(dimuxFtr, "Open-Charm", "l");
			   if(j==3)legend4->AddEntry(dimuxFtr, "Drell-Yan", "l");
			   if (st[j] != NULL){
			   gPad->Update();
			   st[j] = (TPaveStats*)dimuxFtr->FindObject("stats");
			   st[j]->SetTextColor(color);
			   st[j]->SetY1NDC(.695-(j-1)*.200);
			   st[j]->SetY2NDC(.495-(j-1)*.200);
			   }
			   legend4->Draw();
			   if(j==3)c4->SaveAs("trueDimuxF.png");

			   c5->cd();
			   dimuxFreco->SetLineColor(color);
			   dimuxFreco->SetMarkerColor(color);
			   dimuxFreco->SetMarkerStyle(20);
			   dimuxFreco->Draw("e sames");
			   if(j==1)legend5->AddEntry(dimuxFreco, "#psi'-reco", "l");
			   if(j==2)legend5->AddEntry(dimuxFreco, "Open-Charm-reco", "l");
			   if(j==3)legend5->AddEntry(dimuxFreco, "Drell-Yan-reco", "l");
			   if (st[j] != NULL){
			   gPad->Update();
			   st[j] = (TPaveStats*)dimuxFreco->FindObject("stats");
			st[j]->SetTextColor(color);
			st[j]->SetY1NDC(.695-(j-1)*.200);
			st[j]->SetY2NDC(.495-(j-1)*.200);
		}
		legend5->Draw();
		if(j==3)c5->SaveAs("recoDimuxF.png");


		c6->cd();
		dimuxB_tr->SetLineColor(color);
		dimuxB_tr->SetMarkerColor(color);
		dimuxB_tr->SetMarkerStyle(20);
		dimuxB_tr->Draw("e sames");
		if(j==1)legend6->AddEntry(dimuxB_tr, "#psi'", "l");
		if(j==2)legend6->AddEntry(dimuxB_tr, "Open-Charm", "l");
		if(j==3)legend6->AddEntry(dimuxB_tr, "Drell-Yan", "l");
		if (st[j] != NULL){
			gPad->Update();
			st[j] = (TPaveStats*)dimuxB_tr->FindObject("stats");
			st[j]->SetTextColor(color);
			st[j]->SetY1NDC(.695-(j-1)*.200);
			st[j]->SetY2NDC(.495-(j-1)*.200);
		}
		legend6->Draw();
		if(j==3)c6->SaveAs("trueDimuxB.png");


		c7->cd();
		dimuxB_reco->SetLineColor(color);
		dimuxB_reco->SetMarkerColor(color);
		dimuxB_reco->SetMarkerStyle(20);
		dimuxB_reco->Draw("e sames");
		if(j==1)legend7->AddEntry(dimuxB_reco, "#psi'-reco", "l");
		if(j==2)legend7->AddEntry(dimuxB_reco, "Open-Charm-reco", "l");
		if(j==3)legend7->AddEntry(dimuxB_reco, "Drell-Yan-reco", "l");
		if (st[j] != NULL){
			gPad->Update();
			st[j] = (TPaveStats*)dimuxB_reco->FindObject("stats");
			st[j]->SetTextColor(color);
			st[j]->SetY1NDC(.695-(j-1)*.200);
			st[j]->SetY2NDC(.495-(j-1)*.200);
		}
		legend7->Draw();
		if(j==3)c7->SaveAs("recoDimuxB.png");

		c8->cd();
		dimuxT_reco->SetLineColor(color);
		dimuxT_reco->SetMarkerColor(color);
		dimuxT_reco->SetMarkerStyle(20);
		dimuxT_reco->Draw("e sames");
		if(j==1)legend8->AddEntry(dimuxT_reco, "#psi'-reco", "l");
		if(j==2)legend8->AddEntry(dimuxT_reco, "Open-Charm-reco", "l");
		if(j==3)legend8->AddEntry(dimuxT_reco, "Drell-Yan -reco", "l");
		if (st[j] != NULL){
			gPad->Update();
			st[j] = (TPaveStats*)dimuxT_reco->FindObject("stats");
			st[j]->SetTextColor(color);
			st[j]->SetY1NDC(.695-(j-1)*.200);
			st[j]->SetY2NDC(.495-(j-1)*.200);
		}
		legend8->Draw();
		if(j==3)c8->SaveAs("recoDimuxT.png");

		c9->cd();
		dimuxT_tr->SetLineColor(color);
		dimuxT_tr->SetMarkerColor(color);
		dimuxT_tr->SetMarkerStyle(20);
		dimuxT_tr->Draw("e sames");
		if(j==1)legend9->AddEntry(dimuxT_tr, "#psi'", "l");
		if(j==2)legend9->AddEntry(dimuxT_tr, "Open-Charm", "l");
		if(j==3)legend9->AddEntry(dimuxT_tr, "Drell-Yan", "l");
		if (st[j] != NULL){
			gPad->Update();
			st[j] = (TPaveStats*)dimuxT_tr->FindObject("stats");
			st[j]->SetTextColor(color);
			st[j]->SetY1NDC(.695-(j-1)*.200);
			st[j]->SetY2NDC(.495-(j-1)*.200);
		}
		legend9->Draw();
		if(j==3)c9->SaveAs("trueDimuxT.png");


		c10->cd();
		dimuqT->SetLineColor(color);
		dimuqT->SetMarkerColor(color);
		dimuqT->SetMarkerStyle(20);
		dimuqT->Draw("esames");
		if(j==1)legend10->AddEntry(dimuqT, "#psi'", "l");
		if(j==2)legend10->AddEntry(dimuqT, "Open-Charm", "l");
		if(j==3)legend10->AddEntry(dimuqT, "Drell-Yan", "l");
		if (st[j] != NULL){
			gPad->Update();
			st[j] = (TPaveStats*)dimuqT->FindObject("stats");
			st[j]->SetTextColor(color);
			st[j]->SetY1NDC(.695-(j-1)*.200);
			st[j]->SetY2NDC(.495-(j-1)*.200);
		}
		legend10->Draw();
		if(j==3)c10->SaveAs("trueDimuqT.png");


		c11->cd();
		dimuqT_reco->SetLineColor(color);
		dimuqT_reco->SetMarkerColor(color);
		dimuqT_reco->SetMarkerStyle(20);
		dimuqT_reco->Draw("esames");
		if(j==1)legend11->AddEntry(dimuqT_reco, "#psi'-reco", "l");
		if(j==2)legend11->AddEntry(dimuqT_reco, "Open-Charm-reco", "l");
		if(j==3)legend11->AddEntry(dimuqT_reco, "Drell-Yan-reco", "l");
		if (st[j] != NULL){
			gPad->Update();
			st[j] = (TPaveStats*)dimuqT_reco->FindObject("stats");
			st[j]->SetTextColor(color);
			st[j]->SetY1NDC(.695-(j-1)*.200);
			st[j]->SetY2NDC(.495-(j-1)*.200);
		}
		legend11->Draw();
		if(j==3)c11->SaveAs("recoDimuqT.png");


		c12->cd();
		dimuPz->SetLineColor(color);
		dimuPz->SetMarkerColor(color);
		dimuPz->SetMarkerStyle(20);
		dimuPz->Draw("esames");
		if(j==1)legend12->AddEntry(dimuPz, "#psi'", "l");
		if(j==2)legend12->AddEntry(dimuPz, "Open-Charm", "l");
		if(j==3)legend12->AddEntry(dimuPz, "Drell-Yan", "l");
		if (st[j] != NULL){
			gPad->Update();
			st[j] = (TPaveStats*)dimuPz->FindObject("stats");
			st[j]->SetTextColor(color);
			st[j]->SetY1NDC(.695-(j-1)*.200);
			st[j]->SetY2NDC(.495-(j-1)*.200);
		}
		legend12->Draw();
		if(j==3)c12->SaveAs("trueDimuPz.png");


		c13->cd();
		dimuPz_reco->SetLineColor(color);
		dimuPz_reco->SetMarkerColor(color);
		dimuPz_reco->SetMarkerStyle(20);
		dimuPz_reco->Draw("esames");
		if(j==1)legend13->AddEntry(dimuPz_reco, "#psi'-reco", "l");
		if(j==2)legend13->AddEntry(dimuPz_reco, "Open-Charm-reco", "l");
		if(j==3)legend13->AddEntry(dimuPz_reco, "Drell-Yan-reco", "l");
		if (st[j] != NULL){
			gPad->Update();
			st[j] = (TPaveStats*)dimuPz_reco->FindObject("stats");
			st[j]->SetTextColor(color);
			st[j]->SetY1NDC(.695-(j-1)*.200);
			st[j]->SetY2NDC(.495-(j-1)*.200);
		}
		legend13->Draw();
		if(j==3)c13->SaveAs("recoDimuPz.png");

		c14->cd();
		mupPz->SetLineColor(color);
		mupPz->SetMarkerColor(color);
		mupPz->SetMarkerStyle(20);
		mupPz->Draw("esames");

		if(j==1)legend14->AddEntry(mupPz, "#psi'-reco", "l");
		if(j==2)legend14->AddEntry(mupPz, "Open-Charm-reco", "l");
		if(j==3)legend14->AddEntry(mupPz, "Drell-Yan-reco", "l");
		if (st[j] != NULL){
			gPad->Update();
			st[j] = (TPaveStats*)mupPz->FindObject("stats");
			st[j]->SetTextColor(color);
			st[j]->SetY1NDC(.695-(j-1)*.200);
			st[j]->SetY2NDC(.495-(j-1)*.200);
		}
		legend14->Draw();
		if(j==3)c14->SaveAs("truemupPz.png");

		c15->cd();
		mumPz->SetLineColor(color);
		mumPz->SetMarkerColor(color);
		mumPz->SetMarkerStyle(20);
		mumPz->Draw("esames");

		if(j==1)legend15->AddEntry(mumPz, "#psi'-reco", "l");
		if(j==2)legend15->AddEntry(mumPz, "Open-Charm-reco", "l");
		if(j==3)legend15->AddEntry(mumPz, "Drell-Yan-reco", "l");
		if (st[j] != NULL){
			gPad->Update();
			st[j] = (TPaveStats*)mumPz->FindObject("stats");
			st[j]->SetTextColor(color);
			st[j]->SetY1NDC(.695-(j-1)*.200);
			st[j]->SetY2NDC(.495-(j-1)*.200);
		}
		legend15->Draw();
		if(j==3)c15->SaveAs("truemumPz.png");



		c16->cd();
		mupPzReco->SetLineColor(color);
		mupPzReco->SetMarkerColor(color);
		mupPzReco->SetMarkerStyle(20);
		mupPzReco->Draw("esames");

		if(j==1)legend16->AddEntry(mupPzReco, "#psi'-reco", "l");
		if(j==2)legend16->AddEntry(mupPzReco, "Open-Charm-reco", "l");
		if(j==3)legend16->AddEntry(mupPzReco, "Drell-Yan-reco", "l");
		if (st[j] != NULL){
			gPad->Update();
			st[j] = (TPaveStats*)mupPzReco->FindObject("stats");
			st[j]->SetTextColor(color);
			st[j]->SetY1NDC(.695-(j-1)*.200);
			st[j]->SetY2NDC(.495-(j-1)*.200);
		}
		legend16->Draw();
		if(j==3)c16->SaveAs("recomupPz.png");


		c17->cd();
		mumPzReco->SetLineColor(color);
		mumPzReco->SetMarkerColor(color);
		mumPzReco->SetMarkerStyle(20);
		mumPzReco->Draw("esames");

		if(j==1)legend17->AddEntry(mumPzReco, "#psi'-reco", "l");
		if(j==2)legend17->AddEntry(mumPzReco, "Open-Charm-reco", "l");
		if(j==3)legend17->AddEntry(mumPzReco, "Drell-Yan-reco", "l");
		if (st[j] != NULL){
			gPad->Update();
			st[j] = (TPaveStats*)mumPzReco->FindObject("stats");
			st[j]->SetTextColor(color);
			st[j]->SetY1NDC(.695-(j-1)*.200);
			st[j]->SetY2NDC(.495-(j-1)*.200);
		}
		legend17->Draw();
		if(j==3)c17->SaveAs("recomumPz.png");
		*/
		}
	}
	out_File.close();
	OCout_File.close();
}

#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TStyle.h"
#include "TGraphErrors.h"
//if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,0)
R__LOAD_LIBRARY(libana_sim_dst)
//endif
using namespace std;

namespace {
	int nfiles = 1;
	const char* inputs [] = {
		"data/noOffset_june27.root"
	//	"data/offSet2cm_june30.root"
	};
	float config_type[] = {0, 1,2,3,4};
	float config_type_error[] = {0, 0,0,0,0};
}


void ana_beamOffset() {
	gStyle->SetOptFit();

  const int nbin = 14;
  const double xmin = -M_PI;
  const double xmax = M_PI;

    //TCanvas *c0 = new TCanvas("c0","c0"); c0->SetGrid();
  TCanvas *c1 = new TCanvas("c1","c1"); c1->SetGrid();
  TCanvas *c2 = new TCanvas("c2","c2"); c2->SetGrid();
  //TCanvas *c3 = new TCanvas("c3","c3"); c3->SetGrid();
  //TCanvas *c4 = new TCanvas("c4","c4"); c3->SetGrid();

  int color = kBlack;
  for (int j=0; j<nfiles; ++j) {

	  TH1D *phi_diff = new TH1D("phi_diff","phi_diff; (reco-true)Phi; Counts", nbin, 2*xmin, 2*xmax);
	  TH1D *phi_true = new TH1D("phi_true","TruePhi; TruePhi; Counts", nbin, xmin, xmax);
	  TH1D *phi_reco = new TH1D("phi_reco","RecoPhi; RecoPhi; Counts", nbin, xmin, xmax);
	  TH1D *recoM = new TH1D("recoM","RecoM; RecoM; Counts", nbin, 0.0, 8.0);
	  TH1D *trueM = new TH1D("trueM","RecoM; RecoM; Counts", nbin, 0.0, 8.0);
	  TH1D *pos_phi = new TH1D("pos_phi","pos_phi; pos_phi; Counts", nbin, xmax, xmin);
	  TH1D *neg_phi = new TH1D("pos_phi","neg_phi; neg_phi; Counts", nbin, xmax, xmin);
	  TFile *f = TFile::Open(inputs[j],"read");
	  TTree *tr = (TTree*) f->Get("tree");

	  EventData* ed = new EventData();
	  DimuonList* trueDimu = new DimuonList();
	  DimuonList* recoDimu = new DimuonList();  
	  tr->SetBranchAddress("evt", &ed);
	  tr->SetBranchAddress("dim_true", &trueDimu);
	  tr->SetBranchAddress("dim_reco", &recoDimu);

	  int n_ent = tr->GetEntries();
	  //for (int i_ent = 0; i_ent < n_ent; i_ent++) {
	  for (int i_ent = 0; i_ent < 750000; i_ent++) {
		  tr->GetEntry(i_ent);


		  if(ed->rec_stat != 0) continue;
		  if(!(trueDimu->size() > 0)) continue;
		  if(!(recoDimu->at(0).mom.M() > 4.0)) continue;
		  cout << "==================="<<endl;
		  cout <<"true dimu phi " << trueDimu->at(0).mom.Phi()<<endl;
		  cout <<"reco dimu phi " << recoDimu->at(0).mom.Phi()<<endl;
		  phi_diff->Fill(  trueDimu->at(0).mom.Phi() -  recoDimu->at(0).mom.Phi());
		  phi_true->Fill(trueDimu->at(0).mom.Phi());
		  phi_reco->Fill(recoDimu->at(0).mom.Phi());
		  trueM->Fill(trueDimu->at(0).mom.M());
		  recoM->Fill( recoDimu->at(0).mom.M());
	  }


	  if(j==0){
		  gStyle->SetOptStat(111111);
		  color = kBlack;
		  c1->cd();
		  phi_diff->SetLineColor(kBlack);
		  phi_diff->SetMarkerColor(kBlack);
		  phi_diff->SetMarkerStyle(20);
		  phi_diff->Draw("e");

		  TF1 *cosFit = new TF1("cosFit","[0]*([1]*cos(2*x)+1)");
		  c2->cd();
                  phi_reco->SetLineColor(kBlack);
                  phi_reco->SetMarkerColor(kBlack);
                  phi_reco->SetMarkerStyle(20);
                  phi_reco->Draw("e");
		  phi_reco->SetMinimum(0.0);
		  phi_reco->Fit("cosFit");
                  cosFit->SetLineColor(1);
	  }
	  else
	  {
		  gStyle->SetOptStat(111111);
		  color = j+1;
		  c1->cd();
		  phi_diff->SetLineColor(color);
		  phi_diff->SetMarkerColor(color);
		  phi_diff->SetMarkerStyle(20);
		  phi_diff->Draw("e sames");
		  //TF1 *fcos = new TF1("fcos","([0]*cos(2*x)+[1])");
		  //fcos->SetLineColor(color);

		  TF1 *cosFit = new TF1("cosFit","[0]*([1]*cos(2*x)+1)");
		  c1->SaveAs("phidiff.png");



		  c2->cd();
                  phi_reco->SetLineColor(color);
                  phi_reco->SetMarkerColor(color);
                  phi_reco->SetMarkerStyle(20);
                  phi_reco->Draw("esames");
		  phi_reco->Fit("cosFit");
		  phi_reco->SetMinimum(0.0);
		  cosFit->SetLineColor(2);
		  c2->SaveAs("phi_reco.png");
	  }

  }

}


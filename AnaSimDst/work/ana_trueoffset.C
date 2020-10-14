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

double ratio_error(
    const double a,
    const double b,
    const double ea,
    const double eb
    ) {
  double r = a/b;
  double er = r*sqrt(
      (ea/a)*(ea/a) + (eb/b)*(eb/b)
      );
  return er;
}

double binom_error(
    const double a,
    const double b
    ){
  double r=a/b;
  return sqrt(r*(1-r)/b);
}
TH1D * getEffHist(
    const char* hname,
    const TH1D* h1,
    const TH1D* h2
    ) {
  TH1D *h = (TH1D*)h1->Clone(hname);
  int nbin = h->GetNbinsX();
  for(int ibin=1; ibin<=nbin; ++ibin) {
    double a = h1->GetBinContent(ibin);
    double ea = sqrt(h1->GetBinContent(ibin));
    double b = h2->GetBinContent(ibin);
    double eb = sqrt(h2->GetBinContent(ibin));
    double r = a/b;
    //double e = binom_error(a, b);
    double e = ratio_error(a, b, ea, eb);
    h->SetBinContent(ibin, r);
    h->SetBinError(ibin, e);
  }

  return h;
}

namespace {
	int nfiles = 6;
	const char* inputs [] = {
		"data/noOffset.root",
 		"data/Yoffset0.2cmVtxKnown.root",
 		"data/Yoffset0.5cmVtxKnown.root",
		"data/Yoffset1.5cmVtxKnown.root",
		"data/Yoffset1cmVtxKnown.root",
		"data/Yoffset2.0cmVtxKnown.root"
	};
}
void ana_trueoffset() {
	gStyle->SetOptFit();

	const int nbin = 14;
	const double xmin = -M_PI;
	const double xmax = M_PI;
	gStyle->SetStatX(0.65);
	gStyle->SetStatY(0.45);
	int color = kBlack;
	int tot_trk =0;
	int hist_n =6;
        TH1D* true_phi[hist_n];
        char buffer[20];
	std::string Names[6] = {"noOffset", "Yoffset0.2cm", "Yoffset0.5cm", "Yoffset1.0cm", "Yoffset1.5cm","Yoffset2.0cm"};

	TH1D *phi_ratio[hist_n];
	for(int i = 0; i <6; ++i) {
		sprintf(buffer, "%s_true_phi", Names[i].c_str());
		true_phi[i] = new TH1D(buffer, buffer, 20, -M_PI, M_PI);
	}

	TH1D *truePhi = new TH1D("truePhi","truePhi; truePhi; Counts", 20, -M_PI, M_PI);

	for (int j=0; j<nfiles; ++j) {
		TFile *f = TFile::Open(inputs[j],"read");
		TTree *tr = (TTree*) f->Get("tree");

		EventData* ed = new EventData();
		DimuonList* trueDimu = new DimuonList();
		DimuonList* recoDimu = new DimuonList();
		TrackList* true_trk = new TrackList;

		tr->SetBranchAddress("evt", &ed);
		tr->SetBranchAddress("dim_true", &trueDimu);
		tr->SetBranchAddress("trk_true", &true_trk);
		tr->SetBranchAddress("dim_reco", &recoDimu);
		int n_ent = tr->GetEntries();
		for (int i_ent = 0; i_ent < n_ent; i_ent++) {
		//for (int i_ent = 0; i_ent < 100028; i_ent++) {
			tr->GetEntry(i_ent);
			if(ed->rec_stat != 0) continue;
			if(!(trueDimu->size() > 0)) continue;
			if(!(recoDimu->at(0).mom.M() > 4.0)) continue;
			if(!(recoDimu->at(0).mom.M() >= 4.0)) continue;
			
			//hodo-acceptance applied
			if (true_trk->at(0).nhodo <8 )continue;
			if (true_trk->at(0).nhodo >=8 )tot_trk++;
			if (true_trk->at(1).nhodo <8 )continue;
			if (true_trk->at(1).nhodo >=8 )tot_trk++;
			if (!(tot_trk >=2)) continue;
			tot_trk =0; 
			
			truePhi->Fill(trueDimu->at(0).mom.Phi());
			if (j==0)true_phi[j]->Fill(trueDimu->at(0).mom.Phi());
			if (j==1)true_phi[j]->Fill(trueDimu->at(0).mom.Phi());
			if (j==2)true_phi[j]->Fill(trueDimu->at(0).mom.Phi());
			if (j==3)true_phi[j]->Fill(trueDimu->at(0).mom.Phi());
			if (j==4)true_phi[j]->Fill(trueDimu->at(0).mom.Phi());
			if (j==5)true_phi[j]->Fill(trueDimu->at(0).mom.Phi());
		}
	}
	gStyle->SetOptTitle(0);
		
	//scaling the statbox
	gStyle->SetStatW(0.2);
	gStyle->SetStatH(0.1);
	gStyle->SetOptStat(0000000000111);
	TF1 *cosFit = new TF1("cosFit","[0]*([1]*cos(2*x)+1)");
	TCanvas *c4 = new TCanvas("c4","c4",1024,600);
	c4->Divide(3,2,0.0008,0.0008);

	//beamoffset-dimuon-phi
	for (int j=0; j<=5; j++){
		c4->cd(j+1);
		true_phi[j]->Draw("e");			
		true_phi[j]->SetLineColor(kBlack);	
		true_phi[j]->SetMarkerColor(kBlack);
		true_phi[j]->SetMarkerStyle(20);
		true_phi[j]->SetMinimum(0);
		true_phi[j]->SetAxisRange(0, 2000,"Y");
		true_phi[j]->Fit("cosFit");			
	}
	c4->SaveAs("offset.png");
	c4->SaveAs("offset.pdf");

	//ratio plot
	TCanvas *c5 = new TCanvas("c5","c5",1024,600);
	c5->Divide(3,2);
	ostringstream oss;
	for (int j=0; j<=5; j++){
		c5->cd(j+1);
		true_phi[j]->Scale(true_phi[0]->GetEntries() /true_phi[j]->GetEntries());
		true_phi[0]->Sumw2();
		true_phi[j]->Sumw2();
		phi_ratio[j] = getEffHist("phi_ratio[j]",
				true_phi[0],
				true_phi[j]
				);
		oss.str("");
		oss << true_phi[j]->GetName();
		cout << "j:  "<< j << " conditions "<< true_phi[j]->GetName() << endl;
		phi_ratio[j]->SetName(oss.str().c_str());
		phi_ratio[j]->Draw("e");
		phi_ratio[j]->SetMinimum(-.20);
		phi_ratio[j]->SetMaximum(2.0);
		phi_ratio[j]->SetMarkerStyle(20);
		phi_ratio[j]->SetMarkerColor(kBlack);
		phi_ratio[j]->Draw("e");
		phi_ratio[j]->Fit("pol0");
	}
	c5->SaveAs("ratio_offset.pdf");
	}

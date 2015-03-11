//#include <iostream>
//#include <math.h>
//#include "TROOT.h"
//#include "TFile.h"
//#include "TTree.h"
//#include "TH1F.h"
//#include "TCanvas.h"

using namespace std;

double deltaR(double eta1, double phi1, double eta2, double phi2){
  double dEta = eta1 - eta2;

  double dPhi = phi1 - phi2;
  while (dPhi > TMath::Pi()) dPhi -= 2*TMath::Pi();
  while (dPhi <= -TMath::Pi()) dPhi += 2*TMath::Pi();

  double dR = sqrt(dEta*dEta + dPhi*dPhi);
  return dR;
}


TH1F* makeHisto(std::string fname){
  //{
  bool display = 0;

  //TFile *f = new TFile("RadionHH_M800_R12_r15_correctIVFmatch_mc_subjets.root");
  TFile *f = new TFile((fname).c_str());
  TTree *tf = f->Get("analysis");

  //TCanvas *c = new TCanvas("c","SD fatjet higgs matching",800,600);
  TH1F *h = new TH1F(fname.c_str(),fname.c_str(),40,-18,-2);

  const int nEvent = tf->GetEntries();

  double Fj_chi;

  tf->SetBranchAddress("h_sb",&Fj_chi);

  //loop over events
  if(display)cout << endl;
  for (int i =0 ; i<nEvent;i++){
    //for (int i =0 ; i<2;i++){
    tf->GetEntry(i);
    if(Fj_chi>0) h->Fill(log(Fj_chi));


  }//end event loop

  gStyle->SetOptStat("nemrou");
  h->SetLineColor(kBlue);
  h->SetLineWidth(2);
  h->GetXaxis()->SetTitle("log(#chi)");
  h->GetXaxis()->SetTitleOffset(1.3);
  h->GetYaxis()->SetTitle("Entries");
  h->GetYaxis()->SetTitleOffset(1.2);
  h->SetTitle("");

  return h;
}

void Plot_logChi(){
  TH1F *h = (TH1F*) makeHisto("HZgen_HBBmodel_R12_r015_ptCut15_50000events_minHardPt400_useBtag_dR015_2btagcondition.root");
  TH1F *h2 = (TH1F*) makeHisto("TTgen_HBBmodel_R12_r015_ptCut15_50000events_minHardPt400_useBtag_dR015_2btagcondition.root");

  TCanvas *c = new TCanvas("SD - log chi", "SD - log chi",800, 600);
  //c->SetGrid();

  Double_t norm1 = h->GetEntries();
  h->Scale(1/norm1);
  Double_t norm2 = h2->GetEntries();
  h2->Scale(1/norm2);  

  //gStyle->SetOptStat("nemrou");

  c->cd();

  h->SetTitle("HZ");
  h->SetLineColor(kBlue);
  h->Draw();

  h2->SetTitle("ttbar");
  h2->SetLineColor(kRed);
  h2->Draw("SAMES");

  leg = new TLegend(0.15,0.7,.35,0.85);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->AddEntry(h,"HZ","L");
  leg->AddEntry(h2,"ttbar","L");

  leg->Draw("SAME");

  gPad->Update();  

  TPaveStats *tps1 = (TPaveStats*) h->FindObject("stats");
  tps1->SetName("HZ");
  tps1->SetTextColor(kBlue);
  double X1 = tps1->GetX1NDC();
  double Y1 = tps1->GetY1NDC();
  double X2 = tps1->GetX2NDC();
  double Y2 = tps1->GetY2NDC();
  
  TPaveStats *tps2 = (TPaveStats*) h2->FindObject("stats");
  tps2->SetName("tt");
  tps2->SetTextColor(kRed);
  tps2->SetX1NDC(X1);
  tps2->SetX2NDC(X2);
  tps2->SetY1NDC(Y1-(Y2-Y1));
  tps2->SetY2NDC(Y1);

  X1 = tps2->GetX1NDC();
  Y1 = tps2->GetY1NDC();
  X2 = tps2->GetX2NDC();
  Y2 = tps2->GetY2NDC();

  gPad->Update();

  c->SaveAs("SD_Original_HZvsTTbar__HBBmodel_R12_r15_minPtHard400_mjptCut15_50000events_useBtag_dR015_2btagcondition_Plot_logChi_All.eps");
  c->SetLogy();
  c->SaveAs("SD_Original_HZvsTTbar__HBBmodel_R12_r15_minPtHard400_mjptCut15_50000events_useBtag_dR015_2btagcondition_Plot_logChi_All_logyscale.eps");

}


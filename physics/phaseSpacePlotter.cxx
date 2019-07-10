#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <math.h> 
#include "Riostream.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TRandom.h"
#include "TMath.h"
#include <vector>
#include <TAxis.h>
#include <TLorentzRotation.h>
#include<vector>
#include <algorithm>
#include <functional>
#include <TCutG.h>
#include <TStyle.h>

using namespace std;

int phaseSpacePlotter(const char* input, int run, const char* field_config, const char* dataType){

  TFile *fIn = new TFile(input,"");

  double p_min = 1.5;
  double Ebeam = 2.221;

  if(Ebeam > 10) p_min = 1.5;
  if(Ebeam < 10) p_min = 1.0;
  if(Ebeam < 3)  p_min = 0.5;


  TF1 *qcut_clas6 = new TF1("qcut","2*0.938*x*(5.887 - (2*0.938*5.954*x)/(2*0.938*x + 4*5.887*sin((47.0/2.0)*(3.14159/180.0))*sin((47.0/2.0)*(3.14159/180.0))))",0.072,0.624);
  TF1 *wcut_clas6 = new TF1("wcut","(0.938*0.938 - 2*2)*( x/(x-1) )",0.18,0.625);

  TF1 *qcut_clas12 = new TF1("qcut","2*0.938*x*(10.5 - (2*0.938*10.5*x)/(2*0.938*x + 4*10.5*sin((47.0/2.0)*(3.14159/180.0))*sin((47.0/2.0)*(3.14159/180.0))))",0.072,0.824);
  TF1 *wcut_clas12 = new TF1("wcut","(0.938*0.938 - 2*2)*( x/(x-1) )",0.18,0.825);
  

  //TCanvas *c1 = new TCanvas("c1","c1",900,900);
  //c1->cd(0);
  //TH2F *h_q2x = (TH2F*)fIn->Get("kinematics/hist_Q2_x");
  

  TCanvas *c1a = new TCanvas("c1a","c1a",900,900);
  TH2F *h_q2x_gen = (TH2F*)fIn->Get("generated_kin/h_q2xb_gen");
  h_q2x_gen->GetXaxis()->SetTitle("Xb");
  h_q2x_gen->GetYaxis()->SetTitle("Q^2 (GeV^2)");
  h_q2x_gen->GetXaxis()->CenterTitle();
  h_q2x_gen->GetYaxis()->CenterTitle();
  qcut_clas6->SetLineColor(kMagenta);
  wcut_clas6->SetLineColor(kMagenta);
  qcut_clas12->SetLineColor(kRed);
  wcut_clas12->SetLineColor(kRed);
  h_q2x_gen->Draw("colz");
  qcut_clas6->Draw("same");
  wcut_clas6->Draw("same");
  qcut_clas12->Draw("same");
  wcut_clas12->Draw("same");
  c1a->SaveAs(Form("h_q2x_gen_%s_%s.pdf",field_config, dataType));


  TCanvas *c2a = new TCanvas("c2a","c2a",900,900);
  TH2F *h_q2t_gen = (TH2F*)fIn->Get("generated_kin/h_q2t_gen");
  h_q2t_gen->GetXaxis()->SetTitle("t (GeV^2)");
  h_q2t_gen->GetYaxis()->SetTitle("Q^2 (GeV^2)");
  h_q2t_gen->GetXaxis()->CenterTitle();
  h_q2t_gen->GetYaxis()->CenterTitle();
  h_q2t_gen->Draw("colz");
  c2a->SaveAs(Form("h_q2t_gen_%s_%s.pdf",field_config, dataType));


  TCanvas *c3a = new TCanvas("c3a","c3a",900,900);
  TH2F *h_q2w_gen = (TH2F*)fIn->Get("generated_kin/h_q2w_gen");
  h_q2w_gen->GetXaxis()->SetTitle("W (GeV^2)");
  h_q2w_gen->GetYaxis()->SetTitle("Q^2 (GeV^2)");
  h_q2w_gen->GetXaxis()->CenterTitle();
  h_q2w_gen->GetYaxis()->CenterTitle();
  h_q2w_gen->Draw("colz");
  c3a->SaveAs(Form("h_q2w_gen_%s_%s.pdf",field_config, dataType));


  TCanvas *c4a = new TCanvas("c4a","c4a",900,900);
  TH2F *h_tminq2_gen = (TH2F*)fIn->Get("kinematics/hist_tmin_vs_q2");
  h_tminq2_gen->GetXaxis()->SetTitle("t_{min} (GeV^2)");
  h_tminq2_gen->GetYaxis()->SetTitle("Q^2 (GeV^2)");
  h_tminq2_gen->GetXaxis()->CenterTitle();
  h_tminq2_gen->GetYaxis()->CenterTitle();
  h_tminq2_gen->Draw("colz");
  c4a->SaveAs(Form("h_tminq2_gen_%s_%s.pdf",field_config, dataType));

  TCanvas *c5a = new TCanvas("c5a","c5a",900,900);
  TH2F *h_tminxb_gen = (TH2F*)fIn->Get("kinematics/hist_tmin_vs_xb");
  h_tminxb_gen->GetXaxis()->SetTitle("xb");
  h_tminxb_gen->GetYaxis()->SetTitle("Q^2 (GeV^2)");
  h_tminxb_gen->GetXaxis()->CenterTitle();
  h_tminxb_gen->GetYaxis()->CenterTitle();
  h_tminxb_gen->Draw("colz");
  c5a->SaveAs(Form("h_tminxb_gen_%s_%s.pdf",field_config, dataType));

  TCanvas *c6a = new TCanvas("c6a","c6a",900,900);
  c6a->Divide(2,5);
  for( int bb = 0; bb < 10; bb++ ){
    c6a->cd(bb+1);
    TH2F *h_tminq2_xb_bin = (TH2F*)fIn->Get(Form("kinematics/hist_tmin_vs_q2_xb_bin%d",bb));
    h_tminq2_xb_bin->GetXaxis()->SetTitle("Q^{2} (GeV^{2})");
    h_tminq2_xb_bin->GetYaxis()->SetTitle("t_{min} (GeV^{2})");
    double xb_bin_size = 0.1*(bb+1);
    h_tminq2_xb_bin->SetTitle(Form("t_{min} vs Q^{2} xb %f", xb_bin_size) );
    h_tminq2_xb_bin->Draw("colz");
  }
  c6a->SaveAs(Form("h_tminq2_xb_fix_%s_%s.pdf",field_config, dataType));
   

  


  return 0;
}

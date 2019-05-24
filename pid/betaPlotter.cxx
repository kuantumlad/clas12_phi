#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <math.h> 
#include "Riostream.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TFrame.h"
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

int betaPlotter(const char* input, int run){

  TFile *fIn = new TFile(input,"");

  int windowX = 1350;

  TCanvas *c_beta_all = new TCanvas("c_beta_all","c_betal_all",windowX,900);

  TF1 *f_prot = new TF1("f_prot","x/sqrt(x*x + 0.938*0.938)",0.5,6);
  f_prot->SetLineStyle(0);
  f_prot->SetLineWidth(2);
  f_prot->SetLineColor(kRed);

  TF1 *f_kaon = new TF1("f_kaon","x/sqrt(x*x + 0.493*0.493)",0.5,6);
  f_kaon->SetLineStyle(0);
  f_kaon->SetLineWidth(2);
  f_kaon->SetLineColor(kRed);

  TF1 *f_pion = new TF1("f_pion","x/sqrt(x*x + 0.139*0.139)",0.5,6);
  f_pion->SetLineStyle(0);
  f_pion->SetLineWidth(2);
  f_pion->SetLineColor(kRed);

  gStyle->SetPalette(kRainBow);
  /*  TH2F *h2_betap_all = (TH2F*)fIn->Get("FD_PID_hadron_beta_plots/beta_vs_momentum_cut_01");
  gPad->SetLogz();
  gStyle->SetOptStat(0);
  h2_betap_all->SetTitle("#beta vs p All Positives");
  h2_betap_all->GetXaxis()->SetTitle("p [GeV]");
  h2_betap_all->GetXaxis()->CenterTitle();
  h2_betap_all->GetYaxis()->CenterTitle();
  h2_betap_all->Draw("colz");
  f_prot->Draw("same");
  f_kaon->Draw("same");
  f_pion->Draw("same");
  c_beta_all->SaveAs("beta_all_pos.pdf");

  */
  TCanvas *c_beta_p_mle = new TCanvas("c_beta_p_mle","c_beta_p_mle",windowX,900);
  c_beta_p_mle->Divide(3,2);
  
  for( int s = 1; s <=6; s++ ){
    c_beta_p_mle->cd(s);
    TH1D *hframe = new TH1D(Form("hframe%d",s),"",100,-10.0,10.0);
    hframe->SetTitleFont(63,"XYZ");
    hframe->SetLabelFont(63,"XYZ");
    hframe->SetTitleSize(35,"XYZ");
    hframe->SetLabelSize(33,"XYZ");

    TH2F *h2_betap_mle1 = (TH2F*)fIn->Get(Form("FD_PID_hadron_beta_plots/beta_vs_momentum_sec%d_cut_00",s));    
    TH2F *h2_betap_mle2 = (TH2F*)fIn->Get(Form("FD_PID_hadron_beta_plots/beta_vs_momentum_sec%d_cut_20",s));    
    TH2F *h2_betap_mle3 = (TH2F*)fIn->Get(Form("FD_PID_hadron_beta_plots/beta_vs_momentum_sec%d_cut_40",s));    
    h2_betap_mle1->SetTitle(Form("#bf{#beta} vs p Sector %d",s));
    h2_betap_mle1->GetXaxis()->SetTitle("p [GeV]");
    h2_betap_mle1->GetXaxis()->CenterTitle();

    h2_betap_mle1->GetYaxis()->SetTitle("#bf{#beta}");
    h2_betap_mle1->GetYaxis()->CenterTitle();

    gStyle->SetOptStat(0);

    if( s == 5 ){
      TCanvas *c5 = new TCanvas("c5","c5",900,900);
      c5->cd(1);
      h2_betap_mle1->Draw("colz");
      h2_betap_mle2->Draw("col same");
      h2_betap_mle3->Draw("same");
      c5->SaveAs("positives_mle_sector5.pdf");
    }

    
    h2_betap_mle1->Draw("colz");
    h2_betap_mle2->Draw("col same");
    h2_betap_mle3->Draw("same");

      

  }

  c_beta_p_mle->SaveAs(Form("positives_mle_sd_%d.pdf",run));

  TCanvas *c_betap_all_neg = new TCanvas("c_betap_all_neg","c_betap_all_neg",windowX,900);
  c_betap_all_neg->Divide(3,2);

  for( int s = 1; s <=6; s++ ){
    c_betap_all_neg->cd(s);
    TH2F *h2_betap_mle1 = (TH2F*)fIn->Get(Form("FD_PID_hadron_beta_plots/beta_vs_momentum_sec%d_cut_30",s));    
    TH2F *h2_betap_mle2 = (TH2F*)fIn->Get(Form("FD_PID_hadron_beta_plots/beta_vs_momentum_sec%d_cut_50",s));    
    h2_betap_mle1->SetTitle(Form("#bf{#beta} vs p Sector %d",s));
    h2_betap_mle1->GetXaxis()->SetTitle("p [GeV]");
    h2_betap_mle1->GetXaxis()->CenterTitle();

    h2_betap_mle1->GetYaxis()->SetTitle("#bf{#beta}");
    h2_betap_mle1->GetYaxis()->CenterTitle();

    gStyle->SetOptStat(0);

    if( s == 5 ){
      TCanvas *c5 = new TCanvas("c5","c5",900,900);
      c5->cd(1);
      h2_betap_mle1->Draw("colz");
      h2_betap_mle2->Draw("same");
      c5->SaveAs("negatives_mle_sector5.pdf");
    }


    h2_betap_mle1->Draw("colz");
    h2_betap_mle2->Draw("same");
  }

  c_betap_all_neg->SaveAs(Form("negatives_mle_sd_%d.pdf",run));


  TCanvas *c_beta_p_eb = new TCanvas("c_beta_p_eb","c_beta_p_eb",windowX,900);
  c_beta_p_eb->Divide(3,2);
  
  for( int s = 1; s <=6; s++ ){
    c_beta_p_eb->cd(s);
    TH2F *h2_betap_eb = (TH2F*)fIn->Get(Form("FD_PID_hadron_beta_plots/beta_vs_momentum_sec%d_cut_00",s));    
    h2_betap_eb->SetTitle(Form("PROTON #beta vs p Sector %d",s));
    h2_betap_eb->GetXaxis()->SetTitle("p [GeV]");
    h2_betap_eb->GetXaxis()->CenterTitle();

    h2_betap_eb->GetYaxis()->SetTitle("#beta");
    h2_betap_eb->GetYaxis()->CenterTitle();

    gStyle->SetOptStat(0);
    h2_betap_eb->Draw("colz");
    f_prot->Draw("same");
  }

  c_beta_p_eb->SaveAs(Form("proton_eb_sd_%d.pdf",run));

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //
  //
  //
  //    PIONS PLUS
  //
  //
  //
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  TCanvas *c_betap_pip = new TCanvas("c_betap_pip","c_betap_pip",windowX,900);
  c_betap_pip->Divide(3,2);
  
  for( int s = 1; s <=6; s++ ){
    c_betap_pip->cd(s);
    TH2F *h2_temp = (TH2F*)fIn->Get(Form("FD_PID_hadron_beta_plots/beta_vs_momentum_sec%d_cut_29",s));    
    h2_temp->SetTitle(Form("PION PLUS #beta vs p Sector %d",s));
    h2_temp->GetXaxis()->SetTitle("p [GeV]");
    h2_temp->GetXaxis()->CenterTitle();

    h2_temp->GetYaxis()->SetTitle("#beta");
    h2_temp->GetYaxis()->CenterTitle();

    gStyle->SetOptStat(0);
    h2_temp->Draw("colz");
  }

  c_betap_pip->SaveAs(Form("pionP_sd_%d.pdf",run));

  TCanvas *c_betap_pip2 = new TCanvas("c_betap_pip2","c_betap_pip2",windowX,900);
  c_betap_pip2->Divide(3,2);
  
  for( int s = 1; s <=6; s++ ){
    c_betap_pip2->cd(s);
    TH2F *h2_temp = (TH2F*)fIn->Get(Form("FD_PID_hadron_beta_plots/beta_vs_momentum_sec%d_cut_20",s));    
    h2_temp->SetTitle(Form("PION PLUS #beta vs p Sector %d",s));
    h2_temp->GetXaxis()->SetTitle("p [GeV]");
    h2_temp->GetXaxis()->CenterTitle();

    h2_temp->GetYaxis()->SetTitle("#beta");
    h2_temp->GetYaxis()->CenterTitle();

    gStyle->SetOptStat(0);
    h2_temp->Draw("colz");
    f_pion->Draw("same");
  }

  c_betap_pip2->SaveAs(Form("pionP_sd_%d.pdf",run));

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //
  //
  //
  //    PIONS MINUS
  //
  //
  //
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  TCanvas *c_betap_pim = new TCanvas("c_betap_pim","c_betap_pim",windowX,900);
  c_betap_pim->Divide(3,2);
  
  for( int s = 1; s <=6; s++ ){
    c_betap_pim->cd(s);
    TH2F *h2_temp = (TH2F*)fIn->Get(Form("FD_PID_hadron_beta_plots/beta_vs_momentum_sec%d_cut_39",s));    
    h2_temp->SetTitle(Form("PION MINUS #beta vs p Sector %d",s));
    h2_temp->GetXaxis()->SetTitle("p /GeV");
    h2_temp->GetXaxis()->CenterTitle();

    h2_temp->GetYaxis()->SetTitle("#beta");
    h2_temp->GetYaxis()->CenterTitle();

    gStyle->SetOptStat(0);
    h2_temp->Draw("colz");
    f_pion->Draw("same"); 
  }

  c_betap_pim->SaveAs(Form("pionP_sd_%d.pdf",run));

  TCanvas *c_betap_pim2 = new TCanvas("c_betap_pim2","c_betap_pim2",windowX,900);
  c_betap_pim2->Divide(3,2);
  
  for( int s = 1; s <=6; s++ ){
    c_betap_pim2->cd(s);
    TH2F *h2_temp = (TH2F*)fIn->Get(Form("FD_PID_hadron_beta_plots/beta_vs_momentum_sec%d_cut_30",s));    
    h2_temp->SetTitle(Form("PION MINUS #beta vs p Sector %d",s));
    h2_temp->GetXaxis()->SetTitle("p /GeV");
    h2_temp->GetXaxis()->CenterTitle();

    h2_temp->GetYaxis()->SetTitle("#beta");
    h2_temp->GetYaxis()->CenterTitle();

    gStyle->SetOptStat(0);
    h2_temp->Draw("colz");
    f_pion->Draw("colz");
  }

  c_betap_pim2->SaveAs(Form("pionM_sd_%d_eb.pdf",run));



  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //
  //
  //
  //    KAON PLUS
  //
  //
  //
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  TCanvas *c_betap_kp = new TCanvas("c_betap_kp","c_betap_kp",windowX,900);
  c_betap_kp->Divide(3,2);
  
  for( int s = 1; s <=6; s++ ){
    c_betap_kp->cd(s);
    TH2F *h2_temp = (TH2F*)fIn->Get(Form("FD_PID_hadron_beta_plots/beta_vs_momentum_sec%d_cut_49",s));    
    h2_temp->SetTitle(Form("KAON PLUS #beta vs p Sector %d",s));
    h2_temp->GetXaxis()->SetTitle("p /GeV");
    h2_temp->GetXaxis()->CenterTitle();

    h2_temp->GetYaxis()->SetTitle("#beta");
    h2_temp->GetYaxis()->CenterTitle();

    gStyle->SetOptStat(0);
    h2_temp->Draw("colz");
    f_kaon->Draw("same"); 
  }

  c_betap_kp->SaveAs(Form("kaonP_sd_%d.pdf",run));

  TCanvas *c_betap_kp2 = new TCanvas("c_betap_kp2","c_betap_kp2",windowX,900);
  c_betap_kp2->Divide(3,2);
  
  for( int s = 1; s <=6; s++ ){
    c_betap_kp2->cd(s);
    TH2F *h2_temp = (TH2F*)fIn->Get(Form("FD_PID_hadron_beta_plots/beta_vs_momentum_sec%d_cut_40",s));    
    h2_temp->SetTitle(Form("KAON PLUS EB #beta vs p Sector %d",s));
    h2_temp->GetXaxis()->SetTitle("p /GeV");
    h2_temp->GetXaxis()->CenterTitle();

    h2_temp->GetYaxis()->SetTitle("#beta");
    h2_temp->GetYaxis()->CenterTitle();

    gStyle->SetOptStat(0);
    h2_temp->Draw("colz");
    f_kaon->Draw("same");
  }

  c_betap_kp2->SaveAs(Form("kaonP_sd_%d_eb.pdf",run));

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //
  //
  //
  //    KAON MINUS
  //
  //
  //
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  TCanvas *c_betap_km = new TCanvas("c_betap_km","c_betap_km",windowX,900);
  c_betap_km->Divide(3,2);
  
  for( int s = 1; s <=6; s++ ){
    c_betap_km->cd(s);
    TH2F *h2_temp = (TH2F*)fIn->Get(Form("FD_PID_hadron_beta_plots/beta_vs_momentum_sec%d_cut_59",s));    
    h2_temp->SetTitle(Form("KAON MINUS #beta vs p Sector %d",s));
    h2_temp->GetXaxis()->SetTitle("p /GeV");
    h2_temp->GetXaxis()->CenterTitle();

    h2_temp->GetYaxis()->SetTitle("#beta");
    h2_temp->GetYaxis()->CenterTitle();

    gStyle->SetOptStat(0);
    h2_temp->Draw("colz");
    f_kaon->Draw("same"); 
  }

  c_betap_km->SaveAs(Form("kaonM_sd_%d.pdf",run));

  TCanvas *c_betap_km2 = new TCanvas("c_betap_km2","c_betap_km2",windowX,900);
  c_betap_km2->Divide(3,2);
  
  for( int s = 1; s <=6; s++ ){
    c_betap_km2->cd(s);
    TH2F *h2_temp = (TH2F*)fIn->Get(Form("FD_PID_hadron_beta_plots/beta_vs_momentum_sec%d_cut_50",s));    
    h2_temp->SetTitle(Form("KAON MINUS EB #beta vs p Sector %d",s));
    h2_temp->GetXaxis()->SetTitle("p /GeV");
    h2_temp->GetXaxis()->CenterTitle();

    h2_temp->GetYaxis()->SetTitle("#beta");
    h2_temp->GetYaxis()->CenterTitle();

    gStyle->SetOptStat(0);
    h2_temp->Draw("colz");
    f_kaon->Draw("same"); 
  }

  c_betap_km2->SaveAs(Form("kaonM_sd_%d_eb.pdf",run));
  //*/


  //left side are all positive charges and negatives are on right
  TCanvas *c_betap_all_charges = new TCanvas("c_betap_all_charges","c_betap_all_charges",900,900);
  TH2F *h2_betap_all_pos = (TH2F*)fIn->Get("FD_PID_hadron_beta_plots/beta_vs_momentum_cut_01");
  gPad->SetLogz();
  gStyle->SetOptStat(0);
  h2_betap_all_pos->SetTitle("#beta vs p All Positives");
  h2_betap_all_pos->GetXaxis()->SetTitle("p [GeV]");
  h2_betap_all_pos->GetXaxis()->CenterTitle();
  h2_betap_all_pos->GetYaxis()->CenterTitle();
  h2_betap_all_pos->Draw("colz");
  f_prot->Draw("same");
  f_kaon->Draw("same");
  f_pion->Draw("same");
  c_betap_all_charges->Update();
  c_betap_all_charges->GetFrame()->SetFillColor(kGray);
  c_betap_all_charges->SaveAs(Form("h2_betap_all_charges_%d.pdf",run));

  TCanvas *c_betap_all_charges_neg = new TCanvas("c_betap_all_charges_neg","c_betap_all_charges_neg",900,900);
  TH2F *h2_betap_all_neg = (TH2F*)fIn->Get("FD_PID_hadron_beta_plots/beta_vs_momentum_cut_31");
  gPad->SetLogz();
  gStyle->SetOptStat(0);
  h2_betap_all_neg->SetTitle("#beta vs p All Negatives");
  h2_betap_all_neg->GetXaxis()->SetTitle("p [GeV]");
  h2_betap_all_neg->GetXaxis()->CenterTitle();
  h2_betap_all_neg->GetYaxis()->CenterTitle();
  h2_betap_all_neg->Draw("colz");
  f_kaon->Draw("same");
  f_pion->Draw("same");
  c_betap_all_charges_neg->Update();
  c_betap_all_charges_neg->GetFrame()->SetFillColor(kGray);  
  c_betap_all_charges_neg->SaveAs(Form("h2_betap_all_charges_%d.pdf",run));


 

  return 0;
}

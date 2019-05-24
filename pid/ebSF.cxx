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

int ebSF(const char* input, int run){

  TFile *fIn = new TFile(input,"");

  TCanvas *c_el_sf = new TCanvas("c_el_sf", "c_el_sf", 1800, 900);
  c_el_sf->Divide(3,2);

  for( int s = 1; s <= 6; s++ ){
    c_el_sf->cd(s);
    
    TCanvas *c_temp_sf = new TCanvas(Form("c_temp_sf_s%d",s),Form("c_temp_sf_s%d",s),900,600);


    TH2F *h_el_sf_all = (TH2F*)fIn->Get(Form("FD_PID_electron_EC_plots/EC_total_sampling_fraction_sec%d_cut_00",s));
    TH2F *h_el_sf = (TH2F*)fIn->Get(Form("FD_PID_electron_EC_plots/EC_total_sampling_fraction_sec%d_cut_01",s));
    
    h_el_sf->SetTitle(Form("Sampling Fraction vs p Sector %d",s));
    h_el_sf->GetXaxis()->SetTitle("p [GeV]");
    h_el_sf->GetXaxis()->CenterTitle();
    h_el_sf->GetYaxis()->SetTitle("SF");
    h_el_sf->GetYaxis()->CenterTitle();
    gStyle->SetOptStat(0);

    h_el_sf->Draw("colz");
    c_temp_sf->SaveAs(Form("h_el_sf_sector_%d_run%d.pdf",s,run));   
  }

  return 0;
}

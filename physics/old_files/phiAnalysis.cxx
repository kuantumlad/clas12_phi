#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <math.h> 
#include "Riostream.h"
#include "TFile.h"
#include "TLatex.h"
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


Float_t m_Kp = 0.4937;
Float_t m_Km = 0.4937;
double eBeam = 10.594;

int phiAnalysis(const char* input, const char* output, const char* inmodel){


  std::string model(inmodel);

  TFile *inFile = new TFile(input,"");

  TTree *tree1 = (TTree*)inFile->Get("out_tree_epPhi");

  TFile *out = new TFile(output, "RECREATE");

  double W_out;
  double Q2_out;
  double x_out;
  double y_out;
  double nu_out;
  double t1_out;
  double cmphi1_out;
  double cmcostheta1_out;
  double pt1_out;
  double eta1_out;
  double z1_out;
  double fcup;
  double E_ele;
  double px_ele;
  double py_ele;
  double pz_ele;
  double E_prot;
  double px_prot;
  double py_prot;
  double pz_prot;
  double E_kaonP;
  double px_kaonP;
  double py_kaonP;
  double pz_kaonP;
  double E_kaonM;
  double px_kaonM;
  double py_kaonM;
  double pz_kaonM;
  double perp_mntm;
  double missing_e;
  double missing_mm2;
  double epX_mm2;
  double epX_mm;
  int channel;
  
  tree1->SetBranchAddress("W",&W_out);
  tree1->SetBranchAddress("Q2", &Q2_out);
  tree1->SetBranchAddress("x", &x_out);
  tree1->SetBranchAddress("y", &y_out);
  tree1->SetBranchAddress("nu", &nu_out);
  tree1->SetBranchAddress("minus_t", &t1_out);
  tree1->SetBranchAddress("cmphi", &cmphi1_out);
  tree1->SetBranchAddress("cmcostheta", &cmcostheta1_out);
  tree1->SetBranchAddress("pt", &pt1_out);
  tree1->SetBranchAddress("eta", &eta1_out);
  tree1->SetBranchAddress("z", &z1_out);
  tree1->SetBranchAddress("fcup", &fcup);
  //ELECTRON
  tree1->SetBranchAddress("E_ele", &E_ele);
  tree1->SetBranchAddress("px_ele", &px_ele);
  tree1->SetBranchAddress("py_ele", &py_ele);
  tree1->SetBranchAddress("pz_ele", &pz_ele);
  //PROTON
  tree1->SetBranchAddress("E_prot", &E_prot);
  tree1->SetBranchAddress("px_prot", &px_prot);
  tree1->SetBranchAddress("py_prot", &py_prot);
  tree1->SetBranchAddress("pz_prot", &pz_prot);
  //KAON PLUS
  tree1->SetBranchAddress("E_kaonP", &E_kaonP);
  tree1->SetBranchAddress("px_kaonP",&px_kaonP);
  tree1->SetBranchAddress("py_kaonP",&py_kaonP);
  tree1->SetBranchAddress("pz_kaonP",&pz_kaonP);
  //KAON MINUS
  tree1->SetBranchAddress("E_kaonM", &E_kaonM);
  tree1->SetBranchAddress("px_kaonM",&px_kaonM);
  tree1->SetBranchAddress("py_kaonM",&py_kaonM);
  tree1->SetBranchAddress("pz_kaonM",&pz_kaonM);
  tree1->SetBranchAddress("channel",&channel);

  tree1->SetBranchAddress("perp_mntm",&perp_mntm); 
  tree1->SetBranchAddress("missing_e",&missing_e);
  tree1->SetBranchAddress("missing_mm2",&missing_mm2);
  tree1->SetBranchAddress("epX_mm2",&epX_mm2);
  tree1->SetBranchAddress("epX_mm",&epX_mm);
  

  TH1D *h_phi_mass = new TH1D("h_phi_mass","#phi Mass from K^{+}K^{-}",40,0.7,1.8);
  TH1D *h_t = new TH1D("h_t","-t",200,0.0,5.5);
  TH1D *h_mass_epkX = new TH1D("h_mass_epkX","Mass of Kaon Minus from #phi Event",200,0.0, 0.7);
  

  TH1F *hist_epX_mass;
  TH1F *hist_epKpKmX_mass;
  TH1F *hist_epX_mm2;

  TH1F *hist_epKpKmX_mm2;
  TH1F *hist_epKpKmX_me;
  TH1F *hist_epKpKmX_mperpmntm;

  TH1F *hist_epKpKmX_mm2_cute;
  TH1F *hist_epKpKmX_me_cute;
  TH1F *hist_epKpKmX_mperpmntm_cute;

  TH1F *hist_epKpKmX_mm2_cutmm;
  TH1F *hist_epKpKmX_me_cutmm;
  TH1F *hist_epKpKmX_mperpmntm_cutmm;

  TH1F *hist_phi_mass;
  TH1F *hist_phi_mass_noexcl;
  TH1F *hist_phi_mass_cuts;
  TH1F *hist_phi_mass_pangle_cut;

  TH1F *hist_ephiX;
  TH1F *hist_ephiX_mpx;
  TH1F *hist_ephiX_mpy;
  TH1F *hist_ephiX_mpz;
  TH1F *hist_ephiX_me;
  
  TH1F *hist_q2;
  TH1F *hist_w;
  TH1F *hist_xb;
  TH1F *hist_t;
  TH1F *hist_cmcos;
  TH1F *hist_cmphi;
 
  TH2F *hist_q2x;
  TH2F *hist_q2w;
  TH2F *hist_q2t;
  TH2F *hist_wxb;
  TH2F *hist_wt;

  TH1F *hist_q2_final;
  TH1F *hist_w_final;
  TH1F *hist_xb_final;
  TH1F *hist_t_final;
  TH1F *hist_cmcos_final;
  TH1F *hist_cmphi_final;
 
  TH2F *hist_q2x_final;
  TH2F *hist_q2w_final;
  TH2F *hist_q2t_final;
  TH2F *hist_wxb_final;
  TH2F *hist_wt_final;
  TH2F *hist_xt_final;

  TH2F *hist_pr_betap;
  TH2F *hist_kp_betap;
  TH2F *hist_km_betap;

  TH2F *hist_pr_betap_final;
  TH2F *hist_kp_betap_final;
  TH2F *hist_km_betap_final;

  TH2F *hist_el_thetap;
  TH2F *hist_el_phip;
  TH2F *hist_el_thetaphi;

  TH2F *hist_pr_thetap;
  TH2F *hist_pr_phip;
  TH2F *hist_pr_thetaphi;

  TH2F *hist_kp_thetap;
  TH2F *hist_kp_phip;
  TH2F *hist_kp_thetaphi;
  
  TH2F *hist_km_thetap;
  TH2F *hist_km_phip;
  TH2F *hist_km_thetaphi;


  TH1F *h_el_p;
  TH1F *h_pr_p;
  TH1F *h_kp_p;
  TH1F *h_km_p;

  h_el_p = new TH1F("h_el_p","h_el_p",100, 0.0, 10.5);
  h_pr_p = new TH1F("h_pr_p","h_pr_p",100, 0.0, 6.5);
  h_kp_p = new TH1F("h_kp_p","h_kp_p",100, 0.0, 6.5);
  h_km_p = new TH1F("h_km_p","h_km_p",100, 0.0, 6.5);
  

  hist_el_thetap = new TH2F("hist_el_thetap","hist_el_thetap",100,0.0,10.5, 100,0.0, 60.0);
  hist_el_phip = new TH2F("hist_el_phip","hist_el_phip",300, -180.0, 180.0, 300, 0.0, 11.0);
  hist_el_thetaphi = new TH2F("hist_el_thetaphi","hist_el_thetaphi",300, 0.0, 60.0, 300, -180.0, 180.0);

  hist_pr_thetap = new TH2F("hist_pr_thetap","hist_pr_thetap",100,0.0, 6.5, 100,0.0, 60.0);
  hist_pr_phip = new TH2F("hist_pr_phip","hist_pr_phip",300, -180, 180.0, 300,0.0, 6.5);
  hist_pr_thetaphi = new TH2F("hist_pr_thetaphi","hist_pr_thetaphi",300, 0.0, 60.0, 300, -180.0, 180.0);

  hist_kp_thetap = new TH2F("hist_kp_thetap","hist_kp_thetap",100,0.0, 6.5, 100,0.0, 60.0);
  hist_kp_phip = new TH2F("hist_kp_phip","hist_kp_phip",300, -180.0, 180.0, 300, 0.0, 6.50);
  hist_kp_thetaphi = new TH2F("hist_kp_thetaphi","hist_kp_thetaphi",300, 0.0, 60.0, 300, -180.0, 180.0);

  hist_km_thetap = new TH2F("hist_km_thetap","hist_km_thetap",100,0.0, 6.5, 100,0.0, 60.0);
  hist_km_phip = new TH2F("hist_km_phip","hist_km_phip",300, -180.0, 180.0, 300, 0.0, 6.50);
  hist_km_thetaphi = new TH2F("hist_km_thetaphi","hist_km_thetaphi",300, 0.0, 60.0, 300, -180.0, 180.0);
  
  TH2F *hist_pr_phi_pangle;
  TH2F *hist_pr_phi_pangle_final;

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // Binning histograms 
  int xb_bins=2;
  TH1F *hist_xb_bin = new TH1F("hist_xb_bin","hist_xb_bin",xb_bins, 0.0, 0.75);
  vector<TH1F*> hist_t_xbins;
  for( int b = 0; b < xb_bins; b++ ){
    hist_t_xbins.push_back( new TH1F(Form("hist_t_xbin%d",b),Form("hist_t_xbin%d",b), 40, 0.0, 5.0) );
  }


  hist_epX_mass = new TH1F("hist_epX_mass","hist_epX_mass",100, 0.0, 3.0);
  hist_epX_mass->GetXaxis()->SetTitle("Missing Mass of epX [GeV]");
  hist_epX_mass->GetXaxis()->CenterTitle();

  hist_epX_mm2 = new TH1F("hist_epX_mm2","hist_epX_mm2", 100, -8.0, 8.);
  hist_epX_mm2->GetXaxis()->SetTitle("Missing Mass^2 of epX [GeV^2]");
  hist_epX_mm2->GetXaxis()->CenterTitle();

  hist_epKpKmX_mass = new TH1F("hist_epKpKmX_mass","hist_epKpKmX_mass",100, 0.0, 8.0);
  hist_epKpKmX_mass->GetXaxis()->SetTitle("Missing Mass of epKpKmX [GeV]");
  hist_epKpKmX_mass->GetXaxis()->CenterTitle();
  
  hist_epKpKmX_mm2 = new TH1F("hist_epKpKmX_mm2","hist_epKpKmX_mm2", 75, -0.1, 0.1);//-0.204561*5.0, 0.204561*5.0 );
  hist_epKpKmX_mm2->GetXaxis()->SetTitle("MM^2 [GeV^2]");
  hist_epKpKmX_mm2->GetXaxis()->CenterTitle();

  hist_epKpKmX_me = new TH1F("hist_epKpKmX_me","hist_epKpKmX_me", 100, -2.0, 6.0);//-0.35*5.0, 0.35*5.0 );
  hist_epKpKmX_me->GetXaxis()->SetTitle("Missing Energy [GeV]");
  hist_epKpKmX_me->GetXaxis()->CenterTitle();

  hist_epKpKmX_mperpmntm = new TH1F("hist_epKpKmX_mperpmntm","hist_epKpKmX_mperpmntm",100, -3.0, 8.0);
  hist_epKpKmX_mperpmntm->GetXaxis()->SetTitle("Missing Perp. Mntm [GeV]");
  hist_epKpKmX_mperpmntm->GetXaxis()->CenterTitle();

  //cut on missing energy
  hist_epKpKmX_mm2_cute = new TH1F("hist_epKpKmX_mm2_cute","hist_epKpKmX_mm2_cute", 75, -0.1, 0.1);//-8.0, 8.0 );
  hist_epKpKmX_mm2_cute->SetTitle("Missing Mass^2 with Cut on Missing E");
  hist_epKpKmX_mm2_cute->GetXaxis()->SetTitle("MM^2 [GeV^2]");
  hist_epKpKmX_mm2_cute->GetXaxis()->CenterTitle();
  
  hist_epKpKmX_mperpmntm_cute = new TH1F("hist_epKpKmX_mperpmntm_cute","hist_epKpKmX_mperpmntm_cute",100, -3.0, 8.0);
  hist_epKpKmX_mperpmntm_cute->SetTitle("Missing Perp Mntm with Cut on Missing E");
  hist_epKpKmX_mperpmntm_cute->GetXaxis()->SetTitle("Missing Perp. Mntm [GeV]");
  hist_epKpKmX_mperpmntm_cute->GetXaxis()->CenterTitle();

  //cut on missing mass
  hist_epKpKmX_me_cutmm = new TH1F("hist_epKpKmX_me_cutmm","hist_epKpKmX_me_cutmm", 100, -2.0, 6.0);//-8.0, 8.0 );
  hist_epKpKmX_me_cutmm->SetTitle("Missing Energy with Cut on MM^2");
  hist_epKpKmX_me_cutmm->GetXaxis()->SetTitle("Missing Energy [GeV]");
  hist_epKpKmX_me_cutmm->GetXaxis()->CenterTitle();
  
  hist_epKpKmX_mperpmntm_cutmm = new TH1F("hist_epKpKmX_mperpmntm_cutmm","hist_epKpKmX_mperpmntm_cutmm",100, -3.0, 8.0);
  hist_epKpKmX_mperpmntm_cutmm->SetTitle("Missing Perp Mntm with Cut on MM^2");
  hist_epKpKmX_mperpmntm_cutmm->GetXaxis()->SetTitle("Missing Perp. Mntm [GeV]");
  hist_epKpKmX_mperpmntm_cutmm->GetXaxis()->CenterTitle();

  // phi mass plots
  double phi_mass_max = 1.15;
  double phi_mass_min = 0.95;
  if ( model == "sim" ) phi_mass_max = 1.2;
  std::cout << " >> phi mass max " << phi_mass_max << std::endl;
  hist_phi_mass = new TH1F("hist_phi_mass","Invariant Mass K^+ K^-",200, phi_mass_min, phi_mass_max);
  hist_phi_mass->GetXaxis()->SetTitle("Mass [GeV]");
  hist_phi_mass->GetXaxis()->CenterTitle();

  hist_phi_mass_noexcl = new TH1F("hist_phi_mass_noexcl","Invariant Mass K^+ K^- No Excl Cuts",200, phi_mass_min, 2.15);
  hist_phi_mass_noexcl->GetXaxis()->SetTitle("Mass [GeV]");
  hist_phi_mass_noexcl->GetXaxis()->CenterTitle();


  hist_phi_mass_cuts = new TH1F("hist_phi_mass_cuts","Invariant Mass K^+ K^- Final", 200, 0.8, 2.2);//phi_mass_max);
  hist_phi_mass_cuts->GetXaxis()->SetTitle("Mass [GeV]");
  hist_phi_mass_cuts->GetXaxis()->CenterTitle();

  hist_phi_mass_pangle_cut = new TH1F("hist_phi_mass_pangle_cut","Invariant Mass K^+ K^- Pangle Cut",60, 0.90, 2.5);
  hist_phi_mass_pangle_cut->GetXaxis()->SetTitle("Mass [GeV]");
  hist_phi_mass_pangle_cut->GetXaxis()->CenterTitle();


  //only look at when measuring phi's in final state.
  hist_ephiX = new TH1F("hist_ephiX","hist_ephiX",200, -10.0, 10.0);
  hist_ephiX_mpx = new TH1F("hist_ephiX_mpx","hist_ephiX_mpx",200, -3.0, 3.0);
  hist_ephiX_mpy = new TH1F("hist_ephiX_mpy","hist_ephiX_mpy",200, -3.0, 3.0);
  hist_ephiX_mpz = new TH1F("hist_ephiX_mpz","hist_ephiX_mpz",200, -3.0, 3.0);
  hist_ephiX_me = new TH1F("hist_ephiX_me","hist_ephiX_me",200, -1.0, 6.0);


  // phase space plots
  hist_q2 = new TH1F("hist_q2","hist_q2",200,0.0,10.0);
  hist_q2->SetTitle("Q^2");
  hist_q2->GetXaxis()->SetTitle("Q^2 [GeV^2]");
  hist_w = new TH1F("hist_w","hist_w",200,0.0, 10.0);
  hist_w->SetTitle("W");
  hist_w->GetXaxis()->SetTitle("W [GeV]");

  hist_xb = new TH1F("hist_xb","hist_xb",200, 0.0, 1.1);
  hist_xb->SetTitle("Xb");
  hist_xb->GetXaxis()->SetTitle("Xb");

  hist_t = new TH1F("hist_t","hist_t",200, 0.0, 10.0);
  hist_t->SetTitle("t");
  hist_t->GetXaxis()->SetTitle("t [GeV]");

  hist_q2x = new TH2F("hist_q2x","hist_q2x",100, 0.0, 1.1, 100, 0.0, 10.5);
  hist_q2x->SetTitle("Q^2 vx Xb");
  hist_q2x->GetXaxis()->SetTitle("Xb");
  hist_q2x->GetXaxis()->CenterTitle();
  hist_q2x->GetYaxis()->SetTitle("Q^2 [GeV^2]");
  hist_q2x->GetYaxis()->CenterTitle();

  hist_q2w = new TH2F("hist_q2w","hist_q2w",100, 0.0, 5.5, 100, 0.0, 10.5);
  hist_q2w->SetTitle("Q^2 vx W");
  hist_q2w->GetXaxis()->SetTitle("W [GeV]");
  hist_q2w->GetXaxis()->CenterTitle();
  hist_q2w->GetYaxis()->SetTitle("Q^2 [GeV^2]");
  hist_q2w->GetYaxis()->CenterTitle();

  hist_q2t = new TH2F("hist_q2t","hist_q2t",100, 0.0, 5.5, 100, 0.0, 10.5);
  hist_q2t->SetTitle("Q^2 vx -t");
  hist_q2t->GetXaxis()->SetTitle("-t [GeV]");
  hist_q2t->GetXaxis()->CenterTitle();
  hist_q2t->GetYaxis()->SetTitle("Q^2 [GeV^2]");
  hist_q2t->GetYaxis()->CenterTitle();

  hist_wxb = new TH2F("hist_wxb","hist_wxb",100, 0.0, 5.5, 100, 0.0, 1.1);
  hist_wxb->SetTitle("W vs Xb");
  hist_wxb->GetXaxis()->SetTitle("W [GeV]");
  hist_wxb->GetXaxis()->CenterTitle();
  hist_wxb->GetYaxis()->SetTitle("Xb");
  hist_wxb->GetYaxis()->CenterTitle();

  hist_wt = new TH2F("hist_wt","hist_wt",100, 0.0, 5.5, 100, 0.0, 5.5);
  hist_wt->SetTitle("W vs -t");
  hist_wt->GetXaxis()->SetTitle("-t [GeV]");
  hist_wt->GetXaxis()->CenterTitle();
  hist_wt->GetYaxis()->SetTitle("W [GeV]");
  hist_wt->GetYaxis()->CenterTitle();

  hist_cmcos = new TH1F("hist_cmcos","hist_cmcos", 50, -1.0, 1.0);
  hist_cmcos->GetXaxis()->SetTitle("Cos(#theta_{cm})");
  hist_cmcos->GetXaxis()->CenterTitle();
  hist_cmphi = new TH1F("hist_cmphi","hist_cmphi", 80, -180.0, 180.0);
  hist_cmphi->GetXaxis()->SetTitle("#phi_{cm} [deg]");
  hist_cmphi->GetXaxis()->CenterTitle();

  hist_q2_final = new TH1F("hist_q2_final","hist_q2_final",200,0.0,10.0);
  hist_q2_final->SetTitle("Q^2 Final");
  hist_q2_final->GetXaxis()->SetTitle("Q^2 [GeV^2]");
  hist_w_final = new TH1F("hist_w_final","hist_w_final",200,0.0, 10.0);
  hist_w_final->SetTitle("W");
  hist_w_final->GetXaxis()->SetTitle("W [GeV]");

  hist_xb_final = new TH1F("hist_xb_final","hist_xb_final",200, 0.0, 1.1);
  hist_xb_final->SetTitle("Xb Final");
  hist_xb_final->GetXaxis()->SetTitle("Xb");

  hist_t_final = new TH1F("hist_t_final","hist_t_final",200, 0.0, 10.0);
  hist_t_final->SetTitle("-t Final");
  hist_t_final->GetXaxis()->SetTitle("t [GeV]");

  hist_q2x_final = new TH2F("hist_q2x_final","hist_q2x_final",100, 0.0, 1.1, 100, 0.0, 10.5);
  hist_q2x_final->SetTitle("Q^2 vx Xb Final");
  hist_q2x_final->GetXaxis()->SetTitle("Xb");
  hist_q2x_final->GetXaxis()->CenterTitle();
  hist_q2x_final->GetYaxis()->SetTitle("Q^2 [GeV^2]");
  hist_q2x_final->GetYaxis()->CenterTitle();

  hist_q2w_final = new TH2F("hist_q2w_final","hist_q2w_final",100, 0.0, 5.5, 100, 0.0, 10.5);
  hist_q2w_final->SetTitle("Q^2 vx W Final");
  hist_q2w_final->GetXaxis()->SetTitle("W [GeV]");
  hist_q2w_final->GetXaxis()->CenterTitle();
  hist_q2w_final->GetYaxis()->SetTitle("Q^2 [GeV^2]");
  hist_q2w_final->GetYaxis()->CenterTitle();

  hist_q2t_final = new TH2F("hist_q2t_final","hist_q2t_final",100, 0.0, 5.5, 100, 0.0, 10.5);
  hist_q2t_final->SetTitle("Q^2 vx -t Final");
  hist_q2t_final->GetXaxis()->SetTitle("-t [GeV]");
  hist_q2t_final->GetXaxis()->CenterTitle();
  hist_q2t_final->GetYaxis()->SetTitle("Q^2 [GeV^2]");
  hist_q2t_final->GetYaxis()->CenterTitle();

  hist_wxb_final = new TH2F("hist_wxb_final","hist_wxb_final",100, 0.0, 5.5, 100, 0.0, 1.1);
  hist_wxb_final->SetTitle("W vs Xb Final");
  hist_wxb_final->GetXaxis()->SetTitle("W [GeV]");
  hist_wxb_final->GetXaxis()->CenterTitle();
  hist_wxb_final->GetYaxis()->SetTitle("Xb");
  hist_wxb_final->GetYaxis()->CenterTitle();

  hist_xt_final = new TH2F("hist_xt_final","hist_xt_final",100, 0.0, 5.5, 100, 0.0, 1.1);
  hist_xt_final->SetTitle("-t vs Xb Final");
  hist_xt_final->GetXaxis()->SetTitle("-t[GeV]");
  hist_xt_final->GetXaxis()->CenterTitle();
  hist_xt_final->GetYaxis()->SetTitle("Xb");
  hist_xt_final->GetYaxis()->CenterTitle();

  hist_wt_final = new TH2F("hist_wt_final","hist_wt_final",100, 0.0, 5.5, 100, 0.0, 5.5);
  hist_wt_final->SetTitle("W vs -t Final");
  hist_wt_final->GetXaxis()->SetTitle("-t [GeV]");
  hist_wt_final->GetXaxis()->CenterTitle();
  hist_wt_final->GetYaxis()->SetTitle("W [GeV]");
  hist_wt_final->GetYaxis()->CenterTitle();

  hist_cmcos_final = new TH1F("hist_cmcos_final","hist_cmcos_final", 50, -1.0, 1.0);
  hist_cmcos_final->GetXaxis()->SetTitle("Cos(#theta_{cm}) Final");
  hist_cmcos_final->GetXaxis()->CenterTitle();
  hist_cmphi_final = new TH1F("hist_cmphi_final","hist_cmphi_final", 80, -180.0, 180.0);
  hist_cmphi_final->GetXaxis()->SetTitle("#phi_{cm} [deg] Final");
  hist_cmphi_final->GetXaxis()->CenterTitle();

  hist_pr_betap = new TH2F("hist_pr_betap","hist_pr_betap",100,0.0,3.5,100, 0.0, 1.05);
  hist_kp_betap = new TH2F("hist_kp_betap","hist_kp_betap",100,0.0,3.5,100, 0.0, 1.05);
  hist_km_betap = new TH2F("hist_km_betap","hist_km_betap",100,0.0,3.5,100, 0.0, 1.05);

  hist_pr_betap_final = new TH2F("hist_pr_betap_final","hist_pr_betap_final",100,0.0,3.5,100, 0.0, 1.05);
  hist_kp_betap_final = new TH2F("hist_kp_betap_final","hist_kp_betap_final",100,0.0,3.5,100, 0.0, 1.05);
  hist_km_betap_final = new TH2F("hist_km_betap_final","hist_km_betap_final",100,0.0,3.5,100, 0.0, 1.05);

  hist_pr_phi_pangle = new TH2F("hist_pr_phi_pangle","hist_pr_phi_pangle", 75, 0.0, 4.0, 75, 0.0, 180.0);
  hist_pr_phi_pangle_final = new TH2F("hist_pr_phi_pangle_final","hist_pr_phi_pangle_final", 75, 0.0, 4.0, 75, 0.0, 180.0);

  TH1D *hist_missing_var = new TH1D("hist_missing_var","hist_missing_var",100,-1.0,10.0);
  std::vector<TH1D*> v_phi_mass;
  for( int c = 0; c < 7; c++ ){
    v_phi_mass.push_back( new TH1D(Form("h_phi_mass_final_selection_cut_%d",c), Form("h_phi_mass_final_selection_cut_%d",c), 50, 0.8, 1.6) );
  }
			  

  int n_chanel_1 = 0;
  int n_chanel_2 = 0;

  TLorentzVector lv_beam;
  lv_beam.SetPxPyPzE(0,0, eBeam,  eBeam);
  TLorentzVector lv_target;
  lv_target.SetPxPyPzE(0,0,0,0.938);

  double lowMM2ValueFit[1];
  double highMM2ValueFit[1];

  double lowMEValueFit[1];
  double highMEValueFit[1];
  int sector0, sector1;
  double mean0, mean1, sig0, sig1;

  ifstream readFromMM2Cut("mm2_cut_limits_"+model+".txt");
  if( readFromMM2Cut.is_open() ){
    while(readFromMM2Cut >> sector0 ) {//std::getline (readFromWCut, line) ){
	    
      readFromMM2Cut >> mean0 >> sig0;
      std::cout << " >> MM2 CUT PARAMETERS: " << sector0 << " " << mean0 << " " << sig0 << std::endl;
	    
      lowMM2ValueFit[sector0]= mean0 - 3.5*sig0;
      highMM2ValueFit[sector0]= mean0 + 3.5*sig0;
    }
  }

  ifstream readFromMECut("me_cut_limits_"+model+".txt");
  if( readFromMECut.is_open() ){
    while(readFromMECut >> sector1 ) {//std::getline (readFromWCut, line) ){
	    
      readFromMECut >> mean1 >> sig1;
      std::cout << " >> ME CUT PARAMETERS: " << sector1 << " " << mean1 << " " << sig1 << std::endl;
	    
      lowMEValueFit[sector0]= mean1 - 3.5*sig1;
      highMEValueFit[sector0]= mean1 + 3.5*sig1;
    }
  }


  double bot_mm2 = 0;
  double top_mm2 = 0;
  double bot_me = 0.0;
  double top_me = 0.0;
  

  if( model == "data" ){
    bot_mm2 = -0.118;// mean_mm2 + 3.5*sig_mm2;
    top_mm2 = 0.07172;//mean_mm2 - 3.5*sig_mm2;
    
    bot_me = -0.508;//mean_me + 3.5*sig_me;
    top_me = 1.102;//mean_me - 3.5*sig_me;
  }
  else if( model == "sim" ){
    bot_mm2 = -0.006145 - 3*0.02334;
    top_mm2 = -0.006145 + 3*0.02334;
    
    bot_me = 0.0996 - 3*0.8972;
    top_me = 0.0996 + 3*0.8972;
  }
  else if( model == "SIMTESTV5" ){
    top_me = highMEValueFit[0];
    bot_me =  lowMEValueFit[0];
    
    top_mm2 = highMM2ValueFit[0];
    bot_mm2 =  lowMM2ValueFit[0];
  }

  
  


  for( int ev = 0; ev < tree1->GetEntries(); ev++ ){
    tree1->GetEntry(ev);
    
    if( channel == 2 ){
      //std::cout << " >> " << px_ele << std::endl;
      TLorentzVector lv_el(px_ele,py_ele,pz_ele,E_ele);
      TLorentzVector lv_pr(px_prot,py_prot,pz_prot,E_prot);

      double energy_kp = sqrt(px_kaonP*px_kaonP + py_kaonP*py_kaonP + pz_kaonP*pz_kaonP + m_Kp*m_Kp);
      double energy_km = sqrt(px_kaonM*px_kaonP + py_kaonM*py_kaonM + pz_kaonM*pz_kaonM + m_Km*m_Km);

      TLorentzVector lv_kp;
      lv_kp.SetPxPyPzE( px_kaonP, py_kaonP, pz_kaonP, E_kaonP);
      TLorentzVector lv_km;
      lv_km.SetPxPyPzE( px_kaonM, py_kaonM, pz_kaonM, E_kaonM);
      TLorentzVector lv_phi = lv_kp + lv_km;
      
      TLorentzVector lv_epX = lv_beam + lv_target - lv_el - lv_pr;
      TLorentzVector lv_epKpKmX = lv_beam + lv_target - lv_el - lv_pr - lv_kp - lv_km;
      //std::cout<< " >> " << lv_kp.P() << std::endl;

      TLorentzVector lv_ephiX = lv_beam + lv_target - lv_el - lv_kp - lv_km;

      TLorentzVector lv_epkpX = lv_beam + lv_target - lv_el - lv_pr - lv_kp;
      TLorentzVector lv_epkmX = lv_beam + lv_target - lv_el - lv_pr - lv_km;
      TLorentzVector lv_ekpkmX = lv_beam + lv_target - lv_el - lv_km - lv_kp;

      
      //cout << " >> KAON MASS " << lv_kp.M() << endl;    
      //cout<< " >> MASS OF PHI MESON IS " << lv_phi.M() << endl;
      hist_phi_mass_noexcl->Fill(lv_phi.M());
      hist_epX_mass->Fill( lv_epX.M() );

      t1_out = 2.0*0.938*(lv_pr.E() - 0.938);
      //if ( lv_el.Theta() * 180.0/TMath::Pi() < 12.0 ) continue;
      ///std::cout << " >> " << t1_out << std::endl;
      if( Q2_out > 1.0 && W_out > 2.0 ){

	hist_t->Fill(t1_out);	
	hist_q2->Fill(Q2_out);
	hist_w->Fill(W_out);
	hist_xb->Fill(x_out);
	hist_cmcos->Fill(cmcostheta1_out);
	hist_cmphi->Fill(cmphi1_out);
		
	hist_q2x->Fill(x_out, Q2_out);
	hist_q2w->Fill(W_out, Q2_out);
	hist_q2t->Fill(t1_out, Q2_out);
	hist_wxb->Fill(W_out, x_out);
	hist_wt->Fill(t1_out, W_out);

	h_phi_mass->Fill(lv_phi.M());
	hist_epX_mm2->Fill( lv_epX.M2() );

	hist_epKpKmX_mm2->Fill( lv_epKpKmX.M2() );
	hist_epKpKmX_me->Fill(  lv_epKpKmX.E() );
	hist_epKpKmX_mperpmntm->Fill(  lv_epKpKmX.Vect().Perp() );
	
	hist_pr_phi_pangle->Fill( lv_pr.P(), (lv_phi.Vect()).Angle(lv_pr.Vect()) * 180.0/TMath::Pi() );
	

	//mm2 from quickfit -0.0126548 0.00713083
	//double mean_mm2 = -0.0126548;//-0.0663764;
	//double sig_mm2 =  0.00713083; //0.204561; 


	//from fitting ME after mm2 cut using ROOT FIT PANEL
	//double mean_me =  0.147; // 0.240121;
	//double sig_me = 0.251; //0.359356;
	/*
	double bot_mm2 = 0;
	double top_mm2 = 0;
	double bot_me = 0.0;
	double top_me = 0.0;
	*/
	//bot_mm2 = -0.118;// mean_mm2 + 3.5*sig_mm2;
	//top_mm2 = 0.07172;//mean_mm2 - 3.5*sig_mm2;

	//bot_me = -0.508;//mean_me + 3.5*sig_me;
	//top_me = 1.102;//mean_me - 3.5*sig_me;

	
	
	//	if( lv_phi.M() < 1.051 && lv_phi.M() > 0.95 ){
	  
	  hist_ephiX->Fill(lv_ephiX.M());
	  hist_ephiX_mpx->Fill(lv_ephiX.Px());
	  hist_ephiX_mpy->Fill(lv_ephiX.Py());
	  hist_ephiX_mpz->Fill(lv_ephiX.Pz());
	  hist_ephiX_me->Fill(lv_ephiX.E());
	  //}
	
	if( lv_epKpKmX.M2() < top_mm2 && lv_epKpKmX.M2() > bot_mm2  ){
	  hist_epKpKmX_me_cutmm->Fill(  lv_epKpKmX.E() );
	  hist_epKpKmX_mperpmntm_cutmm->Fill(  lv_epKpKmX.Vect().Perp() );
	}
	if( lv_epKpKmX.E() < top_me && lv_epKpKmX.E() > bot_me  ){
	  hist_epKpKmX_mm2_cute->Fill(  lv_epKpKmX.M2() );
	  hist_epKpKmX_mperpmntm_cute->Fill(  lv_epKpKmX.Vect().Perp() );
	}

	//Check if we get more of a signal 
	if( lv_pr.P() > 0.70  && (lv_phi.Vect()).Angle(lv_pr.Vect()) * 180.0/TMath::Pi() < 60.0 ){
	  hist_phi_mass_pangle_cut->Fill(lv_phi.M());	
	}

	double top_ephiX_mass = 0.983 + 3.0*0.223;
	double bot_ephiX_mass = 0.983 - 3.0*0.223;

	//if( lv_epKpKmX.M2() < top_mm2 && lv_epKpKmX.M2() > bot_mm2  ){
	// if( lv_epKpKmX.E() < top_me && lv_epKpKmX.E() > bot_me  ){
	double missingvar  = lv_ephiX.E()*lv_ephiX.E() - lv_ephiX.Vect().Mag()*lv_ephiX.Vect().Mag();
	hist_missing_var->Fill( missingvar );
	
	//if( true ){
	//if( lv_ephiX.M() < top_ephiX_mass && lv_ephiX.M() > bot_ephiX_mass){
	//if ( missingvar < 1.2 && missingvar > 0.61 ){
	double kp_theta = lv_kp.Theta()*180.0/TMath::Pi();
	double km_theta = lv_km.Theta()*180.0/TMath::Pi();

	bool miss_e_cut = lv_epKpKmX.E() < 1.0;//2.5;
	bool miss_tran_p = lv_epKpKmX.Perp() < 0.5;// 0.500;
	bool miss_m2_kp = (lv_kp.E()*lv_kp.E() - lv_kp.P()*lv_kp.P()) < 0.5;//2.25;
	bool miss_m2_km = (lv_km.E()*lv_km.E() - lv_km.P()*lv_km.P()) < 0.5;//2.25;
	bool miss_m2_epkpkmX = lv_epKpKmX.M2() < 0.3;//0.6;
	bool proton_collinearity = lv_pr.Vect().Angle(lv_el.Vect()) * 180.0/TMath::Pi() < 30;//40;
	bool kaon_collinearity = lv_kp.Vect().Angle(lv_km.Vect()) * 180.0/TMath::Pi() < 20;// 30;


	//if( (lv_kp.P() > 1.0 && lv_kp.P() < 2.5 ) && ( lv_km.P() > 1.0 && lv_km.P() < 2.5  ) ){
	double radians_to_deg = 180.0/3.141592658;
	bool low_epkpkm_mm2 = std::pow(lv_epKpKmX.M(),2) < 0.6;

	double colin_pr_ang = TMath::ACos( lv_pr.Vect().Dot(lv_ekpkmX.Vect() ) / (lv_pr.Vect().Mag() * lv_ekpkmX.Vect().Mag() ) ) * radians_to_deg;
	double colin_km_ang = TMath::ACos( lv_km.Vect().Dot(lv_epkpX.Vect() ) / (lv_km.Vect().Mag() * lv_epkpX.Vect().Mag() ) ) * radians_to_deg; 
	double colin_kp_ang = TMath::ACos( lv_kp.Vect().Dot(lv_epkmX.Vect() ) / (lv_kp.Vect().Mag() * lv_epkmX.Vect().Mag() ) ) * radians_to_deg;

	double miss_perp_mntm = sqrt(lv_epKpKmX.Px()*lv_epKpKmX.Px() + lv_epKpKmX.Py()*lv_epKpKmX.Py());

	bool colin_pr = colin_pr_ang < 30;
	bool colin_kp = colin_kp_ang < 20;
	bool colin_km = colin_km_ang < 20;
	bool miss_perp = miss_perp_mntm < 0.5;
	

	if( lv_el.P() > 1.5  && lv_kp.Theta()*radians_to_deg < 35 && lv_km.Theta()*radians_to_deg < 40  && lv_kp.P() < 2.5 && lv_km.P() < 2.5 && lv_pr.P() > 0.5 && lv_pr.P() < 4.0 ) {
	  if( true && miss_tran_p && low_epkpkm_mm2 && colin_pr && colin_km && miss_perp && colin_kp){
	    std::cout << " phi event " << " el momentum " << lv_el.P() << std::endl;
	    std::cout << "           " << " pr momentum " << lv_pr.P() << std::endl;
	    std::cout << "           " << " kp momentum " << lv_kp.P() << std::endl;
	    std::cout << "           " << " km momentum " << lv_km.P() << std::endl;
	    std::cout << "           " << " phi mass    " << lv_phi.M() << std::endl;
	    if( miss_e_cut ){
	      v_phi_mass[0]->Fill(lv_phi.M());
	      if( miss_tran_p ){
		v_phi_mass[1]->Fill(lv_phi.M());
		if( miss_m2_kp ){
		  v_phi_mass[2]->Fill(lv_phi.M());
		  if( miss_m2_km ){
		    v_phi_mass[3]->Fill(lv_phi.M());
		    if( miss_m2_epkpkmX ){
		      v_phi_mass[4]->Fill(lv_phi.M());
		      if( proton_collinearity ){
			v_phi_mass[5]->Fill(lv_phi.M());
			if( kaon_collinearity ){
			  v_phi_mass[6]->Fill(lv_phi.M());
			}
		      }
		    }
		  }
		}
	      }
	    }
	
	
	    //	  if( Q2_out > 1 && W_out > 2) {
	      hist_phi_mass_cuts->Fill(lv_phi.M());
	      
	    hist_t_final->Fill(t1_out);	
	    hist_q2_final->Fill(Q2_out);
	    hist_w_final->Fill(W_out);
	    hist_xb_final->Fill(x_out);
	    hist_cmcos_final->Fill(cmcostheta1_out);
	    hist_cmphi_final->Fill(cmphi1_out);
	      
	    hist_q2x_final->Fill(x_out, Q2_out);
	    hist_q2w_final->Fill(W_out, Q2_out);
	    hist_q2t_final->Fill(t1_out, Q2_out);
	    hist_wxb_final->Fill(W_out, x_out);
	    hist_wt_final->Fill(t1_out, W_out);
	    hist_xt_final->Fill(t1_out, x_out);
	      
	    hist_pr_phi_pangle_final->Fill( lv_pr.P(), (lv_phi.Vect()).Angle(lv_pr.Vect()) * 180.0/TMath::Pi() );
	      
	    //	    int event_x_bin = hist_xb_bin->FindBin(x_out);
	    //	    //std::cout << " EVENT X BIN " << event_x_bin << std::endl;
	    //	    if( event_x_bin > 0 ){
	    //      hist_t_xbins[event_x_bin-1]->Fill(t1_out);	      
	    //  }
	      

	      
	    hist_el_thetap->Fill(lv_el.P(), lv_el.Theta()*180.0/TMath::Pi() );
	    hist_el_phip->Fill(lv_el.Phi()*180.0/TMath::Pi(), lv_el.P() );
	    hist_el_thetaphi->Fill( lv_el.Theta()*180.0/TMath::Pi(), lv_el.Phi()*180.0/TMath::Pi() );

	    hist_pr_thetap->Fill(lv_pr.P(), lv_pr.Theta()*180.0/TMath::Pi() );
	    hist_pr_phip->Fill(lv_pr.Phi()*180.0/TMath::Pi(), lv_pr.P() );
	    hist_pr_thetaphi->Fill( lv_pr.Theta()*180.0/TMath::Pi(), lv_pr.Phi()*180.0/TMath::Pi() );

	    hist_kp_thetap->Fill(lv_kp.P(), lv_kp.Theta()*180.0/TMath::Pi() );
	    hist_kp_phip->Fill(lv_kp.Phi()*180.0/TMath::Pi(), lv_kp.P() );
	    hist_kp_thetaphi->Fill( lv_kp.Theta()*180.0/TMath::Pi(), lv_kp.Phi()*180.0/TMath::Pi() );
	    
	    hist_km_thetap->Fill(lv_km.P(), lv_km.Theta()*180.0/TMath::Pi() );
	    hist_km_phip->Fill(lv_km.Phi()*180.0/TMath::Pi(), lv_km.P() );
	    hist_km_thetaphi->Fill( lv_km.Theta()*180.0/TMath::Pi(), lv_km.Phi()*180.0/TMath::Pi() );
	    
	    h_el_p->Fill(lv_el.P());
	    h_pr_p->Fill(lv_pr.P());
	    h_kp_p->Fill(lv_kp.P());
	    h_km_p->Fill(lv_km.P());
	  	    
	    

	  }
	}
		

	
     		             
     }
      
      //*/
      //
      //
      //FIT PARAMETERS FOR MM2 -0.0663764 STD 0.204561
      //>> FIT PARAMETERS FOR ME 0.240121 STD 0.359356
      //>> FIT PARAMETERS FOR MP 0.295894 STD 0.179074
      //
      //
    }
  }

//std::cout << " channel 2 events " << n_chanel_2 << std::endl;
//  TCanvas *c1 = new TCanvas("c1","c1",800,800);
// h_phi_mass->Draw();
// c1->SaveAs("r3432_phi_mass.pdf");
  
  TF1 *fit_epKpKm_mm2 = new TF1("fit_epKpKm_mm2","gaus",-0.35, 0.35);
  TF1 *fit_epKpKm_me = new TF1("fit_epKpKm_me","gaus",-0.1,0.1);

  //hist_epKpKmX_mm2_cute->Fit("fit_epKpKm_mm2","R+") ;//hist_epKpKmX_mm2->Fit("fit_epKpKm_mm2","R+");
  //hist_epKpKmX_me_cutmm->Fit("fit_epKpKm_me","R+");

  double mean_mm2 = fit_epKpKm_mm2->GetParameter(1);
  double sig_mm2 = fit_epKpKm_mm2->GetParameter(2);

  double mean_me = fit_epKpKm_me->GetParameter(1);
  double sig_me = fit_epKpKm_me->GetParameter(2);
    
  std::cout << " >> FIT PARAMETERS FOR MM2 " << mean_mm2 << " STD " << sig_mm2 << std::endl;
  std::cout << " >> FIT PARAMETERS FOR ME " << mean_me << " STD " << sig_me << std::endl;
  
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // fit the phi mass
  TF1 *fit_phi_mass = new TF1("fit_phi_mass","gaus",1.005, 1.025);
  TF1 *fit_phi_mass2 = new TF1("fit_phi_mass2","gaus",1.005, 1.025);
  TF1 *fit_phi_mass3 = new TF1("fit_phi_mass3","gaus",1.005, 1.025);
  //hist_phi_mass_noexcl->Fit("fit_phi_mass","R+");
  //hist_phi_mass_cuts->Fit("fit_phi_mass2","R+");
  //hist_phi_mass_pangle_cut->Fit("fit_phi_mass3","R+");
  
  double mean_phi_mass = fit_phi_mass->GetParameter(1);
  double sig_phi_mass = fit_phi_mass->GetParameter(2);

  double mean_phi_mass2 = fit_phi_mass2->GetParameter(1);
  double sig_phi_mass2 = fit_phi_mass2->GetParameter(2);

  double mean_phi_mass3 = fit_phi_mass3->GetParameter(1);
  double sig_phi_mass3 = fit_phi_mass3->GetParameter(2);
 
  
  std::cout << "PHI MASS TABLE BELOW " << std::endl;
  std::cout << ">> no excl cuts:" << mean_phi_mass << " sigma: " << sig_phi_mass << std::endl;
  std::cout << ">> all cuts:" << mean_phi_mass2 << " sigma: " << sig_phi_mass2 << std::endl;
  std::cout << ">> Q2>1 & W>2 with pangle cut:" << mean_phi_mass3 << " sigma: " << sig_phi_mass3 << std::endl;


  out->Write();
  out->Close();

  return 0;

}

///  To execute run ROOT and compile macro with .L selector_simple.cxx++
///
///  Then call function:  selector_simple("output/2391_junemap.root", "output/selected_2391_junemap.root")
///   
/// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
#include "TH1D.h"
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
using namespace std;

#include "Math/GenVector/PxPyPzE4D.h"
#include "Math/GenVector/PtEtaPhiE4D.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"

using namespace ROOT::Math;

/// //////////////////////////////////////////////////////////////////////////////////////////////////////
/// settings:

//double Ebeam = 6.42313;
//double Ebeam = 2.22193;
double Ebeam = 10.6;//594;
int process_Events = -1;            // process all events
//int process_Events = 500000;     // process given number of events


/// /////////////////////////////////////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////////////////////////////////////
/// common variables
///

Float_t Pival = 3.14159265359;

Float_t m_e = 0.000511;
Float_t m_p = 0.93827;
Float_t m_n = 0.9396;
Float_t m_pip = 0.1396;
Float_t m_pim = 0.1396;
Float_t m_Kp = 0.4937;
Float_t m_Km = 0.4937;
Float_t m_phi = 1.019;
Float_t c = 299792458;

// vectors with components of the Lorentzvector to fill into the tree

int helicity;
double fcup;
vector<double> *p4_ele_px = 0;
vector<double> *p4_ele_py = 0;
vector<double> *p4_ele_pz = 0;
vector<double> *p4_ele_E = 0;
vector<double> *p4_prot_px = 0;
vector<double> *p4_prot_py = 0;
vector<double> *p4_prot_pz = 0;
vector<double> *p4_prot_E = 0;
vector<double> *p4_neutr_px = 0;
vector<double> *p4_neutr_py = 0;
vector<double> *p4_neutr_pz = 0;
vector<double> *p4_neutr_E = 0;
vector<double> *p4_pip_px = 0;
vector<double> *p4_pip_py = 0;
vector<double> *p4_pip_pz = 0;
vector<double> *p4_pip_E = 0;
vector<double> *p4_pim_px = 0;
vector<double> *p4_pim_py = 0;
vector<double> *p4_pim_pz = 0;
vector<double> *p4_pim_E = 0;
vector<double> *p4_Kp_px = 0;
vector<double> *p4_Kp_py = 0;
vector<double> *p4_Kp_pz = 0;
vector<double> *p4_Kp_E = 0;
vector<double> *p4_Km_px = 0;
vector<double> *p4_Km_py = 0;
vector<double> *p4_Km_pz = 0;
vector<double> *p4_Km_E = 0;
vector<double> *p4_phot_px = 0;
vector<double> *p4_phot_py = 0;
vector<double> *p4_phot_pz = 0;
vector<double> *p4_phot_E = 0;
vector<double> *prot_beta_final = 0;
vector<double> *Kp_beta_final = 0;
vector<double> *Km_beta_final = 0;
vector<double> *el_sector = 0;
vector<double> *prot_sector = 0;
vector<double> *Kp_sector = 0;
vector<double> *Km_sector = 0;

/// varibales for particles:

const static int BUFFER = 20;

TLorentzVector ele[BUFFER]; 
TLorentzVector prot[BUFFER]; 
TLorentzVector neutr[BUFFER]; 
TLorentzVector pip[BUFFER]; 
TLorentzVector pim[BUFFER]; 
TLorentzVector Kp[BUFFER]; 
TLorentzVector Km[BUFFER]; 
TLorentzVector phot[BUFFER];


//added for beta vs p
double prot_beta[BUFFER];
double kaonP_beta[BUFFER];
double kaonM_beta[BUFFER];

int sect_el[BUFFER];
int sect_pr[BUFFER];
int sect_kp[BUFFER];
int sect_km[BUFFER];

///  counting variables:

double ele_count;
double prot_count;
double neutr_count;
double pip_count;
double pim_count;
double Kp_count;
double Km_count;
double phot_count;


// kinematic variables
int hel;
double W, Q2, nu, x, y;
double t1, cmphi1, cmcostheta1, pt1, eta1, z1;
double M_e_p_X_miss, M_e_p_X_miss2;

double W_out, Q2_out, nu_out, x_out, y_out;
double t1_out, cmphi1_out, cmcostheta1_out, pt1_out, eta1_out, z1_out;
double perp_mntm, missing_e, missing_mm2, epX_mm2, epX_mm;
 

double W_2, Q2_2, nu_2, x_2, y_2;
double t1_2, cmphi1_2, cmcostheta1_2, pt1_2, eta1_2, z1_2;


///////////////////////////////////////////////////////////////////////////////////////////////////////
///  selected particles:

Int_t evcat; 
Double_t E_ele; Double_t px_ele; Double_t py_ele; Double_t pz_ele;
Double_t E_prot_1; Double_t px_prot_1; Double_t py_prot_1; Double_t pz_prot_1;

//ADDED 
Double_t E_kaonP_1; Double_t px_kaonP_1; Double_t py_kaonP_1; Double_t pz_kaonP_1;
Double_t E_kaonM_1; Double_t px_kaonM_1; Double_t py_kaonM_1; Double_t pz_kaonM_1;

Double_t prot_beta_out; Double_t kaonP_beta_out; Double_t kaonM_beta_out;
Int_t el_sect; Int_t pr_sect; Int_t kp_sect; Int_t km_sect;

Int_t channel;

/// ////////////////////////////////////////////////////////////////////////////////////////////////////
///  index of selected particle:

int select_ele;
int select_prot_1; 

//NEW
int select_kaonP_1;
int select_kaonM_1;

/// ////////////////////////////////////////////////////////////////////////////////////////////////////
///  input and output files:

TFile *out;

char name[200];
char title[200];


int selector_phi( const char* inFile, const char* outputfile ){

  Char_t *inTree="out_tree";
  cout << "Initalize the input tree ... " << endl;

  /// ///////////////////////////////////////////////////////////////////////////
  ///  get input file:
  /// ///////////////////////////////////////////////////////////////////////////

   TFile *f; 
  Char_t tmpstr[80];
  Double_t fraction;

  f = new TFile(inFile,"");   // Input File

  if(f->IsZombie()){   // Check if TFile exists!
    cout<<"Input file " << inFile << "doesn't exist!" << endl;
    cout<<"Exit program" << endl;
    return 0;
  }

  cout << "Reading from File: " << inFile << endl;

  TTree *anaTree=(TTree *) f->Get(inTree);

  if(anaTree==0){  // Check if TTree exists!
    cout << "Tree " << inTree << " doesn't exist!!!" << endl;
    cout <<"Exit program" << endl;
    return 0;
  }
    
    
  /// /////////////////////////////////////////////////////////////////////////////    
  ///  get branches from input file:
  /// ///////////////////////////////////////////////////////////////////////////// 
    
  anaTree->SetBranchAddress("helicity", &helicity);
  anaTree->SetBranchAddress("fcup", &fcup);
  anaTree->SetBranchAddress("p4_ele_px", &p4_ele_px);
  anaTree->SetBranchAddress("p4_ele_py", &p4_ele_py);
  anaTree->SetBranchAddress("p4_ele_pz", &p4_ele_pz);
  anaTree->SetBranchAddress("p4_ele_E", &p4_ele_E);
  anaTree->SetBranchAddress("p4_prot_px", &p4_prot_px);
  anaTree->SetBranchAddress("p4_prot_py", &p4_prot_py);
  anaTree->SetBranchAddress("p4_prot_pz", &p4_prot_pz);
  anaTree->SetBranchAddress("p4_prot_E", &p4_prot_E);
  anaTree->SetBranchAddress("p4_neutr_px", &p4_neutr_px);
  anaTree->SetBranchAddress("p4_neutr_py", &p4_neutr_py);
  anaTree->SetBranchAddress("p4_neutr_pz", &p4_neutr_pz);
  anaTree->SetBranchAddress("p4_neutr_E", &p4_neutr_E);
  anaTree->SetBranchAddress("p4_pip_px", &p4_pip_px);
  anaTree->SetBranchAddress("p4_pip_py", &p4_pip_py);
  anaTree->SetBranchAddress("p4_pip_pz", &p4_pip_pz);
  anaTree->SetBranchAddress("p4_pip_E", &p4_pip_E);
  anaTree->SetBranchAddress("p4_pim_px", &p4_pim_px);
  anaTree->SetBranchAddress("p4_pim_py", &p4_pim_py);
  anaTree->SetBranchAddress("p4_pim_pz", &p4_pim_pz);
  anaTree->SetBranchAddress("p4_pim_E", &p4_pim_E);
  anaTree->SetBranchAddress("p4_Kp_px", &p4_Kp_px);
  anaTree->SetBranchAddress("p4_Kp_py", &p4_Kp_py);
  anaTree->SetBranchAddress("p4_Kp_pz", &p4_Kp_pz);
  anaTree->SetBranchAddress("p4_Kp_E", &p4_Kp_E);
  anaTree->SetBranchAddress("p4_Km_px", &p4_Km_px);
  anaTree->SetBranchAddress("p4_Km_py", &p4_Km_py);
  anaTree->SetBranchAddress("p4_Km_pz", &p4_Km_pz);
  anaTree->SetBranchAddress("p4_Km_E", &p4_Km_E);
  anaTree->SetBranchAddress("p4_phot_px", &p4_phot_px);
  anaTree->SetBranchAddress("p4_phot_py", &p4_phot_py);
  anaTree->SetBranchAddress("p4_phot_pz", &p4_phot_pz);
  anaTree->SetBranchAddress("p4_phot_E", &p4_phot_E);  

  out = new TFile(outputfile, "RECREATE");

  TTree out_tree1("out_tree_epPhi","out_tree_epPhi");
  out_tree1.Branch("W",&W_out);
  out_tree1.Branch("Q2", &Q2_out);
  out_tree1.Branch("x", &x_out);
  out_tree1.Branch("y", &y_out);
  out_tree1.Branch("nu", &nu_out);
  out_tree1.Branch("minus_t", &t1_out);
  out_tree1.Branch("helicity", &hel);
  out_tree1.Branch("cmphi", &cmphi1_out);
  out_tree1.Branch("cmcostheta", &cmcostheta1_out);
  out_tree1.Branch("pt", &pt1_out);
  out_tree1.Branch("eta", &eta1_out);
  out_tree1.Branch("z", &z1_out);
  out_tree1.Branch("fcup", &fcup);
  out_tree1.Branch("E_ele", &E_ele);
  out_tree1.Branch("px_ele", &px_ele);
  out_tree1.Branch("py_ele", &py_ele);
  out_tree1.Branch("pz_ele", &pz_ele);
  out_tree1.Branch("E_prot", &E_prot_1);
  out_tree1.Branch("px_prot", &px_prot_1);
  out_tree1.Branch("py_prot", &py_prot_1);
  out_tree1.Branch("pz_prot", &pz_prot_1);
  out_tree1.Branch("px_kaonP",&px_kaonP_1);
  out_tree1.Branch("py_kaonP",&py_kaonP_1);
  out_tree1.Branch("pz_kaonP",&pz_kaonP_1);
  out_tree1.Branch("E_kaonP",&E_kaonP_1);
  out_tree1.Branch("px_kaonM",&px_kaonM_1);
  out_tree1.Branch("py_kaonM",&py_kaonM_1);
  out_tree1.Branch("pz_kaonM",&pz_kaonM_1);
  out_tree1.Branch("E_kaonM",&E_kaonM_1);
  out_tree1.Branch("beta_prot",&prot_beta_out);
  out_tree1.Branch("beta_kaonP",&kaonP_beta_out);
  out_tree1.Branch("beta_kaonM",&kaonM_beta_out);
  out_tree1.Branch("perp_mntm",&perp_mntm);
  out_tree1.Branch("missing_e",&missing_e);
  out_tree1.Branch("missing_mm2",&missing_mm2);
  out_tree1.Branch("epX_mm2",&epX_mm2);
  out_tree1.Branch("epX_mm",&epX_mm);
  //out_tree1.Branch("pr_sect",&el_sect);
  //out_tree1.Branch("kp_sect",&kp_sect);
  ///out_tree1.Branch("km_sect",&km_sect);  
  out_tree1.Branch("channel",&channel);


  out->mkdir("particle_histograms_all");				
  out->cd ("particle_histograms_all");

  TH1F *hist_all_electron_p; TH1F *hist_all_electron_theta; TH1F *hist_all_electron_phi;
  TH2F *hist_all_electron_p_vs_theta; TH2F *hist_all_electron_p_vs_phi; TH2F *hist_all_electron_theta_vs_phi;

  TH1F *hist_all_proton_p; TH1F *hist_all_proton_theta; TH1F *hist_all_proton_phi;
  TH2F *hist_all_proton_p_vs_theta; TH2F *hist_all_proton_p_vs_phi; TH2F *hist_all_proton_theta_vs_phi;
  TH1F *hist_all_kp_p; TH1F *hist_all_kp_theta; TH1F *hist_all_kp_phi;
  TH2F *hist_all_kp_p_vs_theta; TH2F *hist_all_kp_p_vs_phi; TH2F *hist_all_kp_theta_vs_phi;
  TH1F *hist_all_km_p; TH1F *hist_all_km_theta; TH1F *hist_all_km_phi;
  TH2F *hist_all_km_p_vs_theta; TH2F *hist_all_km_p_vs_phi; TH2F *hist_all_km_theta_vs_phi;

  hist_all_electron_p = new TH1F("hist_all_electron_p", "electron momentum", 500,0,Ebeam+0.3);   
  hist_all_electron_p->GetXaxis()->SetTitle("p /GeV");
  hist_all_electron_p->GetYaxis()->SetTitle("counts");
  hist_all_electron_theta = new TH1F("hist_all_electron_theta", "electron #Theta", 140,0,140);   
  hist_all_electron_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_all_electron_theta->GetYaxis()->SetTitle("counts");
  hist_all_electron_phi = new TH1F("hist_all_electron_phi", "electron #phi", 180,-180,180);   
  hist_all_electron_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_all_electron_phi->GetYaxis()->SetTitle("counts");
  hist_all_electron_p_vs_theta = new TH2F("hist_all_electron_p_vs_theta", "electron p vs #Theta", 140,0,140,500,0,Ebeam+0.3);   
  hist_all_electron_p_vs_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_all_electron_p_vs_theta->GetYaxis()->SetTitle("p /GeV");
  hist_all_electron_p_vs_phi = new TH2F("hist_all_electron_p_vs_phi", "electron p vs #phi", 180,-180,180, 500,0,Ebeam+0.3);   
  hist_all_electron_p_vs_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_all_electron_p_vs_phi->GetYaxis()->SetTitle("p /GeV");
  hist_all_electron_theta_vs_phi = new TH2F("hist_all_electron_theta_vs_phi", "electron #Theta vs #phi", 180,-180,180, 140,0,140);   
  hist_all_electron_theta_vs_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_all_electron_theta_vs_phi->GetYaxis()->SetTitle("#Theta /deg");

  hist_all_proton_p = new TH1F("hist_all_proton_p", "proton momentum", 500,0,Ebeam+0.3);   
  hist_all_proton_p->GetXaxis()->SetTitle("p /GeV");
  hist_all_proton_p->GetYaxis()->SetTitle("counts");
  hist_all_proton_theta = new TH1F("hist_all_proton_theta", "proton #Theta", 140,0,140);   
  hist_all_proton_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_all_proton_theta->GetYaxis()->SetTitle("counts");
  hist_all_proton_phi = new TH1F("hist_all_proton_phi", "proton #phi", 180,-180,180);   
  hist_all_proton_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_all_proton_phi->GetYaxis()->SetTitle("counts");
  hist_all_proton_p_vs_theta = new TH2F("hist_all_proton_p_vs_theta", "proton p vs #Theta", 140,0,140,500,0,Ebeam+0.3);   
  hist_all_proton_p_vs_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_all_proton_p_vs_theta->GetYaxis()->SetTitle("p /GeV");
  hist_all_proton_p_vs_phi = new TH2F("hist_all_proton_p_vs_phi", "proton p vs #phi", 180,-180,180, 500,0,Ebeam+0.3);   
  hist_all_proton_p_vs_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_all_proton_p_vs_phi->GetYaxis()->SetTitle("p /GeV");
  hist_all_proton_theta_vs_phi = new TH2F("hist_all_proton_theta_vs_phi", "proton #Theta vs #phi", 180,-180,180, 140,0,140);   
  hist_all_proton_theta_vs_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_all_proton_theta_vs_phi->GetYaxis()->SetTitle("#Theta /deg");


  hist_all_kp_p = new TH1F("hist_all_kp_p", "Kp momentum", 500,0,Ebeam+0.3);   
  hist_all_kp_p->GetXaxis()->SetTitle("p /GeV");
  hist_all_kp_p->GetYaxis()->SetTitle("counts");
  hist_all_kp_theta = new TH1F("hist_all_kp_theta", "Kp #Theta", 140,0,140);   
  hist_all_kp_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_all_kp_theta->GetYaxis()->SetTitle("counts");
  hist_all_kp_phi = new TH1F("hist_all_kp_phi", "kp #phi", 180,-180,180);   
  hist_all_kp_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_all_kp_phi->GetYaxis()->SetTitle("counts");
  hist_all_kp_p_vs_theta = new TH2F("hist_all_kp_p_vs_theta", "Kp p vs #Theta", 140,0,140,500,0,Ebeam+0.3);   
  hist_all_kp_p_vs_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_all_kp_p_vs_theta->GetYaxis()->SetTitle("p /GeV");
  hist_all_kp_p_vs_phi = new TH2F("hist_all_kp_p_vs_phi", "Kp p vs #phi", 180,-180,180, 500,0,Ebeam+0.3);   
  hist_all_kp_p_vs_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_all_kp_p_vs_phi->GetYaxis()->SetTitle("p /GeV");
  hist_all_kp_theta_vs_phi = new TH2F("hist_all_kp_theta_vs_phi", "Kp #Theta vs #phi", 180,-180,180, 140,0,140);   
  hist_all_kp_theta_vs_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_all_kp_theta_vs_phi->GetYaxis()->SetTitle("#Theta /deg");

  hist_all_km_p = new TH1F("hist_all_km_p", "km momentum", 500,0,Ebeam+0.3);   
  hist_all_km_p->GetXaxis()->SetTitle("p /GeV");
  hist_all_km_p->GetYaxis()->SetTitle("counts");
  hist_all_km_theta = new TH1F("hist_all_km_theta", "Km #Theta", 140,0,140);   
  hist_all_km_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_all_km_theta->GetYaxis()->SetTitle("counts");
  hist_all_km_phi = new TH1F("hist_all_km_phi", "km #phi", 180,-180,180);   
  hist_all_km_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_all_km_phi->GetYaxis()->SetTitle("counts");
  hist_all_km_p_vs_theta = new TH2F("hist_all_km_p_vs_theta", "Km p vs #Theta", 140,0,140,500,0,Ebeam+0.3);   
  hist_all_km_p_vs_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_all_km_p_vs_theta->GetYaxis()->SetTitle("p /GeV");
  hist_all_km_p_vs_phi = new TH2F("hist_all_km_p_vs_phi", "Km p vs #phi", 180,-180,180, 500,0,Ebeam+0.3);   
  hist_all_km_p_vs_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_all_km_p_vs_phi->GetYaxis()->SetTitle("p /GeV");
  hist_all_km_theta_vs_phi = new TH2F("hist_all_km_theta_vs_phi", "Km #Theta vs #phi", 180,-180,180, 140,0,140);   
  hist_all_km_theta_vs_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_all_km_theta_vs_phi->GetYaxis()->SetTitle("#Theta /deg");

  out->mkdir("particle_histograms_selected");				
  out->cd ("particle_histograms_selected");

  TH1F *hist_electron_p; TH1F *hist_electron_theta; TH1F *hist_electron_phi;
  TH2F *hist_electron_p_vs_theta; TH2F *hist_electron_p_vs_phi; TH2F *hist_electron_theta_vs_phi;

  TH1F *hist_proton_p; TH1F *hist_proton_theta; TH1F *hist_proton_phi;
  TH2F *hist_proton_p_vs_theta; TH2F *hist_proton_p_vs_phi; TH2F *hist_proton_theta_vs_phi;
  TH1F *hist_kp_p; TH1F *hist_kp_theta; TH1F *hist_kp_phi;
  TH2F *hist_kp_p_vs_theta; TH2F *hist_kp_p_vs_phi; TH2F *hist_kp_theta_vs_phi;
  TH1F *hist_km_p; TH1F *hist_km_theta; TH1F *hist_km_phi;
  TH2F *hist_km_p_vs_theta; TH2F *hist_km_p_vs_phi; TH2F *hist_km_theta_vs_phi;

  hist_electron_p = new TH1F("hist_electron_p", "electron momentum", 500,0,Ebeam+0.3);   
  hist_electron_p->GetXaxis()->SetTitle("p /GeV");
  hist_electron_p->GetYaxis()->SetTitle("counts");
  hist_electron_theta = new TH1F("hist_electron_theta", "electron #Theta", 140,0,140);   
  hist_electron_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_electron_theta->GetYaxis()->SetTitle("counts");
  hist_electron_phi = new TH1F("hist_electron_phi", "electron #phi", 180,-180,180);   
  hist_electron_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_electron_phi->GetYaxis()->SetTitle("counts");
  hist_electron_p_vs_theta = new TH2F("hist_electron_p_vs_theta", "electron p vs #Theta", 140,0,140,500,0,Ebeam+0.3);   
  hist_electron_p_vs_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_electron_p_vs_theta->GetYaxis()->SetTitle("p /GeV");
  hist_electron_p_vs_phi = new TH2F("hist_electron_p_vs_phi", "electron p vs #phi", 180,-180,180, 500,0,Ebeam+0.3);   
  hist_electron_p_vs_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_electron_p_vs_phi->GetYaxis()->SetTitle("p /GeV");
  hist_electron_theta_vs_phi = new TH2F("hist_electron_theta_vs_phi", "electron #Theta vs #phi", 180,-180,180, 140,0,140);   
  hist_electron_theta_vs_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_electron_theta_vs_phi->GetYaxis()->SetTitle("#Theta /deg");

  hist_proton_p = new TH1F("hist_proton_p", "proton momentum", 500,0,Ebeam+0.3);   
  hist_proton_p->GetXaxis()->SetTitle("p /GeV");
  hist_proton_p->GetYaxis()->SetTitle("counts");
  hist_proton_theta = new TH1F("hist_proton_theta", "proton #Theta", 140,0,140);   
  hist_proton_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_proton_theta->GetYaxis()->SetTitle("counts");
  hist_proton_phi = new TH1F("hist_proton_phi", "proton #phi", 180,-180,180);   
  hist_proton_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_proton_phi->GetYaxis()->SetTitle("counts");
  hist_proton_p_vs_theta = new TH2F("hist_proton_p_vs_theta", "proton p vs #Theta", 140,0,140,500,0,Ebeam+0.3);   
  hist_proton_p_vs_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_proton_p_vs_theta->GetYaxis()->SetTitle("p /GeV");
  hist_proton_p_vs_phi = new TH2F("hist_proton_p_vs_phi", "proton p vs #phi", 180,-180,180, 500,0,Ebeam+0.3);   
  hist_proton_p_vs_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_proton_p_vs_phi->GetYaxis()->SetTitle("p /GeV");
  hist_proton_theta_vs_phi = new TH2F("hist_proton_theta_vs_phi", "proton #Theta vs #phi", 180,-180,180, 140,0,140);   
  hist_proton_theta_vs_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_proton_theta_vs_phi->GetYaxis()->SetTitle("#Theta /deg");


  hist_kp_p = new TH1F("hist_kp_p", "Kp momentum", 500,0,Ebeam+0.3);   
  hist_kp_p->GetXaxis()->SetTitle("p /GeV");
  hist_kp_p->GetYaxis()->SetTitle("counts");
  hist_kp_theta = new TH1F("hist_kp_theta", "Kp #Theta", 140,0,140);   
  hist_kp_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_kp_theta->GetYaxis()->SetTitle("counts");
  hist_kp_phi = new TH1F("hist_kp_phi", "kp #phi", 180,-180,180);   
  hist_kp_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_kp_phi->GetYaxis()->SetTitle("counts");
  hist_kp_p_vs_theta = new TH2F("hist_kp_p_vs_theta", "Kp p vs #Theta", 140,0,140,500,0,Ebeam+0.3);   
  hist_kp_p_vs_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_kp_p_vs_theta->GetYaxis()->SetTitle("p /GeV");
  hist_kp_p_vs_phi = new TH2F("hist_kp_p_vs_phi", "Kp p vs #phi", 180,-180,180, 500,0,Ebeam+0.3);   
  hist_kp_p_vs_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_kp_p_vs_phi->GetYaxis()->SetTitle("p /GeV");
  hist_kp_theta_vs_phi = new TH2F("hist_kp_theta_vs_phi", "Kp #Theta vs #phi", 180,-180,180, 140,0,140);   
  hist_kp_theta_vs_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_kp_theta_vs_phi->GetYaxis()->SetTitle("#Theta /deg");

  hist_km_p = new TH1F("hist_km_p", "km momentum", 500,0,Ebeam+0.3);   
  hist_km_p->GetXaxis()->SetTitle("p /GeV");
  hist_km_p->GetYaxis()->SetTitle("counts");
  hist_km_theta = new TH1F("hist_km_theta", "Km #Theta", 140,0,140);   
  hist_km_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_km_theta->GetYaxis()->SetTitle("counts");
  hist_km_phi = new TH1F("hist_km_phi", "km #phi", 180,-180,180);   
  hist_km_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_km_phi->GetYaxis()->SetTitle("counts");
  hist_km_p_vs_theta = new TH2F("hist_km_p_vs_theta", "Km p vs #Theta", 140,0,140,500,0,Ebeam+0.3);   
  hist_km_p_vs_theta->GetXaxis()->SetTitle("#Theta /deg");
  hist_km_p_vs_theta->GetYaxis()->SetTitle("p /GeV");
  hist_km_p_vs_phi = new TH2F("hist_km_p_vs_phi", "Km p vs #phi", 180,-180,180, 500,0,Ebeam+0.3);   
  hist_km_p_vs_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_km_p_vs_phi->GetYaxis()->SetTitle("p /GeV");
  hist_km_theta_vs_phi = new TH2F("hist_km_theta_vs_phi", "Km #Theta vs #phi", 180,-180,180, 140,0,140);   
  hist_km_theta_vs_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_km_theta_vs_phi->GetYaxis()->SetTitle("#Theta /deg");

  out->mkdir("kinematics");				
  out->cd ("kinematics");

  // Look at total phase space here as supplied to the detector
  TH1F *hist_W;
  TH1F *hist_Q2;
  TH1F *hist_x;
  TH1F *hist_t;
  TH1F *hist_cmphi;
  TH1F *hist_trentphi;
  TH1F *hist_trentphi_final;
  TH1F *hist_cmcostheta;

  //ADDED 
  TH2F *hist_Q2_x;
  TH2F *hist_x_t;
  TH2F *hist_t_phi;
  TH2F *hist_Q2_phi;
  TH2F *hist_Q2_vs_W;
  TH2F *hist_W_vs_phi;

  TH1F *hist_phi_q2;
  TH1F *hist_phi_w;
  TH1F *hist_phi_x;
  TH1F *hist_phi_t;
  TH1F *hist_phi_pT;
  TH1F *hist_phi_xb;
  TH1F *hist_phi_cmphi;
  TH1F *hist_phi_cmcostheta;

  TH1F *hist_phi_q2_final;
  TH1F *hist_phi_w_final;
  TH1F *hist_phi_x_final;
  TH1F *hist_phi_t_final;
  TH1F *hist_phi_pT_final;
  TH1F *hist_phi_xb_final;
  TH1F *hist_phi_cmphi_final;
  TH1F *hist_phi_cmcostheta_final;

  TH1F *hist_tmin_final;
  TH1F *hist_phi_cos_theta_H_final;
 
  hist_W = new TH1F("W", "W all sectors", 700, 0.5, Ebeam-2);   
  hist_W->GetXaxis()->SetTitle("W /GeV");
  hist_W->GetYaxis()->SetTitle("counts");

  hist_Q2 = new TH1F("Q2", "Q^{2} all sectors", 800, 0, Ebeam+2);   
  hist_Q2->GetXaxis()->SetTitle("Q^{2} /GeV^{2}");
  hist_Q2->GetYaxis()->SetTitle("counts");

  hist_Q2_vs_W = new TH2F("Q2_vs_W", "Q^{2} vs W all sectors", 700, 0.5, Ebeam-2, 800, 0, Ebeam+2);   
  hist_Q2_vs_W->GetXaxis()->SetTitle("W /GeV");
  hist_Q2_vs_W->GetYaxis()->SetTitle("Q^{2}");

  hist_W_vs_phi = new TH2F("W_vs_phi", "W vs phi all sectors", 180, -180, 180, 700, 0.5, Ebeam-2);   
  hist_W_vs_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_W_vs_phi->GetYaxis()->SetTitle("W /GeV");

  hist_t = new TH1F("t", "t #pi^{0} all sectors", 50, 0, 5.0);   
  hist_t->GetXaxis()->SetTitle("t / GeV^{2}");
  hist_t->GetYaxis()->SetTitle("counts");

  hist_x = new TH1F("x", "x all sectors", 500, 0, 1.25);   
  hist_x->GetXaxis()->SetTitle("x");
  hist_x->GetYaxis()->SetTitle("counts");

  hist_cmphi = new TH1F("cmphi", "#phi_{CM} for p of e p --> e p #pi^{0}", 180, -180, 180);   
  hist_cmphi->GetXaxis()->SetTitle("#phi_{CM} /deg");
  hist_cmphi->GetYaxis()->SetTitle("counts");
  hist_cmcostheta = new TH1F("cmcostheta_p", "cos(#Theta_{CM})", 100, -1, 1);   
  hist_cmcostheta->GetXaxis()->SetTitle("cos(#Theta_{CM})");
  hist_cmcostheta->GetYaxis()->SetTitle("counts");
  hist_trentphi = new TH1F("hist_trentphi","hist_trentphi",100,0.0, 360.0);
  
  //ADDED
  hist_Q2_x = new TH2F("hist_Q2_x","Q^{2} vs Bjorken x",200, 0.0, 1.1, 200, 0.0, Ebeam);  
  hist_x_t = new TH2F("hist_x_t","Xb vs -t", 100, 0.0, 4.0, 100, 0.0, 1.1);
  hist_x_t->GetXaxis()->SetTitle("-t [GeV]");
  hist_x_t->GetYaxis()->SetTitle("Xb");
  hist_t_phi = new TH2F("hist_t_phi","#phi vs -t", 200, -180.0, 180.0, 100, 0.0, 4.0);
  hist_t_phi->GetXaxis()->SetTitle("#phi [deg]");
  hist_t_phi->GetYaxis()->SetTitle("-t [GeV]");
  hist_Q2_phi = new TH2F("hist_q2_phi","#phi vs Q^{2}",200, -180.0, 180.0, 100, 0.0, Ebeam);
  hist_Q2_phi->GetXaxis()->SetTitle("#phi [deg]");
  hist_Q2_phi->GetYaxis()->SetTitle("Q^{2} [GeV^{2}]");  

  hist_phi_q2 = new TH1F("hist_phi_q2","hist_phi_q2",100,0.0,10.5);
  hist_phi_xb = new TH1F("hist_phi_xb","hist_phi_xb",100,0.0,1.1);
  hist_phi_w = new TH1F("hist_phi_w","hist_phi_w",100,0.0, 4.5);
  hist_phi_t = new TH1F("hist_phi_t","hist_phi_t", 100, 0.0, 4.0);
  hist_phi_pT = new TH1F("hist_phi_pT","hist_phi_pT", 100, 0.0, 4.0);
  hist_phi_cmphi = new TH1F("hist_phi_cmphi","hist_phi_cmphi", 50, -180.0, 180.0);
  hist_phi_cmcostheta = new TH1F("hist_phi_cmcostheta","hist_phi_cmcostheta", 20, -1.0, 1.0);
  hist_trentphi_final = new TH1F("hist_trentphi_final","hist_trentphi_final",100, 0.0, 360.0);

  hist_phi_q2_final = new TH1F("hist_phi_q2_final","hist_phi_q2_final", 100, 0.0, 10.5);
  hist_phi_xb_final = new TH1F("hist_phi_xb_final","hist_phi_xb_final", 100, 0.0, 1.1);
  hist_phi_w_final = new TH1F("hist_phi_w_final","hist_phi_w_final", 100, 0.0, 4.5);
  hist_phi_t_final = new TH1F("hist_phi_t_final","hist_phi_t_final", 20, 0.0, 4.0);
  hist_phi_pT_final = new TH1F("hist_phi_pT_final","hist_phi_pT_final", 20, 0.0, 4.0);
  hist_phi_cmphi_final = new TH1F("hist_phi_cmphi_final","hist_phi_cmphi_final", 50, -0.0, 180.0);
  hist_phi_cmcostheta_final = new TH1F("hist_phi_cmcostheta_final","hist_phi_cmcostheta_final",20,-1.0,1.0);

  hist_tmin_final = new TH1F("t_{min} final","t_{min} final", 70, -1.0, 5.0);
  hist_tmin_final->GetXaxis()->SetTitle("t_{min} (GeV^{2})");  
  hist_tmin_final->GetYaxis()->SetTitle("counts");

  hist_phi_cos_theta_H_final = new TH1F("cos_theta_H","cos(#theta_{H})",70, -1.0, 5.0);
  hist_phi_cos_theta_H_final->GetXaxis()->SetTitle("cos(#theta_{H})");  
  hist_phi_cos_theta_H_final->GetYaxis()->SetTitle("counts");

  // end phase space plots


  out->mkdir("event_selection_kin");				
  out->cd ("event_selection_kin");

  TH1D *h_miss_perp;
  TH1D *h_colin_pr;
  TH1D *h_colin_kp;
  TH1D *h_colin_km;
  TH1D *h_colin_ang_pr;
  TH1D *h_colin_ang_kp;
  TH1D *h_colin_ang_km;
  TH1D *h_ang_elpr;
  TH1D *h_ang_kpkm;

  TH1D *h_miss_perp_final;
  TH1D *h_colin_pr_final;
  TH1D *h_colin_kp_final;
  TH1D *h_colin_km_final;
  TH1D *h_ang_elpr_final;
  TH1D *h_ang_kpkm_final;
  
  h_miss_perp = new TH1D("h_miss_perp","h_miss_perp",75, -0.5, 1.5); 
  h_colin_pr = new TH1D("h_colin_pr","h_colin_pr",100,-10, 60.0);
  h_colin_kp = new TH1D("h_colin_kp","h_colin_kp",100,-10, 60.0);
  h_colin_km = new TH1D("h_colin_km","h_colin_km",100,-10, 60.0);   
  h_colin_ang_pr = new TH1D("h_colin_ang_pr","h_colin_ang_pr",100,-10, 60.0);
  h_colin_ang_kp = new TH1D("h_colin_ang_kp","h_colin_ang_kp",100,-10, 60.0);
  h_colin_ang_km = new TH1D("h_colin_ang_km","h_colin_ang_km",100,-10, 60.0);   
  h_ang_elpr = new TH1D("h_ang_elpr","h_ang_elpr",100,-10, 60.0);
  h_ang_kpkm = new TH1D("h_ang_kpkm","h_ang_kpkm",100,-10, 60.0);  

  h_miss_perp_final = new TH1D("h_miss_perp_final","h_miss_perp_final",75, -0.5, 1.5); 
  h_colin_pr_final = new TH1D("h_colin_pr_final","h_colin_pr_final",100,-10, 60.0);
  h_colin_kp_final = new TH1D("h_colin_kp_final","h_colin_kp_final",100,-10, 60.0);
  h_colin_km_final = new TH1D("h_colin_km_final","h_colin_km_final",100,-10, 60.0);   
  h_ang_elpr_final = new TH1D("h_ang_elpr_final","h_ang_elpr_final",100,-10, 60.0);
  h_ang_kpkm_final = new TH1D("h_ang_kpkm_final","h_ang_kpkm_final",100,-10, 60.0);  

  out->mkdir("masses");
  out->cd("masses");

  TH1D *h_eX;
  TH1D *h_epX;
  TH1D *h_epkpX;
  TH1D *h_epkmX;
  TH1D *h_epkpkmX;
  TH1D *h_epkpkmXMM2;
  TH1D *h_ekpX;
  TH1D *h_ekpkmX;
  TH1D *h_ekpkmXMM2;
  TH1D *h_ekpkmXPx;
  TH1D *h_ekpkmXPy;
  TH1D *h_ekpkmXPz;
  TH1D *h_inv_kpkm;

  TH1D *h_eX_final;
  TH1D *h_epX_final;
  TH1D *h_epkpX_final;
  TH1D *h_epkmX_final;
  TH1D *h_epkpkmX_final;
  TH1D *h_epkpkmXMM2_final;
  TH1D *h_ekpX_final;
  TH1D *h_ekpkmX_final;
  TH1D *h_ekpkmXMM2_final;
  TH1D *h_ekpkmXPx_final;
  TH1D *h_ekpkmXPy_final;
  TH1D *h_ekpkmXPz_final;
  TH1D *h_inv_kpkm_final;
  TH1D *h_inv_kpkm_final_ekpkm;

  TH2D *h2_epX_eX;
  TH2D *h2_epkpX_eX;
  TH2D *h2_epkmX_eX;
  TH2D *h2_epkpkmX_eX;
  TH2D *h2_ekpkmX_eX;

  TH2D *h2_epkpX_epX;
  TH2D *h2_epkmX_epX;
  TH2D *h2_ekpkmX_epX;
  TH2D *h2_epkpkmX_epX;

  TH2D *h2_epkpX_epkmX;
  TH2D *h2_epkpkmX_epkpX;
  TH2D *h2_epkpkmX_epkmX;

  TH2D *h2_kpkm_pkm;

  h2_epX_eX = new TH2D("h2_epX_eX","h2_epX_eX",100,-6.5,4.0, 100, 3.5, 12.0);
  h2_epkpX_eX = new TH2D("h2_epkpX_eX","h2_epkpX_eX",100,-6.5, 3.0, 100, 3.5, 12.0);
  h2_epkmX_eX = new TH2D("h2_epkmX_eX","h2_epkmX_eX",100,-6.5, 3.0, 100, 3.5, 12.0);
  h2_epkpkmX_eX = new TH2D("h2_epkpkmX_eX","h2_epkpkmX_eX",100, -4.0, 3.0, 100, 3.5, 12.0);
  h2_ekpkmX_eX = new TH2D("h2_ekpkmX_eX","h2_ekpkmX_eX",100, -4.0, 4.0, 100, 3.5, 12.0);

  h2_epkpX_epX = new TH2D("h2_epkpX_epX","h2_epkpX_epX",200, -6.0, 4.0, 200, -6.0, 4.0);
  h2_epkmX_epX = new TH2D("h2_epkmX_epX","h2_epkmX_epX",200, -6.0, 4.0, 200, -6.0, 4.0);
  h2_ekpkmX_epX = new TH2D("h2_ekpkmX_epX","h2_ekpkmX_epX",150, -3.0, 6.0, 150, -3.0, 6.0);
  h2_epkpkmX_epX = new TH2D("h2_epkpkmX_epX","h2_epkpkmX_epX",100, -4.0, 3.0, 100, -6.0, 4.0);

  h2_epkpX_epkmX = new TH2D("h2_epkpX_epkmX","h2_epkpX_epkmX",150, -3.0, 4.0, 150, -3.0, 4.0);
  h2_epkpkmX_epkpX = new TH2D("h2_epkpkmX_epkpX","h2_epkpkmX_epkpX",150, -6.0, 4.0, 150, -6.0, 4.0);
  h2_epkpkmX_epkmX = new TH2D("h2_epkpkmX_epkmX","h2_epkpkmX_epkmX",150, -6.0, 4.0, 150, -6.0, 4.0);
  
  h2_kpkm_pkm = new TH2D("h2_kpkm_pkm","h2_kpkm_pkm", 70, 1.4, 2.0, 70, 0.95, 1.2);

  h_eX = new TH1D("h_eX","h_eX",200, -1.0, 6.0);
  h_epX = new TH1D("h_epX","h_epX",200, -6.5, 4.0);
  h_epkpX = new TH1D("h_epkpX","h_epkpX",200, -6.50, 3.0);
  h_epkmX = new TH1D("h_epkmX","h_epkmX",200, -6.50, 3.0);
  h_epkpkmX = new TH1D("h_epkpkmX","h_epkpkmX",200, -3.0, 3.0);
  h_epkpkmXMM2 = new TH1D("h_epkpkmXMM2","h_epkpkmXMM2",75, -3.0, 3.0);
  h_ekpX = new TH1D("h_ekpX","h_ekpX",200, -6.0, 6.0);
  h_ekpkmX = new TH1D("h_ekpkmX","h_ekpkmX",100, -6.0, 4.0);
  h_ekpkmXMM2 = new TH1D("h_ekpkmXMM2","h_ekpkmXMM2",200, -6.0, 6.0);
  h_ekpkmXPx = new TH1D("h_ekpkmXPx","h_ekpkmXPx",200, -6.0, 6.0);
  h_ekpkmXPy = new TH1D("h_ekpkmXPy","h_ekpkmXPy",200, -6.0, 6.0);
  h_ekpkmXPz = new TH1D("h_ekpkmXPz","h_ekpkmXPz",200, -6.0, 6.0);
  h_inv_kpkm = new TH1D("h_inv_kpkm","h_inv_kpkm",75,0.8,1.5); 

  h_eX_final = new TH1D("h_eX_final","h_eX_final",200, -1.0, 6.0);
  h_epX_final = new TH1D("h_epX_final","h_epX_final",200, -6.5, 4.0);
  h_epkmX_final = new TH1D("h_epkmX_final","h_epkmX_final",200, -6.50, 3.0);
  h_epkpX_final = new TH1D("h_epkpX_final","h_epkpX_final",200, -6.50, 3.0);
  h_epkpkmX_final = new TH1D("h_epkpkmX_final","h_epkpkmX_final",200, -3.0, 3.0);
  h_epkpkmXMM2_final = new TH1D("h_epkpkmXMM2_final","h_epkpkmXMM2_final",75, -3.0, 3.0);
  h_ekpX_final = new TH1D("h_ekpX_final","h_ekpX_final",200, -1.0, 6.0);
  h_ekpkmX_final = new TH1D("h_ekpkmX_final","h_ekpkmX_final",100, -1.0, 4.0);
  h_ekpkmXMM2_final = new TH1D("h_ekpkmXMM2_final","h_ekpkmXMM2_final",200, -1.0, 6.0);
  h_ekpkmXPx_final = new TH1D("h_ekpkmXPx_final","h_ekpkmXPx_final",200, -2.0, 6.0);
  h_ekpkmXPy_final = new TH1D("h_ekpkmXPy_final","h_ekpkmXPy_final",200, -2.0, 6.0);
  h_ekpkmXPz_final = new TH1D("h_ekpkmXPz_final","h_ekpkmXPz_final",200, -2.0, 6.0);
  h_inv_kpkm_final = new TH1D("h_inv_kpkm_final","h_inv_kpkm_final",75,0.8,1.5); 
  h_inv_kpkm_final_ekpkm = new TH1D("h_inv_kpkm_final_ekpkm","h_inv_kpkm_final_ekpkm",75,0.8,1.5); 

  cout << "Analysing Tree: " << inTree << endl;
  cout << "Event Loop starting ... " << endl;
  cout << endl;
  
  int channel_counter = 0;
  double toDeg = 180.0/TMath::Pi();
  TLorentzVector lv_ebeam(0,0,Ebeam, Ebeam);
  
  int nentries = 1000;// anaTree->GetEntriesFast();
  for(Int_t k=0; k < nentries; k++){    

    anaTree->GetEntry(k);

    if(process_Events > 0 && k == process_Events) break;

    /// progress:

    if(k % 100000 == 0){

      double events = anaTree->GetEntriesFast();
      double percent = k/(events/100);
      
      printf("Analysing event number %i of %.00f (%.01f percent)\n", k, events, percent);
    }

    hel = 0;
    hel = helicity;

    // assign number of particles of each type PER EVENT:
    // (SIZE OF ARRAYS)
    ele_count = p4_ele_px->size();
    prot_count = p4_prot_px->size();
    neutr_count = p4_neutr_px->size();
    pip_count = p4_pip_px->size();
    pim_count = p4_pim_px->size();
    Kp_count = p4_Kp_px->size();
    Km_count = p4_Km_px->size();
    phot_count = p4_phot_px->size();
   
    /// initalisatize variables:

    TLorentzVector beam(0,0,Ebeam,Ebeam);
    TLorentzVector target(0,0,0,0.93827);
    double phimass = 1.019;

    // index of selected particle for the different categories:
    // out tree variables
    select_ele = 0;  
    select_prot_1 = 0;
    evcat = 0; 
    E_ele = 0; px_ele = 0; py_ele = 0; pz_ele = 0;
    E_prot_1 = 0; px_prot_1 = 0; py_prot_1 = 0; pz_prot_1 = 0;
    //ADDED
    E_kaonP_1 = 0; px_kaonP_1 = 0; py_kaonP_1 = 0; pz_kaonP_1 = 0;
    E_kaonM_1 = 0; px_kaonM_1 = 0; py_kaonM_1 = 0; pz_kaonM_1 = 0;

    W = 0; Q2 = 0; nu = 0; x = 0; y = 0;
    t1 = 0; cmphi1 = 0; cmcostheta1 = 0; pt1 = 0; eta1 = 0; z1 = 0;

    M_e_p_X_miss = 0; M_e_p_X_miss2 = 0;

    el_sect=0; pr_sect=0; kp_sect=0; km_sect=0;
    
    //  Assign tree components to Lorentz vectors:

    for(Int_t i = 0; i < BUFFER; i++){ 
      if(ele_count > i){ele[i].SetPxPyPzE(p4_ele_px->at(i),p4_ele_py->at(i),p4_ele_pz->at(i),p4_ele_E->at(i));}
      else{ele[i].SetPxPyPzE(0, 0, 0, 0);}
      if(prot_count > i){prot[i].SetPxPyPzE(p4_prot_px->at(i),p4_prot_py->at(i),p4_prot_pz->at(i),p4_prot_E->at(i));}
      else{prot[i].SetPxPyPzE(0, 0, 0, 0);}
      if(neutr_count > i){neutr[i].SetPxPyPzE(p4_neutr_px->at(i),p4_neutr_py->at(i),p4_neutr_pz->at(i),p4_neutr_E->at(i));}
      else{neutr[i].SetPxPyPzE(0, 0, 0, 0);}
      if(pip_count > i){pip[i].SetPxPyPzE(p4_pip_px->at(i),p4_pip_py->at(i),p4_pip_pz->at(i),p4_pip_E->at(i));}
      else{pip[i].SetPxPyPzE(0, 0, 0, 0);}
      if(pim_count > i){pim[i].SetPxPyPzE(p4_pim_px->at(i),p4_pim_py->at(i),p4_pim_pz->at(i),p4_pim_E->at(i));}
      else{pim[i].SetPxPyPzE(0, 0, 0, 0);}
      if(Kp_count > i){Kp[i].SetPxPyPzE(p4_Kp_px->at(i),p4_Kp_py->at(i),p4_Kp_pz->at(i),p4_Kp_E->at(i));}
      else{Kp[i].SetPxPyPzE(0, 0, 0, 0);}
      if(Km_count > i){Km[i].SetPxPyPzE(p4_Km_px->at(i),p4_Km_py->at(i),p4_Km_pz->at(i),p4_Km_E->at(i));}
      else{Km[i].SetPxPyPzE(0, 0, 0, 0);}
      if(phot_count > i){phot[i].SetPxPyPzE(p4_phot_px->at(i),p4_phot_py->at(i),p4_phot_pz->at(i),p4_phot_E->at(i));}
      else{phot[i].SetPxPyPzE(0, 0, 0, 0);}      
    }

    //  Build event:
    TLorentzVector lv_el_max(0,0,0,0);
    TLorentzVector lv_pr_max(0,0,0,0);
    TLorentzVector lv_kp_max(0,0,0,0);
    TLorentzVector lv_km_max(0,0,0,0);
    
    int min_el=0;
    int min_pr=0;
    int min_kp=0;
    int min_km=0;

    //double check that we are working with the most energetic particles
    for( int ii = 0; ii < p4_ele_px->size(); ii++ ){
      TLorentzVector lv_temp = ele[ii];
      if( lv_temp.E() > lv_el_max.E() ){
	lv_el_max = lv_temp;
      }
    }

    for( int ii = 0; ii < p4_prot_px->size(); ii++ ){
      TLorentzVector lv_temp = prot[ii];
      if( lv_temp.E() > lv_pr_max.E() ){
	lv_pr_max = lv_temp;
      }
    }

    for( int ii = 0; ii < p4_Kp_px->size(); ii++ ){
	TLorentzVector lv_temp = Kp[ii];
	if( lv_temp.E() > lv_kp_max.E() ){
	  lv_kp_max = lv_temp;
	}
    }
      
    for( int ii = 0; ii < p4_Km_px->size(); ii++ ){
      TLorentzVector lv_temp = Km[ii];
      if( lv_temp.E() > lv_km_max.E() ){
	lv_km_max = lv_temp;
      }
    }

    //particle kinematics before any selection criteria or cuts
    if( ele_count > min_el ){
      hist_all_electron_p->Fill(lv_el_max.P());
      hist_all_electron_theta->Fill(lv_el_max.Theta()*toDeg);
      hist_all_electron_phi->Fill(lv_el_max.Phi()*toDeg);
      hist_all_electron_p_vs_theta->Fill(lv_el_max.Theta()*toDeg, lv_el_max.P());
      hist_all_electron_p_vs_phi->Fill(lv_el_max.Phi()*toDeg, lv_el_max.P());
      hist_all_electron_theta_vs_phi->Fill(lv_el_max.Phi()*toDeg, lv_el_max.Theta()*toDeg);
    }

    if( prot_count > min_pr ){
      hist_all_proton_p->Fill(lv_pr_max.P());
      hist_all_proton_theta->Fill(lv_pr_max.Theta()*toDeg);
      hist_all_proton_phi->Fill(lv_pr_max.Phi()*toDeg);
      hist_all_proton_p_vs_theta->Fill(lv_pr_max.Theta()*toDeg, lv_pr_max.P());
      hist_all_proton_p_vs_phi->Fill(lv_pr_max.Phi()*toDeg, lv_pr_max.P());
      hist_all_proton_theta_vs_phi->Fill(lv_pr_max.Phi()*toDeg, lv_pr_max.Theta()*toDeg);
    }

    if( Kp_count > min_kp ){
      hist_all_kp_p->Fill(lv_kp_max.P());
      hist_all_kp_theta->Fill(lv_kp_max.Theta()*toDeg);
      hist_all_kp_phi->Fill(lv_kp_max.Phi()*toDeg);
      hist_all_kp_p_vs_theta->Fill(lv_kp_max.Theta()*toDeg, lv_kp_max.P());
      hist_all_kp_p_vs_phi->Fill(lv_kp_max.Phi()*toDeg, lv_kp_max.P());
      hist_all_kp_theta_vs_phi->Fill(lv_kp_max.Phi()*toDeg, lv_kp_max.Theta()*toDeg);
    }

    if( Km_count > min_km ){
      hist_all_km_p->Fill(lv_km_max.P());
      hist_all_km_theta->Fill(lv_km_max.Theta()*toDeg);
      hist_all_km_phi->Fill(lv_km_max.Phi()*toDeg);
      hist_all_km_p_vs_theta->Fill(lv_km_max.Theta()*toDeg, lv_km_max.P());
      hist_all_km_p_vs_phi->Fill(lv_km_max.Phi()*toDeg, lv_km_max.P());
      hist_all_km_theta_vs_phi->Fill(lv_km_max.Phi()*toDeg, lv_km_max.Theta()*toDeg);
    }

    TLorentzVector lv_q = lv_ebeam - lv_el_max;
    TLorentzVector lv_eX = lv_ebeam + target - lv_el_max;
    TLorentzVector lv_ekpX = lv_ebeam + target - lv_el_max - lv_kp_max;
    TLorentzVector lv_epX = lv_ebeam + target - lv_el_max - lv_pr_max;
    TLorentzVector lv_epkpX = lv_ebeam + target - lv_el_max - lv_pr_max - lv_kp_max;
    TLorentzVector lv_epkpkmX = lv_ebeam + target - lv_el_max - lv_pr_max - lv_kp_max - lv_km_max;
    TLorentzVector lv_ekpkmX = lv_ebeam + target - lv_el_max - lv_kp_max - lv_km_max;
    TLorentzVector lv_epkmX = lv_ebeam + target - lv_el_max - lv_pr_max - lv_km_max;





    double el_e = lv_el_max.P();
    double pr_p = lv_pr_max.P();
    double kp_theta = lv_kp_max.Theta()*toDeg;
    double km_theta = lv_km_max.Theta()*toDeg;
    double kp_p = lv_kp_max.P();
    double km_p = lv_km_max.P();
    double miss_perp = sqrt( lv_epkpkmX.Px()*lv_epkpkmX.Px() + lv_epkpkmX.Py()*lv_epkpkmX.Py() );
    double colin_pr = TMath::ACos( lv_pr_max.Vect().Dot(lv_ekpkmX.Vect())/(lv_pr_max.Vect().Mag() * lv_ekpkmX.Vect().Mag() ) ) * toDeg;
    double colin_kp = TMath::ACos( lv_kp_max.Vect().Dot(lv_epkmX.Vect())/(lv_kp_max.Vect().Mag() * lv_epkmX.Vect().Mag() ) ) * toDeg;
    double colin_km = TMath::ACos( lv_km_max.Vect().Dot(lv_epkpX.Vect())/(lv_km_max.Vect().Mag() * lv_epkpX.Vect().Mag() ) ) * toDeg;

    double colin_ang_pr = (lv_pr_max.Theta() - lv_ekpkmX.Theta())*toDeg;
    double colin_ang_kp = (lv_kp_max.Theta() - lv_epkmX.Theta())*toDeg;
    double colin_ang_km = (lv_km_max.Theta() - lv_epkpX.Theta())*toDeg;


    double W = (lv_q + target).M();
    double q2 = 4.0*lv_ebeam.E()*lv_el_max.E()*sin(lv_el_max.Theta()/2.0)*sin(lv_el_max.Theta()/2.0);
    double xb = q2/(2.0 * target.M()*(lv_q.E()));
    double t = 2.0*target.M()*(lv_pr_max.E() - target.M());
    
    //    double W2 = m_p*m_p + q2/xb  - 1;
    /*double E1CM = (W2 + q2 + m_p*m_p)/(2.p*sqrt(W2));
    double P1CM = sqrt(E1CM - m_p*m_p);
    double E3CM = (W2 + m_p*m_p - m_phi*m_phi)/(2*sqrt(W2));
    double P3CM = sqrt(E3CM*E3CM - m_p*m_p);

    double t_min = pow(q2 + m_phi*m_phi,2)/(4.0*W2) - pow(P1CM - P3CM,2);
    */
    double ang_H = lv_kp_max.Vect().Angle(lv_km_max.Vect());
    double cos_angle_H = cos(ang_H);
    std::cout<< " cos angle H (no boost ) " << cos_angle_H << std::endl;




    // phi angle in trento convention
    TVector3 v3L = (lv_ebeam.Vect()).Cross(lv_el_max.Vect());
    TVector3 v3H = (lv_pr_max.Vect()).Cross((lv_ebeam-lv_el_max).Vect());
    float phi = TMath::RadToDeg()* v3L.Angle(v3H);
    if(v3L.Dot(lv_pr_max.Vect())<0)phi = 360.-phi;

    bool high_el_e_cut = el_e > 1.5;
    bool low_p_pr_cut = pr_p > 0.5 && pr_p < 4.0;
    bool low_theta_kp_cut = kp_theta < 35.0;
    bool low_theta_km_cut = km_theta < 40.0;
    bool low_p_kp_cut = kp_p < 2.5;
    bool low_p_km_cut = km_p < 2.5;
    bool w_cut = W > 2.0;
    bool q2_cut = q2 > 1.0;
    bool epkpkm_mm2_cut = fabs(lv_epkpkmX.M2()) <  0.6;
    bool miss_perp_cut = miss_perp < 0.5;
    bool colin_pr_cut = colin_pr < 30.0;
    bool colin_kp_cut = colin_kp < 20.0;
    bool colin_km_cut = colin_km < 20.0;

    if(  ele_count > 0 ){
      h_eX->Fill(lv_eX.M());      
    }

    if( q2_cut && w_cut ){

      hist_W_vs_phi->Fill(lv_el_max.Phi()*toDeg,W);      
      hist_W->Fill(W);
      hist_Q2->Fill(q2);
      hist_Q2_vs_W->Fill(W,q2);
      hist_x->Fill(xb);
      hist_Q2_x->Fill(xb,q2);
      hist_Q2_phi->Fill(lv_el_max.Phi()*toDeg,q2);     
      //hist_phi_cmphi->Fill(phi);
      hist_trentphi->Fill(phi);
	
      if( prot_count > 0 && high_el_e_cut ){
	h_epX->Fill(lv_epX.M2());
	hist_t->Fill(t);
	hist_x_t->Fill(t,xb);
	hist_t_phi->Fill(lv_pr_max.Phi()*toDeg,t);	
	h2_epX_eX->Fill(lv_epX.M2(), lv_eX.M2());
      }

      if( prot_count > 0 && Kp_count > 0 && high_el_e_cut ){
	h_epkpX->Fill(lv_epkpX.M2());	
	h2_epkpX_eX->Fill(lv_epkpX.M2(), lv_eX.M2());
	h2_epkpX_epX->Fill(lv_epkpX.M2(), lv_epX.M2());

      }

      if( prot_count > 0 && Km_count > 0 && high_el_e_cut ){
	h_epkmX->Fill(lv_epkmX.M2());	
	h2_epkmX_eX->Fill(lv_epkmX.M2(), lv_eX.M2());
	h2_epkmX_epX->Fill(lv_epkmX.M2(), lv_epX.M2());
      }

      if( prot_count > 0 && Kp_count > 0 && Km_count > 0 && high_el_e_cut ){
	h_epkpkmX->Fill(lv_epkpkmX.M());
	h_epkpkmXMM2->Fill(lv_epkpkmX.M2());	
	h2_epkpkmX_eX->Fill(lv_epkpkmX.M2(), lv_eX.M2());
	h2_epkpX_epkmX->Fill(lv_epkpX.M2(), lv_epkmX.M2());
	h2_epkpkmX_epkpX->Fill(lv_epkpkmX.M2(), lv_epkpX.M2());
	h2_epkpkmX_epkmX->Fill(lv_epkpkmX.M2(), lv_epkmX.M2());
	h2_epkpkmX_epX->Fill(lv_epkpkmX.M2(), lv_epX.M2());
      }
	
    	         
      if( Kp_count > 0 && Km_count > 0 ){ //high_el_e_cut ){
	h2_ekpkmX_eX->Fill(lv_ekpkmX.M2(), lv_eX.M2());
	h2_ekpkmX_epX->Fill(lv_ekpkmX.M2(), lv_epX.M2());
	
      }

      if( Kp_count > 0 && Km_count > 0 && high_el_e_cut ){
	
	TLorentzVector lv_phi = lv_kp_max + lv_km_max;
	h_ekpX->Fill(lv_ekpX.M2());
	h_ekpkmX->Fill(lv_ekpkmX.M());
	h_ekpkmXMM2->Fill(lv_ekpkmX.M2());
	h_ekpkmXPx->Fill(lv_ekpkmX.Px());
	h_ekpkmXPy->Fill(lv_ekpkmX.Py());
	h_ekpkmXPz->Fill(lv_ekpkmX.Pz());	 
	if( lv_ekpkmX.M2() < 1.15 && lv_ekpkmX.M2() > 0.7 ){
	  h_inv_kpkm_final_ekpkm->Fill(lv_phi.M());
	}

    TVector3 boost = -lv_phi.BoostVector();
    lv_kp_max.Boost(boost);
    lv_km_max.Boost(boost);
    double cos_theta_H = lv_kp_max.Vect().Angle(lv_km_max.Vect());
    std::cout << " cos theta H  (boost) " << cos_theta_H << std::endl;

	
      }

      if( prot_count > 0 && Kp_count > 0 && Km_count > 0 && high_el_e_cut ){
	  h_colin_ang_pr->Fill(colin_ang_pr);
	  h_colin_ang_kp->Fill(colin_ang_kp);
	  h_colin_ang_km->Fill(colin_ang_km);

	  h_colin_pr->Fill(colin_pr);
	  h_colin_kp->Fill(colin_kp);
	  h_colin_km->Fill(colin_km);

	if(low_theta_kp_cut && low_theta_km_cut && low_p_kp_cut && low_p_km_cut && low_p_pr_cut ){
	  
	  hist_phi_q2->Fill(q2);
	  hist_phi_xb->Fill(xb);
	  hist_phi_w->Fill(W);
	  hist_phi_t->Fill(t);
	  hist_phi_pT->Fill(miss_perp);
	  
	  h_miss_perp->Fill(miss_perp);
	  TLorentzVector lv_phi = lv_kp_max + lv_km_max;
	  TLorentzVector lv_pkm = lv_pr_max + lv_km_max;
	  h_inv_kpkm->Fill(lv_phi.M());

	  if( epkpkm_mm2_cut && miss_perp_cut && colin_pr_cut && colin_km_cut ){

	      h_eX_final->Fill(lv_eX.M());
	      h_epX_final->Fill(lv_epX.M2());
	      h_epkpX_final->Fill(lv_epkpX.M2());
	      h_epkpkmX_final->Fill(lv_epkpkmX.M2());
	      h_epkpkmXMM2_final->Fill(lv_epkpkmX.M2());
	      h_ekpX_final->Fill(lv_ekpX.M2());
	      h_ekpkmX_final->Fill(lv_ekpkmX.M2());
	      h_ekpkmXMM2_final->Fill(lv_ekpkmX.M2());
	      h_ekpkmXPx_final->Fill(lv_ekpkmX.Px());
	      h_ekpkmXPy_final->Fill(lv_ekpkmX.Py());
	      h_ekpkmXPz_final->Fill(lv_ekpkmX.Pz());

	      hist_phi_q2_final->Fill(q2);
	      hist_phi_xb_final->Fill(xb);
	      hist_phi_w_final->Fill(W);
	      hist_phi_t_final->Fill(t);
	      hist_phi_pT_final->Fill(miss_perp);
	      hist_trentphi_final->Fill(phi);

	      h_miss_perp_final->Fill(miss_perp);
	      h_colin_pr_final->Fill(colin_pr);
	      h_colin_kp_final->Fill(colin_kp);
	      h_colin_km_final->Fill(colin_km);
	      h_inv_kpkm_final->Fill(lv_phi.M());

	      h2_kpkm_pkm->Fill(lv_pkm.M(), lv_phi.M());
	      
	      //hist_tmin_final->Fill(t_min);
	      hist_phi_cos_theta_H_final->Fill(cos_angle_H);

	      W_out = W;
	      Q2_out=q2;
	      x_out=xb;
	      hel=helicity;
	      cmphi1_out=phi;
	      //cmcostheta1_out= cm_cos_theta;
	      t1_out=t;


	      ////////////////////////////////////////////
	      // added variables to plot from the phi generator
	      /*	      double nu = q2/(2.*target.M()*xb);
	      double y = nu/lv_beam.E();
	      double e1 = TMath::Power(y*xb*target.M(),2)/q2;
	      double EPS = (1.0-y-e1)/(1-y+y*y/2+e1);
	      if(EPS<0.||EPS>1.)return 0.;
	      double W2 = target.M()*target.M() + 2.0*target.M()*nu - q2;
	      if(W2<TMath::Power(ProtonMass+phimass,2.))return 0.;
	      double W = TMath::Sqrt(W2);

	      double Wth = target.M()+phimass;
	      double cT = 400*(1-Wth*Wth/W2)*TMath::Power(W,0.32);
	      double sigmT = cT/TMath::Power(1+Q2/(phimass*phimass),3.0);
	      double R = 0.4*Q2/(phimass*phimass);
	      double sigmL = R*sigmT;
	      
	      double E1cm = target.M()*(target.M() + nu)/w;
	      double P1cm = target.M()*TMath::Sqrt(nu*nu+Q2)/W;
	      double E2cm = (W2 + target.M()*target.M()-phimass*phimass)/(2.*W);
	      if(E2cm<target.M())return 0.;
	      double P2cm = TMath::Sqrt(E2cm*E2cm-target.M()*target.M());
	      double E3cm = (W2-phimass*phimass+target.M()*target.M())/(2*W);
	      double P3cm = TMath::Sqrt(E3cm*E3cm-target.M()*target.M());

	      double tmax = 2.0*(target.M()*target.M() - E1cm*E2cm + P1cm*P2cm);
	      double tmin = 2.0*(target.M()*target.M() - E1cm*E2cm - P1cm*P2cm);
	      if(t<tmin||t>tmax)return 0.;
	      
	      tmin = -(TMath::Power(Q2+phimass*phimass,2)/4./W2-TMath::Power(P1cm-P3cm,2));
	      
	      double T  = TMath::Abs(t);
	      double T0 = TMath::Abs(tmin);
	      double B = 2.2 + 4*0.24*TMath::Log(W);
	      double FF = TMath::Exp(B*(T0-T));
	      TVector3 boost = -lv_phi.BoostVector();
	      lv_kp_max.Boost(boost);
	      lv_km_max.Boost(boost);
	      double cos_theta_H = lv_kp_max.Vect().Angle(lv_km_max.Vect());
	      std::cout << " cos theta H " << cos_theta_H << std::endl;

	      */	      
	      //double r04_00 = EPS*R/(1+EPS*R);
	      ///double CM_F = ( (1.-r04_00) + (3.*r04_00-1.)*cos_theta_H*cos_theta_H )*0.75;
	      //double phi_modul = ( 1 + 0.5*TMath::Cos(PHI_G))/1.5;
	      
	      
	      

	      //std::cout << " helicity " << helicity << " " << hel << std::endl;
	      
	      out_tree1.Fill();
	    }	    
	}
      }
    }	  
  	



      

  }


  out->Write();
  out->Close();

  return 0;
}

/// ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///  ROOT Macro to select physics reactions from filtered clas 6 and clas 12 data  (simple version for inklusive electron scattering)
//
///  Stefan Diehl  (sdiehl@jlab.org)
///
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
double Ebeam = 10.594;

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

/// /////////////////////////////////////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////////////////////////////////////
/// define functions

double kin_W(TLorentzVector ele);
double kin_Q2(TLorentzVector ele);
double kin_x(TLorentzVector ele);
double kin_y(TLorentzVector ele);
double kin_nu(TLorentzVector ele);

double kin_pT(TLorentzVector ele, TLorentzVector hadron);
double kin_eta(TLorentzVector ele, TLorentzVector hadron);
double kin_z(TLorentzVector ele, TLorentzVector hadron);
double kin_t(TLorentzVector ele, TLorentzVector hadron);

double kin_cmphi(TLorentzVector ele, TLorentzVector hadron);
double kin_cmcostheta(TLorentzVector ele, TLorentzVector hadron);

double kin_mismass(TLorentzVector ele, TLorentzVector hadron);
double kin_mismass2(TLorentzVector ele, TLorentzVector hadron);

TLorentzVector miss_X(TLorentzVector ele, TLorentzVector hadron);

double alpha_e_X(TLorentzVector ele, TLorentzVector hadron);
double alpha_p1p2(TLorentzVector particle1, TLorentzVector particle2);

double Phi_Mass(TLorentzVector kaonP, TLorentzVector kaonM);

//ADDED 
double kin_epkXMass(TLorentzVector el, TLorentzVector pr, TLorentzVector kp);
double kin_epKpKmXMass(TLorentzVector el, TLorentzVector pr, TLorentzVector kp, TLorentzVector km);



/// /////////////////////////////////////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// main
///

Int_t selector_simple_v2( Char_t *inFile, Char_t *outputfile)
{
		
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
  /*  anaTree->SetBranchAddress("prot_beta_final",&prot_beta_final);
  anaTree->SetBranchAddress("Kp_beta_final",&Kp_beta_final);
  anaTree->SetBranchAddress("Km_beta_final",&Km_beta_final);
  anaTree->SetBranchAddress("el_sector",&el_sector);
  anaTree->SetBranchAddress("prot_sector",&prot_sector);
  anaTree->SetBranchAddress("Kp_sector",&Kp_sector);
  anaTree->SetBranchAddress("Km_sector",&Km_sector);
  anaTree->SetBranchAddress("el_sector",&el_sector);
  anaTree->SetBranchAddress("prot_sector",&prot_sector);
  anaTree->SetBranchAddress("Kp_sector",&Kp_sector);
  anaTree->SetBranchAddress("Km_sector",&Km_sector);
  */
  /// /////////////////////////////////////////////////////////////////////////////    
  ///  create output tree:
  /// /////////////////////////////////////////////////////////////////////////////

  out = new TFile(outputfile, "RECREATE");

  // categroy 1:  e p -->  e p X

  /*    TTree out_tree1("events_epX","events_epX");

	out_tree1.Branch("W", &W);
	out_tree1.Branch("Q2", &Q2);
	out_tree1.Branch("x", &x);
	out_tree1.Branch("y", &y);
	out_tree1.Branch("nu", &nu);
	out_tree1.Branch("minus_t", &t1);
	out_tree1.Branch("cmphi", &cmphi1);
	out_tree1.Branch("cmcostheta", &cmcostheta1);
	out_tree1.Branch("pt", &pt1);
	out_tree1.Branch("eta", &eta1);
	out_tree1.Branch("z", &z1);
	out_tree1.Branch("missmass", &M_e_p_X_miss);
	out_tree1.Branch("missmass2", &M_e_p_X_miss2);
	out_tree1.Branch("helicity", &helicity);
	out_tree1.Branch("fcup", &fcup);
	out_tree1.Branch("E_ele", &E_ele);
	out_tree1.Branch("px_ele", &px_ele);
	out_tree1.Branch("py_ele", &py_ele);
	out_tree1.Branch("pz_ele", &pz_ele);
	out_tree1.Branch("E_prot", &E_prot_1);
	out_tree1.Branch("px_prot", &px_prot_1);
	out_tree1.Branch("py_prot", &py_prot_1);
	out_tree1.Branch("pz_prot", &pz_prot_1);

  */

  TTree out_tree1("out_tree_epPhi","out_tree_epPhi");
  out_tree1.Branch("W",&W_out);
  out_tree1.Branch("Q2", &Q2_out);
  out_tree1.Branch("x", &x_out);
  out_tree1.Branch("y", &y_out);
  out_tree1.Branch("nu", &nu_out);
  out_tree1.Branch("minus_t", &t1_out);
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
  // out_tree1.Branch("kp_sect",&kp_sect);
  ///out_tree1.Branch("km_sect",&km_sect);
  
  
  out_tree1.Branch("channel",&channel);
  //out_tree1.
  //out_tree1.Branch(

  /// ///////////////////////////////////////////////////////////////
  ///  create histograms:
  /// ///////////////////////////////////////////////////////////////

  out->mkdir("event_information");				
  out->cd ("event_information");

  TH1F *hist_helicity;
  TH1F *hist_faraday_cup;

  hist_helicity = new TH1F("hist_helicity", "helicity", 5, -2.5, 2.5);   
  hist_helicity->GetXaxis()->SetTitle("helicity");
  hist_helicity->GetYaxis()->SetTitle("counts");
  hist_faraday_cup = new TH1F("hist_faraday_cup", "faraday cup", 200, 0, 200);   
  hist_faraday_cup->GetXaxis()->SetTitle("faraday cup");
  hist_faraday_cup->GetYaxis()->SetTitle("counts");

  out->mkdir("particle_histograms_all");				
  out->cd ("particle_histograms_all");

  TH1F *hist_all_electron_p; TH1F *hist_all_electron_theta; TH1F *hist_all_electron_phi;
  TH2F *hist_all_electron_p_vs_theta; TH2F *hist_all_electron_p_vs_phi; TH2F *hist_all_electron_theta_vs_phi;

  TH1F *hist_all_proton_p; TH1F *hist_all_proton_theta; TH1F *hist_all_proton_phi;
  TH2F *hist_all_proton_p_vs_theta; TH2F *hist_all_proton_p_vs_phi; TH2F *hist_all_proton_theta_vs_phi;

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

  out->mkdir("particle_histograms_selected");				
  out->cd ("particle_histograms_selected");

  TH1F *hist_electron_p; TH1F *hist_electron_theta; TH1F *hist_electron_phi; 
  TH2F *hist_electron_p_vs_theta; TH2F *hist_electron_p_vs_phi; TH2F *hist_electron_theta_vs_phi;
  TH1F *hist_proton_p; TH1F *hist_proton_theta; TH1F *hist_proton_phi; 
  TH2F *hist_proton_p_vs_theta; TH2F *hist_proton_p_vs_phi; TH2F *hist_proton_theta_vs_phi;

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
  hist_electron_theta_vs_phi = new TH2F("hist_electron_theta_vs_phi", "electron #Theta vs phi", 180,-180,180, 140,0,140);   
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
  hist_proton_theta_vs_phi = new TH2F("hist_proton_theta_vs_phi", "proton #Theta vs phi", 180,-180,180, 140,0,140);   
  hist_proton_theta_vs_phi->GetXaxis()->SetTitle("#phi /deg");
  hist_proton_theta_vs_phi->GetYaxis()->SetTitle("#Theta /deg");

  TH2F *hist_pr_betap_2b;
  TH2F *hist_kp_betap_2b;
  TH2F *hist_km_betap_2b;

  hist_pr_betap_2b = new TH2F("hist_pr_beta_no_kpkm_anglecut_c2","hist_pr_beta_no_kpkm_anglecut_c2", 100, 0.0, 4.5, 100, 0.0, 1.1);
  hist_kp_betap_2b = new TH2F("hist_kp_beta_no_kpkm_anglecut_c2","hist_kp_beta_no_kpkm_anglecut_c2", 100, 0.0, 4.5, 100, 0.0, 1.1);
  hist_km_betap_2b = new TH2F("hist_km_beta_no_kpkm_anglecut_c2","hist_km_beta_no_kpkm_anglecut_c2", 100, 0.0, 4.5, 100, 0.0, 1.1);

  TH2F *hist_pr_betap_2a;
  TH2F *hist_kp_betap_2a;
  TH2F *hist_km_betap_2a;

  hist_pr_betap_2a = new TH2F("hist_pr_beta_kpkm_anglecut_c2","hist_pr_beta_kpkm_anglecut_c2", 100, 0.0, 4.5, 100, 0.0, 1.1);
  hist_kp_betap_2a = new TH2F("hist_kp_beta_kpkm_anglecut_c2","hist_kp_beta_kpkm_anglecut_c2", 100, 0.0, 4.5, 100, 0.0, 1.1);
  hist_km_betap_2a = new TH2F("hist_km_beta_kpkm_anglecut_c2","hist_km_beta_kpkm_anglecut_c2", 100, 0.0, 4.5, 100, 0.0, 1.1);

  TH2F *hist_pr_betap_1;
  TH2F *hist_kp_betap_1;

  hist_pr_betap_1 = new TH2F("hist_pr_beta_no_kpkm_anglecut_c1","hist_pr_beta_no_kpkm_anglecut_c1", 100, 0.0, 4.5, 100, 0.0, 1.1);
  hist_kp_betap_1 = new TH2F("hist_kp_beta_no_kpkm_anglecut_c1","hist_kp_beta_no_kpkm_anglecut_c1", 100, 0.0, 4.5, 100, 0.0, 1.1);

  TH2F *hist_pr_betap_1_final;
  TH2F *hist_kp_betap_1_final;

  hist_pr_betap_1_final = new TH2F("hist_pr_beta_no_kpkm_anglecut_c1_final","hist_pr_beta_no_kpkm_anglecut_c1_final", 100, 0.0, 4.5, 100, 0.0, 1.1);
  hist_kp_betap_1_final = new TH2F("hist_kp_beta_no_kpkm_anglecut_c1_final","hist_kp_beta_no_kpkm_anglecut_c1_final", 100, 0.0, 4.5, 100, 0.0, 1.1);

  TH2F *hist_pr_betap_final;
  TH2F *hist_kp_betap_final;
  TH2F *hist_km_betap_final;

  hist_pr_betap_final = new TH2F("hist_pr_beta_no_kpkm_anglecut_final","hist_pr_beta_no_kpkm_anglecut_final", 100, 0.0, 4.5, 100, 0.0, 1.1);
  hist_kp_betap_final = new TH2F("hist_kp_beta_no_kpkm_anglecut_final","hist_kp_beta_no_kpkm_anglecut_final", 100, 0.0, 4.5, 100, 0.0, 1.1);
  hist_km_betap_final = new TH2F("hist_km_beta_no_kpkm_anglecut_final","hist_km_beta_no_kpkm_anglecut_final", 100, 0.0, 4.5, 100, 0.0, 1.1);

  
  out->mkdir("kinematics");				
  out->cd ("kinematics");


  // Look at total phase space here as supplied to the detector
  TH1F *hist_W;
  TH1F *hist_Q2;
  TH1F *hist_x;
  TH1F *hist_y;
  TH1F *hist_nu;
  TH1F *hist_t;
  TH1F *hist_cmphi;
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


  std::vector<TH1F*> h_w_cut;
  std::vector<TH1F*> h_q2_cut;
  std::vector<TH1F*> h_x_cut;
  std::vector<TH1F*> h_y_cut;
  std::vector<TH1F*> h_nu_cut;
  std::vector<TH1F*> h_t_cut;
  std::vector<TH1F*> h_cmphi_cut;
  std::vector<TH1F*> h_cmcostheta_cut;

  std::vector<TH2F*> h_q2_x_cut;
  std::vector<TH2F*> h_x_t_cut;
  std::vector<TH2F*> h_t_phi_cut;
  std::vector<TH2F*> h_q2_phi_cut;
  std::vector<TH2F*> h_q2_w_cut;
  std::vector<TH2F*> h_w_phi_cut;
  std::vector<TH2F*> h_q2_t_cut;

  for( int c = 0; c < 12; c++ ){
    h_w_cut.push_back(new TH1F(Form("h_w_cutlvl_%d",c), Form("h_w_cutlvl_%d",c), 100, 0.0, 5.0));
    h_q2_cut.push_back(new TH1F(Form("h_q2_cutlvl_%d",c), Form("h_q2_cutlvl_%d",c), 100, 0.0, 10.0));
    h_x_cut.push_back(new TH1F(Form("h_x_cutlvl_%d",c), Form("h_x_cutlvl_%d",c), 100, 0.0, 1.10));
    h_y_cut.push_back(new TH1F(Form("h_y_cutlvl_%d",c), Form("h_y_cutlvl_%d",c), 100, 0.0, 1.0));
    h_nu_cut.push_back(new TH1F(Form("h_nu_cutlvl_%d",c), Form("h_nu_cutlvl_%d",c), 100, 0.0, 10.0));
    h_t_cut.push_back(new TH1F(Form("h_t_cutlvl_%d",c), Form("h_t_cutlvl_%d",c), 100, 0.0, 5.0));
    h_cmphi_cut.push_back(new TH1F(Form("h_cmphi_cutlvl_%d",c), Form("h_cmphi_cutlvl_%d",c), 200, -180.0, 180.0));
    h_cmcostheta_cut.push_back(new TH1F(Form("h_cmcostheta_cutlvl_%d",c), Form("h_cmcostheta_cutlvl_%d",c), 100, -1.0, 1.0));


    h_q2_x_cut.push_back( new TH2F(Form("h_q2_x_cutlvl_%d",c), Form("h_q2_x_cutlvl_%d",c), 200, 0.0, 1.1, 200, 0.0, 10.5));
    h_x_t_cut.push_back( new TH2F(Form("h_x_t_cutlvl_%d",c), Form("h_x_t_cutlvl_%d",c), 200, 0.0, 1.1, 200, 0.0, 5.0));
    h_t_phi_cut.push_back( new TH2F(Form("h_t_phi_cutlvl_%d",c), Form("h_t_phi_cutlvl_%d",c), 200, -180.0, 180.0, 200, 0.0, 5.0));
    h_q2_phi_cut.push_back( new TH2F(Form("h_q2_phi_cutlvl_%d",c), Form("h_q2_phi_cutlvl_%d",c), 200, -180.0, 180.0, 200, 0.0, 10.5));
    h_q2_w_cut.push_back( new TH2F(Form("h_q2_w_cutlvl_%d",c), Form("h_q2_w_cutlvl_%d",c), 200, 1.0, 5.0, 200, 0.0, 12.0));
    h_w_phi_cut.push_back( new TH2F(Form("h_w_phi_cutlvl_%d",c), Form("h_w_phi_cutlvl_%d",c), 200, -180.0, 180.0, 200, 0.0, 5.0));
    h_q2_t_cut.push_back( new TH2F(Form("h_q2_t_cutlvl_%d",c), Form("h_q2_t_cutlvl_%d",c), 200, 0.0, 12.5, 200, 0.0, 10.0));
    
  }

  hist_phi_q2 = new TH1F("hist_phi_q2","hist_phi_q2",100,0.0,10.5);
  hist_phi_xb = new TH1F("hist_phi_xb","hist_phi_xb",100,0.0,1.1);
  hist_phi_w = new TH1F("hist_phi_w","hist_phi_w",100,0.0, 4.5);
  hist_phi_t = new TH1F("hist_phi_t","hist_phi_t", 100, 0.0, 4.0);
  hist_phi_pT = new TH1F("hist_phi_pT","hist_phi_pT", 100, 0.0, 4.0);
  hist_phi_cmphi = new TH1F("hist_phi_cmphi","hist_phi_cmphi", 50, -180.0, 180.0);
  hist_phi_cmcostheta = new TH1F("hist_phi_cmcostheta","hist_phi_cmcostheta",20,-1.0,1.0);

  hist_phi_q2_final = new TH1F("hist_phi_q2_final","hist_phi_q2_final", 100, 0.0, 10.5);
  hist_phi_xb_final = new TH1F("hist_phi_xb_final","hist_phi_xb_final", 100, 0.0, 1.1);
  hist_phi_w_final = new TH1F("hist_phi_w_final","hist_phi_w_final", 100, 0.0, 4.5);
  hist_phi_t_final = new TH1F("hist_phi_t_final","hist_phi_t_final", 20, 0.0, 4.0);
  hist_phi_pT_final = new TH1F("hist_phi_pT_final","hist_phi_pT_final", 20, 0.0, 4.0);
  hist_phi_cmphi_final = new TH1F("hist_phi_cmphi_final","hist_phi_cmphi_final", 50, -180.0, 180.0);
  hist_phi_cmcostheta_final = new TH1F("hist_phi_cmcostheta_final","hist_phi_cmcostheta_final",20,-1.0,1.0);
  // end phase space plots

  TH1F *h_mass_phi;
  
  TH2F *h_el_ptheta;
  TH2F *h_el_pphi;
  TH2F *h_el_phitheta;

  TH2F *h_pr_ptheta;
  TH2F *h_pr_pphi;
  TH2F *h_pr_phitheta;

  TH2F *h_kp_ptheta;
  TH2F *h_kp_pphi;
  TH2F *h_kp_phitheta;

  TH2F *h_km_ptheta;
  TH2F *h_km_pphi;
  TH2F *h_km_phitheta;

  h_el_ptheta = new TH2F("h_el_ptheta","Electron #theta vs p",75, 0.0, Ebeam, 75, 0.0, 65 );
  h_el_pphi = new TH2F("h_el_pphi","Electron #phi vs p",75, 0.0, Ebeam, 200, -180.0, 180.0 );
  h_el_phitheta = new TH2F("h_el_phitheta","Electron #phi vs #theta",200, -180.0, 180.0, 75, 0.0, 65 );

  h_pr_ptheta = new TH2F("h_pr_ptheta","Proton #theta vs p",75, 0.0, Ebeam, 95, 0.0, 120 );
  h_pr_pphi = new TH2F("h_pr_pphi","Proton #phi vs p",75, 0.0, Ebeam, 150, -180.0, 180.0 );
  h_pr_phitheta = new TH2F("h_pr_phitheta","Proton #phi vs #theta", 200, -180.0, 180.0, 95, 0.0, 120 );

  h_kp_ptheta = new TH2F("h_kp_ptheta","K^{+} #theta vs p", 75, 0.0, Ebeam, 75, 0.0, 120 );
  h_kp_pphi = new TH2F("h_kp_pphi","K^{+} #phi vs p",75, 0.0, Ebeam, 200, -180.0, 180.0 );
  h_kp_phitheta = new TH2F("h_kp_phitheta","K^{+} #phi vs #theta",200, -180.0, 180.0, 95, 0.0, 120 );

  h_km_ptheta = new TH2F("h_km_ptheta","K^{-} #theta vs p",75, 0.0, Ebeam, 75, 0.0, 120 );
  h_km_pphi = new TH2F("h_km_pphi","K^{-} #phi vs p",75, 0.0, Ebeam, 200, -180.0, 180.0 );
  h_km_phitheta = new TH2F("h_km_phitheta","K^{-} #phi vs #theta",200, -180.0, 180.0, 95, 0.0, 120 );

  TH2F *h_el_ptheta_final;
  TH2F *h_el_pphi_final;
  TH2F *h_el_phitheta_final;

  TH2F *h_pr_ptheta_final;
  TH2F *h_pr_pphi_final;
  TH2F *h_pr_phitheta_final;

  TH2F *h_kp_ptheta_final;
  TH2F *h_kp_pphi_final;
  TH2F *h_kp_phitheta_final;

  TH2F *h_km_ptheta_final;
  TH2F *h_km_pphi_final;
  TH2F *h_km_phitheta_final;

  h_el_ptheta_final = new TH2F("h_el_ptheta_final","Electron #theta vs p",75, 0.0, Ebeam, 75, 0.0, 65 );
  h_el_pphi_final = new TH2F("h_el_pphi_final","Electron #phi vs p",75, 0.0, Ebeam, 200, -180.0, 180.0 );
  h_el_phitheta_final = new TH2F("h_el_phitheta_final","Electron #phi vs #theta",200, -180.0, 180.0, 75, 0.0, 65 );

  h_pr_ptheta_final = new TH2F("h_pr_ptheta_final","Proton #theta vs p",75, 0.0, Ebeam, 95, 0.0, 120 );
  h_pr_pphi_final = new TH2F("h_pr_pphi_final","Proton #phi vs p",75, 0.0, Ebeam, 150, -180.0, 180.0 );
  h_pr_phitheta_final = new TH2F("h_pr_phitheta_final","Proton #phi vs #theta", 200, -180.0, 180.0, 95, 0.0, 120 );

  h_kp_ptheta_final = new TH2F("h_kp_ptheta_final","K^{+} #theta vs p", 75, 0.0, Ebeam, 75, 0.0, 120 );
  h_kp_pphi_final = new TH2F("h_kp_pphi_final","K^{+} #phi vs p",75, 0.0, Ebeam, 200, -180.0, 180.0 );
  h_kp_phitheta_final = new TH2F("h_kp_phitheta_final","K^{+} #phi vs #theta",200, -180.0, 180.0, 95, 0.0, 120 );

  h_km_ptheta_final = new TH2F("h_km_ptheta_final","K^{-} #theta vs p",75, 0.0, Ebeam, 75, 0.0, 120 );
  h_km_pphi_final = new TH2F("h_km_pphi_final","K^{-} #phi vs p",75, 0.0, Ebeam, 200, -180.0, 180.0 );
  h_km_phitheta_final = new TH2F("h_km_phitheta_final","K^{-} #phi vs #theta",200, -180.0, 180.0, 95, 0.0, 120 );
 
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

  hist_t = new TH1F("t", "t #pi^{0} all sectors", 50, 0, 2*Ebeam);   
  hist_t->GetXaxis()->SetTitle("t / GeV^{2}");
  hist_t->GetYaxis()->SetTitle("counts");

  hist_x = new TH1F("x", "x all sectors", 500, 0, 1.25);   
  hist_x->GetXaxis()->SetTitle("x");
  hist_x->GetYaxis()->SetTitle("counts");

  hist_y = new TH1F("y", "y all sectors", 440, 0, 1.1);   
  hist_y->GetXaxis()->SetTitle("y");
  hist_y->GetYaxis()->SetTitle("counts");

  hist_nu = new TH1F("nu", "#nu all sectors", 500, 0, Ebeam);   
  hist_nu->GetXaxis()->SetTitle("#nu /GeV");
  hist_nu->GetYaxis()->SetTitle("counts");

  hist_cmphi = new TH1F("cmphi", "#phi_{CM} for p of e p --> e p #pi^{0}", 180, -180, 180);   
  hist_cmphi->GetXaxis()->SetTitle("#phi_{CM} /deg");
  hist_cmphi->GetYaxis()->SetTitle("counts");
  hist_cmcostheta = new TH1F("cmcostheta_p", "cos(#Theta_{CM})", 100, -1, 1);   
  hist_cmcostheta->GetXaxis()->SetTitle("cos(#Theta_{CM})");
  hist_cmcostheta->GetYaxis()->SetTitle("counts");

  //ADDED
  hist_Q2_x = new TH2F("hist_Q2_x","Q^{2} vs Bjorken x",200, 0.0, 1.1, 200, 0.0, Ebeam);  
  h_mass_phi = new TH1F("h_mass_phi","Mass of K^{+}K^{-}",25,0.75,2.8);

  hist_x_t = new TH2F("hist_x_t","Xb vs -t", 100, 0.0, 4.0, 100, 0.0, 1.1);
  hist_x_t->GetXaxis()->SetTitle("-t [GeV]");
  hist_x_t->GetYaxis()->SetTitle("Xb");
  hist_t_phi = new TH2F("hist_t_phi","#phi vs -t", 200, -180.0, 180.0, 100, 0.0, 4.0);
  hist_t_phi->GetXaxis()->SetTitle("#phi [deg]");
  hist_t_phi->GetYaxis()->SetTitle("-t [GeV]");

  hist_Q2_phi = new TH2F("hist_q2_phi","#phi vs Q^{2}",200, -180.0, 180.0, 100, 0.0, Ebeam);
  hist_Q2_phi->GetXaxis()->SetTitle("#phi [deg]");
  hist_Q2_phi->GetYaxis()->SetTitle("Q^{2} [GeV^{2}]");  

  out->mkdir("missing_mass");				
  out->cd ("missing_mass");

  TH1F *hist_e_p_X_mismass;
  TH1F *hist_e_p_X_mismass2;

  hist_e_p_X_mismass = new TH1F("hist_e_p_X_mismass", "e p --> e p X missing mass", 200, 0.0, 2.5);   
  hist_e_p_X_mismass->GetXaxis()->SetTitle("missing mass /GeV");
  hist_e_p_X_mismass->GetYaxis()->SetTitle("counts");
  hist_e_p_X_mismass2 = new TH1F("hist_e_p_X_mismass2", "e p --> e p X missing mass squared", 800, -1, 4);   
  hist_e_p_X_mismass2->GetXaxis()->SetTitle("missing mass /GeV");
  hist_e_p_X_mismass2->GetYaxis()->SetTitle("counts");

  //ADDED
  TH1F *hist_phi_mass;
  TH1F *hist_epkX_mass;
  TH1F *hist_epKpKmX;
  TH1F *hist_epKpKmX2;
  TH2F *hist_kk_pk;

  TH1F *hist_eX_mass;
  TH1F *hist_epX_mass;

  TH2F *hist_epKpKmX_p;
  TH2F *hist_epKpKmX_theta;

  TH1F *hist_prot_mass_c1;
  TH1F *hist_prot_mass_c2;
  TH1F *hist_kp_mass_c1;
  TH1F *hist_kp_mass_c2;
  TH1F *hist_km_mass_c1;
  TH1F *hist_km_mass_c1_phicut;
  TH1F *hist_km_mass_c2;


  //missing mass
  TH1F *hist_epX_mm2_c1;
  TH1F *hist_epKpX_mm2_c1;
  TH1F *hist_eXKpKm_mm2_c1;
  TH1F *hist_eX_mm2_c1;

  TH1F *hist_epX_mm2_c2;
  TH1F *hist_epKpX_mm2_c2;
  TH1F *hist_epKpKmX_mm2_c2;
  TH1F *hist_eXKpKm_mm2_c2;
  TH1F *hist_eX_mm2_c2;

  //missing energy
  TH1F *hist_epX_me_c1;
  TH1F *hist_epKpX_me_c1;
  TH1F *hist_eXKpKm_me_c1;
  TH1F *hist_eX_me_c1;

  TH1F *hist_epX_me_c2;
  TH1F *hist_epKpX_me_c2;
  TH1F *hist_epKpKmX_me_c2;
  TH1F *hist_eXKpKm_me_c2;
  TH1F *hist_eX_me_c2;

  //missing p
  TH1F *hist_epX_mp_c1;
  TH1F *hist_epKpX_mp_c1;
  TH1F *hist_eXKpKm_mp_c1;
  TH1F *hist_eX_mp_c1;

  TH1F *hist_epX_mp_c2;
  TH1F *hist_epKpX_mp_c2;
  TH1F *hist_epKpKmX_mp_c2;
  TH1F *hist_eXKpKm_mp_c2;
  TH1F *hist_eX_mp_c2;

  //missing e and p with mm2 cut
  //missing energy
  TH1F *hist_epX_me_c1_cutmm2;
  TH1F *hist_epKpX_me_c1_cutmm2;
  TH1F *hist_eXKpKm_me_c1_cutmm2;
  TH1F *hist_eX_me_c1_cutmm2;

  TH1F *hist_epX_me_c2_cutmm2;
  TH1F *hist_epKpX_me_c2_cutmm2;
  TH1F *hist_epKpKmX_me_c2_cutmm2;
  TH1F *hist_eXKpKm_me_c2_cutmm2;
  TH1F *hist_eX_me_c2_cutmm2;

  //missing p
  TH1F *hist_epX_mp_c1_cutmm2;
  TH1F *hist_epKpX_mp_c1_cutmm2;
  TH1F *hist_eXKpKm_mp_c1_cutmm2;
  TH1F *hist_eX_mp_c1_cutmm2;

  TH1F *hist_epX_mp_c2_cutmm2;
  TH1F *hist_epKpX_mp_c2_cutmm2;
  TH1F *hist_epKpKmX_mp_c2_cutmm2;
  TH1F *hist_eXKpKm_mp_c2_cutmm2;
  TH1F *hist_eX_mp_c2_cutmm2;

  //missing mm2 and p with cut on the e
  TH1F *hist_epX_mm2_c1_cute;
  TH1F *hist_epKpX_mm2_c1_cute;
  TH1F *hist_eXKpKm_mm2_c1_cute;
  TH1F *hist_eX_mm2_c1_cute;

  TH1F *hist_epX_mm2_c2_cute;
  TH1F *hist_epKpX_mm2_c2_cute;
  TH1F *hist_epKpKmX_mm2_c2_cute;
  TH1F *hist_eXKpKm_mm2_c2_cute;
  TH1F *hist_eX_mm2_c2_cute;

  //missing p
  TH1F *hist_epX_mp_c1_cute;
  TH1F *hist_epKpX_mp_c1_cute;
  TH1F *hist_eXKpKm_mp_c1_cute;
  TH1F *hist_eX_mp_c1_cute;

  TH1F *hist_epX_mp_c2_cute;
  TH1F *hist_epKpX_mp_c2_cute;
  TH1F *hist_epKpKmX_mp_c2_cute;
  TH1F *hist_eXKpKm_mp_c2_cute;
  TH1F *hist_eX_mp_c2_cute;

  //missing mm2 and me for cut on mp
  //missing mass
  TH1F *hist_epX_mm2_c1_cutp;
  TH1F *hist_epKpX_mm2_c1_cutp;
  TH1F *hist_eXKpKm_mm2_c1_cutp;
  TH1F *hist_eX_mm2_c1_cutp;

  TH1F *hist_epX_mm2_c2_cutp;
  TH1F *hist_epKpX_mm2_c2_cutp;
  TH1F *hist_epKpKmX_mm2_c2_cutp;
  TH1F *hist_eXKpKm_mm2_c2_cutp;
  TH1F *hist_eX_mm2_c2_cutp;

  //missing energy
  TH1F *hist_epX_me_c1_cutp;
  TH1F *hist_epKpX_me_c1_cutp;
  TH1F *hist_eXKpKm_me_c1_cutp;
  TH1F *hist_eX_me_c1_cutp;

  TH1F *hist_epX_me_c2_cutp;
  TH1F *hist_epKpX_me_c2_cutp;
  TH1F *hist_epKpKmX_me_c2_cutp;
  TH1F *hist_eXKpKm_me_c2_cutp;
  TH1F *hist_eX_me_c2_cutp;  

  /////////////////
  TH1F *hist_epKpKmX_mm2_c1_cut2;
  TH1F *hist_epKpKmX_me_c1_cut2;
  TH1F *hist_epKpKmX_mp_c1_cut2;

  TH1F *hist_epKpKmX_mm2_c2_cut2;
  TH1F *hist_epKpKmX_me_c2_cut2;
  TH1F *hist_epKpKmX_mp_c2_cut2;


  //////////////
  TH1F *hist_phi_mass_c1;
  TH1F *hist_phi_mass_c2;
  TH1F *hist_phi_mass_c1_missing_cuts;
  TH1F *hist_phi_mass_c2_missing_cuts;
  TH1F *hist_phi_mass_c2_cutmm2;
  TH1F *hist_km_rec_mass;
  TH1F *hist_kp_kmrec_angle;
  TH1F *hist_kp_km_angle_c1;
  TH1F *hist_kp_km_angle_c2;
  TH1F *hist_epphi_mm2_c2;
  TH1F *hist_ekpkm_mm2_c2;
  TH1F *hist_ep_mm2_c2;

  TH1F *hist_epX_c1;
  TH1F *hist_epX_c2;
  TH1F *hist_epX_c1_kmcut;

  TH1F *hist_pr_phi_angle_c1;
  TH1F *hist_pr_phi_angle_c2;
  
  TH2F *hist_pr_phi_pangle_c1;
  TH2F *hist_pr_phi_pangle_c2;

  TH2F *hist_pr_phi_pangle_c1_mm2;
  TH2F *hist_pr_phi_pangle_c2_mm2;

  TH2F *hist_pr_phi_pangle_c1_me;
  TH2F *hist_pr_phi_pangle_c2_me;

  TH2F *hist_pr_phi_pangle_c1_mp;
  TH2F *hist_pr_phi_pangle_c2_mp;
  
  TH1F *hist_phi_mass_c2_cutkaonangles;
  
  hist_epphi_mm2_c2 = new TH1F("hist_epphi_mm2_c2","hist_epphi_mm2_c2",100,-10.0, 10.0);
  hist_ekpkm_mm2_c2 = new TH1F("hist_ekpkm_mm2_c2","hist_ekpkm_mm2_c2",100,-10.0, 10.0 );
  hist_ep_mm2_c2 = new TH1F("hist_ep_mm2_c2","hist_ep_mm2_c2",100,-10.0,5);
  hist_epX_c1 = new TH1F("hist_epX_c1","hist_epX_c1", 100, 0.0, 3.5);
  hist_epX_c2 = new TH1F("hist_epX_c2","hist_epX_c2", 60, 0.0, 3.5);
  hist_epX_c1_kmcut = new TH1F("hist_epX_c1_kmcut","hist_epX_c1_kmcut",100,0.0,3.5);
  hist_phi_mass_c2_cutkaonangles = new TH1F("hist_phi_mass_c2_cutkaonangles","hist_phi_mass_c2_cutkaonangles",100, 0.8, 1.3);

  //missing masss
  hist_epX_mm2_c1 = new TH1F("hist_epX_mm2_c1","hist_epX_mm2_c1",200,-10.0,10.0);
  hist_epKpX_mm2_c1 = new TH1F("hist_epKpX_mm2_c1","hist_epKpX_mm2_c1",200,-10.0,10.0);
  hist_eXKpKm_mm2_c1 = new TH1F("hist_eXKpKm_mm2_c1","hist_eXKpKm_mm2_c1",200,-10.0,10.0);
  hist_eX_mm2_c1 = new TH1F("hist_eX_mm2_c1","hist_eX_mm2_c1",200,-10.0,10.0);

  hist_epX_mm2_c2 = new TH1F("hist_epX_mm2_c2","hist_epX_mm2_c2",200,-10.0,10.0);
  hist_epKpX_mm2_c2 = new TH1F("hist_epKpX_mm2_c2","hist_epKpX_mm2_c2",200,-10.0,10.0);
  hist_epKpKmX_mm2_c2 = new TH1F("hist_epKpKmX_mm2_c2","hist_epKpKmX_mm2_c2",200,-10.0,10.0);
  hist_eXKpKm_mm2_c2 = new TH1F("hist_eXKpKm_mm2_c2","hist_eXKpKm_mm2_c2",200,-10.0,10.0);
  hist_eX_mm2_c2 = new TH1F("hist_eX_mm2_c2","hist_eX_mm2_c2",200,-10.0,10.0);

  //missing energy
  hist_epX_me_c1 = new TH1F("hist_epX_me_c1","hist_epX_me_c1",200,-10.0,10.0);
  hist_epKpX_me_c1 = new TH1F("hist_epKpX_me_c1","hist_epKpX_me_c1",200,-10.0,10.0);
  hist_eXKpKm_me_c1 = new TH1F("hist_eXKpKm_me_c1","hist_eXKpKm_me_c1",200,-10.0,10.0);
  hist_eX_me_c1 = new TH1F("hist_eX_me_c1","hist_eX_me_c1",200,-10.0,10.0);

  hist_epX_me_c2 = new TH1F("hist_epX_me_c2","hist_epX_me_c2",200,-10.0,10.0);
  hist_epKpX_me_c2 = new TH1F("hist_epKpX_me_c2","hist_epKpX_me_c2",200,-10.0,10.0);
  hist_epKpKmX_me_c2 = new TH1F("hist_epKpKmX_me_c2","hist_epKpKmX_me_c2",200,-10.0,10.0);
  hist_eXKpKm_me_c2 = new TH1F("hist_eXKpKm_me_c2","hist_eXKpKm_me_c2",200,-10.0,10.0);
  hist_eX_me_c2 = new TH1F("hist_eX_me_c2","hist_eX_me_c2",200,-10.0,10.0);

  //missing p
  hist_epX_mp_c1 = new TH1F("hist_epX_mp_c1","hist_epX_mp_c1",200,-10.0,10.0);
  hist_epKpX_mp_c1 = new TH1F("hist_epKpX_mp_c1","hist_epKpX_mp_c1",200,-10.0,10.0);
  hist_eXKpKm_mp_c1 = new TH1F("hist_eXKpKm_mp_c1","hist_eXKpKm_mp_c1",200,-10.0,10.0);
  hist_eX_mp_c1 = new TH1F("hist_eX_mp_c1","hist_eX_mp_c1",200,-10.0,10.0);

  hist_epX_mp_c2 = new TH1F("hist_epX_mp_c2","hist_epX_mp_c2",200,-10.0,10.0);
  hist_epKpX_mp_c2 = new TH1F("hist_epKpX_mp_c2","hist_epKpX_mp_c2",200,-10.0,10.0);
  hist_epKpKmX_mp_c2 = new TH1F("hist_epKpKmX_mp_c2","hist_epKpKmX_mp_c2",200,-10.0,10.0);
  hist_eXKpKm_mp_c2 = new TH1F("hist_eXKpKm_mp_c2","hist_eXKpKm_mp_c2",200,-10.0,10.0);
  hist_eX_mp_c2 = new TH1F("hist_eX_mp_c2","hist_eX_mp_c2",200,-10.0,10.0);
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //cut on mm2 plots
  //missing energy
  hist_epX_me_c1_cutmm2 = new TH1F("hist_epX_me_c1_cutmm2","hist_epX_me_c1_cutmm2",200,-10.0,10.0);
  hist_epKpX_me_c1_cutmm2 = new TH1F("hist_epKpX_me_c1_cutmm2","hist_epKpX_me_c1_cutmm2",200,-10.0,10.0);
  hist_eXKpKm_me_c1_cutmm2 = new TH1F("hist_eXKpKm_me_c1_cutmm2","hist_eXKpKm_me_c1_cutmm2",200,-10.0,10.0);
  hist_eX_me_c1_cutmm2 = new TH1F("hist_eX_me_c1_cutmm2","hist_eX_me_c1_cutmm2",200,-10.0,10.0);

  hist_epX_me_c2_cutmm2= new TH1F("hist_epX_me_c2_cutmm2","hist_epX_me_c2_cutmm2",200,-10.0,10.0);
  hist_epKpX_me_c2_cutmm2 = new TH1F("hist_epKpX_me_c2_cutmm2","hist_epKpX_me_c2_cutmm2",200,-10.0,10.0);
  hist_epKpKmX_me_c2_cutmm2 = new TH1F("hist_epKpKmX_me_c2_cutmm2","hist_epKpKmX_me_c2_cutmm2",200,-10.0,10.0);
  hist_eXKpKm_me_c2_cutmm2 = new TH1F("hist_eXKpKm_me_c2_cutmm2","hist_eXKpKm_me_c2_cutmm2",200,-10.0,10.0);
  hist_eX_me_c2_cutmm2 = new TH1F("hist_eX_me_c2_cutmm2","hist_eX_me_c2_cutmm2",200,-10.0,10.0);

  //missing p
  hist_epX_mp_c1_cutmm2 = new TH1F("hist_epX_mp_c1_cutmm2","hist_epX_mp_c1_cutmm2",200,-10.0,10.0);
  hist_epKpX_mp_c1_cutmm2 = new TH1F("hist_epKpX_mp_c1_cutmm2","hist_epKpX_mp_c1_cutmm2",200,-10.0,10.0);
  hist_eXKpKm_mp_c1_cutmm2 = new TH1F("hist_eXKpKm_mp_c1_cutmm2","hist_eXKpKm_mp_c1_cutmm2",200,-10.0,10.0);
  hist_eX_mp_c1_cutmm2 = new TH1F("hist_eX_mp_c1_cutmm2","hist_eX_mp_c1_cutmm2",200,-10.0,10.0);

  hist_epX_mp_c2_cutmm2 = new TH1F("hist_epX_mp_c2_cutmm2","hist_epX_mp_c2_cutmm2",200,-10.0,10.0);
  hist_epKpX_mp_c2_cutmm2 = new TH1F("hist_epKpX_mp_c2_cutmm2","hist_epKpX_mp_c2_cutmm2",200,-10.0,10.0);
  hist_epKpKmX_mp_c2_cutmm2 = new TH1F("hist_epKpKmX_mp_c2_cutmm2","hist_epKpKmX_mp_c2_cutmm2",200,-10.0,10.0);
  hist_eXKpKm_mp_c2_cutmm2 = new TH1F("hist_eXKpKm_mp_c2_cutmm2","hist_eXKpKm_mp_c2_cutmm2",200,-10.0,10.0);
  hist_eX_mp_c2_cutmm2 = new TH1F("hist_eX_mp_c2_cutmm2","hist_eX_mp_c2_cutmm2",200,-10.0,10.0);

  //cut on missing enegy 
  //missing masss
  hist_epX_mm2_c1_cute = new TH1F("hist_epX_mm2_c1_cute","hist_epX_mm2_c1_cute",200,-10.0,10.0);
  hist_epKpX_mm2_c1_cute = new TH1F("hist_epKpX_mm2_c1_cute","hist_epKpX_mm2_c1_cute",200,-10.0,10.0);
  hist_eXKpKm_mm2_c1_cute = new TH1F("hist_eXKpKm_mm2_c1_cute","hist_eXKpKm_mm2_c1_cute",200,-10.0,10.0);
  hist_eX_mm2_c1_cute = new TH1F("hist_eX_mm2_c1_cute","hist_eX_mm2_c1_cute",200,-10.0,10.0);

  hist_epX_mm2_c2_cute = new TH1F("hist_epX_mm2_c2_cute","hist_epX_mm2_c2_cute",200,-10.0,10.0);
  hist_epKpX_mm2_c2_cute = new TH1F("hist_epKpX_mm2_c2_cute","hist_epKpX_mm2_c2_cute",200,-10.0,10.0);
  hist_epKpKmX_mm2_c2_cute = new TH1F("hist_epKpKmX_mm2_c2_cute","hist_epKpKmX_mm2_c2_cute",200,-10.0,10.0);
  hist_eXKpKm_mm2_c2_cute = new TH1F("hist_eXKpKm_mm2_c2_cute","hist_eXKpKm_mm2_c2_cute",200,-10.0,10.0);
  hist_eX_mm2_c2_cute = new TH1F("hist_eX_mm2_c2_cute","hist_eX_mm2_c2_cute",200,-10.0,10.0);
  
  //missing p
  hist_epX_mp_c1_cute = new TH1F("hist_epX_mp_c1_cute","hist_epX_mp_c1_cute",200,-10.0,10.0);
  hist_epKpX_mp_c1_cute = new TH1F("hist_epKpX_mp_c1_cute","hist_epKpX_mp_c1_cute",200,-10.0,10.0);
  hist_eXKpKm_mp_c1_cute = new TH1F("hist_eXKpKm_mp_c1_cute","hist_eXKpKm_mp_c1_cute",200,-10.0,10.0);
  hist_eX_mp_c1_cute = new TH1F("hist_eX_mp_c1_cute","hist_eX_mp_c1_cute",200,-10.0,10.0);

  hist_epX_mp_c2_cute = new TH1F("hist_epX_mp_c2_cute","hist_epX_mp_c2_cute",200,-10.0,10.0);
  hist_epKpX_mp_c2_cute = new TH1F("hist_epKpX_mp_c2_cute","hist_epKpX_mp_c2_cute",200,-10.0,10.0);
  hist_epKpKmX_mp_c2_cute = new TH1F("hist_epKpKmX_mp_c2_cute","hist_epKpKmX_mp_c2_cute",200,-10.0,10.0);
  hist_eXKpKm_mp_c2_cute = new TH1F("hist_eXKpKm_mp_c2_cute","hist_eXKpKm_mp_c2_cute",200,-10.0,10.0);
  hist_eX_mp_c2_cute = new TH1F("hist_eX_mp_c2_cute","hist_eX_mp_c2_cute",200,-10.0,10.0);

  //cute on missing p
  //missing masss
  hist_epX_mm2_c1_cutp = new TH1F("hist_epX_mm2_c1_cutp","hist_epX_mm2_c1_cutp",200,-10.0,10.0);
  hist_epKpX_mm2_c1_cutp = new TH1F("hist_epKpX_mm2_c1_cutp","hist_epKpX_mm2_c1_cutp",200,-10.0,10.0);
  hist_eXKpKm_mm2_c1_cutp = new TH1F("hist_eXKpKm_mm2_c1_cutp","hist_eXKpKm_mm2_c1_cutp",200,-10.0,10.0);
  hist_eX_mm2_c1_cutp = new TH1F("hist_eX_mm2_c1_cutp","hist_eX_mm2_c1_cutp",200,-10.0,10.0);

  hist_epX_mm2_c2_cutp = new TH1F("hist_epX_mm2_c2_cutp","hist_epX_mm2_c2_cutp",200,-10.0,10.0);
  hist_epKpX_mm2_c2_cutp = new TH1F("hist_epKpX_mm2_c2_cutp","hist_epKpX_mm2_c2_cutp",200,-10.0,10.0);
  hist_epKpKmX_mm2_c2_cutp = new TH1F("hist_epKpKmX_mm2_c2_cutp","hist_epKpKmX_mm2_c2_cutp",200,-10.0,10.0);
  hist_eXKpKm_mm2_c2_cutp = new TH1F("hist_eXKpKm_mm2_c2_cutp","hist_eXKpKm_mm2_c2_cutp",200,-10.0,10.0);
  hist_eX_mm2_c2_cutp = new TH1F("hist_eX_mm2_c2_cutp","hist_eX_mm2_c2_cutp",200,-10.0,10.0);

  //missing energy
  hist_epX_me_c1_cutp = new TH1F("hist_epX_me_c1_cutp","hist_epX_me_c1_cutp",200,-10.0,10.0);
  hist_epKpX_me_c1_cutp = new TH1F("hist_epKpX_me_c1_cutp","hist_epKpX_me_c1_cutp",200,-10.0,10.0);
  hist_eXKpKm_me_c1_cutp = new TH1F("hist_eXKpKm_me_c1_cutp","hist_eXKpKm_me_c1_cutp",200,-10.0,10.0);
  hist_eX_me_c1_cutp = new TH1F("hist_eX_me_c1_cutp","hist_eX_me_c1_cutp",200,-10.0,10.0);

  hist_epX_me_c2_cutp = new TH1F("hist_epX_me_c2_cutp","hist_epX_me_c2_cutp",200,-10.0,10.0);
  hist_epKpX_me_c2_cutp = new TH1F("hist_epKpX_me_c2_cutp","hist_epKpX_me_c2_cutp",200,-10.0,10.0);
  hist_epKpKmX_me_c2_cutp = new TH1F("hist_epKpKmX_me_c2_cutp","hist_epKpKmX_me_c2_cutp",200,-10.0,10.0);
  hist_eXKpKm_me_c2_cutp = new TH1F("hist_eXKpKm_me_c2_cutp","hist_eXKpKm_me_c2_cutp",200,-10.0,10.0);
  hist_eX_me_c2_cutp = new TH1F("hist_eX_me_c2_cutp","hist_eX_me_c2_cutp",200,-10.0,10.0);

  // histograms with cuts on two of the missing variables and plotting the third
  hist_epKpKmX_mm2_c1_cut2 = new TH1F("hist_epKpKmX_mm2_c1_cut2","hist_epKpKmX_mm2_c1_cut2",100,-10.0,10.0);
  hist_epKpKmX_me_c1_cut2 = new TH1F("hist_epKpKmX_me_c1_cut2","hist_epKpKmX_me_c1_cut2",100,-10.0,10.0);
  hist_epKpKmX_mp_c1_cut2 = new TH1F("hist_epKpKmX_mp_c1_cut2","hist_epKpKmX_mp_c1_cut2",100,-10.0,10.0);

  hist_epKpKmX_mm2_c2_cut2 = new TH1F("hist_epKpKmX_mm2_c2_cut2","hist_epKpKmX_mm2_c2_cut2",100,-10.0,10.0);
  hist_epKpKmX_me_c2_cut2 = new TH1F("hist_epKpKmX_me_c2_cut2","hist_epKpKmX_me_c2_cut2",100,-10.0,10.0);
  hist_epKpKmX_mp_c2_cut2 = new TH1F("hist_epKpKmX_mp_c2_cut2","hist_epKpKmX_mp_c2_cut2",100,-10.0,10.0);

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  hist_pr_phi_angle_c1 = new TH1F("hist_pr_phi_angle_c1","hist_pr_phi_angle_c1", 200, 0.0, 180.0);
  hist_pr_phi_angle_c2 = new TH1F("hist_pr_phi_angle_c2","hist_pr_phi_angle_c2", 200, 0.0, 180.0);
  hist_pr_phi_pangle_c1 = new TH2F("hist_pr_phi_pvsangle_c1","hist_pr_phi_pvsangle_c1", 200, 0.0, 4.0, 200, 0.0, 180.0);
  hist_pr_phi_pangle_c2 = new TH2F("hist_pr_phi_pvsangle_c2","hist_pr_phi_pvsangle_c2", 200, 0.0, 4.0, 200, 0.0, 180.0);

  hist_pr_phi_pangle_c1_mm2 = new TH2F("hist_pr_phi_pvsangle_c1_mm2","hist_pr_phi_pvsangle_c1_mm2", 200, 0.0, 4.0, 200, 0.0, 180.0);
  hist_pr_phi_pangle_c2_mm2 = new TH2F("hist_pr_phi_pvsangle_c2_mm2","hist_pr_phi_pvsangle_c2_mm2", 200, 0.0, 4.0, 200, 0.0, 180.0);

  hist_pr_phi_pangle_c1_me = new TH2F("hist_pr_phi_pvsangle_c1_me","hist_pr_phi_pvsangle_c1_me", 200, 0.0, 4.0, 200, 0.0, 180.0);
  hist_pr_phi_pangle_c2_me = new TH2F("hist_pr_phi_pvsangle_c2_me","hist_pr_phi_pvsangle_c2_me", 200, 0.0, 4.0, 200, 0.0, 180.0);

  hist_pr_phi_pangle_c1_mp = new TH2F("hist_pr_phi_pvsangle_c1_mp","hist_pr_phi_pvsangle_c1_mp", 200, 0.0, 4.0, 200, 0.0, 180.0);
  hist_pr_phi_pangle_c2_mp = new TH2F("hist_pr_phi_pvsangle_c2_mp","hist_pr_phi_pvsangle_c2_mp", 200, 0.0, 4.0, 200, 0.0, 180.0);

  std::vector<TH1F*> h_eX_mass_cut;
  std::vector<TH1F*> h_epX_mass_cut;
  std::vector<TH1F*> h_epkX_mass_cut;
  std::vector<TH1F*> h_epKpKmX_mass_cut;
  std::vector<TH1F*> h_phi_mass_cut;
  std::vector<TH1F*> h_epKpKmX_missing_energy_cut;

  for( int c = 0; c < 2; c++ ){
    h_eX_mass_cut.push_back( new TH1F(Form("h_eX_mass_cutlvl_%d",c), Form("h_eX_mass_cutlvl_%d",c),2100, 0.0, 4.25));
    h_epX_mass_cut.push_back( new TH1F(Form("h_epX_mass_cutlvl_%d",c), Form("h_epX_mass_cutlvl_%d",c), 200, -10.50, 3.2));
    h_epkX_mass_cut.push_back( new TH1F(Form("h_epkX_mass_cutlvl_%d",c), Form("h_epkX_mass_cutlvl_%d",c), 200, -10.0, 3.0));
    h_epKpKmX_mass_cut.push_back( new TH1F(Form("h_epKpKmX_mass_cutlvl_%d",c), Form("h_epKpKmX_mass_cutlvl_%d",c), 100, -20.0, 20.0));
    h_epKpKmX_missing_energy_cut.push_back( new TH1F(Form("h_epKpKmX_missing_energy_cutlvl_%d",c), Form("h_epKpKmX_missing_energy_cutlvl_%d",c), 100, -10.0, 10.0));
    h_phi_mass_cut.push_back( new TH1F(Form("h_phi_mass_cutlvl_%d",c), Form("h_phi_mass_cutlvl_%d",c), 100, 0.90, 1.85));
  }  

  hist_phi_mass_c1 = new TH1F("hist_phi_mass_c1","hist_phi_mass_c1",200, 0.90, 1.85);
  hist_phi_mass_c2 = new TH1F("hist_phi_mass_c2","hist_phi_mass_c2",300, 0.90, 1.85);
  hist_phi_mass_c1_missing_cuts = new TH1F("hist_phi_mass_c1_missing_cuts","hist_phi_mass_c1_missing_cuts",100, 0.90, 1.85);
  hist_phi_mass_c2_missing_cuts = new TH1F("hist_phi_mass_c2_missing_cuts","hist_phi_mass_c2_missing_cuts",100, 0.90, 1.85);
  hist_phi_mass_c2_cutmm2  = new TH1F("hist_phi_mass_c2_cutmm2","hist_phi_mass_c2_cutmm2",100, 0.90, 1.85);
  hist_km_rec_mass = new TH1F("hist_km_rec_mass","hist_km_rec_mass",200, 0.0, 1.5);
  hist_kp_km_angle_c1 = new TH1F("hist_kp_km_angle_c1","hist_kp_km_angle_c1",100, 0.0, 90.0);
  hist_kp_km_angle_c2 = new TH1F("hist_kp_km_angle_c2","hist_kp_km_angle_c2",100, 0.0, 90.0);

  hist_prot_mass_c1 = new TH1F("hist_prot_mass_c1","hist_prot_mass_c1",100, 0.85, 1.2);
  hist_prot_mass_c2 = new TH1F("hist_prot_mass_c2","hist_prot_mass_c2",100, 0.85, 1.2);

  hist_kp_mass_c1 = new TH1F("hist_kp_mass_c1", "hist_kp_mass_c1", 100, 0.35, 0.85);
  hist_kp_mass_c2 = new TH1F("hist_kp_mass_c2", "hist_kp_mass_c2", 100, 0.35, 0.85);

  hist_km_mass_c1 = new TH1F("hist_km_mass_c1", "hist_km_mass_c1", 100, -1.0, 1.85);
  hist_km_mass_c1_phicut = new TH1F("hist_km_mass_c1_phicut", "hist_km_mass_c1_phicut", 100, -1.0, 1.85);
  hist_km_mass_c2 = new TH1F("hist_km_mass_c2", "hist_km_mass_c2", 100, 0.25, 0.85);
  

  //ADDED
  hist_phi_mass = new TH1F("hist_phi_mass","Mass of #phi Meson from K^{+}K^{-} Pair", 100, -1.0, 5.8);
  hist_phi_mass->GetXaxis()->SetTitle("Mass of K^{+}K^{-} [GeV]");
  hist_phi_mass->GetXaxis()->CenterTitle();

  hist_epkX_mass = new TH1F("hist_epkX_mass","Mass of K^{-} from #phi event", 100, -1.1, 10.6);
  hist_epkX_mass->GetXaxis()->SetTitle("Mass [GeV]");
  hist_epkX_mass->GetXaxis()->CenterTitle();
 
  hist_epKpKmX = new TH1F("hist_epKpKmX","Mass of epK^{+}K^{-}X",200, -25.5, 25.5);
  hist_epKpKmX->GetXaxis()->SetTitle("Mass^{2} [GeV^{2}]");
  hist_epKpKmX->GetXaxis()->CenterTitle();

  hist_epKpKmX2 = new TH1F("hist_epKpKmX2","Mass of epK^{+}K^{-}X",200, -25.5, 25.5);
  hist_epKpKmX2->GetXaxis()->SetTitle("Mass^{2} [GeV^{2}]");
  hist_epKpKmX2->GetXaxis()->CenterTitle();
  
  hist_kk_pk = new TH2F("hist_kk_pk","Invariant Mass of K^{+} K^{-} vs p K^{-}",50, 0.8, 1.5, 50, 1.39, 1.75);
  hist_kk_pk->GetXaxis()->SetTitle("Mass of K^{+} K^{-} [GeV]");
  hist_kk_pk->GetXaxis()->CenterTitle();
  hist_kk_pk->GetYaxis()->SetTitle("Mass of p K^{-} [GeV]");
  hist_kk_pk->GetXaxis()->CenterTitle();

  hist_epKpKmX_p = new TH2F("hist_epKpKmX_p","Missing Mass^{2} vs proton p", 100, 0.0, Ebeam, 100, -25.0, 25);  
  hist_epKpKmX_p->GetXaxis()->SetTitle("p [GeV]");
  hist_epKpKmX_p->GetYaxis()->SetTitle("Mass^{2} [GeV^{2}]");

  hist_epKpKmX_theta = new TH2F("hist_epKpKmX_theta","Missing Mass^{2} vs proton #theta", 100, 0.0, 120, 100, -25.0, 25);  
  hist_epKpKmX_theta->GetXaxis()->SetTitle("#theta [#deg]");
  hist_epKpKmX_theta->GetYaxis()->SetTitle("Mass^{2} [GeV^{2}]");

  hist_eX_mass = new TH1F("hist_eX_mass","hist_eX_mass",100, 0.0, 4.0);
  hist_eX_mass->GetXaxis()->SetTitle("Missing Mass of eX");
  hist_eX_mass->GetXaxis()->CenterTitle();

  hist_epX_mass = new TH1F("hist_epX_mass","hist_epX_mass",100, 0.0, 3.0);
  hist_epX_mass->GetXaxis()->SetTitle("Missing Mass of epX");
  hist_epX_mass->GetXaxis()->CenterTitle();

  out->mkdir("statistics");				
  out->cd ("statistics");

  TH1F *hist_electron_count;
  TH1F *hist_proton_count;
  TH1F *hist_neutron_count;
  TH1F *hist_pip_count;
  TH1F *hist_pim_count;
  TH1F *hist_Kp_count;
  TH1F *hist_Km_count;
  TH1F *hist_photon_count;

  hist_electron_count = new TH1F("hist_electron_count", "electron count per event", 6, -0.5, 5.5);   
  hist_electron_count->GetXaxis()->SetTitle("number of electrons per event");
  hist_electron_count->GetYaxis()->SetTitle("counts");
  hist_proton_count = new TH1F("hist_proton_count", "proton count per event", 11, -0.5, 10.5);   
  hist_proton_count->GetXaxis()->SetTitle("number of protons per event");
  hist_proton_count->GetYaxis()->SetTitle("counts");
  hist_neutron_count = new TH1F("hist_neutron_count", "neutron count per event", 11, -0.5, 10.5);   
  hist_neutron_count->GetXaxis()->SetTitle("number of neutrons per event");
  hist_neutron_count->GetYaxis()->SetTitle("counts");
  hist_pip_count = new TH1F("hist_pip_count", "pip count per event", 11, -0.5, 10.5);   
  hist_pip_count->GetXaxis()->SetTitle("number of pips per event");
  hist_pip_count->GetYaxis()->SetTitle("counts");
  hist_pim_count = new TH1F("hist_pim_count", "pim count per event", 11, -0.5, 10.5);   
  hist_pim_count->GetXaxis()->SetTitle("number of pims per event");
  hist_pim_count->GetYaxis()->SetTitle("counts");
  hist_Kp_count = new TH1F("hist_Kp_count", "Kp count per event", 11, -0.5, 10.5);   
  hist_Kp_count->GetXaxis()->SetTitle("number of Kps per event");
  hist_Kp_count->GetYaxis()->SetTitle("counts");
  hist_Km_count = new TH1F("hist_Km_count", "Km count per event", 11, -0.5, 10.5);   
  hist_Km_count->GetXaxis()->SetTitle("number of Kms per event");
  hist_Km_count->GetYaxis()->SetTitle("counts");
  hist_photon_count = new TH1F("hist_photon_count", "photon count per event", 11, -0.5, 10.5);   
  hist_photon_count->GetXaxis()->SetTitle("number of photons per event");
  hist_photon_count->GetYaxis()->SetTitle("counts");


  /// ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///  start of the event loop     ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  cout << "Analysing Tree: " << inTree << endl;
  cout << "Event Loop starting ... " << endl;
  cout << endl;
  
  int channel_counter = 0;
  for(Int_t k=0; k < anaTree->GetEntriesFast();k++){    

    anaTree->GetEntry(k);

    if(process_Events > 0 && k == process_Events) break;

    /// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// progress:

    if(k % 100000 == 0){

      double events = anaTree->GetEntriesFast();
      double percent = k/(events/100);
      
      printf("Analysing event number %i of %.00f (%.01f percent)\n", k, events, percent);
    }

    //pr_sect=sect_pr[select_prot_1];//
    //kp_sect=sect_kp[select_kaonP_1];//
    //km_sect=sect_kp[select_kaonM_1];

    //std::cout << " pr sector " << pr_sect << std::endl;
	//	std::cout << " kp sector " << kp_sect << std::endl;
    //	std::cout << " km sector " << km_sect << std::endl;


    /// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// assign number of particles of each type PER EVENT:
    // (SIZE OF ARRAYS)
    ele_count = p4_ele_px->size();
    prot_count = p4_prot_px->size();
    neutr_count = p4_neutr_px->size();
    pip_count = p4_pip_px->size();
    pim_count = p4_pim_px->size();
    Kp_count = p4_Kp_px->size();
    Km_count = p4_Km_px->size();
    phot_count = p4_phot_px->size();

    //if ( ele_count==1 && prot_count > 0  && Kp_count > 0 ) std::cout << " kaon present "  << " number of protns " << p4_prot_px->size() << std::endl; 
    

    /// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// initalisatize variables:

    TLorentzVector beam(0,0,Ebeam,Ebeam);
    TLorentzVector target(0,0,0,0.93827);

    ///  index of selected particle for the different categories:
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
    
    /// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //  Assign tree components to Lorentz vectors:
    /// //////////////////////////////////////////////////////////////////

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


    // Assign beta from tree to array
    for( Int_t i = 0; i < BUFFER; i++ ){
      //if( prot_count > i ){ prot_beta[i]=prot_beta_final->at(i); }
      //if( Kp_count > i ){ kaonP_beta[i]=Kp_beta_final->at(i); }
      //if( Km_count > i ){ kaonM_beta[i]=Km_beta_final->at(i); }

      // and assign DC sector for hadrons from tree to array
      //if( prot_count > i ){ sect_pr[i]=prot_sector->at(i); }
      //if( Kp_count > i ){ sect_kp[i]=Kp_sector->at(i); }
      //if( Km_count > i ){ sect_km[i]=Km_sector->at(i); }

    }


    /// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //  Build event:
    /// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    if(ele_count == 1){    // one detected electron is the basic trigger condition for all topologies 

      select_ele = 0;   // if more than one electron is detected (< 0.2 percent of the events), take the one with the hihest momentum
          
      W  = kin_W(ele[select_ele]);
      Q2 = kin_Q2(ele[select_ele]);
      x  = kin_x(ele[select_ele]);
      y  = kin_y(ele[select_ele]);
      nu = kin_nu(ele[select_ele]);
      
      E_ele  = ele[select_ele].E();
      px_ele = ele[select_ele].Px();
      py_ele = ele[select_ele].Py();
      pz_ele = ele[select_ele].Pz();
            
      
      hist_helicity->Fill(helicity);
      hist_faraday_cup->Fill(fcup);
      
      // fill kinematics:
      
      if(W > 0) hist_W->Fill(W);
      if(Q2 > 0) hist_Q2->Fill(Q2);
      //if(W > 0 && Q2 > 0) hist_Q2_vs_W->Fill(W, Q2);
      if(ele[select_ele].Phi() != 0 && W > 0) hist_W_vs_phi->Fill(ele[select_ele].Phi()*180/Pival, W);
      if(x > 0) hist_x->Fill(x);
      if(y > 0) hist_y->Fill(y);
      if(nu > 0) hist_nu->Fill(nu);
      

      if( Q2 > 0 && x > 0 )h_q2_x_cut[0]->Fill(x,Q2);
      if( Q2 > 0 && W > 0 )h_q2_w_cut[0]->Fill(W,Q2);

      if( Q2 > 1 && x > 0 )h_q2_x_cut[1]->Fill(x,Q2);
      if( Q2 > 1 && x > 0 && W > 2 )h_q2_x_cut[2]->Fill(x,Q2);
      if( Q2 > 1 && x > 0 && W > 2 )h_q2_w_cut[1]->Fill(W,Q2);

      if(W > 2 && Q2 > 1) hist_Q2_vs_W->Fill(W, Q2);
      if( Q2 > 1 && x > 0 && W > 2 ) hist_Q2_x->Fill(x,Q2);

      if( Q2 > 1 && W > 0 ) h_w_cut[0]->Fill(W);
      if( Q2 > 1 ) h_q2_cut[0]->Fill(Q2);
      if( Q2 > 1 && x > 0 ) h_x_cut[0]->Fill(x);
      if( Q2 > 1 && y > 0 ) h_y_cut[0]->Fill(y);
      if( Q2 > 1 && nu > 0 ) h_nu_cut[0]->Fill(nu);

      if( W > 2 ) h_w_cut[1]->Fill(W);
      if( W > 2 && Q2 > 0 ) h_q2_cut[1]->Fill(Q2);
      if( W > 2 && x > 0 ) h_x_cut[1]->Fill(x);
      if( W > 2 && y > 0 ) h_y_cut[1]->Fill(y);
      if( W > 2 && nu > 0 ) h_nu_cut[1]->Fill(nu);

      TLorentzVector eX_missing_mass = beam + target - ele[select_ele];
      h_eX_mass_cut[0]->Fill(eX_missing_mass.M());
      if ( Q2 > 1 && W >2 ) h_eX_mass_cut[1]->Fill( eX_missing_mass.M() );

      ///////////////////////////////////////////////////////////
      // category 1: e p X
      
      if(prot_count == 1){	
	evcat = 1;
	select_prot_1 = 0; // FASTEST PROTON AT INDEX 0, SAME FOR OTHER PARTICLES
      
	if(prot_count > select_prot_1){	  
	  
	  TLorentzVector X1 = miss_X(ele[select_ele], prot[select_prot_1]);
	  
	  t1 = -kin_t(ele[select_ele], X1);
	  cmphi1 = kin_cmphi(ele[select_ele], prot[select_prot_1]);
	  cmcostheta1 = kin_cmcostheta(ele[select_ele], prot[select_prot_1]); 
	  pt1 = kin_pT(ele[select_ele], prot[select_prot_1]);
	  eta1 = kin_eta(ele[select_ele], prot[select_prot_1]);
	  z1 = kin_z(ele[select_ele], prot[select_prot_1]); 
	  
	  M_e_p_X_miss = kin_mismass(ele[select_ele], prot[select_prot_1]);
	  M_e_p_X_miss2 = kin_mismass2(ele[select_ele], prot[select_prot_1]);
	  
	  E_prot_1 = prot[select_prot_1].E();
	  px_prot_1 = prot[select_prot_1].Px();
	  py_prot_1 = prot[select_prot_1].Py();
	  pz_prot_1 = prot[select_prot_1].Pz();

	  if( x > 0 && t1 > 0 ) hist_x_t->Fill(t1,x);
	  //if( cmphi1 > 0 && t1 > 0 ) hist_x_t->Fill(cmphi1,t1);
	  if( cmphi1 > 0 && Q2 > 0) hist_Q2_phi->Fill(cmphi1*180/Pival,Q2);
	  
	  
	  if( Q2 > 1 && W > 0 ) h_w_cut[2]->Fill(W);
	  if( Q2 > 1 ) h_q2_cut[2]->Fill(Q2);
	  if( Q2 > 1 && x > 0 ) h_x_cut[2]->Fill(x);
	  if( Q2 > 1 && y > 0 ) h_y_cut[2]->Fill(y);
	  if( Q2 > 1 && nu > 0 ) h_nu_cut[2]->Fill(nu);
	  if( Q2 > 1 && t1 > 0 ) h_t_cut[0]->Fill(t1);
	  if( Q2 > 1 && cmphi1 != 0 ) h_cmphi_cut[0]->Fill(cmphi1*180/Pival);
	  if( Q2 > 1 && cmcostheta1 != 0 ) h_cmcostheta_cut[0]->Fill(cmcostheta1);
	  if( Q2 > 1 && W > 0 ) h_q2_t_cut[0]->Fill(Q2,t1);

	  
	  if( W > 2 ) h_w_cut[3]->Fill(W);
	  if( W > 2 && Q2 > 0 ) h_q2_cut[3]->Fill(Q2);
	  if( W > 2 && x > 0 ) h_x_cut[3]->Fill(x);
	  if( W > 2 && y > 0 ) h_y_cut[3]->Fill(y);
	  if( W > 2 && nu > 0 ) h_nu_cut[3]->Fill(nu);
	  if( W > 2 && t1 > 0 ) h_t_cut[1]->Fill(t1);
	  if( W > 2 && cmphi1 != 0 ) h_cmphi_cut[1]->Fill(cmphi1*180/Pival);
	  if( W > 2 && cmcostheta1 != 0 ) h_cmcostheta_cut[1]->Fill(cmcostheta1);
	  if( W > 2 && t1 > 0 ) h_q2_t_cut[1]->Fill(Q2,t1);

	  
	  if(t1 > 0) hist_t->Fill(t1);
	  if(cmphi1 != 0) hist_cmphi->Fill(cmphi1*180/Pival);
	  if(cmcostheta1 != 0) hist_cmcostheta->Fill(cmcostheta1);
	  
	  h_epX_mass_cut[0]->Fill(X1.M());
	  if ( Q2 > 1 && W > 2 ) h_epX_mass_cut[1]->Fill(X1.M());
	  if( Q2 > 1 && W > 2 ) h_q2_t_cut[2]->Fill(Q2,t1);


	  if( Q2 > 0 && x > 0 )h_q2_x_cut[3]->Fill(x,Q2);
	  if( Q2 > 0 && W > 0 )h_q2_w_cut[2]->Fill(W,Q2);
	  
	  if( Q2 > 1 && x > 0 )h_q2_x_cut[4]->Fill(x,Q2);
	  if( Q2 > 1 && x > 0 && W > 2 )h_q2_x_cut[5]->Fill(x,Q2);
	  if( Q2 > 1 && W > 2 )h_q2_w_cut[3]->Fill(W,Q2);

	  if( x > 0 && t1 > 0 ) h_x_t_cut[0]->Fill(x,t1);
	  if( x > 0 && t1 > 0 && Q2 > 1 ) h_x_t_cut[1]->Fill(x,t1);
	  if( x > 0 && t1 > 0 && W > 2 ) h_x_t_cut[2]->Fill(x,t1);
	  if( x > 0 && t1 > 0 && Q2 > 1 && W > 2 ) h_x_t_cut[3]->Fill(x,t1);

	  if( cmphi1 != 0 && t1 > 0 ) h_t_phi_cut[0]->Fill(cmphi1*180/Pival,t1);
	  if( cmphi1 != 0 && t1 > 0 && Q2 > 1 ) h_t_phi_cut[1]->Fill(cmphi1*180/Pival,t1);
	  if( cmphi1 != 0 && t1 > 0 && W > 2 ) h_t_phi_cut[2]->Fill(cmphi1*180/Pival,t1);
	  if( cmphi1 != 0 && t1 > 0 && Q2 > 1 && W > 2 ) h_t_phi_cut[3]->Fill(cmphi1*180/Pival,t1);

	  if( Q2 > 0 && cmphi1 != 0 ) h_q2_phi_cut[0]->Fill(cmphi1*180/Pival, Q2);
	  if( cmphi1 != 0 && Q2 > 1 ) h_q2_phi_cut[1]->Fill(cmphi1*180/Pival, Q2);
	  if( cmphi1 != 0 && Q2 > 0 && W > 2 ) h_q2_phi_cut[2]->Fill(cmphi1*180/Pival, Q2);
	  if( cmphi1 != 0 && Q2 > 1 && W > 2 ) h_q2_phi_cut[3]->Fill(cmphi1*180/Pival, Q2);

	  if( W > 0 && cmphi1 != 0 ) h_q2_phi_cut[0]->Fill(cmphi1*180/Pival, Q2);
	  if( cmphi1 != 0 && W > 2 ) h_q2_phi_cut[1]->Fill(cmphi1*180/Pival, Q2);
	  if( cmphi1 != 0 && W > 0 && Q2 > 1 ) h_q2_phi_cut[2]->Fill(cmphi1*180/Pival, Q2);
	  if( cmphi1 != 0 && W > 2 && Q2 > 1 ) h_q2_phi_cut[3]->Fill(cmphi1*180/Pival, Q2);
	 
	}
      }
    
    
  
     
      /////////////////////////////////////////////////////////
      // category 2: e p Kp X      
      if(ele_count == 1){    // one detected electron is the basic trigger condition for all topologies 
	select_ele = 0;   // if more than one electron is detected (< 0.2 percent of the events), take the one with the hihest momentum	
	///std::cout << " Electron present " << std::endl;
	
	TLorentzVector lv_pr_max(0,0,0,0);
	for( int ii = 0; ii < p4_prot_px->size(); ii++ ){
	  TLorentzVector lv_temp = prot[ii];
	  if( lv_temp.E() > lv_pr_max.E() ){
	    lv_pr_max = lv_temp;
	  }
	}
	
	

      
	if(prot_count > 0){		  
	  select_prot_1 = 0;
	  if( prot_count > select_prot_1 ){
	    //std::cout << " Proton present " << std::endl;
 	    //std::cout << " Kp count " << Kp_count << std::endl;
	    //TLorentzVector lv_kp = Kp[select_kaonP_1];
	    //std::cout << " kaon plus " << lv_kp.M() << std::endl;

	    if( Kp_count > 0 ){
	      select_kaonP_1 = 0;
	      if( Kp_count > select_kaonP_1 ){
		//std::cout << " Kp  present " << std::endl;
		
		TLorentzVector lv_el = ele[select_ele];
		TLorentzVector lv_pr = lv_pr_max; //prot[select_prot_1];
		TLorentzVector lv_kp = Kp[select_kaonP_1];
		
		double epkX_mass = kin_epkXMass(ele[select_ele], prot[select_prot_1], Kp[select_kaonP_1]);
		hist_epkX_mass->Fill(epkX_mass);
		
		TLorentzVector X1 = miss_X(ele[select_ele], prot[select_prot_1]);
		
		t1 = -kin_t(ele[select_ele], X1);
		cmphi1 = kin_cmphi(ele[select_ele], prot[select_prot_1]);
		cmcostheta1 = kin_cmcostheta(ele[select_ele], prot[select_prot_1]); 
		pt1 = kin_pT(ele[select_ele], prot[select_prot_1]);
		eta1 = kin_eta(ele[select_ele], prot[select_prot_1]);
		z1 = kin_z(ele[select_ele], prot[select_prot_1]); 
		
		M_e_p_X_miss = kin_mismass(ele[select_ele], prot[select_prot_1]);
		M_e_p_X_miss2 = kin_mismass2(ele[select_ele], prot[select_prot_1]);


		if( Q2 > 1 && W > 0 ) h_w_cut[4]->Fill(W);
		if( Q2 > 1 ) h_q2_cut[4]->Fill(Q2);
		if( Q2 > 1 && x > 0 ) h_x_cut[4]->Fill(x);
		if( Q2 > 1 && y > 0 ) h_y_cut[4]->Fill(y);
		if( Q2 > 1 && nu > 0 ) h_nu_cut[4]->Fill(nu);
		if( Q2 > 1 && t1 > 0 ) h_t_cut[2]->Fill(t1);
		if( Q2 > 1 && cmphi1 != 0 ) h_cmphi_cut[2]->Fill(cmphi1*180/Pival);
		if( Q2 > 1 && cmcostheta1 != 0 ) h_cmcostheta_cut[2]->Fill(cmcostheta1);
		if( Q2 > 1 && W > 0 ) h_q2_t_cut[3]->Fill(Q2,t1);
	  
	  
		if( W > 2 ) h_w_cut[5]->Fill(W);
		if( W > 2 && Q2 > 0 ) h_q2_cut[5]->Fill(Q2);
		if( W > 2 && x > 0 ) h_x_cut[5]->Fill(x);
		if( W > 2 && y > 0 ) h_y_cut[5]->Fill(y);
		if( W > 2 && nu > 0 ) h_nu_cut[5]->Fill(nu);
		if( W > 2 && t1 > 0 ) h_t_cut[3]->Fill(t1);
		if( W > 2 && cmphi1 != 0 ) h_cmphi_cut[3]->Fill(cmphi1*180/Pival);
		if( W > 2 && cmcostheta1 != 0 ) h_cmcostheta_cut[3]->Fill(cmcostheta1);	 
		if( W > 2 && Q2 > 0 && t1 > 0) h_q2_t_cut[4]->Fill(Q2,t1);


		h_epkX_mass_cut[0]->Fill(epkX_mass);
		if( Q2 > 1 && W > 2 ) h_epkX_mass_cut[1]->Fill(epkX_mass);
		if( W > 2 && Q2 > 1 && t1 > 0) h_q2_t_cut[5]->Fill(Q2,t1);

		if( Q2 > 0 && x > 0 )h_q2_x_cut[6]->Fill(x,Q2);
		if( Q2 > 0 && W > 0 )h_q2_w_cut[4]->Fill(W,Q2);
	  
		if( Q2 > 1 && x > 0 )h_q2_x_cut[7]->Fill(x,Q2);
		if( Q2 > 1 && x > 0 && W > 2 )h_q2_x_cut[8]->Fill(x,Q2);
		if( Q2 > 1 && W > 2 )h_q2_w_cut[5]->Fill(W,Q2);

		if( x > 0 && t1 > 0 ) h_x_t_cut[4]->Fill(x,t1);
		if( x > 0 && t1 > 0 && Q2 > 1 ) h_x_t_cut[5]->Fill(x,t1);
		if( x > 0 && t1 > 0 && W > 2 ) h_x_t_cut[6]->Fill(x,t1);
		if( x > 0 && t1 > 0 && Q2 > 1 && W > 2 ) h_x_t_cut[7]->Fill(x,t1);

		if( cmphi1 != 0 && t1 > 0 ) h_t_phi_cut[4]->Fill(cmphi1*180/Pival,t1);
		if( cmphi1 != 0 && t1 > 0 && Q2 > 1 ) h_t_phi_cut[5]->Fill(cmphi1*180/Pival,t1);
		if( cmphi1 != 0 && t1 > 0 && W > 2 ) h_t_phi_cut[6]->Fill(cmphi1*180/Pival,t1);
		if( cmphi1 != 0 && t1 > 0 && Q2 > 1 && W > 2 ) h_t_phi_cut[7]->Fill(cmphi1*180/Pival,t1);

		if( Q2 > 0 && cmphi1 != 0 ) h_q2_phi_cut[4]->Fill(cmphi1*180/Pival, Q2);
		if( cmphi1 != 0 && Q2 > 1 ) h_q2_phi_cut[5]->Fill(cmphi1*180/Pival, Q2);
		if( cmphi1 != 0 && Q2 > 0 && W > 2 ) h_q2_phi_cut[6]->Fill(cmphi1*180/Pival, Q2);
		if( cmphi1 != 0 && Q2 > 1 && W > 2 ) h_q2_phi_cut[7]->Fill(cmphi1*180/Pival, Q2);

		if( W > 0 && cmphi1 != 0 ) h_q2_phi_cut[4]->Fill(cmphi1*180/Pival, Q2);
		if( cmphi1 != 0 && W > 2 ) h_q2_phi_cut[5]->Fill(cmphi1*180/Pival, Q2);
		if( cmphi1 != 0 && W > 0 && Q2 > 1 ) h_q2_phi_cut[6]->Fill(cmphi1*180/Pival, Q2);
		if( cmphi1 != 0 && W > 2 && Q2 > 1 ) h_q2_phi_cut[7]->Fill(cmphi1*180/Pival, Q2);

		

	      }
	    }
	  }
	}
      }
    
  
    
      /////////////////////////////////////////////////////////
      // category 3: e p Kp Km X      
      if(ele_count > 0 ){    // one detected electron is the basic trigger condition for all topologies 
	select_ele = 0;   // if more than one electron is detected (< 0.2 percent of the events), take the one with the hihest momentum	
	//std::cout << " Electron present " << std::endl;
	if(prot_count > 0){		  
	  select_prot_1 = 0;
	  
	  TLorentzVector lv_pr_max(0,0,0,0);
	  for( int ii = 0; ii < p4_prot_px->size(); ii++ ){
	    TLorentzVector lv_temp = prot[ii];
	    if( lv_temp.E() > lv_pr_max.E() ){
	      lv_pr_max = lv_temp;
	    }
	  }


	  if( prot_count > select_prot_1 ){
	    ///  std::cout << " Proton present " << std::endl;
	    if( Kp_count > 0 ){
	      select_kaonP_1 = 0;
	      if( Kp_count > select_kaonP_1 ){
		//std::cout << " Kp  present " << std::endl;
		//std::cout << " Km  count  " << Km_count << std::endl;		
		if( Km_count > 0 ){
		  select_kaonM_1 = 0;
		  if( Km_count > select_kaonM_1 ){
		    //std::cout << " here " << std::endl;
		    TLorentzVector lv_el = ele[select_ele];
		    TLorentzVector lv_pr = lv_pr_max;// prot[select_prot_1];
		    TLorentzVector lv_kp = Kp[select_kaonP_1];
		    TLorentzVector lv_km = Km[select_kaonM_1];
		    		    
		    //std::cout << " " << lv_el.M() << " " << lv_pr.M() << " " << lv_kp.M() << " " << lv_km.M() << std::endl;
		    /*std::cout << " phi event " << " el momentum " << lv_el.P() << std::endl;
		    std::cout << "           " << " pr momentum " << lv_pr.P() << std::endl;
		    std::cout << "           " << " kp momentum " << lv_kp.P() << std::endl;
		    std::cout << "           " << " km momentum " << lv_km.P() << std::endl;
		    */

		    TLorentzVector X1 = miss_X(ele[select_ele], prot[select_prot_1]);
		    
		    TLorentzVector epKpKm_missing_mass = beam + target - lv_el - lv_pr - lv_kp - lv_km;
		    h_epKpKmX_mass_cut[0]->Fill(epKpKm_missing_mass.M());
		    h_epKpKmX_missing_energy_cut[0]->Fill(epKpKm_missing_mass.E());
		    if( Q2 > 1 && W > 2 ) h_epKpKmX_mass_cut[1]->Fill(epKpKm_missing_mass.M());
		    if( Q2 > 1 && W > 2 ) h_epKpKmX_missing_energy_cut[1]->Fill(epKpKm_missing_mass.E());
	    
		    t1 = -kin_t(ele[select_ele], (lv_kp+lv_km));
		    cmphi1 = kin_cmphi(ele[select_ele], prot[select_prot_1]);
		    cmcostheta1 = kin_cmcostheta(ele[select_ele], prot[select_prot_1]); 
		    pt1 = kin_pT(ele[select_ele], prot[select_prot_1]);
		    eta1 = kin_eta(ele[select_ele], prot[select_prot_1]);
		    z1 = kin_z(ele[select_ele], prot[select_prot_1]); 
		    
		    M_e_p_X_miss = kin_mismass(ele[select_ele], prot[select_prot_1]);
		    M_e_p_X_miss2 = kin_mismass2(ele[select_ele], prot[select_prot_1]);
		    		    
		    /// phi meson histograms
		    double phi_mass = Phi_Mass( Kp[select_kaonP_1], Km[select_kaonM_1] );
		    //std::cout << " phi mass  " << phi_mass << std::endl;
		    h_mass_phi->Fill(phi_mass);
		    h_phi_mass_cut[0]->Fill(phi_mass);
		    if( Q2 > 1 && W > 2 ) h_phi_mass_cut[1]->Fill(phi_mass);		    
		  }
		}
	      }
	    }
	  }
	}
      }
    }
      
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    // Phi Analysis Here
    // Require an electron
    //      if ( ele_count > 0 ){
    //	W  = kin_W(ele[select_ele]);
    //	Q2 = kin_Q2(ele[select_ele]);
    //	x  = kin_x(ele[select_ele]);
    //	y  = kin_y(ele[select_ele]);
    //	nu = kin_nu(ele[select_ele]);
	
    //E_ele  = ele[select_ele].E();
    //px_ele = ele[select_ele].Px();
    //py_ele = ele[select_ele].Py();
    //pz_ele = ele[select_ele].Pz();
	
    //Apply Exclusivity Cuts
	
    channel = 0;
    if ( Q2 > 0 && W > 0 ){ //change to (true) for only proton analysis
      TLorentzVector lv_el = ele[select_ele];
      
      TLorentzVector lv_pr_max(0,0,0,0);
      for( int ii = 0; ii < p4_prot_px->size(); ii++ ){
	TLorentzVector lv_temp = prot[ii];
	if( lv_temp.E() > lv_pr_max.E() ){
	  lv_pr_max = lv_temp;
	}
      }
      
      // Check the ep -> epKpX channel
      if ( prot_count > 0 && Kp_count > 0 ){

	TLorentzVector lv_kp_max(0,0,0,0);
	for( int ii = 0; ii < Kp_count; ii++ ){
	  TLorentzVector lv_temp = Kp[ii];
	  if( lv_temp.E() > lv_kp_max.E() ){
	    lv_kp_max = lv_temp;
	  }
	}

	
	select_prot_1 = 0; select_kaonP_1 = 0;
	TLorentzVector lv_pr_c1 = lv_pr_max;//prot[select_prot_1];
	TLorentzVector lv_kp_c1 = lv_kp_max;//Kp[select_kaonP_1];
	TLorentzVector lv_km_c1 = beam + target - lv_el - lv_pr_c1 - lv_kp_c1;
	TLorentzVector lv_eXKpKm_c1 = beam + target - lv_el - lv_kp_c1 - lv_km_c1;
	TLorentzVector lv_eX_c1 = beam + target - lv_el;

	double km_c1_mass = lv_km_c1.M();

	hist_km_rec_mass->Fill(km_c1_mass);

	TLorentzVector lv_epX_c1 = beam + target - lv_el - lv_pr_c1;	    
	hist_epX_c1->Fill(lv_epX_c1.M());

	// cut on the epX mass around phi mass about +- 0.4 on each side as the data is still not properly calibrated
	if( lv_epX_c1.M() >= (1.018 - 0.1) && lv_epX_c1.M() <= (1.018 + 0.1) ){
	  hist_km_mass_c1_phicut->Fill(lv_km_c1.M());
	}
	    
	// cut on the km mass (smaller acceptance for field) to then see what the epX spectrum looks like. Choose +- 0.15
	if( lv_km_c1.M() >= (0.493 - 0.1) && lv_km_c1.M() <= (0.493 + 0.1) ){
	  hist_epX_c1_kmcut->Fill(lv_epX_c1.M());
	}

	// check angle between the Kp and recKm
	double angle_kp_kmrec = (lv_kp_c1.Vect()).Angle(lv_km_c1.Vect()) * 180.0/TMath::Pi();
	hist_kp_km_angle_c1->Fill(angle_kp_kmrec);
	    
	// check inv mass of phi candidate
	double inv_mass_kp_kmrec = (lv_kp_c1 + lv_km_c1).M();
	hist_phi_mass_c1->Fill(inv_mass_kp_kmrec);
	    
	hist_prot_mass_c1->Fill(lv_pr_c1.M());
	hist_kp_mass_c1->Fill(lv_kp_c1.M());
	    
	TLorentzVector lv_phi_c1 = lv_kp_c1 + lv_km_c1;

	hist_pr_phi_angle_c1->Fill( (lv_phi_c1.Vect()).Angle(lv_pr_c1.Vect()) * 180.0/TMath::Pi() );
	hist_pr_phi_pangle_c1->Fill( lv_pr_c1.P(), (lv_phi_c1.Vect()).Angle(lv_pr_c1.Vect()) * 180.0/TMath::Pi() );
	TLorentzVector epKpKm_missing_mass = beam + target - lv_el - lv_pr_c1 - lv_kp_c1 - lv_km_c1;

	//from fits in channel 2
	double mm2_sig = 2.5*0.19404;
	double mm2_min =-10.0;// -0.0673325  - mm2_sig;
	double mm2_max =10.0;// -0.0673325  + mm2_sig;

	double me_sig = 2.5*0.307256;
	double me_min = -10.0;//0.246309 - me_sig;
	double me_max = 10.0;//0.246309 + me_sig;

	double mp_sig = 2.5*0.181323;
	double mp_min = -10.0;//0.297822  - mp_sig;
	double mp_max = 10.0;//0.297822  + mp_sig;
	    
	// WORK ON THIS HERE OKAY!!!!! 
	//hist_pr_phi_pangle_c2->Fill( lv_pr_c2.P(), (lv_phi_c2.Vect()).Angle(lv_pr_c2.Vect()) * 180.0/TMath::Pi() );
	//if ( epKpKm_missing_mass.M2 
	    
	if( epKpKm_missing_mass.M2() <= mm2_max && epKpKm_missing_mass.M2() >= mm2_min ){
	  hist_pr_phi_pangle_c1_mm2->Fill(lv_pr_c1.P(), (lv_phi_c1.Vect()).Angle(lv_pr_c1.Vect()) * 180.0/TMath::Pi());
	  hist_epX_me_c1_cutmm2->Fill( epKpKm_missing_mass.E() );
	  hist_epX_mp_c1_cutmm2->Fill( epKpKm_missing_mass.P() );	      
	}
	if( epKpKm_missing_mass.E() <= me_max && epKpKm_missing_mass.E() >= me_min ){
	  hist_pr_phi_pangle_c1_me->Fill(lv_pr_c1.P(), (lv_phi_c1.Vect()).Angle(lv_pr_c1.Vect()) * 180.0/TMath::Pi());
	  hist_epX_mm2_c1_cute->Fill( epKpKm_missing_mass.M2() );
	  hist_epX_mp_c1_cute->Fill( epKpKm_missing_mass.P() );
	}
	if( epKpKm_missing_mass.P() <= mp_max && epKpKm_missing_mass.P() >= mp_min ){
	  hist_pr_phi_pangle_c1_mp->Fill(lv_pr_c1.P(), (lv_phi_c1.Vect()).Angle(lv_pr_c1.Vect()) * 180.0/TMath::Pi());
	  hist_epX_mm2_c1_cutp->Fill( epKpKm_missing_mass.M2() );
	  hist_epX_me_c1_cutp->Fill( epKpKm_missing_mass.E() );
	}

	if( (epKpKm_missing_mass.E() <= me_max && epKpKm_missing_mass.E() >= me_min ) &&
	    (epKpKm_missing_mass.P() <= mp_max && epKpKm_missing_mass.P() >= mp_min ) ){
	  hist_epKpKmX_mm2_c1_cut2->Fill(epKpKm_missing_mass.M2());
	}

	if( (epKpKm_missing_mass.M2() <= mm2_max && epKpKm_missing_mass.M2() >= mm2_min ) &&
	    (epKpKm_missing_mass.P() <= mp_max && epKpKm_missing_mass.P() >= mp_min ) ){
	  hist_epKpKmX_me_c1_cut2->Fill(epKpKm_missing_mass.E());
	}

	if( (epKpKm_missing_mass.M2() <= mm2_max && epKpKm_missing_mass.M2() >= mm2_min ) &&
	    (epKpKm_missing_mass.E() <= me_max && epKpKm_missing_mass.E() >= me_min ) ){
	  hist_epKpKmX_mp_c1_cut2->Fill(epKpKm_missing_mass.P());
	}

	if((epKpKm_missing_mass.M2() <= mm2_max && epKpKm_missing_mass.M2() >= mm2_min) &&
	   (epKpKm_missing_mass.E() <= me_max && epKpKm_missing_mass.E() >= me_min) &&
	   (epKpKm_missing_mass.P() <= mp_max && epKpKm_missing_mass.P() >= mp_min) )
	  {	      
	    hist_phi_mass_c1_missing_cuts->Fill(inv_mass_kp_kmrec);
	  }
	    
	double me_epphi = epKpKm_missing_mass.E();
	double me_ekpkm = (beam + target - lv_el - lv_kp_c1 - lv_km_c1).E();
	double me_ep = (beam + target - lv_el - lv_pr_c1).E();
	double me_epKpX = (beam+target - lv_el - lv_pr_c1 - lv_kp_c1).E();
	double me_eX = (beam + target - lv_el).E();

	double mp_epphi = epKpKm_missing_mass.P();
	double mp_ekpkm = (beam + target - lv_el - lv_kp_c1 - lv_km_c1).P();
	double mp_ep = (beam + target - lv_el - lv_pr_c1).P();
	double mp_epKpX = (beam+target - lv_el - lv_pr_c1 - lv_kp_c1).P();
	double mp_eX = (beam + target - lv_el).P();

	//missing mass
	hist_epX_mm2_c1->Fill(lv_epX_c1.M2());
	hist_epKpX_mm2_c1->Fill(lv_km_c1.M2());
	hist_eXKpKm_mm2_c1->Fill(lv_eXKpKm_c1.M2());
	hist_eX_mm2_c1->Fill(lv_eX_c1.M2());

	//missing energy
	hist_epX_me_c1->Fill(me_ep);
	hist_epKpX_me_c1->Fill(me_epKpX);
	hist_eXKpKm_me_c1->Fill(me_ekpkm);
	hist_eX_me_c1->Fill(me_eX);
	//missing p
	hist_epX_mp_c1->Fill(mp_ep);
	hist_epKpX_mp_c1->Fill(mp_epKpX);
	hist_eXKpKm_mp_c1->Fill(mp_ekpkm);
	hist_eX_mp_c1->Fill(mp_eX);

	double pr_beta = prot_beta[select_prot_1];
	double kp_beta = kaonP_beta[select_kaonP_1];
	double km_beta = kaonM_beta[select_kaonM_1];

	prot_beta_out = pr_beta;
	kaonP_beta_out = kp_beta;
	kaonM_beta_out = km_beta;

	hist_pr_betap_1->Fill(lv_pr_c1.P(), pr_beta);
	hist_kp_betap_1->Fill(lv_kp_c1.P(), kp_beta);
	    
	if((epKpKm_missing_mass.M2() <= mm2_max && epKpKm_missing_mass.M2() >= mm2_min) &&
	   (epKpKm_missing_mass.E() <= me_max && epKpKm_missing_mass.E() >= me_min) &&
	   (epKpKm_missing_mass.P() <= mp_max && epKpKm_missing_mass.P() >= mp_min) )
	  {		
	    hist_pr_betap_1_final->Fill(lv_pr_c1.P(), pr_beta);
	    hist_kp_betap_1_final->Fill(lv_kp_c1.P(), kp_beta);	   
	  }

	W_out  = kin_W(ele[select_ele]);
	Q2_out = kin_Q2(ele[select_ele]);
	x_out  = kin_x(ele[select_ele]);
	y_out  = kin_y(ele[select_ele]);
	nu_out = kin_nu(ele[select_ele]);
	    
	E_ele  = ele[select_ele].E();
	px_ele = ele[select_ele].Px();
	py_ele = ele[select_ele].Py();
	pz_ele = ele[select_ele].Pz();
	    
	E_prot_1 = lv_pr_c1.E();
	px_prot_1 = lv_pr_c1.Px();
	py_prot_1 = lv_pr_c1.Py();
	pz_prot_1 = lv_pr_c1.Pz();	    

	E_kaonP_1 = lv_kp_c1.E();
	px_kaonP_1 = lv_kp_c1.Px();
	py_kaonP_1 = lv_kp_c1.Py();
	pz_kaonP_1 = lv_kp_c1.Pz();

	E_kaonM_1 = lv_km_c1.E();
	px_kaonM_1 = lv_km_c1.Px();
	py_kaonM_1 = lv_km_c1.Py();
	pz_kaonM_1 = lv_km_c1.Pz();
	channel = 1;
	    
	perp_mntm = epKpKm_missing_mass.Vect().Perp();
	missing_e = epKpKm_missing_mass.E();
	missing_mm2 = epKpKm_missing_mass.M2();
	epX_mm2 = lv_epX_c1.M2();
	epX_mm = lv_epX_c1.M();

	    
	double pT0 = -kin_pT(lv_el,lv_pr_c1);
	double t0 = -kin_t(lv_el,lv_phi_c1);
	double cmphi0 = kin_cmphi(lv_el, lv_phi_c1) * 180.0/TMath::Pi();
	double cmcostheta0 = kin_cmcostheta(ele[select_ele], lv_phi_c1);
	    
	t1_out = t0;
	cmphi1_out  = cmphi0;
	cmcostheta1_out = cmcostheta0;
	pt1_out = pT0;


	//end channel 1 for event selection
      }

      // Check the ep -> epKpKmX channel
      //boolean = 
      
      //std::cout<< " >> " <<  prot_count << " " << Kp_count << " " << Km_count << std::endl;
      if(  prot_count == 1 && Kp_count >= 1 && Km_count >= 1 ){ //prot_count == 1
      // get events with Kp and Km - dont care about el or proton until later
      //if( Kp_count >= 1 && Km_count >= 1 ){
	select_ele = 0 ; select_prot_1 = 0; select_kaonP_1 = 0; select_kaonM_1 = 0;	    
	TLorentzVector lv_pr_c2 = prot[select_prot_1];
	TLorentzVector lv_kp_c2 = Kp[select_kaonP_1];
	TLorentzVector lv_km_c2 = Km[select_kaonM_1];
	TLorentzVector lv_phi_c2 = lv_kp_c2 + lv_km_c2;
	    	    
	TLorentzVector epKpKm_missing_mass = beam + target - lv_el - lv_pr_c2 - lv_kp_c2 - lv_km_c2;

	double mm2_epphi = epKpKm_missing_mass.M2();
	double mm2_ekpkm = (beam + target - lv_el - lv_kp_c2 - lv_km_c2).M2();
	double mm2_ep = (beam + target - lv_el - lv_pr_c2).M2();
	double mm2_epKpX = (beam+target - lv_el - lv_pr_c2 - lv_kp_c2).M2();
	double mm2_eX = (beam + target - lv_el).M2();

	double me_epphi = epKpKm_missing_mass.E();
	double me_ekpkm = (beam + target - lv_el - lv_kp_c2 - lv_km_c2).E();
	double me_ep = (beam + target - lv_el - lv_pr_c2).E();
	double me_epKpX = (beam+target - lv_el - lv_pr_c2 - lv_kp_c2).E();
	double me_eX = (beam + target - lv_el).E();

	double mp_epphi = epKpKm_missing_mass.P();
	double mp_ekpkm = (beam + target - lv_el - lv_kp_c2 - lv_km_c2).P();
	double mp_ep = (beam + target - lv_el - lv_pr_c2).P();
	double mp_epKpX = (beam+target - lv_el - lv_pr_c2 - lv_kp_c2).P();
	double mp_eX = (beam + target - lv_el).P();

	TLorentzVector lv_epX_c2 = beam + target - lv_el - lv_pr_c2;
	hist_epX_c2->Fill(lv_epX_c2.M());

	hist_epphi_mm2_c2->Fill(mm2_epphi);
	hist_ekpkm_mm2_c2->Fill(mm2_ekpkm);
	hist_ep_mm2_c2->Fill(mm2_ep);

	hist_epKpKmX_mm2_c2->Fill(epKpKm_missing_mass.M2());
	hist_epKpKmX_me_c2->Fill(epKpKm_missing_mass.E());
	hist_epKpKmX_mp_c2->Fill(epKpKm_missing_mass.P());

	//missing mass
	hist_epX_mm2_c2->Fill(mm2_ep);
	hist_epKpX_mm2_c2->Fill(mm2_epKpX);
	hist_eXKpKm_mm2_c2->Fill(mm2_ekpkm);
	hist_eX_mm2_c2->Fill(mm2_eX);
	//missing energy
	hist_epX_me_c2->Fill(me_ep);
	hist_epKpX_me_c2->Fill(me_epKpX);
	hist_eXKpKm_me_c2->Fill(me_ekpkm);
	hist_eX_me_c2->Fill(me_eX);
	//missing p
	hist_epX_mp_c2->Fill(mp_ep);
	hist_epKpX_mp_c2->Fill(mp_epKpX);
	hist_eXKpKm_mp_c2->Fill(mp_ekpkm);
	hist_eX_mp_c2->Fill(mp_eX);

	//>> FIT PARAMETERS FOR MM2 -0.0857026 STD 0.218661
	//>> FIT PARAMETERS FOR ME 0.237102 STD 0.256921
	//>> FIT PARAMETERS FOR MP 0.292407 STD 0.1589
	double mm2_sig = 2.5*0.218661;
	double me_sig= 2.5*0.256921;
	double mp_sig = 2.5*0.1589;

	double mm2_min = -0.0857026 - mm2_sig;
	double mm2_max = -0.0857026 + mm2_sig;

	double me_min = -0.237102 - me_sig;
	double me_max = -0.237102 + me_sig;

	double mp_min = -0.292407 - mp_sig;
	double mp_max = -0.292407 + mp_sig;
	    
	// WORK ON THIS HERE OKAY!!!!! 
	//hist_pr_phi_pangle_c2->Fill( lv_pr_c2.P(), (lv_phi_c2.Vect()).Angle(lv_pr_c2.Vect()) * 180.0/TMath::Pi() );
	//if ( epKpKm_missing_mass.M2 
	    
	if( epKpKm_missing_mass.M2() <= mm2_max && epKpKm_missing_mass.M2() >= mm2_min ){
	  hist_pr_phi_pangle_c2_mm2->Fill(lv_pr_c2.P(), (lv_phi_c2.Vect()).Angle(lv_pr_c2.Vect()) * 180.0/TMath::Pi());
	  hist_epKpKmX_me_c2_cutmm2->Fill( epKpKm_missing_mass.E() );
	  //std::cout << " here in mm2 cut " << epKpKm_missing_mass.E() << std::endl;
	  hist_epKpKmX_mp_c2_cutmm2->Fill( epKpKm_missing_mass.P() );	      
	}
	if( epKpKm_missing_mass.E() <= me_max && epKpKm_missing_mass.E() >= me_min ){	      
	  hist_pr_phi_pangle_c2_me->Fill(lv_pr_c2.P(), (lv_phi_c2.Vect()).Angle(lv_pr_c2.Vect()) * 180.0/TMath::Pi());
	  hist_epKpKmX_mm2_c2_cute->Fill( epKpKm_missing_mass.M2() );
	  hist_epKpKmX_mp_c2_cute->Fill( epKpKm_missing_mass.P() );	      
	}
	if( epKpKm_missing_mass.P() <= mp_max && epKpKm_missing_mass.P() >= mp_min ){
	  hist_pr_phi_pangle_c2_mp->Fill(lv_pr_c2.P(), (lv_phi_c2.Vect()).Angle(lv_pr_c2.Vect()) * 180.0/TMath::Pi());
	  hist_epKpKmX_mm2_c2_cutp->Fill( epKpKm_missing_mass.M2() );
	  hist_epKpKmX_me_c2_cutp->Fill( epKpKm_missing_mass.E() );	      
	}


	if( (epKpKm_missing_mass.E() <= me_max && epKpKm_missing_mass.E() >= me_min ) &&
	    (epKpKm_missing_mass.P() <= mp_max && epKpKm_missing_mass.P() >= mp_min ) ){
	  hist_epKpKmX_mm2_c2_cut2->Fill(epKpKm_missing_mass.M2());
	}

	if( (epKpKm_missing_mass.M2() <= mm2_max && epKpKm_missing_mass.M2() >= mm2_min ) &&
	    (epKpKm_missing_mass.P() <= mp_max && epKpKm_missing_mass.P() >= mp_min ) ){
	  hist_epKpKmX_me_c2_cut2->Fill(epKpKm_missing_mass.E());
	}

	if( (epKpKm_missing_mass.M2() <= mm2_max && epKpKm_missing_mass.M2() >= mm2_min ) &&
	    (epKpKm_missing_mass.E() <= me_max && epKpKm_missing_mass.E() >= me_min ) ){
	  hist_epKpKmX_mp_c2_cut2->Fill(epKpKm_missing_mass.P());
	}


	TLorentzVector lv_pr_rec = beam + target - lv_el - lv_kp_c2 - lv_km_c2;
	TLorentzVector lv_kp_rec = beam + target - lv_el - lv_pr_c2 - lv_km_c2;
	TLorentzVector lv_km_rec = beam + target - lv_el - lv_pr_c2 - lv_kp_c2;

	hist_prot_mass_c2->Fill( lv_pr_rec.M() );
	hist_kp_mass_c2->Fill( lv_kp_rec.M() );
	hist_km_mass_c2->Fill( lv_km_rec.M() );

	double angle_kp_kmrec = (lv_kp_c2.Vect()).Angle(lv_km_c2.Vect()) * 180.0/TMath::Pi();
	hist_kp_km_angle_c2->Fill(angle_kp_kmrec);

	double inv_mass_kp_km_c2 = (lv_kp_c2 + lv_km_c2).M();
	hist_phi_mass_c2->Fill(inv_mass_kp_km_c2);

	if((epKpKm_missing_mass.M2() <= mm2_max && epKpKm_missing_mass.M2() >= mm2_min) &&
	   (epKpKm_missing_mass.E() <= me_max && epKpKm_missing_mass.E() >= me_min) &&
	   (epKpKm_missing_mass.P() <= mp_max && epKpKm_missing_mass.P() >= mp_min) )
	  {
	    hist_phi_mass_c2_missing_cuts->Fill(inv_mass_kp_km_c2);
	  }

	if( mm2_epphi >= -2.0 && mm2_epphi <= 2.0 ){
	  hist_phi_mass_c2_cutmm2->Fill(inv_mass_kp_km_c2);
	      
	  if( angle_kp_kmrec >= 15.0 && angle_kp_kmrec <= 25.0 ){
	    hist_phi_mass_c2_cutkaonangles->Fill(inv_mass_kp_km_c2);
	  }
	}

	//pr_sect=sect_pr[select_prot_1];
	//kp_sect=sect_kp[select_kaonP_1];
	//km_sect=sect_kp[select_kaonM_1];

	//	std::cout << " pr sector " << pr_sect << std::endl;
	//std::cout << " kp sector " << kp_sect << std::endl;
	//std::cout << " km sector " << km_sect << std::endl;

	///%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	    
	W_out  = kin_W(ele[select_ele]);
	Q2_out = kin_Q2(ele[select_ele]);
	x_out = kin_x(ele[select_ele]);
	y_out  = kin_y(ele[select_ele]);
	nu_out = kin_nu(ele[select_ele]);
	    
	E_ele  = ele[select_ele].E();
	px_ele = ele[select_ele].Px();
	py_ele = ele[select_ele].Py();
	pz_ele = ele[select_ele].Pz();
	    
	E_prot_1 = lv_pr_c2.E();
	px_prot_1 = lv_pr_c2.Px();
	py_prot_1 = lv_pr_c2.Py();
	pz_prot_1 = lv_pr_c2.Pz();	    

	E_kaonP_1 = lv_kp_c2.E();
	px_kaonP_1 = lv_kp_c2.Px();
	py_kaonP_1 = lv_kp_c2.Py();
	pz_kaonP_1 = lv_kp_c2.Pz();

	E_kaonM_1 = lv_km_c2.E();
	px_kaonM_1 = lv_km_c2.Px();
	py_kaonM_1 = lv_km_c2.Py();
	pz_kaonM_1 = lv_km_c2.Pz();
	channel = 2;

	channel_counter++;

	// Look at the beta vs p cut before and after angle cut between kp and km
	//	double pr_beta = prot_beta[select_prot_1];
	//double kp_beta = kaonP_beta[select_kaonP_1];
	//double km_beta = kaonM_beta[select_kaonM_1];

	//prot_beta_out = pr_beta;
	///kaonP_beta_out = kp_beta;
	//kaonM_beta_out = km_beta;
	    
	perp_mntm = epKpKm_missing_mass.Vect().Perp();
	missing_e = epKpKm_missing_mass.E();
	missing_mm2 = epKpKm_missing_mass.M2();
	epX_mm2 = lv_epX_c2.M2();
	epX_mm = lv_epX_c2.M();

	/*hist_pr_betap_2b->Fill(lv_pr_c2.P(), pr_beta);
	hist_kp_betap_2b->Fill(lv_kp_c2.P(), kp_beta);
	hist_km_betap_2b->Fill(lv_km_c2.P(), km_beta);

	if( angle_kp_kmrec >= 15.0 && angle_kp_kmrec <= 25.0 ){
	  hist_pr_betap_2a->Fill(lv_pr_c2.P(), pr_beta);
	  hist_kp_betap_2a->Fill(lv_kp_c2.P(), kp_beta);
	  hist_km_betap_2a->Fill(lv_km_c2.P(), km_beta);
	}
	*/
	////%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// define the kinematic variables to bin in
	double q2 = kin_Q2(lv_el);
	double xb = kin_x(lv_el);
	double w = kin_W(lv_el);
	double pT0 =-kin_pT(lv_el,lv_pr_c2);
	double t0 =  -2.0*target.M()*( lv_pr_c2.E() - target.M() ); //-kin_t(lv_el,lv_phi_c2);
	double cmphi0 = kin_cmphi(lv_el, lv_phi_c2) * 180.0/TMath::Pi();
	double cmcostheta0 = kin_cmcostheta(ele[select_ele], lv_phi_c2);
	    
	t1_out = t0;
	cmphi1_out  = cmphi0;
	cmcostheta1_out = cmcostheta0;
	pt1_out = pT0;

	////////////////////////////////////////
	// only care about the fully exclusive events for now
	out_tree1.Fill();


	hist_phi_q2->Fill(q2);
	hist_phi_xb->Fill(xb);
	hist_phi_w->Fill(w);
	hist_phi_t->Fill(t0);
	hist_phi_pT->Fill(pT0);
	hist_phi_cmphi->Fill(cmphi0* 180/Pival);
	hist_phi_cmcostheta->Fill(cmcostheta0);


	if((epKpKm_missing_mass.M2() <= mm2_max && epKpKm_missing_mass.M2() >= mm2_min) &&
	   (epKpKm_missing_mass.E() <= me_max && epKpKm_missing_mass.E() >= me_min) &&
	   (epKpKm_missing_mass.P() <= mp_max && epKpKm_missing_mass.P() >= mp_min) )
	  {
	    //hist_pr_betap_final->Fill(lv_pr_c2.P(), pr_beta);
	    ///hist_kp_betap_final->Fill(lv_kp_c2.P(), kp_beta);
	    ///hist_km_betap_final->Fill(lv_km_c2.P(), km_beta);

	    hist_phi_q2_final->Fill(q2);
	    hist_phi_xb_final->Fill(xb);
	    hist_phi_w_final->Fill(w);		
	    hist_phi_t_final->Fill(t0);
	    hist_phi_pT_final->Fill(pT0);
	    hist_phi_cmphi_final->Fill(cmphi0* 180/Pival);
	    hist_phi_cmcostheta_final->Fill(cmcostheta0);

	    h_el_ptheta_final->Fill(lv_el.P(),lv_el.Theta() * 180.0/TMath::Pi());
	    h_el_pphi_final->Fill(lv_el.P(),lv_el.Phi() * 180.0/TMath::Pi());
	    h_el_phitheta_final->Fill(lv_el.Phi() * 180.0/TMath::Pi(),lv_el.Theta() * 180.0/TMath::Pi());
	    
	    h_pr_ptheta_final->Fill(lv_pr_c2.P(),lv_pr_c2.Theta() * 180.0/TMath::Pi());
	    h_pr_pphi_final->Fill(lv_pr_c2.P(),lv_pr_c2.Phi() * 180.0/TMath::Pi());
	    h_pr_phitheta_final->Fill(lv_pr_c2.Phi() * 180.0/TMath::Pi(),lv_pr_c2.Theta() * 180.0/TMath::Pi());

	    h_kp_ptheta_final->Fill(lv_kp_c2.P(),lv_kp_c2.Theta() * 180.0/TMath::Pi());
	    h_kp_pphi_final->Fill(lv_kp_c2.P(),lv_kp_c2.Phi() * 180.0/TMath::Pi());
	    h_kp_phitheta_final->Fill(lv_kp_c2.Phi() * 180.0/TMath::Pi(),lv_kp_c2.Theta() * 180.0/TMath::Pi());
	   	  
	    h_km_ptheta_final->Fill(lv_km_c2.P(),lv_km_c2.Theta() * 180.0/TMath::Pi());
	    h_km_pphi_final->Fill(lv_km_c2.P(),lv_km_c2.Phi() * 180.0/TMath::Pi());
	    h_km_phitheta_final->Fill(lv_km_c2.Phi() * 180.0/TMath::Pi(),lv_km_c2.Theta() * 180.0/TMath::Pi());

	  }

	hist_pr_phi_angle_c2->Fill( (lv_phi_c2.Vect()).Angle(lv_pr_c2.Vect()) * 180.0/TMath::Pi() );
	hist_pr_phi_pangle_c2->Fill( lv_pr_c2.P(), (lv_phi_c2.Vect()).Angle(lv_pr_c2.Vect()) * 180.0/TMath::Pi() );
	   
	// particle kinematic ranges
	h_el_ptheta->Fill(lv_el.P(),lv_el.Theta() * 180.0/TMath::Pi());
	h_el_pphi->Fill(lv_el.P(),lv_el.Phi() * 180.0/TMath::Pi());
	h_el_phitheta->Fill(lv_el.Phi() * 180.0/TMath::Pi(),lv_el.Theta() * 180.0/TMath::Pi());
	    
	h_pr_ptheta->Fill(lv_pr_c2.P(),lv_pr_c2.Theta() * 180.0/TMath::Pi());
	h_pr_pphi->Fill(lv_pr_c2.P(),lv_pr_c2.Phi() * 180.0/TMath::Pi());
	h_pr_phitheta->Fill(lv_pr_c2.Phi() * 180.0/TMath::Pi(),lv_pr_c2.Theta() * 180.0/TMath::Pi());

	h_kp_ptheta->Fill(lv_kp_c2.P(),lv_kp_c2.Theta() * 180.0/TMath::Pi());
	h_kp_pphi->Fill(lv_kp_c2.P(),lv_kp_c2.Phi() * 180.0/TMath::Pi());
	h_kp_phitheta->Fill(lv_kp_c2.Phi() * 180.0/TMath::Pi(),lv_kp_c2.Theta() * 180.0/TMath::Pi());
	   	  
	h_km_ptheta->Fill(lv_km_c2.P(),lv_km_c2.Theta() * 180.0/TMath::Pi());
	h_km_pphi->Fill(lv_km_c2.P(),lv_km_c2.Phi() * 180.0/TMath::Pi());
	h_km_phitheta->Fill(lv_km_c2.Phi() * 180.0/TMath::Pi(),lv_km_c2.Theta() * 180.0/TMath::Pi());
	    
      }
      

    }
  
	//	if( channel > 0 ){
	//	}


    
      /*
      if( ele_count > 0 ){
	TLorentzVector lv_kp = Kp[select_kaonP_1];
	
      }

      */
      /*
      TLorentzVector eX_miss   = beam + target - ele[select_ele];

      hist_eX_mass->Fill(eX_miss.M());

      if( prot_count == 1 ){
	TLorentzVector lv_km = Km[select_kaonM_1];

	hist_epX_mass->Fill( (beam + target - ele[select_ele] - prot[select_prot_1]).M() );
      }

      //if( prot_count > 0 ) std::cout << " particle multiplicity " << prot_count << " " << Kp_count << " " << Km_count << std::endl;

      if( prot_count == 1 && Kp_count >= 1 ){

	TLorentzVector lv_el = ele[select_ele];
	TLorentzVector lv_pr = prot[select_prot_1];
	TLorentzVector lv_kp = Kp[select_kaonP_1];

	//if( Q2 > 1 && x > 0 && W > 2 ) hist_Q2_x->Fill(x,Q2);
 



 
	double phi_mass = Phi_Mass( Kp[select_kaonP_1], Km[select_kaonM_1] );

	h_mass_phi->Fill(phi_mass);

	t1 = -kin_t( ele[select_ele], prot[select_prot_1]);
	
	
	hist_phi_mass->Fill(phi_mass);

	double epkX_mass = kin_epkXMass(ele[select_ele], prot[select_prot_1], Kp[select_kaonP_1]);
	double epKpKmX_mass = kin_epKpKmXMass( ele[select_ele], prot[select_prot_1], Kp[select_kaonP_1], Km[select_kaonM_1]);
	TLorentzVector pKm = prot[select_prot_1] + Km[select_kaonM_1];

	hist_epkX_mass->Fill(epkX_mass);
	hist_epKpKmX->Fill(epKpKmX_mass);

	hist_epKpKmX_p->Fill(prot[select_prot_1].P() * 180.0/TMath::Pi(), epKpKmX_mass);
	hist_epKpKmX_theta->Fill(prot[select_prot_1].Theta() * 180.0/TMath::Pi(), epKpKmX_mass);

	hist_kk_pk->Fill(phi_mass, pKm.M() );

	double phi_mean = 1.069;
	double sig_3 = 0.01752*3.0;
	
	double phi_top = phi_mean + sig_3;
	double phi_bot = phi_mean - sig_3;


	if( phi_mass <= phi_top && phi_mass >= phi_bot   ){	 
	  hist_epKpKmX2->Fill(epKpKmX_mass);

	  //FILL PARTICLE KINEMATICS HERE OF FINAL STATE PARTICLES

	  h_el_ptheta->Fill(lv_el.P(),lv_el.Theta() * 180.0/TMath::Pi());
	  h_el_pphi->Fill(lv_el.P(),lv_el.Phi() * 180.0/TMath::Pi());
	  h_el_phitheta->Fill(lv_el.Phi() * 180.0/TMath::Pi(),lv_el.Theta() * 180.0/TMath::Pi());

	  h_pr_ptheta->Fill(lv_pr.P(),lv_pr.Theta() * 180.0/TMath::Pi());
	  h_pr_pphi->Fill(lv_pr.P(),lv_pr.Phi() * 180.0/TMath::Pi());
	  h_pr_phitheta->Fill(lv_pr.Phi() * 180.0/TMath::Pi(),lv_pr.Theta() * 180.0/TMath::Pi());

	  h_kp_ptheta->Fill(lv_kp.P(),lv_kp.Theta() * 180.0/TMath::Pi());
	  h_kp_pphi->Fill(lv_kp.P(),lv_kp.Phi() * 180.0/TMath::Pi());
	  h_kp_phitheta->Fill(lv_kp.Phi() * 180.0/TMath::Pi(),lv_kp.Theta() * 180.0/TMath::Pi());





	}


	
 


	///hist_epkX_mass->Fill();
	
	px_prot_1 = prot[select_prot_1].Px();
	py_prot_1 = prot[select_prot_1].Py();
	pz_prot_1 = prot[select_prot_1].Pz();
	E_prot_1 = prot[select_prot_1].E();
       
	px_kaonP_1 = Kp[select_kaonP_1].Px();
	py_kaonP_1 = Kp[select_kaonP_1].Py();
	pz_kaonP_1 = Kp[select_kaonP_1].Pz();
	E_kaonP_1 = Kp[select_kaonP_1].E();

	px_kaonM_1 = Km[select_kaonM_1].Px();
	py_kaonM_1 = Km[select_kaonM_1].Py();
	pz_kaonM_1 = Km[select_kaonM_1].Pz();
	E_kaonM_1 = Km[select_kaonM_1].E();

	out_tree1.Fill();
	
      }
      
      */
    
    // }   // end of event builder condition
    /// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    
    /// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///  fill the particle historgrams:
    /// ///////////////////////////////////////////////////////////////

    // fill event properties


    for(Int_t i = 0; i < BUFFER; i++){ 

      if(ele[i].P() > 0) hist_all_electron_p->Fill(ele[i].P());
      if(ele[i].Theta() > 0) hist_all_electron_theta->Fill(ele[i].Theta()*180/Pival);
      if(ele[i].Phi() != 0) hist_all_electron_phi->Fill(ele[i].Phi()*180/Pival);
      if(ele[i].P() > 0) hist_all_electron_p_vs_theta->Fill(ele[i].Theta()*180/Pival, ele[i].P());
      if(ele[i].P() > 0) hist_all_electron_p_vs_phi->Fill(ele[i].Phi()*180/Pival, ele[i].P());
      if(ele[i].Theta() > 0 && ele[i].Phi() != 0) hist_all_electron_theta_vs_phi->Fill(ele[i].Phi()*180/Pival, ele[i].Theta()*180/Pival);

      if(prot[i].P() > 0) hist_all_proton_p->Fill(prot[i].P());
      if(prot[i].Theta() > 0) hist_all_proton_theta->Fill(prot[i].Theta()*180/Pival);
      if(prot[i].Phi() != 0) hist_all_proton_phi->Fill(prot[i].Phi()*180/Pival);
      if(prot[i].P() > 0) hist_all_proton_p_vs_theta->Fill(prot[i].Theta()*180/Pival, prot[i].P());
      if(prot[i].P() > 0) hist_all_proton_p_vs_phi->Fill(prot[i].Phi()*180/Pival, prot[i].P());
      if(prot[i].Theta() > 0 && prot[i].Phi() != 0) hist_all_proton_theta_vs_phi->Fill(prot[i].Phi()*180/Pival, prot[i].Theta()*180/Pival);

    }


    if(ele[select_ele].P() > 0) hist_electron_p->Fill(ele[select_ele].P());
    if(ele[select_ele].Theta() > 0) hist_electron_theta->Fill(ele[select_ele].Theta()*180/Pival);
    if(ele[select_ele].Phi() != 0) hist_electron_phi->Fill(ele[select_ele].Phi()*180/Pival);
    if(ele[select_ele].P() > 0) hist_electron_p_vs_theta->Fill(ele[select_ele].Theta()*180/Pival, ele[select_ele].P());
    if(ele[select_ele].P() > 0) hist_electron_p_vs_phi->Fill(ele[select_ele].Phi()*180/Pival, ele[select_ele].P());
    if(ele[select_ele].Theta() > 0 && ele[select_ele].Phi() != 0) hist_electron_theta_vs_phi->Fill(ele[select_ele].Phi()*180/Pival, ele[select_ele].Theta()*180/Pival);

    if(prot[select_prot_1].P() > 0) hist_proton_p->Fill(prot[select_prot_1].P());
    if(prot[select_prot_1].Theta() > 0) hist_proton_theta->Fill(prot[select_prot_1].Theta()*180/Pival);
    if(prot[select_prot_1].Phi() != 0) hist_proton_phi->Fill(prot[select_prot_1].Phi()*180/Pival);
    if(prot[select_prot_1].P() > 0) hist_proton_p_vs_theta->Fill(prot[select_prot_1].Theta()*180/Pival, prot[select_prot_1].P());
    if(prot[select_prot_1].P() > 0) hist_proton_p_vs_phi->Fill(prot[select_prot_1].Phi()*180/Pival, prot[select_prot_1].P());
    if(prot[select_prot_1].Theta() > 0 && prot[select_prot_1].Phi() != 0) hist_proton_theta_vs_phi->Fill(prot[select_prot_1].Phi()*180/Pival, prot[select_prot_1].Theta()*180/Pival);


    // fill missing mass and missing energy plots:

    if(M_e_p_X_miss > 0 && prot_count > 0) hist_e_p_X_mismass->Fill(M_e_p_X_miss);
    if(M_e_p_X_miss2 != 0 && prot_count > 0) hist_e_p_X_mismass2->Fill(M_e_p_X_miss2);



    /// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///  fill the statistics histograms:
    /// ///////////////////////////////////////////////////////////////

    hist_electron_count->Fill(ele_count);
    hist_proton_count->Fill(prot_count);
    hist_neutron_count->Fill(neutr_count);
    hist_pip_count->Fill(pip_count);
    hist_pim_count->Fill(pim_count);
    hist_Kp_count->Fill(Kp_count);
    hist_Km_count->Fill(Km_count);
    hist_photon_count->Fill(phot_count);


    /// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  
  }  // end of event loop
  /// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // FIT THE MISSING MASS2 and Missing energy and missing momentum plots here.

  TF1 *fit_epKpKm_mm2 = new TF1("fit_epKpKm_mm2","gaus",-0.5,0.5);
  TF1 *fit_epKpKm_me = new TF1("fit_epKpKm_me","gaus",-0.5,0.5);
  TF1 *fit_epKpKm_mp = new TF1("fit_epKpKm_mp","gaus",-0.5,0.5);

  TCanvas *c1 = new TCanvas("c1","",1800,600);
  c1->Divide(3,1);
  c1->cd(1);
  hist_epKpKmX_mm2_c2->Fit("fit_epKpKm_mm2","R+");
  hist_epKpKmX_mm2_c2->Draw("same");
  fit_epKpKm_mm2->Draw("same");
  c1->cd(2);
  hist_epKpKmX_me_c2->Fit("fit_epKpKm_me","R+");
  hist_epKpKmX_me_c2->Draw("same");
  fit_epKpKm_me->Draw("same");
  c1->cd(3);
  hist_epKpKmX_mp_c2->Fit("fit_epKpKm_mp","R+");
  hist_epKpKmX_mp_c2->Draw("same");
  fit_epKpKm_mp->Draw("same");

  c1->SaveAs("h_missing_variables_fit.pdf");

  double mean_mm2 = fit_epKpKm_mm2->GetParameter(1);
  double sig_mm2 = fit_epKpKm_mm2->GetParameter(2);

  double mean_me = fit_epKpKm_me->GetParameter(1);
  double sig_me = fit_epKpKm_me->GetParameter(2);
  
  double mean_mp = fit_epKpKm_mp->GetParameter(1);
  double sig_mp = fit_epKpKm_mp->GetParameter(2);
  
  std::cout << " >> FIT PARAMETERS FOR MM2 " << mean_mm2 << " STD " << sig_mm2 << std::endl;
  std::cout << " >> FIT PARAMETERS FOR ME " << mean_me << " STD " << sig_me << std::endl;
  std::cout << " >> FIT PARAMETERS FOR MP " << mean_mp << " STD " << sig_mp << std::endl;

  std::cout << " EVENTS IN CHANNEL 2 " << channel_counter << std::endl;

  cout << endl;
  cout << "Tree successfully analysed!" << endl;
  cout << "Writing the output file ... " << endl;
  out->Write(); // Saving Histograms
  cout << "Histograms saved in File: " << outputfile << endl;
  out->Close(); // Closing Output File
  f->Close();   // Closing Input File
  cout << "... Completed!" << endl;


  
  return 1;
    

  /// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
}   /// end of main
/// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////
/// angles between particles

double alpha_e_X(TLorentzVector ele, TLorentzVector hadron){

  TLorentzVector beam(0,0,Ebeam,Ebeam);
  TLorentzVector target(0,0,0,0.93827);
  TLorentzVector miss   = beam + target - ele - hadron;

  return ele.Vect().Angle(miss.Vect());
}

double alpha_p1p2(TLorentzVector particle1, TLorentzVector particle2){

  return particle1.Vect().Angle(particle2.Vect());
}


///////////////////////////////////////////////////////////////
///  kinematics

double kin_W(TLorentzVector ele){
  TLorentzVector beam(0,0,Ebeam,Ebeam);
  TLorentzVector target(0,0,0,0.93827);
  TLorentzVector fGamma = beam - ele;
  TLorentzVector fCM = fGamma + target;
  return fCM.M();
}

double kin_Q2(TLorentzVector ele){
  TLorentzVector beam(0,0,Ebeam,Ebeam);
  TLorentzVector target(0,0,0,0.93827);
  TLorentzVector fGamma = beam - ele;
  return -fGamma.M2();
}

double kin_x(TLorentzVector ele){
  TLorentzVector beam(0,0,Ebeam,Ebeam);
  TLorentzVector target(0,0,0,0.93827);
  TLorentzVector fGamma = beam - ele;
  return -fGamma.M2()/(2*0.938272*fGamma.E());
}


double kin_y(TLorentzVector ele){
  TLorentzVector beam(0,0,Ebeam,Ebeam);
  TLorentzVector target(0,0,0,0.93827);
  TLorentzVector fGamma = beam - ele;
  return fGamma.E()/Ebeam;
}

double kin_nu(TLorentzVector ele){
  TLorentzVector beam(0,0,Ebeam,Ebeam);
  TLorentzVector target(0,0,0,0.93827);
  TLorentzVector fGamma = beam - ele;
  return fGamma.E();
}

double kin_pT(TLorentzVector ele, TLorentzVector hadron){
  TLorentzVector beam(0,0,Ebeam,Ebeam);
  TLorentzVector target(0,0,0,0.93827);
  TLorentzVector fGamma = beam - ele;
  TLorentzVector fCM = fGamma + target;
  TVector3 boost = -fCM.BoostVector();
  TLorentzVector hadronClone(hadron); 
  hadronClone.Boost(boost);
  return hadronClone.Perp(); 
}

double kin_eta(TLorentzVector ele, TLorentzVector hadron){ 
  TLorentzVector beam(0,0,Ebeam,Ebeam);
  TLorentzVector target(0,0,0,0.93827);
  TLorentzVector fGamma = beam - ele;
  TLorentzVector fCM = fGamma + target;
  TVector3 boost = -fCM.BoostVector();
  TLorentzVector hadronClone(hadron); 
  hadronClone.Boost(boost);
  return hadronClone.Rapidity(); 
}

double kin_z(TLorentzVector ele, TLorentzVector hadron){
  TLorentzVector beam(0,0,Ebeam,Ebeam);
  TLorentzVector target(0,0,0,0.93827);
  TLorentzVector fGamma = beam - ele;  
  return hadron.T()/fGamma.T();
}

double kin_t(TLorentzVector ele, TLorentzVector hadron){

  TLorentzVector beam(0,0,Ebeam,Ebeam);
  TLorentzVector target(0,0,0,0.93827);
  TLorentzVector fGamma = beam - ele;
  TLorentzVector transfer = fGamma - hadron;
  return transfer.M2();
}

double kin_cmphi(TLorentzVector ele, TLorentzVector hadron){

  TLorentzVector beam(0,0,Ebeam,Ebeam);
  TLorentzVector target(0,0,0,0.93827);

  TLorentzVector fGamma = beam - ele;
  TLorentzVector fCM = fGamma + target;

  TVector3 boost = -fCM.BoostVector();  // get boost vector to CM frame

  TLorentzVector gammaClone(fGamma);    // clone particles
  TLorentzVector eleClone(ele); 
  TLorentzVector hadronClone(hadron); 

  gammaClone.Boost(boost);     // boost particles to CM frame
  eleClone.Boost(boost);
  hadronClone.Boost(boost);

  TVector3 v0, v1;
  v0 = gammaClone.Vect().Cross(eleClone.Vect());
  v1 = gammaClone.Vect().Cross(hadronClone.Vect());    // hadronClone.Vect().Cross(gammaClone.Vect());

  Double_t c0, c1, c2, c3;
  c0 = v0.Dot(hadronClone.Vect());
  c1 = v0.Dot(v1);
  c2 = v0.Mag();
  c3 = v1.Mag();

  return (c0/TMath::Abs(c0)) * TMath::ACos(c1 /(c2*c3));
}

double kin_cmcostheta(TLorentzVector ele, TLorentzVector hadron){

  TLorentzVector beam(0,0,Ebeam,Ebeam);
  TLorentzVector target(0,0,0,0.93827);

  TLorentzVector fGamma = beam - ele;
  TLorentzVector fCM = fGamma + target;

  TVector3 boost = -fCM.BoostVector();

  TLorentzVector hadronClone(hadron); 
  hadronClone.Boost(boost);

  return TMath::Cos(hadronClone.Theta()); 
}


/////////////////////////////////////////////////////////////////////
///  missing mass and energy

TLorentzVector miss_X(TLorentzVector ele, TLorentzVector hadron){

  TLorentzVector beam(0,0,Ebeam,Ebeam);
  TLorentzVector target(0,0,0,0.93827);
  TLorentzVector miss   = beam + target - ele - hadron;
  return miss;
}

double kin_mismass(TLorentzVector ele, TLorentzVector hadron){
  TLorentzVector beam(0,0,Ebeam,Ebeam);
  TLorentzVector target(0,0,0,0.93827);
  TLorentzVector miss   = beam + target - ele - hadron;
  return sqrt(miss.E()*miss.E() - miss.Px()*miss.Px() - miss.Py()*miss.Py() - miss.Pz()*miss.Pz());
}

double kin_mismass2(TLorentzVector ele, TLorentzVector hadron){
  TLorentzVector beam(0,0,Ebeam,Ebeam);
  TLorentzVector target(0,0,0,0.93827);
  TLorentzVector miss   = beam + target - ele - hadron;
  return (miss.E()*miss.E() - miss.Px()*miss.Px() - miss.Py()*miss.Py() - miss.Pz()*miss.Pz());
}

double Phi_Mass( TLorentzVector kaonP, TLorentzVector kaonM ){

  TLorentzVector lv_phi(0,0,0,0);
  lv_phi =  kaonP + kaonM;
   
  return lv_phi.M();
}


double kin_epkXMass( TLorentzVector el, TLorentzVector pr, TLorentzVector kp){

  TLorentzVector beam(0,0,Ebeam,Ebeam);
  TLorentzVector target(0,0,0,0.93827);

  TLorentzVector temp = beam + target - el - pr - kp;

  return temp.M2();

}

double kin_epKpKmXMass( TLorentzVector el, TLorentzVector pr, TLorentzVector kp, TLorentzVector km){

  TLorentzVector beam(0,0,Ebeam,Ebeam);
  TLorentzVector target(0,0,0,0.93827);

  //cout << " >> " <<  beam.E() << " " << beam.Pz() << " " << el.E() << " " << pr.E() << " " << kp.E() << " " << km.E() << endl;

  TLorentzVector temp = beam + target - el - pr - kp - km; 

  return temp.M2();

}


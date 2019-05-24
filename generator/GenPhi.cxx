#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <TROOT.h>
#include <TFile.h>
#include <TChain.h>
#include <TNtuple.h>
#include <TLorentzVector.h>
#include <TRandom3.h>
#include <TF1.h>
#include <time.h>
#include <TSystem.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TStyle.h>

using namespace std;

const Double_t alpha_em = 1./137.036;
const Double_t EBeam = 10.6;
const Double_t ProtonMass = 0.93827;
const Double_t phimass = 1.01946;
const Double_t Kmass = 0.49368;
const Double_t me2 = TMath::Power(0.000511,2.);

TFile *outfile;
TTree *outtree;

TRandom3 *ThisMCrandom;
TLorentzVector v4Beam, v4Target, v4Elec, v4Prot, v4Kplus, v4Kminus, v4phi, v4Q;
Double_t gen_Q2, gen_xB, gen_t, gen_phi, gen_phie;
Double_t gen_elec_p, gen_elec_t, gen_elec_f;
Double_t gen_prot_p, gen_prot_t, gen_prot_f;
Double_t gen_Kplus_P, gen_Kplus_theta, gen_Kplus_phi, gen_Kplus_el_ang;
Double_t gen_Kminus_P, gen_Kminus_theta, gen_Kminus_phi, gen_Kminus_el_ang;
Double_t gen_phi_p, gen_phi_t, gen_phi_f, gen_gg_angle;
Double_t gen_phi_CMt, gen_phi_CMf;
Double_t gen_vx, gen_vy, gen_vz;
Double_t phi_xsec;

void Kinematics(Double_t x, Double_t q, Double_t t, Double_t phi);
Bool_t GenAcc(Double_t x, Double_t q, Double_t t);
void MakeTree(void);
void FillTree(void);

Double_t GetMinQ2(Double_t Xbj);
Double_t GetMaxQ2(Double_t Xbj);
Double_t GetQ2_thetaX(Double_t Xbj, Double_t theta);
Double_t GetQ2_WX(Double_t Xbj, Double_t W);
Double_t GetQ2_EpX(Double_t Xbj, Double_t Ep);

Bool_t AboveTMin(Double_t DVCS_t_Pr, Double_t DVCS_Xbj, Double_t DVCS_Q2);
double dvmpx(double t, double xb, double Q2, double PHI_G , double cos_theta_H, double E);

int main(Int_t argc, Char_t *argv[]){
cout << "---------------------------------------------------" << endl;
cout << "|               phi CLAS12 generator               |" << endl;
cout << "---------------------------------------------------" << endl;
cout << "!!!!!!!!!!!!!! EBeam = ";cout<<Form("%1.2f",EBeam);cout<<" GeV !!!!!!!!!!!!!" << endl;

gROOT->SetStyle("Plain");
gStyle->SetOptTitle(0);
gStyle->SetPalette(1);
gStyle->SetOptStat(0);
gStyle->SetOptFit(0);
gStyle->SetOptDate(0);

int expo = 4;
if(argc>1)expo=atoi(argv[1]);

ThisMCrandom = new TRandom3(0);
for(int MCinit=0;MCinit<100;MCinit++)ThisMCrandom->Uniform();
cout << "My seed " << ThisMCrandom->Uniform() << endl;
time_t prevtimeout, currenttimeout;
cout << endl << "The time is now ";cout.flush(); gSystem->Exec("date");cout<<endl;

Double_t minxB, minQ2, mint, minphi;
Double_t maxxB, maxQ2, maxt, maxphi;

minxB = 0.05;
maxxB = 0.8;
minQ2 = 0.8;
maxQ2 = 14.;
mint = 0.01;
maxt = 5.;

minphi = 0.5*TMath::DegToRad();
maxphi = 359.5*TMath::DegToRad();

Long_t NeventsTot = (Long_t)TMath::FloorNint(TMath::Power(10.,expo));
time ( &prevtimeout);

Double_t thisX, thisQ, thisT, thisPhi;

MakeTree();

ofstream outlund;
outlund.open("rename_phi_lund.txt");
outlund.precision(4);
outlund << right;
for(int ie=0;ie<NeventsTot;){
	Long_t evtsPrint = (Long_t)TMath::FloorNint(ie/(NeventsTot/100));
	evtsPrint = evtsPrint*TMath::FloorNint(NeventsTot/100);
	if(ie==evtsPrint){
		evtsPrint = TMath::FloorNint(ie/(NeventsTot/100));
		time ( &currenttimeout );
		double dif = difftime (currenttimeout,prevtimeout);
		cout << "\r" << evtsPrint << "\% ("<<ie<<"/"
			<<NeventsTot<<") last time to process 1\% of events : "<< dif << "s ... ("<<difftime(currenttimeout,prevtimeout)
		        <<" , "<<currenttimeout<<" , "<<prevtimeout<<")         ";cout.flush();
		time ( &prevtimeout);
	}
	thisX = minxB+(maxxB-minxB)*ThisMCrandom->Uniform();
	thisQ = minQ2+(maxQ2-minQ2)*ThisMCrandom->Uniform();
	thisT = mint+(maxt-mint)*ThisMCrandom->Uniform();
	thisPhi = minphi+(maxphi-minphi)*ThisMCrandom->Uniform();
	if(GenAcc(thisX,thisQ,thisT)){
		Kinematics(thisX,thisQ,thisT,thisPhi*TMath::DegToRad());
		phi_xsec = dvmpx(-thisT,thisX,thisQ,thisPhi,TMath::Cos(TMath::DegToRad()*gen_phi_CMt),EBeam);
		if(phi_xsec>ThisMCrandom->Uniform()){
			FillTree();
			ie++;
                 outlund << "4 1 1 0 -1 " 
                         << gen_xB << " " << gen_Q2/(2*ProtonMass*gen_xB*EBeam) << " " << TMath::Sqrt(ProtonMass*ProtonMass+gen_Q2*(1/gen_xB-1)) << " " << gen_Q2 << " " << gen_Q2/(2*ProtonMass*gen_xB) << endl;
                 outlund << "1 -1 1   11 0 0 " 
                         << setw(8) << v4Elec.Px() << " " << setw(8) << v4Elec.Py() << " " << setw(8) << v4Elec.Pz() << " " << setw(8) << v4Elec.E() << " 0.0005 " << gen_vx << " " << gen_vy << " " << gen_vz << endl;
                 outlund << "2  1 1 2212 0 0 " 
                         << setw(8) << v4Prot.Px() << " " << setw(8) << v4Prot.Py() << " " << setw(8) << v4Prot.Pz() << " " << setw(8) << v4Prot.E() << " 0.9383 " << gen_vx << " " << gen_vy << " " << gen_vz << endl;
                 outlund << "3  1 1  321 0 0 " 
                         << setw(8) << v4Kplus.Px() << " " << setw(8) << v4Kplus.Py() << " " << setw(8) << v4Kplus.Pz() << " " << setw(8) << v4Kplus.E() << " 0.4937 " << gen_vx << " " << gen_vy << " " << gen_vz << endl;
                 outlund << "4 -1 1 -321 0 0 " 
                         << setw(8) << v4Kminus.Px() << " " << setw(8) << v4Kminus.Py() << " " << setw(8) << v4Kminus.Pz() << " " << setw(8) << v4Kminus.E() << " 0.4937 " << gen_vx << " " << gen_vy << " " << gen_vz << endl;

		}
	}
}
 cout << endl << "Over" << endl;
 cout << endl << "The time is now ";cout.flush(); gSystem->Exec("date");cout<<endl;
 outtree->Write();
 outfile->Close();
 outlund.close();
 return 0;
}

void MakeTree(void){
 outfile = new TFile("phi.root","recreate");
 outtree = new TTree("phi","phi");
 outtree->Branch("gen_Q2",&gen_Q2,"gen_Q2/D");
 outtree->Branch("gen_xB",&gen_xB,"gen_xB/D");
 outtree->Branch("gen_t",&gen_t,"gen_t/D");
 outtree->Branch("gen_phi",&gen_phi,"gen_phi/D");
 outtree->Branch("gen_phie",&gen_phie,"gen_phie/D");
 outtree->Branch("gen_elec_p",&gen_elec_p,"gen_elec_p/D");
 outtree->Branch("gen_elec_t",&gen_elec_t,"gen_elec_t/D");
 outtree->Branch("gen_elec_f",&gen_elec_f,"gen_elec_f/D");
 outtree->Branch("gen_prot_p",&gen_prot_p,"gen_prot_p/D");
 outtree->Branch("gen_prot_t",&gen_prot_t,"gen_prot_t/D");
 outtree->Branch("gen_prot_f",&gen_prot_f,"gen_prot_f/D");
 outtree->Branch("gen_Kplus_P",&gen_Kplus_P,"gen_Kplus_P/D");
 outtree->Branch("gen_Kplus_theta",&gen_Kplus_theta,"gen_Kplus_theta/D");
 outtree->Branch("gen_Kplus_phi",&gen_Kplus_phi,"gen_Kplus_phi/D");
 outtree->Branch("gen_Kminus_P",&gen_Kminus_P,"gen_Kminus_P/D");
 outtree->Branch("gen_Kminus_theta",&gen_Kminus_theta,"gen_Kminus_theta/D");
 outtree->Branch("gen_Kminus_phi",&gen_Kminus_phi,"gen_Kminus_phi/D");
 outtree->Branch("gen_phi_p",&gen_phi_p,"gen_phi_p/D");
 outtree->Branch("gen_phi_t",&gen_phi_t,"gen_phi_t/D");
 outtree->Branch("gen_phi_f",&gen_phi_f,"gen_phi_f/D");
 outtree->Branch("gen_gg_angle",&gen_gg_angle,"gen_gg_angle/D");
 outtree->Branch("gen_phi_CMt",&gen_phi_CMt,"gen_phi_CMt/D");
 outtree->Branch("gen_phi_CMf",&gen_phi_CMf,"gen_phi_CMf/D");
 outtree->Branch("gen_vx",&gen_vx,"gen_vx/D");
 outtree->Branch("gen_vy",&gen_vy,"gen_vy/D");
 outtree->Branch("gen_vz",&gen_vz,"gen_vz/D");
 outtree->Branch("phi_xsec",&phi_xsec,"phi_xsec/D");
}

void FillTree(void){
 outtree->Fill();
}

void Kinematics(Double_t x, Double_t q, Double_t t, Double_t phi){
 //kinematic variables in the lab

 v4Beam.SetXYZT(0,0,EBeam,TMath::Sqrt( EBeam*EBeam+me2));
 v4Target.SetXYZT(0.,0.,0.,ProtonMass);
 Double_t Px, Py, Pz, P, E, Theta, Phi;
 E = EBeam - q/(2*ProtonMass*x);
 P = TMath::Sqrt(E*E-me2);
 Theta = 2.*TMath::ASin( TMath::Sqrt(q/(4*EBeam*E)) );
 Phi = ThisMCrandom->Uniform()*2.*TMath::Pi();
 Px = P*TMath::Sin(Theta)*TMath::Cos(Phi);
 Py = P*TMath::Sin(Theta)*TMath::Sin(Phi);
 Pz = P*TMath::Cos(Theta);
 v4Elec.SetXYZT(Px,Py,Pz,E);
 v4Q = v4Beam-v4Elec;

 gen_Q2 = q;
 gen_xB = x;
 gen_t = t;
 gen_phi = TMath::RadToDeg()*phi;
 gen_phie = TMath::RadToDeg()*Phi;
 float W2 = ProtonMass*ProtonMass+gen_Q2*(1./gen_xB-1.);
 float E1CM = (W2+gen_Q2+ProtonMass*ProtonMass)/(2*TMath::Sqrt(W2));
 float P1CM = TMath::Sqrt(E1CM*E1CM-ProtonMass*ProtonMass);
 float E3CM = (W2+ProtonMass*ProtonMass-phimass*phimass)/(2*TMath::Sqrt(W2));
 float P3CM = TMath::Sqrt(E3CM*E3CM-ProtonMass*ProtonMass);
 float TMIN = TMath::Power(gen_Q2+phimass*phimass,2)/(4.*W2)-TMath::Power(P1CM-P3CM,2);
 float thetacm = 2.*TMath::ASin(TMath::Sqrt( (gen_t-TMath::Abs(TMIN))/(4*P1CM*P3CM) ));
 TLorentzVector v4ProtCM, v4phiCM, v4KplusCM, v4KminusCM, v4QCM, v4TargetCM ;
 v4QCM = v4Q;
 v4TargetCM = v4Target;
 TVector3 BETA, CMboost;
 BETA = (v4Q+v4Target).BoostVector();
 TVector3 LABboost = BETA;
 CMboost = - BETA;
 v4QCM.Boost(CMboost);
 Double_t QThetaCM = (v4QCM.Vect()).Theta();
 Double_t QPhiCM = (v4QCM.Vect()).Phi();
 v4TargetCM.Boost(CMboost);
 v4TargetCM.RotateZ(-QPhiCM);
 v4QCM.RotateZ(-QPhiCM);
 v4TargetCM.RotateY(-QThetaCM);
 v4QCM.RotateY(-QThetaCM);

 E = TMath::Sqrt(phimass*phimass+P3CM*P3CM);
 P = P3CM;
 Px = P*TMath::Sin(thetacm)*TMath::Cos(phi+TMath::Pi());
 Py = P*TMath::Sin(thetacm)*TMath::Sin(phi+TMath::Pi());
 Pz = P*TMath::Cos(thetacm);
 v4phiCM.SetXYZT(Px,Py,Pz,E);

 thetacm = TMath::Pi()-thetacm;
 E = TMath::Sqrt(P*P+ProtonMass*ProtonMass);
 Px = P*TMath::Sin(thetacm)*TMath::Cos(phi);
 Py = P*TMath::Sin(thetacm)*TMath::Sin(phi);
 Pz = P*TMath::Cos(thetacm);
 v4ProtCM.SetXYZT(Px,Py,Pz,E);

 v4ProtCM.RotateY(QThetaCM);
 v4phiCM.RotateY(QThetaCM);
 v4ProtCM.RotateZ(QPhiCM);
 v4phiCM.RotateZ(QPhiCM);
 v4ProtCM.Boost(LABboost);
 v4phiCM.Boost(LABboost);
 v4Prot = v4ProtCM;
 v4phi = v4phiCM;

 //generate Kaons here
 BETA = v4phi.BoostVector();
 Theta = 2*ThisMCrandom->Uniform()-1.;//cos(theta) is uniform
 Phi = 2*TMath::Pi()*ThisMCrandom->Uniform();//phi is uniform
 E = phimass/2.;
 P = TMath::Sqrt(E*E-Kmass*Kmass);
 Px = P*TMath::Sqrt(1-Theta*Theta)*TMath::Cos(Phi);
 Py = P*TMath::Sqrt(1-Theta*Theta)*TMath::Sin(Phi);
 Pz = P*Theta;
 v4Kplus.SetXYZT(Px,Py,Pz,E);
 v4Kminus.SetXYZT(-Px,-Py,-Pz,E);
 v4Kplus.Boost(BETA);
 v4Kminus.Boost(BETA);

 gen_elec_p = v4Elec.P();
 gen_elec_t = v4Elec.Theta()*TMath::RadToDeg();
 gen_elec_f = v4Elec.Phi()*TMath::RadToDeg();
 gen_prot_p = v4Prot.P();
 gen_prot_t = v4Prot.Theta()*TMath::RadToDeg();
 gen_prot_f = v4Prot.Phi()*TMath::RadToDeg();
 gen_Kplus_P = v4Kplus.P();
 gen_Kplus_theta = TMath::RadToDeg()*v4Kplus.Theta();
 gen_Kplus_phi = TMath::RadToDeg()*v4Kplus.Phi();
 gen_Kplus_el_ang = TMath::RadToDeg()*(v4Kplus.Vect()).Angle(v4Elec.Vect());
 gen_Kminus_P = v4Kminus.P();
 gen_Kminus_theta = TMath::RadToDeg()*v4Kminus.Theta();
 gen_Kminus_phi = TMath::RadToDeg()*v4Kminus.Phi();
 gen_Kminus_el_ang = TMath::RadToDeg()*(v4Kminus.Vect()).Angle(v4Elec.Vect());
 gen_phi_p = v4phi.P();
 gen_phi_t = TMath::RadToDeg()*v4phi.Theta();
 gen_phi_f = TMath::RadToDeg()*v4phi.Phi();
 gen_gg_angle = TMath::RadToDeg() * (v4Kplus.Vect()).Angle(v4Kminus.Vect()) ;
 gen_phi_CMt = TMath::RadToDeg()*TMath::ACos(Theta);
 gen_phi_CMf = TMath::RadToDeg()*Phi;

 gen_vx = 0;//0.04*ThisMCrandom->Gaus();
 gen_vy = 0;//0.04*ThisMCrandom->Gaus();
 gen_vz = -2.5 + 5.0*ThisMCrandom->Uniform();

}

Bool_t GenAcc(Double_t x, Double_t q, Double_t t){
 Double_t Eprime = EBeam-q/(2*ProtonMass*x);
 Double_t El_Theta = 2*TMath::RadToDeg()*TMath::ASin(TMath::Sqrt(q/(4*EBeam*Eprime) ));
 Double_t W = TMath::Sqrt(ProtonMass*ProtonMass+q*(1./x-1.));
 if( q>0.8 && x>0.05 && El_Theta>4. && W>1.5 && Eprime>0.5 && AboveTMin(t,x,q) )return kTRUE;
 return kFALSE;
}

Bool_t AboveTMin(Double_t DVCS_t_Pr, Double_t DVCS_Xbj, Double_t DVCS_Q2){
 Double_t DVCS_W2 = ProtonMass*ProtonMass+DVCS_Q2*(1./DVCS_Xbj-1.); 
 Double_t E1CM = (DVCS_W2+DVCS_Q2+ProtonMass*ProtonMass)/(2*TMath::Sqrt(DVCS_W2));
 Double_t P1CM = TMath::Sqrt(E1CM*E1CM-ProtonMass*ProtonMass);
 Double_t E3CM = (DVCS_W2-phimass*phimass+ProtonMass*ProtonMass)/(2*TMath::Sqrt(DVCS_W2));
 Double_t P3CM = TMath::Sqrt(E3CM*E3CM-ProtonMass*ProtonMass);
 Double_t TMIN = TMath::Power(DVCS_Q2+phimass*phimass,2)/(4.*DVCS_W2)-TMath::Power(P1CM-P3CM,2);
 Double_t TMAX = TMath::Power(DVCS_Q2+phimass*phimass,2)/(4.*DVCS_W2)-TMath::Power(P1CM+P3CM,2);
 return (TMath::Abs(DVCS_t_Pr)>TMath::Abs(TMIN)&&TMath::Abs(DVCS_t_Pr)<TMath::Abs(TMAX));
}

Double_t GetTmin(Double_t DVCS_Xbj, Double_t DVCS_Q2){
 Double_t DVCS_W2 = ProtonMass*ProtonMass+DVCS_Q2*(1./DVCS_Xbj-1.); 
 Double_t E1CM = (DVCS_W2+DVCS_Q2+ProtonMass*ProtonMass)/(2*TMath::Sqrt(DVCS_W2));
 Double_t P1CM = TMath::Sqrt(E1CM*E1CM-ProtonMass*ProtonMass);
 Double_t E3CM = (DVCS_W2-phimass*phimass+ProtonMass*ProtonMass)/(2*TMath::Sqrt(DVCS_W2));
 Double_t P3CM = TMath::Sqrt(E3CM*E3CM-ProtonMass*ProtonMass);
 Double_t TMIN = TMath::Power(DVCS_Q2+phimass*phimass,2)/(4.*DVCS_W2)-TMath::Power(P1CM-P3CM,2);
 return TMIN;
}

Double_t GetTmax(Double_t DVCS_Xbj, Double_t DVCS_Q2){
 Double_t DVCS_W2 = ProtonMass*ProtonMass+DVCS_Q2*(1./DVCS_Xbj-1.);
 Double_t E1CM = (DVCS_W2+DVCS_Q2+ProtonMass*ProtonMass)/(2*TMath::Sqrt(DVCS_W2));
 Double_t P1CM = TMath::Sqrt(E1CM*E1CM-ProtonMass*ProtonMass);
 Double_t E3CM = (DVCS_W2-phimass*phimass+ProtonMass*ProtonMass)/(2*TMath::Sqrt(DVCS_W2));
 Double_t P3CM = TMath::Sqrt(E3CM*E3CM-ProtonMass*ProtonMass);
 Double_t TMAX = TMath::Power(DVCS_Q2+phimass*phimass,2)/(4.*DVCS_W2)-TMath::Power(P1CM+P3CM,2);
 return TMAX;
}

Double_t GetT(Double_t DVCS_Xbj, Double_t DVCS_Q2, Double_t theta_CM){
 Double_t DVCS_W2 = ProtonMass*ProtonMass+DVCS_Q2*(1./DVCS_Xbj-1.);
 Double_t E1CM = (DVCS_W2+DVCS_Q2+ProtonMass*ProtonMass)/(2*TMath::Sqrt(DVCS_W2));
 Double_t P1CM = TMath::Sqrt(E1CM*E1CM-ProtonMass*ProtonMass);
 Double_t E3CM = (DVCS_W2-phimass*phimass+ProtonMass*ProtonMass)/(2*TMath::Sqrt(DVCS_W2));
 Double_t P3CM = TMath::Sqrt(E3CM*E3CM-ProtonMass*ProtonMass);
 return GetTmin(DVCS_Xbj,DVCS_Q2)-4.*P1CM*P3CM*TMath::Power(TMath::Sin(TMath::DegToRad()*theta_CM/2.),2.);
}

Double_t Fun_T(Double_t *X, Double_t *PAR){
 return GetT(PAR[0],PAR[1],X[0]);
}

Double_t GetMinQ2(Double_t Xbj){
 Double_t result = 1.;
 Double_t thetalim = GetQ2_thetaX(Xbj,21.);
 Double_t Wlim = GetQ2_WX(Xbj,2.);
 //cout << "thetalim = " << thetalim << endl;
 //cout << "Wlim = " << Wlim << endl;
 if(TMath::Max(thetalim,Wlim)>1.)result = TMath::Max(thetalim,Wlim);
 return result;
}

Double_t GetMaxQ2(Double_t Xbj){
 Double_t result = 1.;
 Double_t eplim = GetQ2_EpX(Xbj,0.8);
 Double_t thetalim = GetQ2_thetaX(Xbj,45.);
 //cout << "thetalim = " << thetalim << endl;
 //cout << "elim = " << eplim << endl;
 if(TMath::Min(eplim,thetalim)>1.)result = TMath::Min(eplim,thetalim);
 return result;
}

Double_t GetQ2_thetaX(Double_t Xbj, Double_t theta){
 Double_t THETA  = 2*EBeam*TMath::Power(TMath::Sin(TMath::DegToRad()*theta/2.),2.);
 return 2.*EBeam*THETA*ProtonMass*Xbj/(THETA+ProtonMass*Xbj);
}

Double_t GetQ2_WX(Double_t Xbj, Double_t W){
 return (W*W-ProtonMass*ProtonMass)/(1./Xbj-1.);
}

Double_t GetQ2_EpX(Double_t Xbj, Double_t Ep){
 return 2.*ProtonMass*(EBeam-Ep)*Xbj;
}

/**********************************************************************************/
/*                             phi Model begins here                              */
/**********************************************************************************/

//const Double_t ProtonMass = 0.93827;
//const Double_t phimass = 1.01946;
//const Double_t Kmass = 0.49368;
double dvmpx(double t, double xb, double Q2, double PHI_G, double cos_theta_H, double E){
 //  dsigma/dQ2 dX dt dphi for ep-->ep phi
 if(xb<Q2/(2.0*ProtonMass*E)||xb>1.0)return 0.;
 double nu  = Q2/(2.*ProtonMass*xb);
 double y = nu/E;
 double e1 = TMath::Power(y*xb*ProtonMass,2)/Q2;
 double EPS = (1.0-y-e1)/(1-y+y*y/2+e1);
 if(EPS<0.||EPS>1.)return 0.;
 double W2  = ProtonMass*ProtonMass + 2.0*ProtonMass*nu - Q2;
 if(W2<TMath::Power(ProtonMass+phimass,2.))return 0.;
 double W = TMath::Sqrt(W2);

 double Wth = ProtonMass+phimass;
 double cT = 400*(1-Wth*Wth/W2)*TMath::Power(W,0.32);
 double sigmT = cT/TMath::Power(1+Q2/(phimass*phimass),3.0);
 double R = 0.4*Q2/(phimass*phimass);
 double sigmL = R*sigmT;

 double E1cm = ProtonMass*(ProtonMass + nu)/W;
 double P1cm = ProtonMass*TMath::Sqrt(nu*nu+Q2)/W;
 double E2cm = (W2 + ProtonMass*ProtonMass-phimass*phimass)/(2.*W);
 if(E2cm<ProtonMass)return 0.;
 double P2cm = TMath::Sqrt(E2cm*E2cm-ProtonMass*ProtonMass);
 double E3cm = (W2-phimass*phimass+ProtonMass*ProtonMass)/(2*W);
 double P3cm = TMath::Sqrt(E3cm*E3cm-ProtonMass*ProtonMass);

 double tmax = 2.0*(ProtonMass*ProtonMass - E1cm*E2cm + P1cm*P2cm);
 double tmin = 2.0*(ProtonMass*ProtonMass - E1cm*E2cm - P1cm*P2cm);
 if(t<tmin||t>tmax)return 0.;

 tmin = -(TMath::Power(Q2+phimass*phimass,2)/4./W2-TMath::Power(P1cm-P3cm,2));

 double T  = TMath::Abs(t);
 double T0 = TMath::Abs(tmin);
 double B = 2.2 + 4*0.24*TMath::Log(W);
 double FF = TMath::Exp(B*(T0-T));

 double r04_00 = EPS*R/(1+EPS*R);
 double CM_F = ( (1.-r04_00) + (3.*r04_00-1.)*cos_theta_H*cos_theta_H )*0.75;
 double phi_modul = ( 1 + 0.5*TMath::Cos(PHI_G))/1.5;

 phi_modul=phi_modul*phi_modul;


 double FLUX = alpha_em /(8*TMath::Pi()) * Q2/(ProtonMass*ProtonMass*E*E) * (1-xb)/(xb*xb*xb) * 1/(1-EPS);
 double DVMPX = FLUX * ( sigmT + EPS*sigmL ) * FF * CM_F ;
 if(DVMPX<0.) DVMPX=0.;
 return DVMPX;
}


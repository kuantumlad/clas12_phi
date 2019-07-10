#include <iostream>
#include <vector>
#include <map>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TMath.h>
#include <string>
#include <sstream>
#include <TLorentzVector.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TVector3.h>

int calcThetaCM(){

  double el_mass = 0.000511;
  double pr_mass = 0.93827;
  double kp_mass = 0.49368;
  double km_mass = 0.49368;
  double phi_mass = 1.01946;
  double beam_energy = 10.6;
  TLorentzVector lv_ebeam(0, 0, beam_energy, beam_energy);
  TLorentzVector lv_target(0,0,0,pr_mass);
  
  int helicity;
  double el_px, el_py, el_pz;
  double pr_px, pr_py, pr_pz;
  double kp_px, kp_py, kp_pz;
  double km_px, km_py, km_pz;

  TFile *fOut = new TFile("out.root","RECREATE");
  
  ifstream infile;
  infile.open("/Users/bclary/Documents/work/clas12/phi_analyis/data/all_final_state_particles.txt");
  //infile.open("/Users/bclary/Documents/work/clas12/phi_analyis/theory_notebooks/test.txt");

  TH1D *h_cos_theta_cm_kp = new TH1D("h_cos_theta_cm_kp","h_cos_theta_cm_kp",12,-1.1, 1.1);
  TH1D *h_w_cos_theta_cm_kp  = new TH1D("h_w_cos_theta_cm_kp","h_w_cos_theta_cm_kp",12,-1.0 , 1.0 );


  // using test.txt to test calculation methods:

    for(int ip=0;ip<3387;ip++){
    infile >> helicity >> el_px >> el_py >> el_pz
	   >> pr_px >> pr_py >> pr_pz
	   >> kp_px >> kp_py >> kp_pz
	   >> km_px >> km_py >> km_pz;
      

    TLorentzVector lv_el(el_px,el_py,el_pz,sqrt( pow(el_px,2) + pow(el_py,2) + pow(el_pz,2) + el_mass*el_mass));
    TLorentzVector lv_pr(pr_px,pr_py,pr_pz,sqrt( pow(pr_px,2) + pow(pr_py,2) + pow(pr_pz,2) + pr_mass*pr_mass));
    TLorentzVector lv_kp(kp_px,kp_py,kp_pz,sqrt( pow(kp_px,2) + pow(kp_py,2) + pow(kp_pz,2) + kp_mass*kp_mass));
    TLorentzVector lv_km(km_px,km_py,km_pz,sqrt( pow(km_px,2) + pow(km_py,2) + pow(km_pz,2) + km_mass*km_mass));
    TLorentzVector lv_phi;   
    lv_phi = lv_kp + lv_km;
    //if ( lv_phi.M() > 1.04 ) continue;


    ////// Calculate phi angle
    TLorentzVector lv_b_el(lv_el), lv_b_pr(lv_pr), lv_b_kp(lv_kp), lv_b_km(lv_km), lv_b_phi(lv_phi);
    TLorentzVector lv_vph = lv_ebeam - lv_el;

    TVector3 v_boost_vph = lv_vph.BoostVector();
    TVector3 v_boost_pr = lv_pr.BoostVector();
    TVector3 v_boost_kp = lv_kp.BoostVector();
    TVector3 v_boost_phi = (lv_kp + lv_km).BoostVector();
    
    TVector3 v_lepton = (lv_vph.Vect()).Cross(lv_b_el.Vect());
    TVector3 v_hadron = (lv_vph.Vect()).Cross(lv_b_phi.Vect());
    
    TVector3 n_lepton = v_lepton.Unit();
    TVector3 n_hadron = v_hadron.Unit();
    double c0 = v_lepton.Dot(lv_b_phi.Vect());
    double n0 = n_lepton.Dot(n_hadron);
    double n1 = n_lepton.Mag();
    double n2 = n_hadron.Mag();
    double n_ang = (180.0/3.14159265)*TMath::ACos(n0/(n1*n2)); // (c0/TMath::Abs(c0))* TMath::ACos(n0/(n1*n2));
    
    //kin quantities
    double q2 = 4*lv_ebeam.E()*lv_el.E()*sin(lv_el.Theta()/2.0)*sin(lv_el.Theta()/2.0);
    double xb = q2 / (2.0*pr_mass*(lv_ebeam.E() - lv_el.E()));
    double w = sqrt(-q2 + pr_mass*pr_mass + 2*pr_mass*(lv_ebeam.E() - lv_el.E()));
    double t = 2.0*pr_mass*(lv_pr.E() - pr_mass);    

    TVector3 v3L = (lv_ebeam.Vect()).Cross(lv_el.Vect());
    TVector3 v3H = (lv_pr.Vect()).Cross((lv_ebeam-lv_el).Vect());
    //phi in trento convention
    
    float phi = TMath::RadToDeg()* v3L.Angle(v3H);
    if(v3L.Dot(lv_vph.Vect())<0)phi = 360.-phi;
    std::cout << " phi " << phi << std::endl;

    /// first lets calculate the vector meson production angle (phi meson decays to KPKM)
    double cos_phi_meson_production_angle = lv_vph.Dot(lv_phi)/(lv_vph.Vect().Mag() * lv_phi.Vect().Mag());    
    std::cout << " cos phi meson production angle " << cos_phi_meson_production_angle << std::endl;

    // Now let us calculate the decay distribution of the charged positive Kaon based on the helicity coordinate system.
    //define the helicity frame based on Wolfe & Schiller (sp?)
    TLorentzVector lv_q = (lv_ebeam - lv_el);    

    ////////////////////////////////////////////////////////
    // setting up boosts using convention from FX
    TVector3 beta_q2 = (lv_q + lv_target).BoostVector();

    TVector3 Z_hel_unit = -(lv_pr.Vect()).Unit(); // do these need a boost? 
    TVector3 Y_hel_unit = ((lv_q).Vect().Cross(lv_phi.Vect())).Unit(); // need a boost?
    TVector3 X_hel_unit = (Y_hel_unit.Cross(Z_hel_unit));
    
    //set up boost vectors to boost measured Kp into phi COM frame
    TVector3 boost_phiCM;       
    TLorentzVector lv_kpCMPHI;
    TLorentzVector lv_kmCMPHI;
    lv_kpCMPHI=lv_kp;
    lv_kmCMPHI=lv_km;   
    boost_phiCM = lv_phi.BoostVector();   
    std::cout << " phi boost components " << boost_phiCM.Px() << " " << boost_phiCM.Py() << " " << boost_phiCM.Pz() << std::endl;
   
    lv_kpCMPHI.Boost(-boost_phiCM);
    lv_kmCMPHI.Boost(-boost_phiCM);

    // sanity check
    // this should give us kaon pairs with an angle between of pi, physically going back to back
    // must set boost_phiCM = -lv_phi.BoostVector()
    //double Theta_CM = lv_kpCMPHI.Vect().Angle( lv_kmCMPHI.Vect() ); 

    std::cout << " kaon plus components " << lv_kpCMPHI.Px() << " " << lv_kpCMPHI.Py() << " " << lv_kpCMPHI.Pz() << " " << std::endl;
    std::cout << " kaon minus components " << lv_kmCMPHI.Px() << " " << lv_kmCMPHI.Py() << " " << lv_kmCMPHI.Pz() << " " << std::endl;
    
    double cos_theta_cm_kp = TMath::Cos(lv_kpCMPHI.Theta());  // this is called Theta in FX's generator, technically it is cos(theta) variables
    std::cout << " cos_theta_cm_kp " << cos_theta_cm_kp << std::endl;
    
    double gg_angle = TMath::RadToDeg() * (lv_kp.Vect().Angle( lv_km.Vect() ) ); /// in the lab frame    
    std::cout << " gg_angle " << gg_angle << std::endl;

    double theta = TMath::RadToDeg()*TMath::ACos(cos_theta_cm_kp);  
    std::cout << " theta " << theta << std::endl;
    
    h_cos_theta_cm_kp->Fill(cos_theta_cm_kp);

    // now we can continue to determine other event related variables
    double nu = q2/(2.0*pr_mass*xb);
    double y = nu/lv_ebeam.E();
    double e1 = TMath::Power(y*xb*pr_mass,2)/q2;

    // photon polarization parameter 
    double EPS = (1.0 - y - e1)/(1 - y + y*y/2 +e1);

    //W2 and W
    double W2 = pr_mass*pr_mass + 2.0*pr_mass*nu - q2;
    double W = sqrt(W2);
    double Wth = pr_mass + phi_mass;

    //lets calculate tmin
    double e1cm = pr_mass*(pr_mass + nu)/W;
    double p1cm = pr_mass*TMath::Sqrt(nu*nu + q2)/W;
    double e2cm = (W2 + pr_mass*pr_mass - phi_mass*phi_mass)/(2.0*W);

    double p2cm = TMath::Sqrt(e2cm*e2cm - pr_mass*pr_mass);
    double e3cm = (W2-phi_mass*phi_mass + pr_mass*pr_mass)/(2.0*W);
    double p3cm = TMath::Sqrt(e3cm*e3cm - pr_mass*pr_mass);

    // first method to get the tmin and tmax values
    double tmin = 2.0*(pr_mass*pr_mass - e1cm*e2cm - p1cm*p2cm);    
    //double tmin2 = TMath::Power(q2 + phi_mass*phi_mass,2)/(4.0*W2) - TMath::Power(p1cm-p3cm,2);
    double tmax = 2.0*(pr_mass*pr_mass - e1cm*e2cm + p1cm*p2cm);
    //double tmax2 = TMath::Power(q2 + phi_mass*phi_mass,2)/(4.0*W2) - TMath::Power(p1cm+p3cm,2);
    
    std::cout <<" tmin " << tmin << std::endl;
    //std::cout <<" tmax " << tmax << " tmax2 " << tmax2 << std::endl;

    tmin = -(TMath::Power(q2+phi_mass*phi_mass,2)/4./W2-TMath::Power(p1cm-p3cm,2));
    
    // parameterization of the model
    double cT = 400*(1-Wth*Wth/W2)*TMath::Power(W,0.32);    
    double sigmaT = cT/TMath::Power(1+q2/(phi_mass*phi_mass),3.0);
    double R = 0.4*q2/(phi_mass*phi_mass); // assume indep of W
    double sigmaL = R*sigmaT;

    //exponential t-dep of cross section
    double T = TMath::Abs(t);
    double t0 = TMath::Abs(tmin);
    double B = 2.2 + 4.0*0.24*TMath::Log(W);
    double FF = TMath::Exp(B*(t0-T));

    double r04_00 = EPS*R/(1+EPS*R);
    double cm_f = ( (1.-r04_00) + (3.*r04_00-1.)*cos_theta_cm_kp*cos_theta_cm_kp)*0.75;

    std::cout << " event details " << std::endl;
    std::cout << "q2 " << q2 << std::endl;
    std::cout << "xb " << xb << std::endl;
    std::cout << "w " << w << std::endl;
    std::cout << "t " << t << std::endl;
    std::cout << "nu " << nu << std::endl;
    std::cout << "trento phi " << phi << std::endl;
    std::cout << "nu " << nu << std::endl;
    std::cout << "y " << y << std::endl;
    std::cout << "e1 " << e1 << std::endl;
    std::cout << "EPS " << EPS << std::endl;
    std::cout << "Wth " << Wth << std::endl;
    std::cout << "tmax " << tmax << std::endl;
    std::cout << "tmin " << tmin << std::endl;
    std::cout << "cT " << cT << std::endl;
    std::cout << "sigmaT " << sigmaT << std::endl;
    std::cout << "R " << R << std::endl;
    std::cout << "sigmaL " << sigmaL << std::endl;
    std::cout << "T " << T << std::endl;
    std::cout << "t0 " << t0 << std::endl;
    std::cout << "B " << B << std::endl;
    std::cout << "FF " << FF << std::endl;
    std::cout << "r04_00 " << r04_00 << std::endl;
    std::cout <<"cm_f " << cm_f << std::endl;
    



    // com function 
    
     

    
    }

    TF1 *f_cos_theta_cm_kp = new TF1("f_cos_theta_cm_kp","( (1.-[0]) + (3.*[0]-1.)*cos(x)*cos(x))*0.75", -2*3.14, 2*3.14); //-1.0, 1.0 );
    //const Double_t *buffer = h_cos_theta_cm_kp->GetBuffer();
    // number of entry is first entry in the buffer
    //int n = buffer[0];
    // when creating the data object it is important to create with the size of the data
    //ROOT::Fit::UnBinData data(n);
    //for (int i = 0; i < n; ++i){
    //  data.(buffer[2*i+1]);
    // }
    
    TCanvas *c1 = new TCanvas("c1","c1",900,900);
    c1->cd(0);
    h_cos_theta_cm_kp->Draw();
    h_cos_theta_cm_kp->Fit(f_cos_theta_cm_kp);
    return 0;
}

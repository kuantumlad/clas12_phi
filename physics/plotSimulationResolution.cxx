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

  for( int ss = 0; ss < 6; ss++ ){
    TProfile *temp_prof_epkpX = (TProfile*)fIn->Get((Form("prof_epkpX_pkp_sect%d",ss)));
    
						    

  }

  return 0;
}

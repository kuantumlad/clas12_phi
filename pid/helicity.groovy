package org.jlab.detector.helicity;
import org.jlab.clas.physics.Particle;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.fitter.DataFitter;
import org.jlab.groot.group.DataGroup;
import org.jlab.groot.math.F1D;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import org.jlab.detector.calib.utils.ConstantsManager;
import org.jlab.utils.groups.IndexedTable;
import java.io.IOException;
import java.io.PrintWriter;
import org.jlab.clas.physics.LorentzVector;
import org.jlab.clas.physics.Vector3;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.io.hipo.HipoDataSource;
import org.jlab.io.hipo.HipoDataSync;
import org.jlab.detector.decode.CodaEventDecoder;
import org.jlab.detector.decode.DetectorEventDecoder;
import org.jlab.detector.view.DetectorListener;
import org.jlab.groot.base.GStyle;
import org.jlab.groot.data.H2F;
import org.jlab.groot.data.TDirectory;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.io.base.DataEventType;
import org.jlab.io.task.DataSourceProcessorPane;
import org.jlab.io.task.IDataEventListener;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.groot.data.GraphErrors;
import org.jlab.groot.data.TDirectory;

import java.util.List;
import java.util.ArrayList;
import java.util.Arrays;

import org.jlab.jnp.hipo4.io.HipoReader;
import org.jlab.jnp.hipo4.data.Event;
import org.jlab.jnp.hipo4.data.Bank;
import org.jlab.jnp.hipo4.data.SchemaFactory;

/**
 * An example of reading the helicity flips, analyzing the sequence, and getting
 * the state for any event.
 * 
 * @author baltzell
 */
public class HelicityAnalysis {

         public static double eBeam = 10.646;
         public static double mTarget = 0.9327;
	 public static double plusSinusPredicted = 0;
	 public static double minusSinusPredicted = 0;
	 public static double plusSinusCounterPredicted = 0;
         public static double minusSinusCounterPredicted = 0;

  	 public static double plusSinusLevel3 = 0;
         public static double minusSinusLevel3 = 0;
         public static double plusSinusCounterLevel3 = 0;
         public static double minusSinusCounterLevel3 = 0;

  	 public static double plusSinusMeasured = 0;
         public static double minusSinusMeasured = 0;
         public static double plusSinusCounterMeasured = 0;
         public static double minusSinusCounterMeasured = 0;


	 public static double eventCounter = 0;

	 public static H1F wH;
         public static H1F xH;
	 public static H2F wq2H;
         public static H2F xq2H;
	 public static H1F q2H;


         public static H1F pxeH;
         public static H1F pyeH;
         public static H1F pzeH;

         public static H1F pxpipH;
         public static H1F pypipH;
         public static H1F pzpipH;

         public static H1F mm2H;
         public static H1F mm2HN;
         public static H1F mm2HNWithCut;

         public static H1F mm2HNFromThree;

         public static H1F mm2CutH;
         public static H1F thetaEH;
         public static H1F thetaPiPH;
         public static H1F thetaPiCMSH;
         public static H1F phiPiCMSH;
         public static H1F phiPiCMSHSmallBins;

 	 public static H1F thetaPiCMSHMinusMeasured;
         public static H1F thetaPiCMSHPlusMeasured;
         public static H1F phiPiCMSHMinusMeasured;
         public static H1F phiPiCMSHPlusMeasured;
         public static H1F phiPiCMSHRatioMeasured;


	 public static H1F thetaPiCMSHMinusPredicted;	
         public static H1F thetaPiCMSHPlusPredicted;
         public static H1F phiPiCMSHMinusPredicted;
	 public static H1F phiPiCMSHPlusPredicted;
         public static H1F phiPiCMSHRatioPredicted;

  	 public static H1F thetaPiCMSHMinusLevel3;
         public static H1F thetaPiCMSHPlusLevel3;
         public static H1F phiPiCMSHMinusLevel3;
	 public static H1F phiPiCMSHPlusLevel3;

	 

    /**
     * This reads tag=1 events for HEL::flip banks, and initializes and returns
     * a {@link HelicitySequenceDelayed} with delay set to zero.  The delay can
     * be changed later before a user tries to access the sequence.
     * 
     * @param filenames list of names of HIPO file to read
     * @return  unanalyzed sequence
     */
    public static HelicitySequenceDelayed readSequence(List<String> filenames) {
        
        HelicitySequenceDelayed seq=new HelicitySequenceDelayed(0);
        seq.setVerbosity(1);
       
        for (String filename : filenames) {

            HipoReader reader = new HipoReader();
            reader.setTags(1);
            reader.open(filename);
        
            SchemaFactory schema = reader.getSchemaFactory();
        
            while (reader.hasNext()) {
            
                Event event=new Event();
                Bank flipBank=new Bank(schema.getSchema("HEL::flip"));
            
                reader.nextEvent(event);
                event.read(flipBank);
            
                if (flipBank.getRows()<1) continue;
        
                seq.addState(HelicityState.createFromFlipBank(flipBank));
            }

            reader.close();
        }
        
        return seq;
    }

public static void ana(){
      for (int t = 0; t < 20; t++){
      System.out.println("ratio " + phiPiCMSHPlusMeasured.getBinContent(t + 1) + phiPiCMSHMinusMeasured.getBinContent(t + 1));
      if (phiPiCMSHPlusMeasured.getBinContent(t + 1) + phiPiCMSHMinusMeasured.getBinContent(t + 1) > 0)
      	  phiPiCMSHRatioMeasured.setBinContent(t + 1, (phiPiCMSHPlusMeasured.getBinContent(t + 1) - phiPiCMSHMinusMeasured.getBinContent(t + 1))/(phiPiCMSHPlusMeasured.getBinContent(t + 1) + phiPiCMSHMinusMeasured.getBinContent(t + 1))); 
     if (phiPiCMSHPlusPredicted.getBinContent(t + 1) + phiPiCMSHMinusPredicted.getBinContent(t + 1) > 0)
          phiPiCMSHRatioPredicted.setBinContent(t + 1, (phiPiCMSHPlusPredicted.getBinContent(t + 1) - phiPiCMSHMinusPredicted.getBinContent(t + 1))/(phiPiCMSHPlusPredicted.getBinContent(t + 1) + phiPiCMSHMinusPredicted.getBinContent(t + 1)));
      }
}


public static void plot(){


int runNumber= 5038;	    
       EmbeddedCanvas momCanvas = new EmbeddedCanvas();
        momCanvas.setSize(2400,600);
        momCanvas.divide(3, 1);

        momCanvas.setAxisTitleSize(14);
        momCanvas.setAxisFontSize(14);
        momCanvas.setTitleSize(14);
        momCanvas.cd(0);
	momCanvas.draw(pxeH);
	momCanvas.cd(1);
	momCanvas.draw(pyeH);
        momCanvas.cd(2);
	momCanvas.draw(pzeH);
	momCanvas.save("momCanvasE" + runNumber + ".png");


	momCanvas.cd(0);
	momCanvas.draw(wH);
        momCanvas.cd(1);
        momCanvas.draw(wq2H);
	momCanvas.cd(2);
	momCanvas.draw(q2H);
        momCanvas.save("wQ2CanvasE" + runNumber + ".png");


	momCanvas.cd(0);
        momCanvas.draw(xH);
        momCanvas.cd(1);
        momCanvas.draw(xq2H);
        momCanvas.cd(2);
	momCanvas.save("xQ2CanvasE" + runNumber + ".png");


	momCanvas.cd(0);
        momCanvas.draw(pxpipH);
        momCanvas.cd(1);
	momCanvas.draw(pypipH);
        momCanvas.cd(2);
	momCanvas.draw(pzpipH);
	momCanvas.save("momCanvasPIP" + runNumber + ".png");

        momCanvas.cd(0);
        momCanvas.draw(thetaEH);
        momCanvas.cd(1);
	momCanvas.draw(thetaPiPH);
	momCanvas.cd(2);
	momCanvas.draw(mm2H);
        momCanvas.draw(mm2HN, "same");
        momCanvas.draw(mm2HNFromThree, "same");
        momCanvas.draw(mm2HNWithCut, "same");


	momCanvas.save("momCanvasMM2" + runNumber + ".png");


	momCanvas.cd(0);
	momCanvas.draw(thetaPiCMSH);
        momCanvas.cd(1);
        momCanvas.draw(phiPiCMSH);
        phiPiCMSHSmallBins.setLineColor(5);
	momCanvas.draw(phiPiCMSHSmallBins, "same");
        momCanvas.cd(2);
	momCanvas.draw(mm2H);
        momCanvas.save("momCanvasCMS" + runNumber + ".png");



	momCanvas.cd(0);
	momCanvas.draw(thetaPiCMSH);
        momCanvas.draw(thetaPiCMSHMinusPredicted, "same");
        momCanvas.draw(thetaPiCMSHPlusPredicted, "same");

	momCanvas.cd(1);
	momCanvas.draw(phiPiCMSH);
	momCanvas.draw(phiPiCMSHMinusPredicted, "same");
	momCanvas.draw(phiPiCMSHPlusPredicted, "same");

	momCanvas.cd(2);
	momCanvas.draw(mm2CutH);
	momCanvas.save("momCanvasPiCMSPredicted" + runNumber + ".png");

 	momCanvas.cd(0);
        momCanvas.draw(thetaPiCMSH);
        momCanvas.draw(thetaPiCMSHMinusMeasured, "same");
        momCanvas.draw(thetaPiCMSHPlusMeasured, "same");



        momCanvas.cd(1);
        momCanvas.draw(phiPiCMSH);
        momCanvas.draw(phiPiCMSHMinusMeasured, "same");
        momCanvas.draw(phiPiCMSHPlusMeasured, "same");

        momCanvas.cd(2);
        momCanvas.draw(mm2CutH);
        momCanvas.save("momCanvasPiCMSMeasured" + runNumber + ".png");

  	momCanvas.cd(0);
        momCanvas.draw(thetaPiCMSH);
        momCanvas.draw(thetaPiCMSHMinusLevel3, "same");
        momCanvas.draw(thetaPiCMSHPlusLevel3, "same");

        momCanvas.cd(1);
        momCanvas.draw(phiPiCMSH);
        momCanvas.draw(phiPiCMSHMinusLevel3, "same");
        momCanvas.draw(phiPiCMSHPlusLevel3, "same");

        momCanvas.cd(2);
        momCanvas.draw(mm2CutH);
        momCanvas.save("momCanvasPiCMSLevel3" + runNumber + ".png");



        momCanvas.cd(0);
        thetaPiCMSHMinusMeasured.setLineColor(3);
        momCanvas.draw(thetaPiCMSHMinusMeasured);
	thetaPiCMSHMinusPredicted.setLineColor(2);
        momCanvas.draw(thetaPiCMSHMinusPredicted, "same");
        momCanvas.draw(thetaPiCMSHMinusLevel3, "same");

	momCanvas.cd(1);
	thetaPiCMSHPlusMeasured.setLineColor(3);
   	momCanvas.draw(thetaPiCMSHPlusMeasured);
	thetaPiCMSHPlusPredicted.setLineColor(2);
   	momCanvas.draw(thetaPiCMSHPlusPredicted, "same");
   	momCanvas.draw(thetaPiCMSHPlusLevel3, "same");

        momCanvas.save("momCanvasPiCMSHelicityComparedTheta" + runNumber + ".png");

    	momCanvas.cd(0);
	phiPiCMSHMinusMeasured.setLineColor(3);
    	momCanvas.draw(phiPiCMSHMinusMeasured);
	phiPiCMSHMinusPredicted.setLineColor(2);
    	momCanvas.draw(phiPiCMSHMinusPredicted, "same");
    	momCanvas.draw(phiPiCMSHMinusLevel3, "same");

        momCanvas.cd(1);
	phiPiCMSHPlusMeasured.setLineColor(3);
        momCanvas.draw(phiPiCMSHPlusMeasured);
	phiPiCMSHPlusPredicted.setLineColor(2);
        momCanvas.draw(phiPiCMSHPlusPredicted, "same");
        momCanvas.draw(phiPiCMSHPlusLevel3, "same");

        momCanvas.save("momCanvasPiCMSHelicityComparedPhi" + runNumber + ".png");

 	momCanvas.cd(0);
        momCanvas.draw(phiPiCMSHRatioMeasured);
   
        momCanvas.cd(1);
        momCanvas.draw(phiPiCMSHRatioPredicted);


        momCanvas.save("momCanvasPiCMSRatioPhi" + runNumber + ".png");






}

public static double returnX(Particle e, double eBeam, double mTarget, double q2) {
       double pEl = Math.sqrt(e.vector().px() * e.vector().px() + e.vector().py() * e.vector().py() + e.vector().pz() * e.vector().pz());
       double x = q2 / (2 * (eBeam - pEl));
       return x;
}


public static double returnW(double mTarget, double q2, double x) {
       double w = Math.sqrt(mTarget * mTarget - q2 + q2 / x);
       return w;
}

public static double returnQ2(Particle e, double eBeam, double mTarget) {
       double pEl = Math.sqrt(e.vector().px() * e.vector().px() + e.vector().py() * e.vector().py() + e.vector().pz() * e.vector().pz());
       double q2 = (eBeam - pEl) * (eBeam - pEl) - e.vector().px() * e.vector().px() - e.vector().py() * e.vector().py() - (eBeam - e.vector().pz()) * (eBeam - e.vector().pz());
       return -q2;
}

    public static LorentzVector returnMissing(Particle e, Particle pip, double eBeam, double mTarget){
       LorentzVector vecE = new LorentzVector(e.vector().px(), e.vector().py(), e.vector().pz(), Math.sqrt(e.p()*e.p() + 0.000511*0.000511));
       LorentzVector vecPip = new LorentzVector(pip.vector().px(), pip.vector().py(), pip.vector().pz(),  Math.sqrt(pip.p()*pip.p() + 0.140*0.140));
       LorentzVector beam = new LorentzVector(0.0,0.0,eBeam,eBeam);
       LorentzVector target = new LorentzVector(0.0,0.0,0.0,mTarget);
       LorentzVector miss = new LorentzVector(0,0,0,0);
       miss.add(beam);
       miss.add(target);
       miss.sub(vecE);
       miss.sub(vecPip);
       return miss;
}

    public static LorentzVector returnMissingFromThree(Particle e, Particle pip, Particle n, double eBeam, double mTarget){
       LorentzVector vecE = new LorentzVector(e.vector().px(), e.vector().py(), e.vector().pz(), Math.sqrt(e.p()*e.p() + 0.000511*0.000511));
       LorentzVector vecPip = new LorentzVector(pip.vector().px(), pip.vector().py(), pip.vector().pz(),  Math.sqrt(pip.p()*pip.p() + 0.140*0.140));
       LorentzVector vecN = new LorentzVector(n.vector().px(), n.vector().py(), n.vector().pz(),  Math.sqrt(n.p()*n.p() + 0.9396*0.9396));

       LorentzVector beam = new LorentzVector(0.0,0.0,eBeam,eBeam);
       LorentzVector target = new LorentzVector(0.0,0.0,0.0,mTarget);
       LorentzVector miss = new LorentzVector(0,0,0,0);
       miss.add(beam);
       miss.add(target);
       miss.sub(vecE);
       miss.sub(vecPip);
       miss.sub(vecN);
       return miss;
}



   public static LorentzVector returnResonance(Particle e, Particle pip, double eBeam, double mTarget, LorentzVector missingVector){
       LorentzVector vecE = new LorentzVector(e.vector().px(), e.vector().py(), e.vector().pz(), Math.sqrt(e.p()*e.p() + 0.000511*0.000511));
       LorentzVector vecPip = new LorentzVector(pip.vector().px(), pip.vector().py(), pip.vector().pz(),  Math.sqrt(pip.p()*pip.p() + 0.140*0.140));
       LorentzVector beam = new LorentzVector(0.0,0.0,eBeam,eBeam);
       LorentzVector target = new LorentzVector(0.0,0.0,0.0,mTarget);
       LorentzVector resonance = new LorentzVector(0,0,0,0);
       resonance.add(beam);
       resonance.add(target);
       resonance.sub(vecE);
       return resonance;
}
   public static LorentzVector labToRes(LorentzVector h, Particle  pip){
       LorentzVector pion = new LorentzVector(pip.vector().px(), pip.vector().py(), pip.vector().pz(),  Math.sqrt(pip.p()*pip.p() + 0.140*0.140));
       Vector3 boostVector =  h.boostVector();
//       Vector3 boostVector = new Vector3(0, 0, h.vector().p()/h.vector.e());
//       pion.rotateZ(h.phi() + 3.1415926);
//       pion.rotateY(-h.theta());
       pion.boost(boostVector);
       return pion;
}


    public static void processEvent(Event event, SchemaFactory schema, HelicityBit level3, HelicityBit predicted, HelicityBit measured){
    	   Bank recBankPart = new Bank(schema.getSchema("REC::Particle"));
	   event.read(recBankPart);
	   eventCounter++;
	        for (int loopE = 0; loopE < 1; loopE++) {
	     	  if (recBankPart.getInt("pid", loopE) == 11) {
			    Particle electron = null;
                            electron = new Particle (
			    	11,
                                0.000511,
				(byte) -1,
				recBankPart.getFloat("px", loopE),
				recBankPart.getFloat("py", loopE),
                                recBankPart.getFloat("pz", loopE),
                                recBankPart.getFloat("vx", loopE),
                                recBankPart.getFloat("vy", loopE),
				recBankPart.getFloat("vz", loopE));
                            	pxeH.fill(electron.vector().px());
                                pyeH.fill(electron.vector().py());
                            	pzeH.fill(electron.vector().pz());
                      		for (int loopP = 0; loopP < recBankPart.getRows(); loopP++) {
                              	if (recBankPart.getInt("pid", loopP) == 211) {
                               Particle pion = null;
                               pion = new Particle (
                                211,
                                0.140,
                                (byte) 1,
                                recBankPart.getFloat("px", loopP),
                                recBankPart.getFloat("py", loopP),
                                recBankPart.getFloat("pz", loopP),
                                recBankPart.getFloat("vx", loopP),
                                recBankPart.getFloat("vy", loopP),
                                recBankPart.getFloat("vz", loopP));
                               pxpipH.fill(pion.vector().px());
                               pypipH.fill(pion.vector().py());
                               pzpipH.fill(pion.vector().pz());
			       for (int loopN = 0; loopN < recBankPart.getRows(); loopN++) {
			       	   	      if (recBankPart.getInt("pid", loopN) == 2112 && electron.p() > 2.0 && pion.p() > 2.0) {
					       LorentzVector missingVectorN = returnMissing(electron, pion, eBeam, mTarget);
                                               mm2HN.fill(missingVectorN.mass2());
					           Particle neutron = null;
                              			    neutron = new Particle (
			       	      	  	    2112,
                                		    0.9396,
                                		    (byte) 1,
                                		    recBankPart.getFloat("px", loopN),
                                		    recBankPart.getFloat("py", loopN),
                                		    recBankPart.getFloat("pz", loopN),
						    recBankPart.getFloat("vx", loopN),
                                		    recBankPart.getFloat("vy", loopN),
                                		    recBankPart.getFloat("vz", loopN));
                                                    LorentzVector missingVectorFromThreeN = returnMissingFromThree(electron, pion, neutron, eBeam, mTarget);
					            mm2HNFromThree.fill(missingVectorFromThreeN.mass2());
						    if (missingVectorFromThreeN.mass2() > -0.2 && missingVectorFromThreeN.mass2() <0.2){
						    mm2HNWithCut.fill(missingVectorN.mass2());
						    }
					       }
				}
                               if (electron.p() > 2.0 && pion.p() > 2.0){
                                  thetaEH.fill(Math.toDegrees(electron.theta()));
                                  thetaPiPH.fill(Math.toDegrees(pion.theta()));
                                  if (Math.toDegrees(electron.theta()) > 10 && Math.toDegrees(electron.theta()) < 35 && Math.toDegrees(pion.theta()) < 35){
                                     LorentzVector missingVector = returnMissing(electron, pion, eBeam, mTarget);
                                     mm2H.fill(missingVector.mass2());
                                     if (missingVector.mass2() < 1.1 && missingVector.mass2() > 0.8){
                                     LorentzVector resonanceVector = returnResonance(electron, pion, eBeam, mTarget, missingVector);
                                     LorentzVector pionCMS = labToRes(resonanceVector, pion);
				     mm2CutH.fill(missingVector.mass2());
                                     thetaPiCMSH.fill(Math.cos(pionCMS.theta()));
                                     phiPiCMSH.fill(Math.toDegrees(pionCMS.phi()));
                                     phiPiCMSHSmallBins.fill(Math.toDegrees(pionCMS.phi()));

				     double Q2 = returnQ2(electron, eBeam, mTarget);
                                     double X = returnX(electron, eBeam, mTarget, Q2);
                                     double W = returnW(mTarget, Q2, X);

				     wH.fill(W);
				     q2H.fill(Q2);
				     xH.fill(X);
				     xq2H.fill(X, Q2);
				     wq2H.fill(W, Q2);
                                     if (predicted.value==1) {
                                        plusSinusPredicted = plusSinusPredicted + Math.sin(pionCMS.phi());
                                        plusSinusCounterPredicted++;
                                        phiPiCMSHPlusPredicted.fill(Math.toDegrees(pionCMS.phi()));
                                        thetaPiCMSHPlusPredicted.fill(Math.cos(pionCMS.theta()));
                                     
				     }
                                     if (predicted.value==-1){
                                      phiPiCMSHMinusPredicted.fill(Math.toDegrees(pionCMS.phi()));
                                      thetaPiCMSHMinusPredicted.fill(Math.cos(pionCMS.theta()));
                                      minusSinusPredicted = minusSinusPredicted + Math.sin(pionCMS.phi());
                                      minusSinusCounterPredicted++;
                                     }


                                     if (level3.value==1) {
                                        plusSinusLevel3 = plusSinusLevel3 + Math.sin(pionCMS.phi());
                                        plusSinusCounterLevel3++;
                                        phiPiCMSHPlusLevel3.fill(Math.toDegrees(pionCMS.phi()));
                                        thetaPiCMSHPlusLevel3.fill(Math.cos(pionCMS.theta()));
                                     }
                                     if (level3.value==-1){
                                      phiPiCMSHMinusLevel3.fill(Math.toDegrees(pionCMS.phi()));
                                      thetaPiCMSHMinusLevel3.fill(Math.cos(pionCMS.theta()));
                                      minusSinusLevel3 = minusSinusLevel3 + Math.sin(pionCMS.phi());
                                      minusSinusCounterLevel3++;
                                     }


				     if (measured.value==1) {
				     	plusSinusMeasured = plusSinusMeasured + Math.sin(pionCMS.phi());
					plusSinusCounterMeasured++;
					phiPiCMSHPlusMeasured.fill(Math.toDegrees(pionCMS.phi()));
					thetaPiCMSHPlusMeasured.fill(Math.cos(pionCMS.theta()));
				     }
                                     if (measured.value==-1){
				      phiPiCMSHMinusMeasured.fill(Math.toDegrees(pionCMS.phi()));
                                      thetaPiCMSHMinusMeasured.fill(Math.cos(pionCMS.theta()));
				      minusSinusMeasured = minusSinusMeasured + Math.sin(pionCMS.phi());
				      minusSinusCounterMeasured++;
                                     }

                                  }
                               }
                            }

                }
        }
        }
        }
}
    public static void main(String[] args) {
        

wH = new H1F("wH", 250, 0, 5);
wH.setTitleX("W, GeV")
xH = new H1F("xH", 250, 0, 1);
xH.setTitleX("X");
q2H = new H1F("q2H", 150, 0, 10);
q2H.setTitleX("Q2, GeV2");
wq2H = new H2F("wq2H", 250, 0, 5, 250, 0, 10);
wq2H.setTitleY("Q2, GeV2");
wq2H.setTitleX("W, GeV");

xq2H = new H2F("xq2H", 250, 0, 1, 250, 0, 10);
xq2H.setTitleX("X");
xq2H.setTitleY("Q2, GeV2");


	
pxeH = new H1F("pxeH", 1000, -2, 2);
pyeH = new H1F("pyeH", 1000, -2, 2);
pzeH = new H1F("pzeH", 1000, 0, 10);
pxeH.setOptStat(1111);
pxpipH = new H1F("pxpipH", 1000, -2, 2);
pypipH = new H1F("pypipH", 1000, -2, 2);
pzpipH = new H1F("pzpipH", 1000, 0, 10);
pxpipH.setOptStat(1111);
mm2H = new H1F("mm2H", 780, -1, 12);
mm2H.setOptStat(1111);

mm2HN = new H1F("mm2HN", 780, -1, 12);
mm2HN.setOptStat(1111);
mm2HN.setLineColor(2);

mm2HNWithCut = new H1F("mm2HNWithCut", 780, -1, 12);
mm2HNWithCut.setOptStat(1111);
mm2HNWithCut.setLineColor(7);


mm2HNFromThree = new H1F("mm2HN", 780, -1, 12);
mm2HNFromThree.setOptStat(1111);
mm2HNFromThree.setLineColor(4);



mm2CutH = new H1F("mm2CutH", 500, -0.5, 2);
mm2CutH.setOptStat(1111);

thetaEH = new H1F("thetaEH", 100, 0, 40);
thetaPiPH = new H1F("thetaPiPH", 100, 0, 70);


thetaPiCMSH = new H1F("thetaPiCMSH", 250, 0, 1.1);
phiPiCMSH = new H1F("phiPiCMSH", 20, -200, 200);

phiPiCMSHSmallBins = new H1F("phiPiCMSHSmallBins", 100, -200, 200);

thetaPiCMSHMinusPredicted = new H1F("thetaPiCMSHMinusPredicted", 250, 0.5, 1.1);
thetaPiCMSHPlusPredicted = new H1F("thetaPiCMSHPlusPredicted", 250, 0.5, 1.1);
thetaPiCMSHPlusPredicted.setLineColor(2);
thetaPiCMSHMinusPredicted.setLineColor(4);

phiPiCMSHMinusPredicted = new H1F("phiPiCMSHMinusPredicted", 20, -200, 200);
phiPiCMSHPlusPredicted = new H1F("phiPiCMSHPlusPredicted", 20, -200, 200);
phiPiCMSHPlusPredicted.setLineColor(2);
phiPiCMSHMinusPredicted.setLineColor(4);

phiPiCMSHRatioPredicted = new H1F("phiPiCMSHRatioPredicted", 20, -200, 200);


thetaPiCMSHMinusMeasured = new H1F("thetaPiCMSHMinusMeasured", 250, 0.5, 1.1);
thetaPiCMSHPlusMeasured = new H1F("thetaPiCMSHPlusMeasured", 250, 0.5, 1.1);
thetaPiCMSHPlusMeasured.setLineColor(2);
thetaPiCMSHMinusMeasured.setLineColor(4);

phiPiCMSHMinusMeasured = new H1F("phiPiCMSHMinusMeasured", 20, -200, 200);
phiPiCMSHPlusMeasured = new H1F("phiPiCMSHPlusMeasured", 20, -200, 200);
phiPiCMSHPlusMeasured.setLineColor(2);
phiPiCMSHMinusMeasured.setLineColor(4);
phiPiCMSHRatioMeasured = new H1F("phiPiCMSHRatioMeasured", 20, -200, 200);



thetaPiCMSHMinusLevel3 = new H1F("thetaPiCMSHMinusLevel3", 250, 0.5, 1.1);
thetaPiCMSHPlusLevel3 = new H1F("thetaPiCMSHPlusLevel3", 250, 0.5, 1.1);
thetaPiCMSHPlusLevel3.setLineColor(2);
thetaPiCMSHMinusLevel3.setLineColor(4);

phiPiCMSHMinusLevel3 = new H1F("phiPiCMSHMinusLevel3", 20, -200, 200);
phiPiCMSHPlusLevel3 = new H1F("phiPiCMSHPlusLevel3", 20, -200, 200);
phiPiCMSHPlusLevel3.setLineColor(2);
phiPiCMSHMinusLevel3.setLineColor(4);



//cooked
//        final String dir="/volatile/clas12/rg-a/production/recon/calib/v5/005030/";
//        final String file="calib_clas_005030.evio.00240-00244.hipo";
//decoded
	final String dir="/cache/clas12/rg-a/production/decoded/6b.2.0/005030/";
	final String file="clas_005030.evio.00240-00244.hipo";
        List<String> filenames=new ArrayList<>();
        if (args.length>0) filenames.addAll(Arrays.asList(args));
        else               filenames.add(dir+file);
       
        // initialize a sequence from tag=1 events:
        HelicitySequenceDelayed seq = HelicityAnalysis.readSequence(filenames);
        final boolean integrity = seq.analyze();
        if (!integrity) {
            System.err.println("\n\n######### OOPS\n\n");
            // We may want to investigate further, or discard events.
        }

        // print the sequence:
        seq.show();

        // set the appropriate delay for this data:
        seq.setDelay(8);

        // now read the full events, e.g. during a normal physics analysis: 
        int nevents=0;
        int nflips=0;
       
        for (String filename : filenames) {

            HipoReader reader = new HipoReader();
            reader.setTags(0);
            reader.open(filename);
            
            SchemaFactory schema = reader.getSchemaFactory();
        
            while (reader.hasNext()) {

                nevents++;
                Bank flipBank=new Bank(schema.getSchema("HEL::flip"));
                Bank rcfgBank=new Bank(schema.getSchema("RUN::config"));
                Bank onliBank=new Bank(schema.getSchema("HEL::online"));
               
                Event event=new Event();
                reader.nextEvent(event);
                event.read(flipBank);
                event.read(rcfgBank);
                event.read(onliBank);

                // just to curtail printouts:
                if (flipBank.getRows()>0) nflips++;
                if (nflips<240) continue;
                if (nevents==5000000) break;
         
                long timestamp = -1;
                HelicityBit level3 = HelicityBit.UDF;
                HelicityBit predicted = HelicityBit.UDF;
                HelicityBit measured = HelicityBit.UDF;
   
                // Get RUN::config.timestamp for this event:
                if (rcfgBank.getRows()>0) 
                    timestamp = rcfgBank.getLong("timestamp",0);

                // Get HEL::online.rawHelicity, the online, delay-corrected
                // helicity for this event (if available):
                if (onliBank.getRows()>0)
                    level3 = HelicityBit.create(onliBank.getByte("helicity",0));

                // Get the offline, delay-corrected helicity for this event based
                // on the measured sequence.  If this timestamp is outside the
                // measured range, the bit will be null. 
                if (seq.find(timestamp)!=null)
                    measured = seq.find(timestamp);

                // Same as previous, except use the pseudo-random generator's
                // prediction which provides for bits later than the measured range.
                // For example, the last N windows in a given file are measured in
                // the next file (or not at all if it's the last file in a run),
                // so will only be accessible with the generator.  If you try to 
                // use a timestamp before the measured sequence, the generator will
                // return null.
                if (seq.findPrediction(timestamp)!=null)
                   predicted = seq.findPrediction(timestamp);


//		HelicityAnalysis.processEvent(event, schema, level3, predicted, measured);
                System.out.println(String.format("%d %5d L3/Predict/Measured = %6s%6s%6s",
                        timestamp,nflips,level3,predicted,measured));
            }

            reader.close();
            HelicityAnalysis.ana();
	    HelicityAnalysis.plot();


	    System.out.println("positive predicted: " + plusSinusPredicted/plusSinusCounterPredicted + " counter " + plusSinusCounterPredicted + " total plus sinus " + plusSinusPredicted);
            System.out.println("negative predicted: " + minusSinusPredicted/minusSinusCounterPredicted  + " counter " + minusSinusCounterPredicted + " total minus sinus " + minusSinusPredicted) ;

	    System.out.println("positive measured: " + plusSinusMeasured/plusSinusCounterMeasured + " counter " + plusSinusCounterMeasured + " total plus sinus " + plusSinusMeasured);
            System.out.println("negative measured: " + minusSinusMeasured/minusSinusCounterMeasured  + " counter " + minusSinusCounterMeasured + " total minus sinus " + minusSinusMeasured) ;


	    System.out.println("positive level3: " + plusSinusLevel3/plusSinusCounterLevel3 + " counter " + plusSinusCounterLevel3 + " total plus sinus " + plusSinusLevel3);
            System.out.println("negative level3: " + minusSinusLevel3/minusSinusCounterLevel3  + " counter " + minusSinusCounterLevel3 + " total minus sinus " + minusSinusLevel3) ;


        }
	
    }
}
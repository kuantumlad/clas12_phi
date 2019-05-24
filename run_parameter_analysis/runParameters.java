import org.jlab.clas.analysis.clary.*;

import java.io.*;
import java.util.*;

import com.google.gson.*;
import org.json.*;
import org.jlab.jnp.utils.json.*;

import org.jlab.io.hipo.HipoDataSource;
import org.jlab.jnp.hipo.io.HipoWriter;
import org.jlab.clas.physics.LorentzVector;
import org.jlab.clas.physics.Particle;

import org.jlab.io.base.DataEvent;
import org.jlab.io.base.DataBank;

import org.jlab.groot.group.*;
import org.jlab.groot.ui.TBrowser;
import org.jlab.groot.tree.*;
import org.jlab.groot.data.TDirectory;
import org.jlab.groot.data.*; //H1F H2F GraphErrors
import org.jlab.groot.math.*;
import org.jlab.groot.fitter.*;
import org.jlab.jnp.hipo.io.*;
import org.jlab.groot.graphics.EmbeddedCanvas;



public class runParameters{
    
   

    public static void main( String[] args){

	int run = 3050;
	double beam_energy = 6.4;
	String s_run = Integer.toString(run);
	String s_run_number = s_run; //TOO LAZY TO CHANGE IT

	TDirectory d=new TDirectory();
	String inFile = args[0];
	d.readFile(inFile);	
	d.cd();

	H2F h2_temp = (H2F)d.getObject("/cutlvls/h2_el_phip/h_3050_el_phip_cutlvl0");
	
	EmbeddedCanvas c_temp = new EmbeddedCanvas();
	c_temp.setSize(800,800);
	c_temp.draw(h2_temp,"colz");
	c_temp.save("testing.png");

	System.out.println(">> DONE " );
	Vector< GraphErrors > g_sf_sect_means = new Vector< GraphErrors >();
	Vector< GraphErrors > g_sf_sect_sigmas = new Vector< GraphErrors >();
	HashMap<Integer, ArrayList<H1F> > m_el_sect_slices_ectotp = new HashMap<Integer, ArrayList<H1F> >();
	

	for( int i = 0; i < 6; i++ ){
	    g_sf_sect_means.add( new GraphErrors() );
	    g_sf_sect_sigmas.add( new GraphErrors() );	    
	}

	for( int i = 0; i < 6; i++ ){
	    m_el_sect_slices_ectotp.put(i, new ArrayList<H1F>() );	   
	}

	
	Vector< Vector<H2F> > h2_el_sect_ectotp = new Vector< Vector<H2F> >();
	
           
	for( int s = 0; s < 6; s++ ){
	    Vector<H2F> v_htemp = new Vector<H2F>();
	    for( int c = 0; c <= 10; c++ ){
		H2F sf_temp = (H2F)d.getObject("/cutlvls/h2_el_sect_ectotp/h2_"+s_run+"_el_"+Integer.toString(s)+"_ectotp_cutlvl"+Integer.toString(c) );
		v_htemp.add(sf_temp);
	    }	    
	    h2_el_sect_ectotp.add( v_htemp );
	}

	
	

	SliceNFit( h2_el_sect_ectotp, g_sf_sect_means, g_sf_sect_sigmas, m_el_sect_slices_ectotp );

	Vector<GraphErrors> g_plus_threesig = new Vector<GraphErrors>();
	Vector<GraphErrors> g_minus_threesig = new Vector<GraphErrors>();
	
	for( int s = 0 ; s < 6; s++ ){
	    GraphErrors g_temp = g_sf_sect_means.get(s);	    
	    GraphErrors g_sig = g_sf_sect_sigmas.get(s); 
	    g_plus_threesig.add( new GraphErrors() );
	    g_minus_threesig.add( new GraphErrors() );
	    for( int b = 0; b < g_temp.getVectorX().size(); b++){		
		//System.out.println(" >> X " + g_temp.getDataX(b) + " Y " + g_temp.getDataY(b) + " ERR " + 3*g_sig.getDataY(b) );
		double y_sig = 2.75*(Math.abs(g_sig.getDataY(b)));
		
		g_plus_threesig.get(s).addPoint(g_temp.getDataX(b), g_temp.getDataY(b) + y_sig, g_temp.getDataEX(b), g_temp.getDataEY(b) );
		g_minus_threesig.get(s).addPoint(g_temp.getDataX(b),  g_temp.getDataY(b) - y_sig, g_sig.getDataEX(b), g_sig.getDataEY(b) );		    
	    }		
	}
	
       	EmbeddedCanvas c_g_meanfits = new EmbeddedCanvas();
	c_g_meanfits.setSize(1200,900);
	c_g_meanfits.divide(2,3);
	System.out.println(">> CHECKING G_SF SIZE BCLARY " + g_sf_sect_means.size() );

	///////////////////////////////////////////////////////
	//CREATE MAP FOR STORING FIT PARAMETERS TO PUSH TO JSON
	Map<Integer, ArrayList<Double> > m_maxfit_para = new HashMap<Integer, ArrayList<Double> >();	
	Map<Integer, ArrayList<Double> > m_minfit_para = new HashMap<Integer, ArrayList<Double> >();	
	ArrayList<Double> temp_maxpara = new ArrayList<Double>();
	ArrayList<Double> temp_minpara = new ArrayList<Double>();

	for( int s = 0; s < 6; s++ ){
	    c_g_meanfits.cd(s);
	    
	    Vector<H2F> v_temp = h2_el_sect_ectotp.get(s);
	    H2F h_temp = v_temp.get(4); //GET THE HISTOGRAM AFTER FIDUCIAL CUTS ARE APPLIED AND FIT THAT
	    GraphErrors g_mean2 = g_sf_sect_means.get(s);
	    GraphErrors g_plus = g_plus_threesig.get(s);
	    GraphErrors g_minus= g_minus_threesig.get(s);

	    for( int b = 0; b < g_mean2.getVectorX().size(); b++ ){
 		System.out.println(" >> " + g_mean2.getDataX(b) + " " + g_mean2.getDataY(b) );
	    }
	    
	    F1D f_sf_mean2 = new F1D("f_sf_mean2","[a] + [b]/(x-0.1)", 0.8, beam_energy);                                                                                                        
	    DataFitter.fit(f_sf_mean2, g_mean2,"REQ");

	    F1D f_sf_top = new F1D("f_sf_top","[a] + [b]/x + [c]/(x*x)", 0.8, beam_energy);
	    DataFitter.fit(f_sf_top, g_plus,"REQ");

	    F1D f_sf_bot = new F1D("f_sf_bot","[a] + [b]/x + [c]/(x*x)", 0.8, beam_energy);	    
	    DataFitter.fit(f_sf_bot, g_minus,"REQ");

	   
	    Double ma = f_sf_mean2.getParameter(0);
	    Double mb = f_sf_mean2.getParameter(1);
	    Double mc = 0.0;//f_sf_mean2.getParameter(2);
	    Double md = 0.0;//f_sf_mean2.getParameter(3);
	    
	    Double topa = f_sf_top.getParameter(0);
	    Double topb = f_sf_top.getParameter(1);
	    Double topc = f_sf_top.getParameter(2);
	    
	    Double bota = f_sf_bot.getParameter(0);
	    Double botb = f_sf_bot.getParameter(1);
	    Double botc = f_sf_bot.getParameter(2);
	    
	    temp_maxpara.clear();
	    temp_minpara.clear();

 	    temp_maxpara.add(topa);
 	    temp_maxpara.add(topb);
 	    temp_maxpara.add(topc);

 	    temp_minpara.add(bota);
 	    temp_minpara.add(botb);
 	    temp_minpara.add(botc);

	    m_maxfit_para.put(s,temp_maxpara);
	    m_minfit_para.put(s,temp_minpara);

	    f_sf_mean2.setLineWidth(5);
	    f_sf_mean2.setLineStyle(0);
	    f_sf_mean2.setLineColor(2);
	    
	    f_sf_top.setLineWidth(5);
	    f_sf_top.setLineStyle(2);
	    f_sf_top.setLineColor(2);

	    f_sf_bot.setLineWidth(5);
	    f_sf_bot.setLineStyle(2);
	    f_sf_bot.setLineColor(2);
	    
	    g_mean2.setTitle("Mean SF with #sigma vs p");
	    g_mean2.setTitleX("p [GeV]");
	    g_mean2.setTitleY("Mean SF");
	    g_mean2.setMarkerSize(2);
	    g_mean2.setMarkerStyle(0);
	    g_mean2.setMarkerColor(2);
	    	    
	    h_temp.setTitle("Sampling Fraction for Sector " + Integer.toString(s));
	    h_temp.setTitleX("p [GeV]");
	    h_temp.setTitleY("Etot/p");

	    c_g_meanfits.draw(h_temp,"colz");       
	    c_g_meanfits.draw(f_sf_top,"same");
	    c_g_meanfits.draw(f_sf_bot,"same");
	    c_g_meanfits.draw(f_sf_mean2,"same");
	    	
	    c_g_meanfits.save("g_meansf_sector_"+s_run+".png");

	}

	System.out.println(">> NOW FIXING EC SF CUT PROPERTIES IN JSON FILE FOR RUN " + s_run );
	RunPropertiesLoader run_properties = new RunPropertiesLoader();
	run_properties.loadRunProperties(run);
	
	System.out.println(" >> MAX FIT FOR SECTOR 1 "  + m_maxfit_para.get(0) );

	run_properties.getRunInfoClass().getRunParametersClass("r"+s_run_number).getCutTypeClass("cut_nom").getECCutParametersClass("ec_sf_cut").setMaxFitParametersSector1( m_maxfit_para.get(0) );
	run_properties.getRunInfoClass().getRunParametersClass("r"+s_run_number).getCutTypeClass("cut_nom").getECCutParametersClass("ec_sf_cut").setMaxFitParametersSector2( m_maxfit_para.get(1) );
	run_properties.getRunInfoClass().getRunParametersClass("r"+s_run_number).getCutTypeClass("cut_nom").getECCutParametersClass("ec_sf_cut").setMaxFitParametersSector3( m_maxfit_para.get(2) );
	run_properties.getRunInfoClass().getRunParametersClass("r"+s_run_number).getCutTypeClass("cut_nom").getECCutParametersClass("ec_sf_cut").setMaxFitParametersSector4( m_maxfit_para.get(3) );
	run_properties.getRunInfoClass().getRunParametersClass("r"+s_run_number).getCutTypeClass("cut_nom").getECCutParametersClass("ec_sf_cut").setMaxFitParametersSector5( m_maxfit_para.get(4) );
	run_properties.getRunInfoClass().getRunParametersClass("r"+s_run_number).getCutTypeClass("cut_nom").getECCutParametersClass("ec_sf_cut").setMaxFitParametersSector6( m_maxfit_para.get(5) );
	
	run_properties.getRunInfoClass().getRunParametersClass("r"+s_run_number).getCutTypeClass("cut_nom").getECCutParametersClass("ec_sf_cut").setMinFitParametersSector1( m_minfit_para.get(0) );
	run_properties.getRunInfoClass().getRunParametersClass("r"+s_run_number).getCutTypeClass("cut_nom").getECCutParametersClass("ec_sf_cut").setMinFitParametersSector2( m_minfit_para.get(1) );
	run_properties.getRunInfoClass().getRunParametersClass("r"+s_run_number).getCutTypeClass("cut_nom").getECCutParametersClass("ec_sf_cut").setMinFitParametersSector3( m_minfit_para.get(2) );
	run_properties.getRunInfoClass().getRunParametersClass("r"+s_run_number).getCutTypeClass("cut_nom").getECCutParametersClass("ec_sf_cut").setMinFitParametersSector4( m_minfit_para.get(3) );
	run_properties.getRunInfoClass().getRunParametersClass("r"+s_run_number).getCutTypeClass("cut_nom").getECCutParametersClass("ec_sf_cut").setMinFitParametersSector5( m_minfit_para.get(4) );
	run_properties.getRunInfoClass().getRunParametersClass("r"+s_run_number).getCutTypeClass("cut_nom").getECCutParametersClass("ec_sf_cut").setMinFitParametersSector6( m_minfit_para.get(5) );
	run_properties.writeRunProperties();
	
	System.out.println(">> FINISHED " );

    }


    public static void SliceNFit(Vector< Vector<H2F> > h_temp_in, Vector< GraphErrors > g_sf_sect_means, Vector< GraphErrors > g_sf_sect_sigmas, HashMap<Integer, ArrayList<H1F> > m_el_sect_slices_ectotp ){

	int counter = 0;
	int cutlvl = 5;
	HashMap<Integer, Vector<Double> > h2_temp_binX = new HashMap<Integer, Vector<Double> >();
		
	for( Vector<H2F> h2_temp : h_temp_in ){
	    //CUTLLVL GETS THE CORRECT SF GRAPH TO SLICE N FIT
	    System.out.println(">> PREPARING TO SLICE 'N' FIT " + h2_temp.get(cutlvl).getTitle() );
	    H2F h2_temp_rebinX = h2_temp.get(cutlvl).rebinX(5); //was 5
	    H2F h2_temp_rebinXY = h2_temp_rebinX.rebinY(5);

	    ArrayList <H1F> h_temp_rebinXY_sliceX = new ArrayList<H1F>();
	    h_temp_rebinXY_sliceX = h2_temp_rebinXY.getSlicesX();
	    
 	    Vector<Double> bin_center_temp = new Vector<Double>();
	    for( int bin = 1; bin < h2_temp_rebinXY.getXAxis().getNBins(); bin++ ){ /// was bin = 1 to start 
		double bin_center = h2_temp_rebinXY.getXAxis().getBinCenter(bin);
		System.out.println(" >> SECTOR " + counter + " BIN CENTER " + bin_center ); 
		bin_center_temp.add(bin_center);
		h2_temp_binX.put(counter,bin_center_temp);
		
	    }
		
	    m_el_sect_slices_ectotp.put(counter, h_temp_rebinXY_sliceX);
	    
	    //ParallelSliceFitter fit_temp = new ParallelSliceFitter(h2_temp_rebinXY);
	    //fit_temp.setMinBin(1);
	    //fit_temp.setMaxBin(18);
 	    //fit_temp.fitSlicesX();
	    //GraphErrors temp_mean = fit_temp.getMeanSlices();
	    //GraphErrors temp_sigma = fit_temp.getSigmaSlices();	   
	    //h_bpid.g_sf_sect_meansfits.add(temp_mean);
	    //h_bpid.g_sf_sect_sigmasfits.add(temp_sigma);
	    
	    counter=counter+1; //SECTOR NUMBER
	}
	
	System.out.println(" >> FITTING HISTOGRAMS ACROSS ALL SECTORS ");
	for( int sector = 0; sector < 6; sector++ ){
	    System.out.println(" FITTING SECTOR " + sector );
	    ArrayList<H1F> al_h_ectotp = m_el_sect_slices_ectotp.get(sector);
	    System.out.println( " >> NUMBER OF SLICES " + al_h_ectotp.size());

	    GraphErrors g_temp = g_sf_sect_means.get(sector);
	    GraphErrors g_temp_sigmas = g_sf_sect_sigmas.get(sector);
	    g_temp.setTitle("Mean Sector " + Integer.toString(sector));
	    g_temp_sigmas.setTitle(" SIGMAS SECTOR " + Integer.toString(sector));
	    
	    for( int n_htemp = 1; n_htemp < al_h_ectotp.size() -1; n_htemp++){
		System.out.println(" >> n_htemp " + n_htemp );
		H1F h_temp = al_h_ectotp.get(n_htemp);
		F1D fit_temp = Calculator.fitHistogram(h_temp);
		double fit_mean = fit_temp.getParameter(1);
		double fit_sigma = fit_temp.getParameter(2);
		double fit_mean_err = fit_temp.parameter(1).error();
		double fit_sigma_err = fit_temp.parameter(2).error();
		//System.out.println(" >> size of bin vector " + h2_temp_binX.get(sector).size() );
		//System.out.println(" >> ADDING POINT " + h2_temp_binX.get(sector).get(n_htemp) +  " " + fit_mean + " ERROR " + fit_mean_err + " SIG ERROR " + fit_sigma_err);
		if( fit_mean > 1.0 || fit_mean < 0.0 || Math.abs(fit_sigma) > 1.0 ){ continue; }
		g_temp.addPoint(h2_temp_binX.get(sector).get(n_htemp), fit_mean, fit_mean_err, fit_sigma_err);
		g_temp_sigmas.addPoint(h2_temp_binX.get(sector).get(n_htemp), fit_sigma, fit_mean_err, fit_sigma_err);
	    }

	}

    }

}

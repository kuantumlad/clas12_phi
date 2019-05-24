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



public class runCalFiducial{
    

    ////////////////////////////////////////
    //MAP FIDUCIAL CUTS ONTO CALORIMETERS
    /* PCAL
     d_left = y_rot > 1.86*x_rot + 51.0;
    d_right = y_rot > -1.876*x_rot + 49.0;
    d_up = y_rot < 330.0;
    d_down = y_rot > 52.0;
    */
     
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


    }

}

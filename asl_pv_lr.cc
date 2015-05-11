/*   asl_pv_lr.cc Partial volume correction on single TI ASL data using linear regression

    Moss Zhao - IBME Quantitative Biomedical Inference (QuBIc) Group

    Copyright (C) 2015 University of Oxford  */

#include <iostream>
#include <math.h>
#include <string>
#include "newmatap.h"
#include "newmatio.h"
#include "newimage/newimageall.h"
#include "miscmaths/miscmaths.h"
#include "utils/tracer_plus.h"
#include "stdlib.h"

#include "readoptions.h"
#include "asl_pv_lr_functions.h"

using namespace Utilities;
using namespace NEWMAT;
using namespace NEWIMAGE;
using namespace MISCMATHS;

using namespace OXASL;

int main(int argc, char *argv[])
{
  try {

    cout << "ASL_PV_LR (1.0)" << endl;

    //parse command line (puts these into the log file)
    ReadOptions& opts = ReadOptions::getInstance();
    opts.parse_command_line(argc,argv);

    cout << "Loading data" << endl;

    // data file name
    string data_filename;
    data_filename = opts.datafile.value();

    // load data (3D matrix)
    volume<float> data;
    read_volume(data, opts.datafile.value());
    //volume4D<float> data;
    //read_volume4D(data,opts.datafile.value());

    // load mask
    volume<float> mask;
    read_volume(mask, opts.maskfile.value());

    // load pvgm
    volume<float> pvgm;
    read_volume(pvgm, opts.pvgmfile.value());

    // load pvwm
    volume<float> pvwm;
    read_volume(pvwm, opts.pvwmfile.value());

    // load output file name
    string outname;
    outname = opts.outname.value();

    // load kernel size
    float kernel;
    kernel = opts.kernel.value();

    // GM and WM corrected single TI file
    string gm_out_filename;
    gm_out_filename = gm_out_filename + data_filename + "_gm";
    //strcat(gm_out_filename, data_filename);
    //strcat(gm_out_filename, "_gm");

    string wm_out_filename;
    wm_out_filename = wm_out_filename + data_filename + "_wm";
    //strcat(wm_out_filename, data_filename);
    //strcat(wm_out_filename, "_wm");

    // Empty 3D matrix to save  PV corrected results
    volume<float> data_gm;
    volume<float> data_wm;

    pv_correct(data, mask, pvgm, kernel, data_gm);
    pv_correct(data, mask, pvwm, kernel, data_wm);

    save_volume(data_gm, gm_out_filename);
    save_volume(data_wm, wm_out_filename);

    /*


    // load mask
    volume<float> mask(data.xsize(),data.ysize(),data.zsize());
    mask.setdims(data.xdim(),data.ydim(),data.zdim());
    read_volume(mask,opts.maskfile.value());

    Matrix asldata;
    asldata = data.matrix(mask); // creates a 2D matrix of all voxels and their time series
    //data.setmatrix(asldata,mask);
    int nvox=asldata.Ncols();
    
    //load aif
    volume4D<float> aif;
    read_volume4D(aif,opts.aif.value());

    
    // select AIF based on metric
    // load metric image if it exists
    volume<float> metric;
    if (opts.metric.set()) {
      cout << "Preparing AIFs" << endl;
      read_volume(metric,opts.metric.value());

      Prepare_AIF(aif, metric, mask, opts.mthresh.value());
      
      //volume4D<float> aifout;
      //aifout.setmatrix(aif,mask);
      save_volume4D(aif,opts.outname.value()+"_aifs");
    }

    Matrix aifmtx;
    aifmtx = aif.matrix(mask);

    // do deconvolution
    cout << "Performing deconvolution" << endl;
    ColumnVector mag; // Column Vector to save CBF of each voxel (taken as the largest value of each column of scaled residue matrix)
    Matrix resid;
    Deconv(asldata,aifmtx,opts.dt.value(),mag,resid);

    // estimate BAT (of tissue)
    ColumnVector batt;
    if (opts.batout.set() | (opts.tcorrect.set() & !opts.batt.set())) {
      cout << "Estimating BAT" << endl;
      Estimate_onset(asldata,batt,opts.dt.value());

      if (opts.batout.set()) {
	//output the BAT image (from the tissue)
        volume4D<float> batout;
        batout.setmatrix(batt.AsMatrix(1,nvox),mask);
        save_volume4D(batout,opts.outname.value()+"_bat");
      }
    }

    // correct aif magntiude for differences in arrival time between aif and tissue 
    ColumnVector batd;
    if(opts.tcorrect.set()) {
      cout << "Performing timing correction" << endl;
      

      if (opts.bata.set()) {
	//calculate difference between tissue and AIF curves using suppled BAT images
        volume4D<float> bat_art;
        read_volume4D(bat_art,opts.bata.value());
        if (opts.batt.set()) {
	  //load supplied tissue BAT
          volume4D<float> bat_tiss;
          read_volume4D(bat_tiss,opts.batt.value());
          batt = (bat_tiss.matrix(mask)).AsColumn();
        }

        if (opts.metric.set()) {
	  // correct the AIF bat image to match the AIFs where a metric image has been supplied
          Prepare_AIF(bat_art,metric,mask,opts.mthresh.value());
        }
        ColumnVector bata;
        bata = (bat_art.matrix(mask)).AsColumn();
        batd = batt-bata;
      }
      else {
	//otherwise estimate BAT difference using the peak in the residue function
	//Estimate_BAT_difference(resid,batd,opts.dt.value());

	// Estiamte BAT difference using edge detection
        ColumnVector bata;
        Estimate_onset(aifmtx,bata,opts.dt.value());
        batd = batt-bata;
      }

      for (int i=1; i<=batd.Nrows(); i++) {
        if (batd(i)<0.0)
          batd(i)=0.0;
      }
      Correct_magnitude(mag,batd,opts.T1.value(),opts.dt.value(),opts.fa.value());
    }

    if(opts.std.set()) {
      // do wild boostrapping std dev estimation for cbf
      cout << "Performing wild bootstrapping for precision estimation" << endl;
      ColumnVector magstd;
      BootStrap(aifmtx, asldata, opts.dt.value(), mag, resid, opts.nwb.value(), magstd);
      if (opts.tcorrect.set()) {
	// if needed we should correct the std dev for timing discrpancies between aif and ctc
        Correct_magnitude(magstd,batd,opts.T1.value(),opts.dt.value(),opts.fa.value());
      }

      // save it
      volume4D<float> stdoutVol; //stdout is a reserved name - can cause weird compile errors
      stdoutVol.setmatrix(magstd.AsMatrix(1,nvox),mask);
      save_volume4D(stdoutVol,opts.outname.value()+"_stddev");
    }
    
    cout << "Saving results" << endl;
    //output 
    volume4D<float> residout;
    residout.setmatrix(resid,mask);
    save_volume4D(residout,opts.outname.value()+"_residuals");

    volume4D<float> magout;
    magout.setmatrix(mag.AsMatrix(1,nvox),mask);
    save_volume4D(magout,opts.outname.value()+"_magntiude");

    */

    cout << "ASL_PV_LR - Done!" << endl;

  }

catch(Exception& e)
    {
      cerr << endl << e.what() << endl;
      return 1;
    }
  catch(X_OptionError& e)
    {
      cerr << endl << e.what() << endl;
      return 1;
    }

  cout << "Done." << endl;

  return 0;


}

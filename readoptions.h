/*  readoptions.h

    Moss Zhao - IBME Quantitative Biomedical Inference (QuBIc) Group

    Copyright (C) 2015 University of Oxford  */

#if !defined(ReadOptions_h)
#define ReadOptions_h

#include <string>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include "utils/options.h"
#include "utils/log.h"

using namespace Utilities;

namespace OXASL {

class ReadOptions {
public:
  static ReadOptions& getInstance();
  ~ReadOptions() { delete ropt; }
  
  Option<bool> help;

  Option<string> datafile;
  Option<string> maskfile;
  Option<string> pvgmfile;
  Option<string> pvwmfile;
  Option<string> outname;
  //Option<string> aif;
  Option<int> kernel;

  //Option<string> metric;
  //Option<float> mthresh;

  // timing correction
  //Option<bool> tcorrect;
  //Option<string> bata;
  //Option<string> batt;
  //Option<bool> batout;
  //Option<float> T1;
  //Option<float> fa;

  // std deviation estimation
  //Option<bool> std;
  //Option<int> nwb;

  void parse_command_line(int argc, char** argv);

private:
  ReadOptions();  
  const ReadOptions& operator=(ReadOptions&);
  ReadOptions(ReadOptions&);

  OptionParser options; 
      
  static ReadOptions* ropt;
  
};

inline ReadOptions& ReadOptions::getInstance() {
  if(ropt == NULL)
    ropt = new ReadOptions();

  return *ropt;
}

inline ReadOptions::ReadOptions():

/* Option class:
  @param k Comma seperated list of key aliases
  @param v Default value for this option
  @param ht Help text to be printed when outputting usage
  @param c If true then this option is compulsory
  @param f This options argument requirements
*/

help(string("-h,--help"), false,
  string("display this message"),
  false, no_argument),

//input files
/*
datafile(string("--data,--datafile"), string("data"),
  string("ASL data file"),
  true, requires_argument),  
*/

// Compulsory arguments
datafile(string("--data,--datafile"), string("data"),
  string("Single_TI data file"),
  true, requires_argument),

maskfile(string("--mask"), string("mask"),
  string("mask"), true, requires_argument),

pvgmfile(string("--pvgm"), string("pvgm"),
  string("GM PV map"), true, requires_argument),

pvwmfile(string("--pvwm"), string("pvwm"),
  string("WM PV map"), true, requires_argument),

outname(string("--out"), string("outname"),
  string("Output directory name"), true, requires_argument),

kernel(string("--kenrel"), 5,
  string("Kernel size (Default: 5)"), true, requires_argument),

/*
outname(string("--out"), string("asl_mfree"),
  string("Output directory name"), true, requires_argument),
*/

/*
aif(string("--aif"),string(""),
  string("Arterial input functions for deconvolution (4D volume, one aif for each voxel within mask)"),
  true,requires_argument),

dt(string("--dt"),1.0,
  string("Temporal spacing in data (s)\n"),
  true,requires_argument),

metric(string("--metric"),string(""),
  string("Metric image"), false,requires_argument),

mthresh(string("--mthresh"),0.0,
  string("Metric threshold\n"), false,requires_argument),

tcorrect(string("--tcorrect"),false,
  string("Apply correction for timing difference between AIF and tissue curve"),
  false,no_argument),

bata(string("--bata"),string(""), string("arterial BAT image"), false,requires_argument),

batt(string("--batt"),string(""), string("tissue BAT image"), false,requires_argument),

batout(string("--bat"),false, 
  string("Estimate tissue BAT from data (and save this image)"), false,no_argument),

T1(string("--t1"),1.6, string("T1 (of blood) value"), false,requires_argument),

fa(string("--fa"),0.0, string("Flip anlge (is using LL readout)"), false,requires_argument),

std(string("--std"),false,
  string("Calculate standard deviations on perfusion values using wild bootstrapping"),
  false,no_argument),

nwb(string("--nwb"),1000,
  string("Number of permuatations for wild bootstrapping"),
  false,requires_argument),
*/
//options("asl_mfree","asl_mfree --verbose\n") {
options("asl_pv_lr","asl_pv_lr --verbose\n") {
  try {
    options.add(help);

    options.add(datafile);
    options.add(maskfile);
    options.add(pvgmfile);
    options.add(pvwmfile);
    options.add(outname);
    options.add(kernel);
    //options.add(aif);
    //options.add(dt);

       //options.add(metric);
       //options.add(mthresh);

       //options.add(tcorrect);
       //options.add(bata);
       //options.add(batt);
       //options.add(batout);
       //options.add(T1);
       //options.add(fa);
       
       //options.add(std);
       //options.add(nwb);

        }
     catch(X_OptionError& e) {
       options.usage();
       cerr << endl << e.what() << endl;
     } 
     catch(std::exception &e) {
       cerr << e.what() << endl;
     }    
     
   }

}
#endif




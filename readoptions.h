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
  Option<string> pvfile;
  Option<string> outname;
  Option<int> kernel;

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

// Compulsory arguments
datafile(string("--data,--datafile"), string("data"),
  string("Single_TI data file"),
  true, requires_argument),

maskfile(string("--mask"), string("mask"),
  string("mask"), true, requires_argument),

pvfile(string("--pv"), string("pvmap"),
  string("GM PV map"), true, requires_argument),

outname(string("--out"), string("outname"),
  string("Output directory name"), true, requires_argument),

kernel(string("--kernel"), 5,
  string("Kernel size (Default: 5)"), true, requires_argument),

options("asl_pv_lr","asl_pv_lr --verbose\n") {
  try {
    options.add(help);

    options.add(datafile);
    options.add(maskfile);
    options.add(pvfile);
    options.add(outname);
    options.add(kernel);

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




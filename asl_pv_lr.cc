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

    // load mask
    volume<float> mask;
    read_volume(mask, opts.maskfile.value());

    // load pvgm
    volume<float> pv_map;
    read_volume(pv_map, opts.pvfile.value());

    // load output file name
    string outname;
    outname = opts.outname.value();

    // load kernel size
    int kernel;
    kernel = opts.kernel.value();

    // Empty 3D matrix to save PV corrected results
    volume<float> data_pv_corr;

    pv_correct(data, mask, pv_map, kernel, data_pv_corr);

    save_volume(data_pv_corr, outname);

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

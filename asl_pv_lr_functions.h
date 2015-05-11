/*   asl_pv_lr_functions.h Partial volume correction on single TI ASL data using linear regression

    Moss Zhao - IBME Quantitative Biomedical Inference (QuBIc) Group

    Copyright (C) 2015 University of Oxford  */

#if !defined(asl_pv_lr_functions)
#define asl_pv_lr_functions

#include <string>
#include "newimage/newimageall.h"
#include "miscmaths/miscmaths.h"

using namespace MISCMATHS;
using namespace NEWIMAGE;

namespace OXASL {


  void pv_correct(const volume<float>& data_in, const volume<float>& mask, const volume<float>& pv_map, float kernel, volume<float>& data_out);

  // PV correction using linear regression (Asllani's method)
  ReturnMatrix correct_pv_lr(const volume<float>& data_in, const volume<float>& mask, const volume<float>& pv_map, float kernel);

  /*

  //output the result of deconvolution with an AIF
  void Deconv(const Matrix& data, const Matrix& aif, float dt,ColumnVector& mag, Matrix& resid);

  // prepare the AIFs
  void Prepare_AIF(volume4D<float>& aif, const volume<float>& metric, const volume<float>& mask, const float mthresh);
  // do SVD convoloution
  ReturnMatrix SVDdeconv(const Matrix& data, const Matrix& aif, float dt);
  // create a (simple) convolution matrix
  ReturnMatrix convmtx(const ColumnVector& invec);
  // create a cicular deconvolution
  ReturnMatrix SVDdeconv_circular(const Matrix& data, const Matrix& aif, float dt);
  ReturnMatrix convmtx_circular(const ColumnVector& invec);
  ReturnMatrix SVDdeconv_wu(const Matrix& data, const Matrix& aif, float dt);

  void Estimate_BAT_difference(const Matrix& resid, ColumnVector& batd, const float dt);
  void Correct_magnitude(ColumnVector& mag, const ColumnVector& batd, const float T1, const float dt, const float fa);

  void Estimate_onset(const Matrix& curves, ColumnVector& bate, const float dt);

  // estimation of magntiude precision
  void BootStrap(const Matrix& aif, const Matrix& data, float dt, const ColumnVector& mag, const Matrix& resid, int Nwb, ColumnVector& magsd);

  */
}

#endif

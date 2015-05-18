/*   asl_pv_lr_functions.cc Partial volume correction on single TI ASL data using linear regression

    Moss Zhao - IBME Quantitative Biomedical Inference (QuBIc) Group

    Copyright (C) 2015 University of Oxford  */

#include "asl_pv_lr_functions.h"

namespace OXASL
{

  void pv_correct(const volume<float>& data_in, const volume<float>& mask, const volume<float>& pv_map, int kernel, volume<float>& data_out)
  {
    // Define variables to store NaN corrected data
    volume<float> data_in_corr(data_in.xsize(), data_in.ysize(), data_in.zsize());
    volume<float> mask_in_corr(mask.xsize(), mask.ysize(), mask.zsize());
    volume<float> pv_map_in_corr(pv_map.xsize(), pv_map.ysize(), pv_map.zsize());

    data_in_corr   = correct_NaN(data_in);
    mask_in_corr   = correct_NaN(mask);
    pv_map_in_corr = correct_NaN(pv_map);

    // Linear regression partial volume correction
    data_out = correct_pv_lr(data_in_corr, mask_in_corr, pv_map_in_corr, kernel);

  }

  // Function to correct PV using LR method
  volume<float> correct_pv_lr(const volume<float>& data_in, const volume<float>& mask, const volume<float>& pv_map, int kernel)
  {

    volume<float> submask;
    volume<float> data_roi;
    volume<float> pv_roi;
    RowVector pseudo_inv;
    Matrix pv_corr_result;

    // Variables to store the boundary index of submask (ROI)
    int x_0;
    int x_1;
    int y_0;
    int y_1;
    int z_0;
    int z_1;

    int count;
    float pv_ave = 0.0f;

    // Get x y z dimension
    int x = data_in.xsize();
    int y = data_in.ysize();
    int z = data_in.zsize();

    volume<float> corr_data(x, y, z); // result matrix

    // Linear regression to correct (smooth) the data
    for (int i = 0; i < x; i++) {
      for (int j = 0; j < y; j++) {
        for (int k = 0; k < z; k++) {
          // Only work with positive voxels
          if(mask.value(i, j, k) > 0) {

            // Determine ROI boundary index
            x_0 = max(i - kernel, 0);
            x_1 = min(i + kernel, x - 1);
            y_0 = max(j - kernel, 0);
            y_1 = min(j + kernel, y - 1);
            z_0 = max(k - kernel, 0);
            z_1 = min(k + kernel, z - 1);

            // create a submask here
            mask.setROIlimits(x_0, x_1, y_0, y_1, z_0, z_1);
            mask.activateROI();
            submask = mask.ROI();

            // calculate the sum of all elements in submask
            // proceed if sum is greater than 5 (arbitrary threshold)
            if(submask.sum() > 5) {
              /* Create an ROI (sub volume of data and PV map),
                then mask it with submask to create sub data and PV map */

              // Obtain ROI volume (must set limits and activate first)
              data_in.setROIlimits(x_0, x_1, y_0, y_1, z_0, z_1);
              pv_map.setROIlimits(x_0, x_1, y_0, y_1, z_0, z_1);
              data_in.activateROI();
              pv_map.activateROI();
              data_roi = data_in.ROI();
              pv_roi = pv_map.ROI();
              
              
              ColumnVector data_roi_v_t = ColumnVector(data_roi.xsize() * data_roi.ysize() * data_roi.zsize());
              ColumnVector pv_roi_v_t = ColumnVector(pv_roi.xsize() * pv_roi.ysize() * pv_roi.zsize());

              count = 0;
              // Apply a mask on data_roi and pv_roi
              // Extract values from data_roi and pv_roi whose mask values are non-zero
              for(int a = 0; a < data_roi.xsize(); a++) {
                for(int b = 0; b < data_roi.ysize(); b++) {
                  for(int c = 0; c < data_roi.zsize(); c++) {
                    if(submask.value(a, b, c) > 0) {
                      data_roi_v_t.element(count) = data_roi.value(a, b, c);
                      pv_roi_v_t.element(count) = pv_roi.value(a, b, c);
                      count++;

                    }
                    else {
                      continue;
                    }
                  }
                }
              }
              
              ColumnVector data_roi_v = ColumnVector(count);
              ColumnVector pv_roi_v = ColumnVector(count);
              
              
              for(int a = 0; a < count; a++) {
                data_roi_v.element(a) = data_roi_v_t.element(a);
                pv_roi_v.element(a) = pv_roi_v_t.element(a);
              }
              
              // Unstable function
              //data_roi_v = apply_mask(data_roi, submask);
              //pv_roi_v = apply_mask(pv_roi, submask);

              
              // Now data_roi_v and pv_roi_v stores the non-zero elemnts within this submask

              // Deactivate ROI
              data_in.deactivateROI();
              pv_map.deactivateROI();

              // If pv_roi is all zeros, then the pseudo inversion matrix will be singular
              // This will cause run time error
              // So we assign the corrected result to zero in such cases
              if(pv_roi_v.IsZero()) {
                corr_data.value(i, j, k) = 0.0f;
              }
              else {
                // Compute pseudo inversion matrix of PV map
                // ((P^t * P) ^ -1) * (P^t)
                pseudo_inv = ( (pv_roi_v.t() * pv_roi_v).i() ) * (pv_roi_v.t());

                //cout << "  hh" << endl;

                // Get average PV value of the current kernel
                pv_ave = (float) pv_roi_v.Sum() / pv_roi_v.Nrows();

                // Calculate PV corrected data only if there is some PV compoment
                // If there is little PV small, make it zero
                if(pv_ave >= 0.01) {
                  pv_corr_result = pseudo_inv * data_roi_v;
                  corr_data.value(i, j, k) = pv_corr_result.element(0, 0);
                }
                else {
                  corr_data.value(i, j, k) = 0.0f;
                }
              }
            }

            else {
              // do nothing at the moment
            } // end submask

            // Discard current submask (ROI)
            mask.deactivateROI();

          }

          else{
            // do nothing at the moment
          } // end mask

        }
      }
    }

    return corr_data;

  } // End function correct_pv_lr

  // Function to apply mask on 3D matrix and output to a column vector
  volume<float> correct_NaN(const volume<float>& data_in) {

    // Clone the input data to output data
    volume<float> data_out = data_in;

    for(int i = 0; i < data_in.xsize(); i++) {
      for(int j = 0; j < data_in.ysize(); j++) {
        for(int k = 0; k < data_in.zsize(); k++) {
          // IEEE standard: comparison between NaN values is always false
          // i.e. NaN == NaN is false
          // In this case, we set it to zero
          if(data_in.value(i, j, k) != data_in.value(i, j, k)) {
            data_out.value(i, j, k) = 0.0f;
          }
          else {
            continue;
          }
        }
      }
    }

    return data_out;

  }

  // This function extracts all non-zero elemtns of input matrix to a column vector
  // Current unstable
  /*
  ColumnVector apply_mask(const volume<float>& data_in, const volume<float>& mask) {
    
    ColumnVector temp = ColumnVector(mask.xsize() * mask.ysize() * mask.zsize());
    int count = 0;

    for(int i = 0; i < mask.xsize(); i++) {
      for(int j = 0; j < mask.ysize(); j++) {
        for(int k = 0; k < mask.zsize(); k++) {
          if(mask.value(i, j, k) > 0) {
           temp.element(count) = data_in.value(i, j, k);
           count++;
          }
          else {
            continue;
          }
        }
      }
    }

    ColumnVector vector_out = ColumnVector(count);

    for(int a = 0; a < count; a++) {
      vector_out.element(a) = temp.element(a);
    }

    return vector_out;
  }
  */
  
  /*
  ReturnMatrix SVDdeconv(const Matrix& data, const Matrix& aif, float dt) {
    // do a singular value deconvolution of the data to get residue function

    int nti = data.Nrows();
    int nvox = data.Ncols();
    float truncfac = 0.2;

    // voxelwise SVD deconvolution
    Matrix aifconv; Matrix residue(nti,nvox);
    DiagonalMatrix S; DiagonalMatrix D; Matrix U; Matrix V;
    for (int v=1; v<=nvox; v++) {
      //make convolution matrix
      aifconv = dt * convmtx(aif.Column(v));
      //SVD
      SVD(aifconv,S,U,V);
      // invert the singular values
      D = S.i();
      // truncate (zero all singular values below threshold)
      for (int i=2; i<=D.Nrows(); i++) {
	if (S(i,i) < truncfac*S(1,1)) D(i,i)=0;
      }
      // calculate resdiue
      residue.Column(v) = V*D*U.t()*data.Column(v);
    }
      
    return residue; 
  }

  ReturnMatrix convmtx(const ColumnVector& invec){
    // create a (simple) convolution matrix

    int nentry = invec.Nrows();
    Matrix cmat(nentry,nentry);
    cmat=0.0;
    for (int i=1; i<=nentry; i++) {
      cmat.SubMatrix(i,i,1,i) = ((invec.Rows(1,i)).Reverse()).AsRow();
    }

    return cmat;
  }

  ReturnMatrix SVDdeconv_circular(const Matrix& data, const Matrix& aif, float dt) {
    // do a singular value deconvolution of the data to get residue function

    int nti = data.Nrows();
    int nvox = data.Ncols();
    float truncfac = 0.2;

    int nextra = floor(nti*1.2);
    ColumnVector padding(nextra);
    padding = 0.0;

    // voxelwise SVD deconvolution
    Matrix aifconv; Matrix residue(nti+nextra,nvox);
    DiagonalMatrix S; DiagonalMatrix D; Matrix U; Matrix V;
    for (int v=1; v<=nvox; v++) {
      //make convolution matrix
      aifconv = dt * convmtx_circular(aif.Column(v) & padding);
      //SVD
      SVD(aifconv,S,U,V);
      // invert the singular values
      D = S.i();
      // truncate (zero all singular values below threshold)
      for (int i=2; i<=D.Nrows(); i++) {
	if (S(i,i) < truncfac*S(1,1)) D(i,i)=0;
      }
      // calculate resdiue
      residue.Column(v) = V*D*U.t()*(data.Column(v) & padding);
    }
      
    residue = residue.Rows(1,nti); // only keep the number of timepoints that we have in the input.

    return residue; 
  }

 ReturnMatrix SVDdeconv_wu(const Matrix& data, const Matrix& aif, float dt) {
    // do a singular value deconvolution of the data to get residue function
   // circular convolution matrix
   // OI index
   // after Wu MRM 2003.

    int nti = data.Nrows();
    int nvox = data.Ncols();
    float oi_thresh = 0.1;

    int nextra = floor(nti*1.2);
    ColumnVector padding(nextra);
    padding = 0.0;

    // voxelwise SVD deconvolution
    Matrix aifconv; Matrix residue(nti+nextra,nvox);
    ColumnVector resid(nti+nextra);// Length (L) of vector must satisfy L >= 2N to avoid time aliasing
    DiagonalMatrix S; DiagonalMatrix D; Matrix U; Matrix V;
    for (int v=1; v<=nvox; v++) {
      //make convolution matrix
      aifconv = dt * convmtx_circular(aif.Column(v) & padding); // & operator is to stack matrix (columns) vertically
      //SVD
      SVD(aifconv,S,U,V);
      // invert the singular values
      D = S.i();

      //first try with all singular values
      resid = V*D*U.t()*(data.Column(v) & padding);
      float oi;
      oi = 1/(nti*resid.Maximum())* (resid.Rows(3,nti+nextra)-2*resid.Rows(2,nti+nextra-1)+resid.Rows(1,nti+nextra-2)).SumAbsoluteValue();

      // start removing singular values one by one until we get the OI we want
      int i=nti+nextra;
      while (oi>oi_thresh & i>1) {
        D(i,i) = 0;
        resid = V*D*U.t()*(data.Column(v) & padding);
        oi = 1/((nti+nextra)*resid.Maximum())* (resid.Rows(3,nti+nextra)-2*resid.Rows(2,nti+nextra-1)+resid.Rows(1,nti+nextra-2)).SumAbsoluteValue();
        i--;
      }

      residue.Column(v) = resid;
    }


    return residue; 
  }

  // Create matrix a block-circulant matrix (matrix C in Wu's paper)
  ReturnMatrix convmtx_circular(const ColumnVector& invec){
    // create a (simple) convolution matrix

    int nentry = invec.Nrows();
    Matrix cmat(nentry,nentry);
    cmat=0.0;
    for (int i=1; i<=nentry; i++) {
      cmat.SubMatrix(i,i,1,i) = ((invec.Rows(1,i)).Reverse()).AsRow();
      cmat.SubMatrix(i,i,i+1,nentry) = ((invec.Rows(i+1,nentry)).Reverse()).AsRow();
    }

    return cmat;
  }

  void Deconv(const Matrix& data, const Matrix& aif, float dt, ColumnVector& mag, Matrix& resid) {
    //perform deconvolution and output the magnitude and residue function

    //do the deconvolution
    resid = SVDdeconv_wu(data,aif,dt);
    // extract magntiude and residue separately
    int nvox = data.Ncols();
    mag.ReSize(nvox);
    for (int v=1; v<=nvox; v++) {
      mag(v) = (resid.Column(v)).Maximum();
      resid.Column(v) /= mag(v);
    }

  }

  void BootStrap(const Matrix& aif, const Matrix& data, float dt, const ColumnVector& mag, const Matrix& resid, int Nwb, ColumnVector& magsd) {
    // do wild boot strapping estimation of CBF precision - returns the s.d. estimate
    int nvox = data.Ncols();
    int ntpts = data.Nrows();

    // calculate the residuals
    Matrix modelfit(data);
    Matrix aifconv; ColumnVector mfittemp;
    int nextra = resid.Nrows() - aif.Nrows();
    ColumnVector padding(nextra);
    padding = 0.0;
    for (int v=1; v<=nvox; v++) {
      //make convolution matrix
      aifconv = dt * convmtx_circular(aif.Column(v) & padding);
      // calculate resdiue
      mfittemp = aifconv*resid.Column(v);
      modelfit.Column(v) = mfittemp.Rows(1,ntpts);
    }
    Matrix residuals(data);
    residuals = data - modelfit;

    // wild bootstrapping
    ColumnVector Radevec(ntpts);
    Matrix WBdata(data);
    Matrix estresid(resid);
    Matrix WBresiduals(residuals);
    Matrix magdist(Nwb,nvox);
    cout << "WB step (of " << Nwb << "): ";
    for (int b=1; b<=Nwb; b++) {
      cout << b << " " << flush;
      // sample from the Rademacher distrbution
      Radevec=1.0;
      for (int i=1; i<=ntpts; i++) {
	if (rand() > RAND_MAX/2) { Radevec(i)=-1.0; }
      }

      // apply this to the residuals to get WB residuals
      for (int v=1; v<=nvox; v++) {
	WBresiduals.Column(v) = SP(residuals.Column(v),Radevec);
      }
      // make the WB data
      WBdata = modelfit + WBresiduals;

      // do the deconvolution to get magntiude estimates
      estresid = SVDdeconv_wu(WBdata,aif,dt);
      // extract magntiude and residue separately

      for (int v=1; v<=nvox; v++) {
	magdist(b,v) = (estresid.Column(v)).Maximum();
      }
    }
    

    magsd.ReSize(nvox);
    magsd = (stdev(magdist)).t();
    cout << endl << "WB done" << endl;
  }

  void Prepare_AIF(volume4D<float>& aif, const volume<float>& metric, const volume<float>& mask, const float mthresh) {
    // Prepeare the AIF by retaining only the AIF that meet our metric and filling in with the nearest good AIF if necessary
    int nx = aif.xsize();
    int ny = aif.ysize();
    int nz = aif.zsize();

    float bestdist; float dist;
    for (int x=0; x<=nx; x++) { for (int y=0; y<=ny; y++) { for (int z=0; z<=nz; z++) {
      // if metric is not above threshold then search for nearest voxel where it is satisfied
	  //cout << "Voxel: " << x << " " << y << " " << z << ": " << metric(x,y,z) << endl;
	  if ( mask(x,y,z)>0 ) {
      if (metric(x,y,z) < mthresh) {
        bestdist = 1e12;
        Matrix aifcand; // candidate AIF
        for (int sx=0; sx<=nx; sx++) {
          for (int sy=0;sy<=ny; sy++) {
            for (int sz=0; sz<=nz; sz++) {
	      // cout << sx << " " << sy << " " << sz << endl;
              if ( (mask(x,y,z)>0) && (metric(sx,sy,sz) >= mthresh) ) {
                dist = pow(sx-x,2.0) + pow(sy-y,2.0) + pow(sz-z,2.0); //strictly this is distance squared
		//cout << dist << endl;
                if (dist < bestdist) {
                  aifcand = aif.voxelts(sx,sy,sz); 
	      //cout << "New cand: " << aifcand.t() << endl;
                  bestdist=dist;
                }
                else if (dist == bestdist) {
                  aifcand |= aif.voxelts(sx,sy,sz); // Concatinate matrix horizontally
	      //cout << aifcand.t() << endl;
                }
              }
            }
          }
        }
	//cout << "AIF: " << mean(aifcand,2).t() << endl;
        aif.setvoxelts(mean(aifcand,2),x,y,z);
      }
    }
    } } }
  }

  void Correct_magnitude(ColumnVector& mag, const ColumnVector& batd, const float T1, const float dt=0.0, const float fa=0.0) {
    //magnitude correction to take acocutn fop differences in BAT between aif and tissue

    mag = SP(mag,exp( batd/T1 ));
    if (fa>0) {
      for (int v=1; v<=mag.Nrows(); v++) {
        mag(v) *= 1/pow( cos(fa/180*M_PI),floor((batd(v)-1e-3)/dt) ); //the 1e-3 deals with the case where batd is a integer mul;tiple of dt
      }
    }
    cout << endl;
  }

  void Estimate_BAT_difference(const Matrix& resid, ColumnVector& batd, const float dt) {
    //Estimate the BAT difference between AIF and tissue curve from residue function

    batd.ReSize(resid.Ncols());

    float magtemp;
    int battemp;
    for (int v=1; v<=resid.Ncols(); v++) {
      magtemp = (resid.Column(v)).Maximum1(battemp);
      if (battemp > resid.Nrows()/2) {
	battemp -= resid.Nrows();
      }
      batd(v) = dt*(battemp-1);
    }

  }

  void Estimate_onset(const Matrix& curves, ColumnVector& bate, const float dt) {
    //Estimate the time of onset of the curve using edge detection

    int ntpts=curves.Nrows();
    bate.ReSize(curves.Ncols());

    ColumnVector kernel(ntpts);
    kernel = 0.0;
    ColumnVector kern(7);
    kern << 0.006 << 0.061 << 0.242 << 0.383 << 0.242 << 0.061 << 0.006;
    kernel.Rows(1,kern.Nrows()) = kern;
    
    Matrix kernmtx;
    kernmtx = convmtx(kernel);
    
    ColumnVector smoothdata(curves.Column(1));
    ColumnVector dgrad(ntpts-1);
    for (int v=1; v<=curves.Ncols(); v++) {
      //apply smoothing kernel
      smoothdata = kernmtx*curves.Column(v);
      // take gradient (forward difference)
      dgrad = smoothdata.Rows(2,ntpts) - smoothdata.Rows(1,ntpts-1);
      float gthresh = 0.2*dgrad.Maximum();
      int i=1; bool cont=true;
      while ((i<ntpts) & cont) {
        if (dgrad(i)>gthresh)
          cont=false;
        else
          i++;
	       }
      //cout << i << " " ;
      bate.Row(v) = i*dt;
    }
    cout << endl;
  }

  */

}

/*****************************************************************************/
// File: dft.cpp
// Author: David Taubman
// Last Revised: 28 August, 2007
/*****************************************************************************/
// Copyright 2007, David Taubman, The University of New South Wales (UNSW)
/*****************************************************************************/

#define _USE_MATH_DEFINES  // Makes `M_PI' available as a constant.

#include <stdlib.h>
#include <math.h>
#include "dft.h"

/*****************************************************************************/
/*                            my_direct_dft::init                            */
/*****************************************************************************/

void my_direct_dft::init(int N, bool is_forward)
{
  cleanup(); // Delete any pre-existing buffers.
  this->N = N;
  real_buf = new float[N];
  imag_buf = new float[N];
  real_trig = new double[N];
  imag_trig = new double[N];
  for (int n=0; n < N; n++)
    {
      real_trig[n] = cos(n * 2.0 * M_PI / N);
      imag_trig[n] = sin(n * 2.0 * M_PI / N);
      if (is_forward)
        imag_trig[n] = -imag_trig[n];
    }
}

/*****************************************************************************/
/*                   my_direct_dft::perform_transform                        */
/*****************************************************************************/

void my_direct_dft::perform_transform(float *real, float *imag, int stride)
{
  // First copy the input values to the real and imaginary temporary buffers.
  int n, k;
  float *rp, *ip;
  for (rp=real, ip=imag, n=0; n < N; n++, rp+=stride, ip+=stride)
    {real_buf[n] = *rp;  imag_buf[n] = *ip; }

  // Now compute each output coefficient in turn(vertical)
  for(int i=0; i < stride; i++)
  for (rp=real+i, ip=imag+i, k=0; k < N; k++, rp+=stride, ip+=stride)
    {
      int index = 0; // This holds n*k mod N; it indexes the trig tables
      double real_sum=0.0, imag_sum=0.0;
      for (n=0; n < N; n++, index+=k)
        {
          if (index >= N)
            index -= N;
          real_sum += real_buf[n]*real_trig[index]
                    - imag_buf[n]*imag_trig[index];
          imag_sum += real_buf[n]*imag_trig[index]
                    + imag_buf[n]*real_trig[index];
        }
      *rp = (float) real_sum;
      *ip = (float) imag_sum;
    }

  //then to horizontal
   for(int i=0; i < stride; i++)
  for (rp=real+i*stride, ip=imag+i*stride, k=0; k < N; k++, rp+=1, ip+=1)
    {
      int index = 0; // This holds n*k mod N; it indexes the trig tables
      double real_sum=0.0, imag_sum=0.0;
      for (n=0; n < N; n++, index+=k)
        {
          if (index >= N)
            index -= N;
          real_sum += real_buf[n]*real_trig[index]
                    - imag_buf[n]*imag_trig[index];
          imag_sum += real_buf[n]*imag_trig[index]
                    + imag_buf[n]*real_trig[index];
        }
      *rp = (float) real_sum;
      *ip = (float) imag_sum;
    }
}

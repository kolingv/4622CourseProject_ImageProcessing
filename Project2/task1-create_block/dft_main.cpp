/*****************************************************************************/
// File: dft_main.cpp
// Author: David Taubman
// Last Revised: 28 August, 2007
/*****************************************************************************/
// Copyright 2007, David Taubman, The University of New South Wales (UNSW)
/*****************************************************************************/

#include <math.h>
#include "io_bmp.h"
#include "image_comps.h"
#include "dft.h"
#define X_EXCEED -6
#define Y_EXCEED -7
#define NOT_AN_EVEN_NUMBER -8

/*****************************************************************************/
/*                  my_image_comp::perform_boundary_extension                */
/*****************************************************************************/

void my_image_comp::perform_boundary_extension()
{
  int r, c;
  
  //symetric expansion

  // First extend upwards
  float *first_line = buf;
  for (r=1; r <= border; r++)
    for (c=0; c < width; c++)
      first_line[-r*stride + c] = first_line[(r-1)*stride+c];

  // Now extend downwards
  float *last_line = buf+(height-1)*stride;
  for (r=1; r <= border; r++)
    for (c=0; c < width; c++)
      last_line[r*stride+c] = last_line[(1-r)*stride+c];

  // Now extend all rows to the left and to the right
  float *left_edge = buf - border*stride;
  float *right_edge = left_edge + width - 1;
  for (r=height+2*border; r > 0; r--, left_edge+=stride, right_edge+=stride)
    for (c=1; c <= border; c++)
      {
        left_edge[-c] = left_edge[c-1];
        right_edge[c] = right_edge[-c+1];
      }
}

/* ========================================================================= */
/*                              Global Functions                             */
/* ========================================================================= */
void get_block(my_image_comp *in, my_image_comp *out, int p1, int p2)
{
	int r, c;
	float *ip = in->buf + (p1 - out->width/2) + (p2 - out->width/2)*in->stride;
	//move coordinate to the top-right point of the block
	float *op = out->buf;
	
	for(r=0; r < out->height; r++ )
	{
	   for(c=0; c < out->width; c++)
	   {
		    
			op[c + r*out->stride] = ip[c + r*in->stride];
			//*op = *ip;
		   // printf("op\n");
	     
	   }
	}

}


/*****************************************************************************/
/*                                    main                                   */
/*****************************************************************************/

int
  main(int argc, char *argv[])
{
  if (argc != 6)
    {
      fprintf(stderr,"Usage: %s <in bmp file> <out bmp file> N p1 p2 \n",argv[0]);
      return -1;
    }

  int err_code=0;
  try {
      // Read the input image
      bmp_in in;
      if ((err_code = bmp_in__open(&in,argv[1])) != 0)
        throw err_code;

	  int N, p1, p2;
	  N = atoi(argv[3]);
	  p1 = atoi(argv[4]);
	  p2 = atoi(argv[5]);
      int width = in.cols, height = in.rows;
      int n, num_comps = in.num_components;

		if (N % 2 == 1)
		{
			err_code = NOT_AN_EVEN_NUMBER;
			throw err_code;
		}

		if (width-1 < p1) {
			err_code = X_EXCEED;
			throw err_code;
		}

		if (height-1 < p2) {
			err_code = Y_EXCEED;
			throw err_code;
		}


 
      my_image_comp *input_comps = new my_image_comp[num_comps];
      for (n=0; n < num_comps; n++)
        input_comps[n].init(height,width,N/2); // N/2 extension required
      
      int r, c; // Declare row index
      io_byte *line = new io_byte[width*num_comps];
      for (r=height-1; r >= 0; r--)
        { // "r" holds the true row index we are reading, since the image is
          // stored upside down in the BMP file.
          if ((err_code = bmp_in__get_line(&in,line)) != 0)
            throw err_code;
          for (n=0; n < num_comps; n++)
            {
              io_byte *src = line+n; // Points to first sample of component n
              float *dst = input_comps[n].buf + r * input_comps[n].stride;
              for (c=0; c < width; c++, src+=num_comps)
                dst[c] = (float) *src; // The cast to type "float" is not
                      // strictly required here, since bytes can always be
                      // converted to floats without any loss of information.
            }
        }
      bmp_in__close(&in);

	  int out_height = N;
	  int out_width = N;

      // Allocate storage for the output image(N*N)
      my_image_comp *output_comps = new my_image_comp[num_comps];
      for (n=0; n < num_comps; n++)
        output_comps[n].init(out_height,out_width,0); // No extension required
	  
	  //symetric extension
	   for (n=0; n < num_comps; n++)
        input_comps[n].perform_boundary_extension();

	  //get the value for the N*N block
	   for(n=0; n < num_comps; n++)
	   get_block(input_comps+n,output_comps+n, p1, p2);
	 
	   /*

      // Allocate storage for DFT buffers
      int max_dim = height;
      if (max_dim < width)
        max_dim = width;
      float *dft_real = new float[height*width];
      float *dft_imag = new float[height*width];

      // Process the image, plane by plane.
      for (n=0; n < num_comps; n++)
        {
          // First copy all samples to the `dft_real' buffer
          int stride = input_comps[n].stride;
          for (r=0; r < height; r++)
            for (c=0; c < width; c++)
              {
                dft_real[r*width+c] = input_comps[n].buf[r*stride+c];
                dft_imag[r*width+c] = 0.0F;
              }

          // Next, perform the 2D DFT

               // Put your code here


          // Write DFT magnitudes to the output image, possibly using a log
          // function to see the values more clearly.

               // Put your code here


          // Normalize the output image so that the maximum value is 255
          // and clip to avoid negative values.
          float max_val = 0.0F;
          for (r=0; r < height; r++)
            for (c=0; c < width; c++)
              {
                float val = output_comps[n].buf[r*stride+c];
                if (val > max_val)
                  max_val = val;
              }
          float scale = 1.0F;
          if (max_val > 0.0F)
            scale = 255.0F / max_val;
          for (r=0; r < height; r++)
            for (c=0; c < width; c++)
              output_comps[n].buf[r*stride+c] *= scale;
        }
		*/

      // Write the image back out again
      bmp_out out;
	  io_byte *new_line = new io_byte[out_width*num_comps];
      if ((err_code = bmp_out__open(&out,argv[2],out_width,out_height,num_comps)) != 0)
        throw err_code;
      for (r=out_height-1; r >= 0; r--)
        { // "r" holds the true row index we are writing, since the image is
          // written upside down in BMP files.
          for (n=0; n < num_comps; n++)
            {
              io_byte *dst = new_line+n; // Points to first sample of component n
              float *src = output_comps[n].buf + r * output_comps[n].stride;
              for (c=0; c < out_width; c++, dst+=num_comps)
                *dst = (io_byte) src[c];
            }
          bmp_out__put_line(&out,new_line);
        }
      bmp_out__close(&out);
      delete[] line;
      delete[] input_comps;
      delete[] output_comps;
	  delete[] new_line;
      //delete[] dft_real;
      //delete[] dft_imag;
    }
  catch (int exc) {
      if (exc == IO_ERR_NO_FILE)
        fprintf(stderr,"Cannot open supplied input or output file.\n");
      else if (exc == IO_ERR_FILE_HEADER)
        fprintf(stderr,"Error encountered while parsing BMP file header.\n");
      else if (exc == IO_ERR_UNSUPPORTED)
        fprintf(stderr,"Input uses an unsupported BMP file format.\n  Current "
                "simple example supports only 8-bit and 24-bit data.\n");
      else if (exc == IO_ERR_FILE_TRUNC)
        fprintf(stderr,"Input or output file truncated unexpectedly.\n");
      else if (exc == IO_ERR_FILE_NOT_OPEN)
        fprintf(stderr,"Trying to access a file which is not open!(?)\n");
	  else if (exc == X_EXCEED)
			fprintf(stderr, "The input location exceeds x dimension!\n");
		else if (exc == Y_EXCEED)
			fprintf(stderr, "The given location exceeds y dimension!\n");
		else if (exc == NOT_AN_EVEN_NUMBER)
			fprintf(stderr, "The given output dimension is not an even number!\n");
      return -1;
    }
  return 0;
}

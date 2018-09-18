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
#define pi 3.1415926
#define forward 1
#define inverse 0

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
	float sum=0, avg=0;
	float *ip = in->buf + (p1 - out->width/2) + (p2 - out->width/2)*in->stride;
	//move coordinate to the top-right point of the block
	float *op = out->buf;
	
	for(r=0; r < out->height; r++ )
	{
	   for(c=0; c < out->width; c++)
	   {
		    
			op[c + r*out->stride] = ip[c + r*in->stride];
			sum += op[c + r*out->stride];
			//*op = *ip;
		   // printf("op\n");
	    }
	}

	//task2
	
	avg = sum/(out->width * out->height);
	for(r=0; r < out->height; r++ )
	{
	  for(c=0; c < out->width; c++)
	  {
		op[c + r*out->stride] -= avg;
	  }
	}
	for(r=0; r < out->height; r++ )
	{
	  for(c=0; c < out->width; c++)
	  {
	      float hanning =  0.25 * (1 + cos((pi*(c - out->width/2))/(out->width/2))) 
			  * (1 + cos((pi*(r - out->width/2))/(out->width/2)));

		  op[c + r*out->stride] = op[c + r*out->stride] * hanning;
		   //printf("%f ", op[c + r*out->stride]);
	  }
	}
}


/*****************************************************************************/
/*                                    main                                   */
/*****************************************************************************/

int
  main(int argc, char *argv[])
{
  if (argc != 7)
    {
      fprintf(stderr,"Usage: %s <in bmp file> <out bmp file> N p1 p2 a \n",argv[0]);
      return -1;
    }

  int err_code=0;
  try {
      // Read the input image
      bmp_in in;
      if ((err_code = bmp_in__open(&in,argv[1])) != 0)
        throw err_code;

	  int N, p1, p2;
	  float a;
	  N = atoi(argv[3]);
	  p1 = atoi(argv[4]);
	  p2 = atoi(argv[5]);
	  a = atof(argv[6]);
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
	  int gamma_height = N+1;
	  int gamma_width = N+1;

      // Allocate storage for the output image(N*N)
      my_image_comp *output_comps = new my_image_comp[num_comps];
      for (n=0; n < num_comps; n++)
        output_comps[n].init(out_height,out_width,0); // No extension required

	  // Allocate storage for the gamma output ((N+1)*(N+1))
      my_image_comp *gamma_comps = new my_image_comp[num_comps];
      for (n=0; n < num_comps; n++)
        gamma_comps[n].init(gamma_height,gamma_width,0); // No extension required
	  
	  //symetric extension
	   for (n=0; n < num_comps; n++)
        input_comps[n].perform_boundary_extension();

	  //get the value for the N*N block
	   for(n=0; n < num_comps; n++)
	   get_block(input_comps+n,output_comps+n, p1, p2);
	 
	   

      // Allocate storage for DFT buffers(block)
      int max_dim = out_height;
      if (max_dim < out_width)
        max_dim = out_width;
      float *dft_real = new float[out_height * out_width];
      float *dft_imag = new float[out_height * out_width];

	  my_direct_dft DFT;
	  DFT.init(max_dim, forward);

      // Process the image, plane by plane.
      for (n=0; n < num_comps; n++)
        {
          // First copy all samples to the `dft_real' buffer
          int stride = output_comps[n].stride;
          for (r=0; r < out_height; r++)
            for (c=0; c < out_width; c++)
              {
                dft_real[r*out_width+c] = output_comps[n].buf[r*stride+c];
                dft_imag[r*out_width+c] = 0.0F;
              }

          // Next, perform the 2D DFT

              //first vetical 
			for(c=0; c < out_width; c++)
			{
			  float *rp, *ip;
			  rp = dft_real + c; ip = dft_imag + c;
			  DFT.perform_transform(rp, ip, stride);
			}

			//then horizontal
			for(r=0; r < out_height; r++)
			{
			  float *rp, *ip;
			  rp = dft_real + r*stride; ip = dft_imag + r*stride;
			  DFT.perform_transform(rp, ip, 1);
			}

			printf("DFT finished!\n\n");

          // square the DFT magnitude, divide by N*N

               for(r=0; r < out_height; r++)
			  {
			    for(c=0; c < out_width; c++)
				{
				   //output_comps[n].buf[r*stride+c] = log10(1 + pow(dft_real[r*width+c],2));
					output_comps[n].buf[r*stride+c] = 
						(pow(dft_real[r*out_width+c],2) + pow(dft_imag[r*out_width+c],2))/(N*N);
				  // printf("%f ",output_comps[n].buf[r*stride+c]);
				}
			  }
				printf("square finished!\n\n");
			  printf("\n");
			  printf("\n");

			  
				 
          // Normalize the output image so that the maximum value is 255
          // and clip to avoid negative values.
          float max_val = 0.0F;
          for (r=0; r < out_height; r++)
            for (c=0; c < out_width; c++)
              {
                float val = output_comps[n].buf[r*stride+c];
                if (val > max_val)
                  max_val = val;
              }
          float scale = 1.0F;
          if (max_val > 0.0F)
            scale = 255.0F / max_val;
          for (r=0; r < out_height; r++)
            for (c=0; c < out_width; c++)
              output_comps[n].buf[r*stride+c] *= scale;
		  printf("normolize finished!\n\n");
		  //re-arrange the quadrants			  
			  //float *gamma_buf = new float[(N+1)*(N+1)];
			  float *gamma = gamma_comps[n].buf + N/2*(N+1) + N/2;  // point to origin
			  //Q1
			  for(r=0; r <= N/2; r++)
				 for(c=0; c <= N/2; c++)
				 {
				   gamma[r * (N+1) + c] = a * output_comps[n].buf[r*stride+c];				   
				 }
			  //Q2
			  for(r=0; r <= N/2; r++)
				 for(c=0; c < N/2; c++)
				 {	
					gamma[r * (N+1) + (-N/2) + c] = a * 
					output_comps[n].buf[r*stride + N/2 + c];
				 }
			 //Q3
			 for(r=0; r < N/2; r++)
				 for(c=0; c <= N/2; c++)
				 {
					gamma[(r-N/2) * (N+1) + c] = a * 
					output_comps[n].buf[(r + N/2)*stride + c];
				 }
			//Q4
			 for(r=0; r < N/2; r++)
				 for(c=0; c < N/2; c++)
				 {
				   gamma[(r - N/2) * (N+1) + (-N/2) + c] = a * 
				   output_comps[n].buf[(r + N/2)*stride+ N/2 + c];				  
				 }	
			
				//calculate R,G,B
						  float R=0.0f, G=0.0f, B=0.0f;
						  float sumR1=0.0f,sumR2=0.0f, sumG1=0.0f,
							  sumG2=0.0f, sumB1=0.0f,sumB2=0.0f;
						  float avgR=0.0f, avgG=0.0f, avgB=0.0f;

						  //R
						 for(r=-N/4; r <= N/4; r++)
							 for(c=-N/4; c <= N/4; c++)					
							   sumR1 += gamma[r*gamma_width + c];					 
						 for(r=-N/8; r <= N/8; r++)
							 for(c=-N/8; c <= N/8; c++)
							   sumR2 += gamma[r*gamma_width + c];					 
						avgR = (sumR1 - sumR2)/((N/4)*(N/4)-(N/8)*(N/8));
					//	RGBoutput_comps[0].buf[r*RGBoutput_comps[0].stride + c]=avgR;
						 printf("R= %f\t", avgR);
						//G
						for(r=-N/2; r <= N/2; r++)
							for(c=-N/2; c <= N/2; c++)
								sumG1 += gamma[r*gamma_width + c];
						for(r=-N/4; r <= N/4; r++)
							for(c= -N/4; c <= N/4; c++)
								sumG2 += gamma[r*gamma_width + c];
						avgG = (sumG1 - sumG2)/((N/2)*(N/2) - (N/4)*(N/4));
						//printf("g= %.2f\t", avgG);
						//RGBoutput_comps[1].buf[r*RGBoutput_comps[1].stride + c]=avgG;
						 printf("G= %f\t", avgG);
						//B
						for(r=-N/8; r <= N/8; r++)
							for(c=-N/8; c <= N/8; c++)
							sumB1 += gamma[r*gamma_width + c];
						for(r=-N/16; r <= N/16; r++)
							for(c= -N/16; c <= N/16; c++)
							sumB2 += gamma[r*gamma_width + c];
						avgB = (sumB1 - sumB2)/((N/8)*(N/8) - (N/16)*(N/16));
						//RGBoutput_comps[2].buf[r*RGBoutput_comps[2].stride + c]=avgB;
						 printf("B= %f\n", avgB);
        }
	   
		printf("re-arrange finished!\n\n");

      // Write the image back out again
      bmp_out out;
	  io_byte *new_line = new io_byte[gamma_width * num_comps];
      if ((err_code = bmp_out__open(&out,argv[2],gamma_width,gamma_height,num_comps)) != 0)
        throw err_code;
      for (r=gamma_height-1; r >= 0; r--)
        { // "r" holds the true row index we are writing, since the image is
          // written upside down in BMP files.
          for (n=0; n < num_comps; n++)
            {
              io_byte *dst = new_line+n; // Points to first sample of component n
              float *src = gamma_comps[n].buf + r * gamma_comps[n].stride;
              for (c=0; c < gamma_width; c++, dst+=num_comps)
			  {
                if(src[c]>255) src[c]=255;
				 else if(src[c] < 0) src[c]=0;				  
				 *dst = (io_byte) src[c];
			  }
            }
          bmp_out__put_line(&out,new_line);
        }
      bmp_out__close(&out);
      delete[] line;
      delete[] input_comps;
      delete[] output_comps;
	  delete[] gamma_comps;
	  delete[] new_line;
      delete[] dft_real;
      delete[] dft_imag;
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

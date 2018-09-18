/*****************************************************************************/
// File: filtering_main.cpp
// Author: David Taubman
// Last Revised: 13 August, 2007
/*****************************************************************************/
// Copyright 2007, David Taubman, The University of New South Wales (UNSW)
/*****************************************************************************/
#include <iostream>
#include "io_bmp.h"
#include "image_comps.h"
#include <math.h>
#include <string.h>
#define pi 3.1415926



/* ========================================================================= */
/*                 Implementation of `my_image_comp' functions               */
/* ========================================================================= */

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
      first_line[-1*stride + c] = 0;

  // Now extend downwards
  float *last_line = buf+(height-1)*stride;
  for (r=1; r <= border; r++)
    for (c=0; c < width; c++)
      last_line[stride+c] = 0;

  // Now extend all rows to the left and to the right
  float *left_edge = buf - border*stride;
  float *right_edge = left_edge + width - 1;
  for (r=height+2*border; r > 0; r--, left_edge+=stride, right_edge+=stride)
    for (c=1; c <= border; c++)
      {
        left_edge[-c] = 0;
        right_edge[c] = 0;
      }
}


/* ========================================================================= */
/*                              Global Functions                             */
/* ========================================================================= */

/*****************************************************************************/
/*                                apply_filter                               */
/*****************************************************************************/

void apply_filter(my_image_comp *in, my_image_comp *out)
{
	#define FILTER_EXTENT 1
	/*
#define FILTER_TAPS (2*FILTER_EXTENT+1)

    float DCgain1=0.0F, DCgain2=0.0F;
	float *SincWindow1= new float[FILTER_TAPS];
	float *SincWindow2= new float[FILTER_TAPS];
	float tau=float(FILTER_EXTENT)+0.5;
	//delete[] SincWindow1
	for (int t=-FILTER_EXTENT;t<=FILTER_EXTENT;t++)
	{
	  SincWindow1[t+FILTER_EXTENT]=sin(pi*0.6*(float(t) - 5.0/3.0))/(pi*0.6*(float(t) - 5.0/3.0));
	 
	  SincWindow2[t+FILTER_EXTENT]=sin(pi*0.6*((float(t) - 10.0/3.0)))/(pi*0.6*((float(t) - 10.0/3.0)));
	  
	  	  float Hanning1 = 0.5*(1+cos(pi*0.6*(float(t) - 5.0/3.0)/tau));
	  float Hanning2 = 0.5*(1+cos(pi*0.6*(float(t) - 5.0/3.0)/tau));

	  SincWindow1[t+FILTER_EXTENT]=SincWindow1[t+FILTER_EXTENT]*Hanning1;
	  SincWindow2[t+FILTER_EXTENT]=SincWindow2[t+FILTER_EXTENT]*Hanning2;

	  DCgain1 += SincWindow1[t+FILTER_EXTENT];
	  DCgain2 += SincWindow2[t+FILTER_EXTENT];
	}
	printf("sincwindow\n");
  // Create the vertical filter PSF as a local array on the stack.
  //float filter_buf[FILTER_TAPS];

  

  float *filter_buf1 = new float[FILTER_TAPS];
  float *filter_buf2 = new float[FILTER_TAPS];

 float *mirror_SincWindow1 ;
   float *mirror_SincWindow2 ;

  mirror_SincWindow1 = filter_buf1 + FILTER_EXTENT;
   mirror_SincWindow2 = filter_buf2 + FILTER_EXTENT;

   // `mirror_SincWindow' points to the central tap in the filter
   for (int t=-FILTER_EXTENT; t <= FILTER_EXTENT; t++)
  {
    //float xy=0.0f;
    mirror_SincWindow1[-t]= SincWindow1[FILTER_EXTENT + t]/DCgain1;
	mirror_SincWindow2[-t]= SincWindow2[FILTER_EXTENT + t]/DCgain2;
	//xy=mirror_SincWindow1[t + FILTER_EXTENT];
  }
   printf("sincwindowmirror\n");

   */

     printf("ready for convolution!\n");
   // Check for consistent dimensions
  // Check for consistent dimensions
  assert(in->border >= FILTER_EXTENT);
  assert((out->height >= in->height) && (out->width >= in->width));

    // Perform the convolution
    int r, c;
	float *Temp_buf = new float[out->height * out->width];
	float *Temp;
	Temp = Temp_buf;  

	printf("horizontal convolution starts!\n");
	//first horizontal
	for(r=0; r < in->height; r++)
		for( c=0; 3*c <= in->width; c++)
		{		   
			float *ip = in->buf + r*in->stride + 3*c;
		  //float *op = out->buf + r*out->stride + 3*c;
		  //float sum1 = 0.0F, sum2= 0.0F;
		    Temp = Temp_buf + r*out->width + 5*c;
			*Temp=*ip;
			if(c==((in->width - 1)%3))
			{
			  int extra_b = out->width -5*c;
			  if(extra_b==1) Temp[1] = 0.4 * float(ip[0]) + 0.6 * float(ip[1]);
			  if(extra_b==2) 
			  {
			    Temp[1] = 0.4 * float(ip[0]) + 0.6 * float(ip[1]);
				Temp[2] = 0.8 * float(ip[1]) + 0.2 * float(ip[2]);
			  }
			  if(extra_b==3)
			  {
			    Temp[1] = 0.4 * float(ip[0]) + 0.6 * float(ip[1]);
				Temp[2] = 0.8 * float(ip[1]) + 0.2 * float(ip[2]);
			    Temp[3] = 0.2 * float(ip[1]) + 0.8 * float(ip[2]);
			  }
			  if(extra_b==4)
			  {
			    Temp[1] = 0.4 * float(ip[0]) + 0.6 * float(ip[1]);
				Temp[2] = 0.8 * float(ip[1]) + 0.2 * float(ip[2]);
			    Temp[3] = 0.2 * float(ip[1]) + 0.8 * float(ip[2]);
			  }
			}
		    
		    Temp[1] = 0.4 * float(ip[0]) + 0.6 * float(ip[1]);
		    Temp[2] = 0.8 * float(ip[1]) + 0.2 * float(ip[2]);
			Temp[3] = 0.2 * float(ip[1]) + 0.8 * float(ip[2]);
			Temp[4] = 0.6 * float(ip[2]) + 0.4 * float(ip[3]);
			
	            /*
		   for(int y = -FILTER_EXTENT; y <= FILTER_EXTENT; y++)
		 {
		    sum1 += ip[y]*mirror_SincWindow1[y];	//3m+1			 
			sum2 += ip[y]*mirror_SincWindow2[y];	//3m+2	
		 }
		  Temp[r*out->width + 3*c + 1] = sum1; 
		  Temp[r*out->width + 3*c + 2] = sum2;
		  */
		}

       printf("horizontal convolution finished! Then starts Temp expansion!\n");
		
		 
	 printf("Temp expansion finished! Then starts vertical convolution\n");
	//then vertical
	 for(r=0; r < out->width; r++)
		 for(c=0; 3*c < in->height-2; c++)
		 {
		   float *ip = Temp_buf + 3*c*out->width + r;
		   float *op = out->buf + 5*c*out->width +r;
		   //float sum1 = 0.0F, sum2=0.0F;

		   *op = *ip;

			if(c==((in->height - 1)%3))
			{
			  int extra_b = out->height -1 -5*c;
			  if(extra_b==1) op[out->width] = 0.4F * float(ip[out->width]) + 0.6F * float(ip[0]);;
			  if(extra_b==2) 
			  {
			     op[out->width] = 0.4F * float(ip[out->width]) + 0.6F * float(ip[0]);
		    op[2*out->width] = 0.8F * float(ip[2*out->width]) + 0.2F * float(ip[out->width]);
			  }
			  if(extra_b==3)
			  {
			 op[out->width] = 0.4F * float(ip[out->width]) + 0.6F * float(ip[0]);
		    op[2*out->width] = 0.8F * float(ip[2*out->width]) + 0.2F * float(ip[out->width]);
			op[3*out->width] = 0.2F * float(ip[2*out->width]) + 0.8F* float(ip[out->width]);
			  }
			  if(extra_b==4)
			  {
			    op[out->width] = 0.4F * float(ip[out->width]) + 0.6F * float(ip[0]);
		    op[2*out->width] = 0.8F * float(ip[2*out->width]) + 0.2F * float(ip[out->width]);
			op[3*out->width] = 0.2F * float(ip[2*out->width]) + 0.8F* float(ip[out->width]);
			op[4*out->width] = 0.6F * float(ip[3*out->width]) + 0.4F * float(ip[2*out->width]);
			  }
			}
		    //op += out->width;
		    op[out->width] = 0.4F * float(ip[out->width]) + 0.6F * float(ip[0]);
		    op[2*out->width] = 0.8F * float(ip[2*out->width]) + 0.2F * float(ip[out->width]);
			op[3*out->width] = 0.2F * float(ip[2*out->width]) + 0.8F* float(ip[out->width]);
			op[4*out->width] = 0.6F * float(ip[3*out->width]) + 0.4F * float(ip[2*out->width]);
				
			
		   /*
		   for(int y= -FILTER_EXTENT; y <= FILTER_EXTENT; y++)  
		   {
		    sum1 += ip[y*out->width]*mirror_SincWindow1[y];	//3m+1			 
			sum2 += ip[y*out->width]*mirror_SincWindow2[y];	//3m+2
		   }
		   *op = sum1; op += out->width; 	
		   *op = sum2;   
		   */
		 }	

		printf("vertical convolution finished\n");
	  	  delete[] Temp_buf;
	}

/*****************************************************************************/
/*                                    main                                   */
/*****************************************************************************/

int
  main(int argc, char *argv[])
{
  if (argc != 3)
    {
      fprintf(stderr,"Usage: %s <in bmp file> <out bmp file>\n",argv[0]);
      return -1;
    }

  
  int err_code=0;
  try {
      // Read the input image
      bmp_in in;

      if ((err_code = bmp_in__open(&in,argv[1])) != 0)
        throw err_code;

	 // int H;
	//  H = atoi(argv[4]);

      int width = in.cols, height = in.rows;
      int n, num_comps = in.num_components;
      my_image_comp *input_comps = new my_image_comp[num_comps];
      for (n=0; n < num_comps; n++)
        input_comps[n].init(height,width,1); // Leave a border of 1
      
      int r; // Declare row index
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
              for (int c=0; c < width; c++, src+=num_comps)
                dst[c] = (float) *src; // The cast to type "float" is not
                      // strictly required here, since bytes can always be
                      // converted to floats without any loss of information.
			  //sampled data stored in dst from line
            }
        }
      bmp_in__close(&in);
	    delete[] line;

	  int out_width = (int)ceilf(width * 5.0F/3.0F);
	  int out_height = (int)ceilf(height * 5.0F/3.0F);
	  	 
      // Allocate storage for the filtered output
      my_image_comp *output_comps = new my_image_comp[num_comps];
      for (n=0; n < num_comps; n++)
        output_comps[n].init(out_height,out_width,0); // Don't need a border for output

      // Process the image, all in floating point (easy)
      for (n=0; n < num_comps; n++)
        input_comps[n].perform_boundary_extension();

	  printf("Strat filtering!\n");
      for (n=0; n < num_comps; n++){
        apply_filter(input_comps+n,output_comps+n);
		printf("comp_num= %d\n", num_comps);
	  }
	  printf("Filtering finished!\n");
	  //filtered data stored in op

      // Write the image back out again
      bmp_out out;
      if ((err_code = bmp_out__open(&out,argv[2],out_width,out_height,num_comps)) != 0)
        throw err_code;
	  printf("bmp_out\n");
	   io_byte *new_line = new io_byte[out_width*num_comps];
      for (r=out_height-1; r >= 0; r--)
        { // "r" holds the true row index we are writing, since the image is
          // written upside down in BMP files.
			
          for (n=0; n < num_comps; n++)
            {

              io_byte *dst = new_line+n; // Points to first sample of component n
              float *src = output_comps[n].buf + r * output_comps[n].stride;
              for (int c=0; c < out_width; c++, dst+=num_comps)
			  {
				 if(src[c]>255) src[c]=255;
				 else if(src[c] < 0) src[c]=0;
                *dst = (io_byte) src[c];
				
			  }// The cast to type "io_byte" is
                      // required here, since floats cannot generally be
                      // converted to bytes without loss of information.  The
                      // compiler will warn you of this if you remove the cast.
                      // There is in fact not the best way to do the
                      // conversion.  You should fix it up in the lab.
            }
		  
          bmp_out__put_line(&out,new_line);
        }
	  printf("dst\n");
      bmp_out__close(&out);
	  
	
	  delete[] new_line;
      delete[] input_comps;
      delete[] output_comps;
	
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
      return -1;
    }
  return 0;
}

/*****************************************************************************/
// File: dft_main.cpp
// Author: David Taubman
// Last Revised: 28 August, 2007
/*****************************************************************************/
// Copyright 2007, David Taubman, The University of New South Wales (UNSW)
/*****************************************************************************/

#include <math.h>
#include<vector>
#include "io_bmp.h"
#include "image_comps.h"
using namespace std;
/*****************************************************************************/
/*                  my_image_comp::perform_boundary_extension                */
/*****************************************************************************/

void my_image_comp::perform_boundary_extension()
{
  int r, c;
  
  //odd symetric expansion

  // First extend upwards
  float *first_line = buf;
  for (r=1; r <= border; r++)
    for (c=0; c < width; c++)
      first_line[-r*stride + c] = first_line[r*stride+c];

  // Now extend downwards
  float *last_line = buf+(height-1)*stride;
  for (r=1; r <= border; r++)
    for (c=0; c < width; c++)
      last_line[r*stride+c] = last_line[-r*stride+c];

  // Now extend all rows to the left and to the right
  float *left_edge = buf - border*stride;
  float *right_edge = left_edge + width - 1;
  for (r=height+2*border; r > 0; r--, left_edge+=stride, right_edge+=stride)
    for (c=1; c <= border; c++)
      {
        left_edge[-c] = left_edge[c];
        right_edge[c] = right_edge[-c];
      }
}


/* ========================================================================= */
/*                              Global Functions                             */
/* ========================================================================= */
int max(int *in, int size)
{
	int i;
	int temp=0;
	for(i=0; i < size; i++)
		 if(in[i] >= temp) temp = in[i];
	return temp;
}
/*******************************************************************************/
int min(int *in, int size)
{
	int i;
	int temp=0;
	for(i=0; i < size; i++)
		 if(in[i] <= temp) temp = in[i];
	return temp;
}

/*****************************************************************************/
/*                                    main                                   */
/*****************************************************************************/

int
  main(int argc, char *argv[])
{
	
	 if (argc != 4) // pairs of elements(vectors) 
    {
      fprintf(stderr,"Usage: %s <in bmp file> <out bmp file> radius\n",argv[0]);
      return -1;
    }
	
  printf("ready\n");
  int err_code=0;
  try {
      // Read the input image
      bmp_in in;
      if ((err_code = bmp_in__open(&in,argv[1])) != 0)
        throw err_code;
		
	  printf("get radius A\n");
	  int R = atof(argv[3]) + 0.5;
	  vector<int>A_x1;
	  vector<int>A_x2;
	  int x1, x2;
	  for(x2=0; x2 < 2*R+1; x2++)  // square(A_x1)+square(A_x2) < square(R+0.5)
		  for(x1=0; x1 < 2*R+1; x1++)
		  {
		    float r_square = (x1 - R)*(x1 - R) + (x2 - R)*(x2 - R);
			if(r_square <= (atof(argv[3]) + 0.5)*(atof(argv[3]) + 0.5)) 
			{
			  A_x1.push_back(x1 - R);
			  A_x2.push_back(x2 - R);
			}
		  }
	   printf("radius A got\n");
	  int aoff_size = A_x1.size(); 
	  int *aoff = new int[aoff_size];
	  /*
	  int aoff_size = (argc - 3)/2; 
	  int *aoff = new int[aoff_size];
	  int *A_x1 = new int[aoff_size];
	  int *A_x2 = new int[aoff_size];
		printf("ready for separate \n");
	  int i, j=0;
	  for(i=0; i < aoff_size; i++, j += 2)
		  A_x1[i] = atoi(argv[3+i]);
	  j=0;
	  for(i=0; i < aoff_size; i++, j += 2)
		  A_x2[i] = atoi(argv[4+i]);
	  printf("separate ready\n");
	  */


	  //calculate the border
		/*  
	  int upper_border, lower_border, left_border, right_border, border;
	   upper_border = abs(min(A_x1, aoff_size));
	   lower_border = max(A_x1, aoff_size);
	   left_border  = abs(min(A_x2, aoff_size));
	   right_border = max(A_x2, aoff_size);
	   int temp[4]={upper_border, lower_border, left_border, right_border};	  
	   border = max(temp, 4);
	   printf("border calculation finish\n");
	   */

      int width = in.cols, height = in.rows;
      int n, num_comps = in.num_components;	  
      my_image_comp *input_comps = new my_image_comp[num_comps];
      for (n=0; n < num_comps; n++)
        input_comps[n].init(height,width,R); // border extension required(R) 
      	 
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
	 
	  int i;
	  //get aoff
	  printf("get aoff\n");
	  for (n=0; n < num_comps; n++)
		  for(i=0; i < aoff_size; i++)
		  {
			aoff[i] = A_x1[i] * input_comps[n].stride + A_x2[i];
			printf("aoff= %d\t", aoff[i]);
		  }
		  

	  //symetric extension
	   for (n=0; n < num_comps; n++)
        input_comps[n].perform_boundary_extension();
	  	
	    // Allocate storage for the output image
	   int out_height = height;
	  int out_width = width;
      my_image_comp *output_comps = new my_image_comp[num_comps];
      for (n=0; n < num_comps; n++)
        output_comps[n].init(out_height,out_width,0); // No extension required
		printf("ready for erosion\n");
	   for (n=0; n < num_comps; n++)
		   for(r=0; r < height; r++)
			   for(c=0; c < width; c++)
			   {
			     float *pn = input_comps[n].buf + r * input_comps[n].stride + c;
				 int val = 255;
				 for(i=0; i < aoff_size; i++)
				 {	
					//printf("in= %d", pn[aoff[i]]);
				   val &= int(pn[aoff[i]]);
				 }
				 output_comps[n].buf[r * output_comps[n].stride + c] = val;
			   }
	  printf("erosion finish\n");




      // Write the image back out again
      bmp_out out;
      if ((err_code = bmp_out__open(&out,argv[2],out_width,out_height,num_comps)) != 0)
        throw err_code;
      for (r=out_height-1; r >= 0; r--)
        { // "r" holds the true row index we are writing, since the image is
          // written upside down in BMP files.
          for (n=0; n < num_comps; n++)
            {
              io_byte *dst = line+n; // Points to first sample of component n
              float *src = output_comps[n].buf + r * output_comps[n].stride;
              for (c=0; c < width; c++, dst+=num_comps)
			  {				 
				 *dst = (io_byte) src[c];
			  }
                
            }
          bmp_out__put_line(&out,line);
        }
      bmp_out__close(&out);
      delete[] line;
      delete[] input_comps;
	  delete[] output_comps;	  
      delete[] aoff;
      //A_x1.clear;
	 // A_x2.clear;
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

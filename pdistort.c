/* The basic idea is to take a source image (PNG for now) and apply a distortion
   field that comes from random Perlin noise that twists and bends it in many
   directions.
   Each pixel of the destination image is computed like this:
   x1 = destpixel.x; y1 = destpixel.y
   k1, k2 exist as parameters (any floats, really)
   distort_amt exists as a parameter (float from 0 to inf but probably <1)
   destPixel.color =
     srcImage[x1+perlin1D(x1/width + k1)*distort_amt]
             [y1+perlin1D(y1/height + k2)*distort_amt]
*/

#include <stdio.h>
#include <stdlib.h>
#include <png.h>
#include <math.h>
#include <stdbool.h>
#include <perlin.h>

typedef unsigned int uint;
typedef unsigned char uchar;

// avoid compiler error with setjmp
#define PNG_SETJMP_NOT_SUPPORTED

/* Given an (x,y) coordinate as floats, return a RGBA pixel whose values are
   linearly interpolated from the surrounding pixels.
   Pixels outside the bounds of the image are considered colorless (not
   contributing to the rest of the colors' summation) and transparent.
   Expects a uchar[4] to be passed in to hold the resulting pixel values.
*/
void getPixel(float x, float y, png_bytep* rows, png_uint_32 width,
	      png_uint_32 height, png_byte* results) {
  if((x <= -1) || (x >= width) || (y <= -1) || (y >= height)) {
    //results[0], [1], [2] are irrelevant
    results[3] = 0x00; // transparent
    return;
  }
  int xlow = (int) floor(x);
  int ylow = (int) floor(y);
  float fracx = x-xlow;
  float fracy = y-ylow;
  float oppx = 1-fracx;
  float oppy = 1-fracy;
  uint rowbytes = width*4;
  float tmp[4]; //results get stored here before being rounded and cast to uchar
  bool corner = false;
  
  /* oppx and fracx are complements to each other.
     +---------+
     |      |  |
     |      |  |
     |      |  |
     |      |  |   fracx = distance from left = 0.7
     |      |  |   oppx = distance from right = 0.3
     |      |  |   fracy = distance from top = 0.8
     |      |  |   oppy = distance from bottom = 0.2
     |------+--|
     |      |  |
     +---------+
  */

  // The only cases where results should be directly modified are the literal
  // corner cases, where the colors just take the value of one pixel.
  // All other calculations should set tmp.
  
  if(ylow == -1) { // the top pixels are off the edge
    if(xlow == -1) {
      // top left corner:
      // just the color of the one pixel mixed with transparency
      // transparent mixing: oppx for horizontal, oppy for vertical
      corner = true;
      results[0] = rows[0][0];
      results[1] = rows[0][1];
      results[2] = rows[0][2];
      tmp[3] = rows[0][3] * oppx * oppy;
    }
    else if(xlow == width-1) {
      // top right corner
      // transparent mixing: fracx * oppy
      corner = true;
      results[0] = rows[0][(width*4)-4];
      results[1] = rows[0][(width*4)-3];
      results[2] = rows[0][(width*4)-2];
      tmp[3] = rows[0][(width*4)-1] * fracx * oppy;
    }
    else {
      // x is normal but ylow is still off the image
      tmp[0] = rows[0][xlow*4]*fracx + rows[0][xlow*4+4]*oppx;
      tmp[1] = rows[0][xlow*4+1]*fracx + rows[0][xlow*4+5]*oppx;
      tmp[2] = rows[0][xlow*4+2]*fracx + rows[0][xlow*4+6]*oppx;
      tmp[3] = (rows[0][xlow*4+3]*fracx + rows[0][xlow*4+7]*oppx)*oppy;
    }
  } else if (ylow == height-1) { // the bottom pixels are off the edge
    if(xlow == -1) {
      // bottom left corner
      corner = true;
      results[0] = rows[height-1][0];
      results[1] = rows[height-1][1];
      results[2] = rows[height-1][2];
      tmp[3] = rows[height-1][3] * oppx * fracy;
    }
    else if(xlow == width-1) {
      // bottom right corner
      corner = true;
      results[0] = rows[height-1][(width*4)-4];
      results[1] = rows[height-1][(width*4)-3];
      results[2] = rows[height-1][(width*4)-2];
      tmp[3] = rows[height-1][(width*4)-1] * fracx * fracy;
    }
    else {
      // x is normal but ylow+1 is still off the image
      tmp[0] = rows[height-1][xlow*4]*fracx + rows[height-1][xlow*4+4]*oppx;
      tmp[1] = rows[height-1][xlow*4+1]*fracx + rows[height-1][xlow*4+5]*oppx;
      tmp[2] = rows[height-1][xlow*4+2]*fracx + rows[height-1][xlow*4+6]*oppx;
      tmp[3] = (rows[height-1][xlow*4+3]*fracx + rows[height-1][xlow*4+7]*oppx)
	* fracy;
    }
  } else { // y is normal
    if(xlow == -1) {
      // y is normal but the left pixels are off the edge
      tmp[0] = rows[ylow][0]*fracy + rows[ylow+1][0]*oppy;
      tmp[1] = rows[ylow][1]*fracy + rows[ylow+1][1]*oppy;
      tmp[2] = rows[ylow][2]*fracy + rows[ylow+1][2]*oppy;
      tmp[3] = (rows[ylow][3]*fracy + rows[ylow+1][3]*oppy) * oppx;
    }
    else if(xlow == width-1) {
      // y is normal but the right pixels are off the edge
      tmp[0] = (rows[ylow][rowbytes-4]*fracy +
			    rows[ylow+1][rowbytes-4]*oppy);
      tmp[1] = (rows[ylow][rowbytes-3]*fracy +
			    rows[ylow+1][rowbytes-3]*oppy);
      tmp[2] = (rows[ylow][rowbytes-2]*fracy +
			    rows[ylow+1][rowbytes-2]*oppy);
      tmp[3] = (rows[ylow][rowbytes-1]*fracy +
			    rows[ylow+1][rowbytes-1]*oppy) * fracx;
    }
    else {
      // normal interpolation between normal pixels
      tmp[0] = (rows[ylow][4*xlow]*fracx + rows[ylow][4*xlow+4]*oppx)*fracy +
	(rows[ylow+1][4*xlow]*fracx + rows[ylow+1][4*xlow+4]*oppx) * oppy;
      
      tmp[1] = (rows[ylow][4*xlow+1]*fracx + rows[ylow][4*xlow+5]*oppx)*fracy +
	(rows[ylow+1][4*xlow+1]*fracx + rows[ylow+1][4*xlow+5]*oppx) * oppy;
      
      tmp[2] = (rows[ylow][4*xlow+2]*fracx + rows[ylow][4*xlow+6]*oppx)*fracy +
	(rows[ylow+1][4*xlow+2]*fracx + rows[ylow+1][4*xlow+6]*oppx) * oppy;
      
      tmp[3] = (rows[ylow][4*xlow+3]*fracx + rows[ylow][4*xlow+7]*oppx)*fracy +
	(rows[ylow+1][4*xlow+3]*fracx + rows[ylow+1][4*xlow+7]*oppx) * oppy;
    }
  }
  if(!corner) {
    results[0] = (uchar) lroundf(tmp[0]);
    results[1] = (uchar) lroundf(tmp[1]);
    results[2] = (uchar) lroundf(tmp[2]);
  }
  results[3] = (uchar) lroundf(tmp[3]);
}

void usage(char* argv0) {
  /* usage: argv[0] srcImg dstImg kDistort [-P obegin oend persistence] [-k k1 k2]
     defaults: obegin = 0, oend = 3, persistence = 0.5, k1 = 0.13, k2 = 0.47 */
  fprintf(stderr, "Usage: %s sourceImg destImg distort-coefficient\n", argv0);
  fprintf(stderr, "          [-P obegin oend persistence]\n");
  fprintf(stderr, "          [-k constant1 constant2]\n");
  fprintf(stderr, "All values are floats except for the image filename (string) and obegin and oend (ints).\n");
}

int main(int argc, char* argv[]) {
  char* srcImg = NULL;
  char* destImg = NULL;
  float kDistort = 0;
  uchar obegin = 0;
  uchar oend = 4;
  float persistence = 0.3;
  float k1 = 0.13;
  float k2 = -0.47;
  float k3 = -2.63;
  float k4 = 1.89;
  char* check;
  //set up perlin seed
  mytime = time(0);
  
  if(argc < 4) {
    usage(argv[0]);
    return 1;
  } else {
    // convert image and distort coefficient
    srcImg = argv[1];
    destImg = argv[2];
    kDistort = strtof(argv[3], &check);
    if(*check != '\0') {
      fprintf(stderr, "distort-coefficient must be a valid float\n");
      return 1;
    }

    // check for extra parameters
    if(argc > 4) {
      if(!strcmp(argv[4], "-P")) {
	if(argc < 8) {
	  usage(argv[0]);
	  return 1;
	}
	obegin = (uchar) strtoul(argv[5], &check, 10);
	if(*check != '\0') {
	  fprintf(stderr, "obegin must be a valid int\n");
	  return 1;
	}
	oend = (uchar) strtoul(argv[6], &check, 10);
	if(*check != '\0') {
	  fprintf(stderr, "oend must be a valid int\n");
	  return 1;
	}
	persistence = strtof(argv[7], &check);
	if(*check != '\0') {
	  fprintf(stderr, "persistence must be a valid float\n");
	  return 1;
	}
      }
    }	  
  }
  // all parameters have been set

  // open source image
  FILE* srcFile = fopen(srcImg, "rb");
  if(!srcFile) {
    perror(srcImg);
    return 1;
  }
  // verify that source image is actually a PNG
  png_byte header[8];
  fread(header, 1, 8, srcFile);
  if (png_sig_cmp(header, 0, 8)) {
    fprintf(stderr, "error: %s is not a PNG.\n", srcImg);
    fclose(srcFile);
    return 1;
  }

  // create necessary PNG structures for reading
  png_structp png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
  if (!png_ptr) {
    fprintf(stderr, "error: png_create_read_struct returned 0.\n");
    fclose(srcFile);
    return 1;
  }
  png_infop info_ptr = png_create_info_struct(png_ptr);
  if (!info_ptr) {
    fprintf(stderr, "error: png_create_info_struct returned 0.\n");
    png_destroy_read_struct(&png_ptr, NULL, NULL);
    fclose(srcFile);
    return 1;
  }
  // the code in this if statement gets called if libpng encounters an error
  if (setjmp(png_jmpbuf(png_ptr))) {
    fprintf(stderr, "error from libpng\n");
    png_destroy_info_struct(png_ptr, &info_ptr);
    png_destroy_read_struct(&png_ptr, NULL, NULL);
    fclose(srcFile);
    return 1;
  }

  // vars for image info
  png_uint_32 width, height;
  int bit_depth, colortype;
  
  // start reading image metadata
  png_init_io(png_ptr, srcFile);
  // let libpng know you already read the first 8 bytes
  png_set_sig_bytes(png_ptr, 8);
  // read all the info up to the image data
  png_read_info(png_ptr, info_ptr);
  // transfer that into the variables
  width = png_get_image_width(png_ptr, info_ptr);
  height = png_get_image_height(png_ptr, info_ptr);
  colortype = png_get_color_type(png_ptr, info_ptr);
  bit_depth = png_get_bit_depth(png_ptr, info_ptr);

  // Output will be RGBA, so force any image to be read as RGBA.
  // Source: https://gist.github.com/niw/5963798
  if(bit_depth == 16)
    png_set_strip_16(png_ptr);
  if(colortype == PNG_COLOR_TYPE_PALETTE)
    png_set_palette_to_rgb(png_ptr);
  if(colortype == PNG_COLOR_TYPE_GRAY && bit_depth < 8)
    png_set_expand_gray_1_2_4_to_8(png_ptr);
  if(png_get_valid(png_ptr, info_ptr, PNG_INFO_tRNS))
    png_set_tRNS_to_alpha(png_ptr);
  if(colortype == PNG_COLOR_TYPE_RGB ||
     colortype == PNG_COLOR_TYPE_GRAY ||
     colortype == PNG_COLOR_TYPE_PALETTE)
    png_set_filler(png_ptr, 0xFF, PNG_FILLER_AFTER);
  if(colortype == PNG_COLOR_TYPE_GRAY ||
     colortype == PNG_COLOR_TYPE_GRAY_ALPHA)
    png_set_gray_to_rgb(png_ptr);
  png_read_update_info(png_ptr, info_ptr);

  png_bytep* row_pointers = (png_bytep*) malloc(sizeof(png_bytep) * height);
  uint x, y;
  uint rowsize = png_get_rowbytes(png_ptr, info_ptr);
  for(y=0; y<height; ++y) {
    row_pointers[y] = (png_byte*) malloc(rowsize);
  }
  png_read_image(png_ptr, row_pointers);
  png_destroy_info_struct(png_ptr, &info_ptr);
  png_destroy_read_struct(&png_ptr,NULL,NULL);
  fclose(srcFile);

  //Now row_pointers[y][x*4] = (x,y) red, [y][x*4+1] = (x,y) green,
  //[y][x*4+2] = (x,y) blue, [y][x*4+3] = (x,y) alpha

  png_bytep* newimg = (png_bytep*) malloc(sizeof(png_bytep) * height);
  for(y=0; y<height; ++y) {
    newimg[y] = (png_byte*) malloc(rowsize);
  }
  uchar tmp_pixel[4];
  float hashx, hashy;
  for(y=0; y<height; ++y) {
    for(x=0; x<width; ++x) {
      float invX = (float)x/width;
      float invY = (float)y/height;
      hashx = (invX +
	       perlin_2d(invX+k1, invY+k2, persistence, obegin, oend)*kDistort)
	* width;
      hashy = (invY +
	       perlin_2d(invX+k3, invY+k4, persistence, obegin, oend)*kDistort)
	* height;
      getPixel(hashx, hashy, row_pointers, width, height, tmp_pixel);
      newimg[y][x*4] = tmp_pixel[0];
      newimg[y][x*4+1] = tmp_pixel[1];
      newimg[y][x*4+2] = tmp_pixel[2];
      newimg[y][x*4+3] = tmp_pixel[3];
      //printf("%f * %d = %f\n", invX, width, hashx);
      //printf("%d %d -> %f %f\n", x, y, hashx, hashy);
    }
  }
  
  // open destination image
  FILE* destFile = fopen(destImg, "wb");
  if(!destFile) {
    perror(destImg);
    return 1;
  }

  // create necessary PNG structures for writing
  png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
  if (!png_ptr) {
    fprintf(stderr, "error: png_create_write_struct returned 0.\n");
    fclose(destFile);
    return 1;
  }
  info_ptr = png_create_info_struct(png_ptr);
  if (!png_ptr) {
    fprintf(stderr, "error: png_create_info_struct returned 0.\n");
    png_destroy_read_struct(&png_ptr, NULL, NULL);
    fclose(destFile);
    return 1;
  }
  if (setjmp(png_jmpbuf(png_ptr))) {
    fprintf(stderr, "error from libpng\n");
    png_destroy_info_struct(png_ptr, &info_ptr);
    png_destroy_read_struct(&png_ptr, NULL, NULL);
    fclose(destFile);
    return 1;
  }
  // start writing metadata
  png_init_io(png_ptr, destFile);
  png_set_IHDR(png_ptr, info_ptr, width, height, 8 /* bit depth */,
	       PNG_COLOR_TYPE_RGBA, PNG_INTERLACE_NONE,
	       PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
  png_write_info(png_ptr, info_ptr);

  // write in image
  png_write_image(png_ptr, newimg);
  png_write_end(png_ptr, NULL);

  // clean up png structs
  png_destroy_info_struct(png_ptr, &info_ptr);
  //png_destroy_read_struct(&png_ptr,NULL,NULL); //why throws memory error?

  // deallocate memory and close file
  for(y=0; y<height; ++y) {
    free(row_pointers[y]);
    free(newimg[y]);
  }
  free(row_pointers);
  free(newimg);
  fclose(destFile);
  
  
  /*
  uchar stuff[4];
  getPixel(-1.5, 3.5, row_pointers, width, height, stuff);
  printf("%d %d %d %d\n", stuff[0], stuff[1], stuff[2], stuff[3]);
  */
  
  return 0;
}

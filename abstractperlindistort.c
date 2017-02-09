#include <perlin.h>
#include <quickpng.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <limits.h>

#define SIZE 400
float height[SIZE][SIZE];
float q[SIZE][SIZE];

typedef unsigned char uchar;

typedef struct {
  uchar r;
  uchar g;
  uchar b;
} color;

//between should fall in [0,1]
color mix (color c1, color c2, float between) {
  color out;
  out.r = ((c2.r - c1.r) * between) + c1.r;
  out.g = ((c2.g - c1.g) * between) + c1.g;
  out.b = ((c2.b - c1.b) * between) + c1.b;
  return out;
}

float ridgenoise(float x, float y, int octaves) {
  float tmp = perlin_2d(x, y, 0.5, 0, octaves);
  if(tmp < 0) tmp = -tmp;
  return 1 - tmp;
}

typedef union {
  unsigned int i;
  float f;
} uni;
typedef unsigned int uint;

int main(int argc, char* argv[]) {
  mytime = time(0);
  float min = INFINITY;
  float max = -INFINITY;
  float tmp;
  float op;

  float qmax = 0;
  float x1, x2, y1, y2, xfin, yfin;
  const float scalecontrib = 0.25;
  quickpng_init(SIZE, SIZE);
  int i,j;
  for(i=0; i<SIZE; ++i) {
    for(j=0; j<SIZE; ++j) {
      x1 = (float)j/SIZE;
      y1 = (float)i/SIZE;
      x2 = x1 + perlin_2d(x1+0.46, y1+0.04, 0.5, 0, 5)*scalecontrib;
      y2 = y1 + perlin_2d(x1+0.02, y1+0.31, 0.5, 0, 5)*scalecontrib;
      //can add further layers of distortion but it is pretty costly
      xfin = x2*8;
      yfin = y2*8;
      tmp = perlin_2d(xfin, yfin, 0.5, 0, 5);
      q[i][j] = sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
      if(q[i][j] > qmax)
	qmax = q[i][j];
      min = (tmp < min ? tmp : min);
      max = (tmp > max ? tmp : max);
      height[i][j] = tmp;
    }
  }
  float b, a;
  color base1 = {0,0,0};//{40, 40, 100};
  color base2 = {210, 210, 255};
  color alt1 = {255,0,0};
  color tmpc;
  for(i=0; i<SIZE; ++i) {
    for(j=0; j<SIZE; ++j) {
      imgA[i][j] = 255;
      b = (height[i][j] - min) / (max - min);
      tmpc = mix(base1, base2, b);
      //a = q[i][j]/qmax;
      //b = sinf(b * a);
      //tmpc = mix(tmpc, alt1, b);
      imgR[i][j] = tmpc.r;
      imgG[i][j] = tmpc.g;
      imgB[i][j] = tmpc.b;
    }
  }
  quickpng_write("hmap.png");
  
  quickpng_destroy();
  
}

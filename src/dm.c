#include<stdio.h>
#include<math.h>
#include"global.h"

double dm(double z)
{
  int i;
  double t = z/1000.0;
  double sum = 0;
  double w;

  for(i=0;i<1000;i++){
      w = (i+0.5)*t;
      sum += 1.0*t/sqrt(OM*(1.0+w)*(1.0+w)*(1.0+w)+OA);
                     }
  sum = sum*LIGHT_SPEED/(H0*(1.0+z));
  return sum;
}


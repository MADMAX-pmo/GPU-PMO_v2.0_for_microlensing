#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include"global.h"

void lens_salpetermf(int nlens,Ls *Lens,struct Ip ImagePlane,struct Mf MassFunction,double *real_massave)
{
 int i,j;
 srand(time(NULL));
 
 double Ixmin = ImagePlane.xmin;
 double Ixmax = ImagePlane.xmax;
 double Iymin = ImagePlane.ymin;
 double Iymax = ImagePlane.ymax;
 double power1 = MassFunction.power+1.0;
 double a2 = -pow(MassFunction.maxmass,power1)/power1;
 double a1 = -pow(MassFunction.minmass,power1)/power1 - a2;
 double base;
 double massall=0;
 
 for(i=0;i<nlens;i++){
     Lens[i].lensx = (rand()/(double)RAND_MAX)*(Ixmax-Ixmin)+Ixmin;
     Lens[i].lensy = (rand()/(double)RAND_MAX)*(Iymax-Iymin)+Iymin;        
     base = -power1 * (a1*rand()/(RAND_MAX+1.0) + a2);     
     //Lens[i].lensx = (double)(drand48()*2.0-1.0) * Ixmax;
     //Lens[i].lensy = (double)(drand48()*2.0-1.0) * Iymax;
     //base = -power1 * (a1*drand48() + a2);
 
     Lens[i].lensm = pow(base,(1.0/power1));
     massall += Lens[i].lensm;
 }

 double real_ave = massall/nlens;
 *real_massave = real_ave;
}

double massave(struct Mf MassFunction)
{
 double power1 = MassFunction.power+1.0;
 double power2 = MassFunction.power+2.0;
 double b2 = (pow(MassFunction.maxmass,power2)-pow(MassFunction.minmass,power2))/power2;
 double b1 = (pow(MassFunction.maxmass,power1)-pow(MassFunction.minmass,power1))/power1;
 double Ave = b2/b1;
 return Ave;
 }




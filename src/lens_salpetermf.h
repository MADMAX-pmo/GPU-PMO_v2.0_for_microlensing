#ifndef LENSPOS_SALPETERMF_H
#define LENSPOS_SALPETERMF_H
#include"global.h"
void lens_salpetermf(int nlens,Ls *Lens,struct Ip ImagePlane,struct Mf MassFunction,double *real_massave);
double massave(struct Mf MassFunction);
#endif

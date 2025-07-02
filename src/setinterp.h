#ifndef __SETINTERP_H
#define __SETINTERP_H
#include"global.h"

#ifdef __cplusplus
extern "C"{
#endif
    void setinterp(Double8 *far_lenses,int nlens,double mass_ave,Ls *Lens,struct Pa Parameters,struct Ip ImagePlane,int *hoc,int *ll,int dev);
#ifdef __cplusplus
    }
#endif
#endif

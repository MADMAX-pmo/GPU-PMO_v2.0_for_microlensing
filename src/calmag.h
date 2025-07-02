#ifndef __CALMAG_H
#define __CALMAG_H
#include"global.h"

#ifdef __cplusplus
extern "C"{
#endif
   
    void calmag(int *num_map,int *hoc,int *ll,Double8 *far_lenses,int nlens,Ls *Lens,double mass_ave,struct Pa Parameters,struct Sp SourcePlane,struct Ip ImagePlane,int dev); 
#ifdef __cplusplus
    }
#endif
#endif

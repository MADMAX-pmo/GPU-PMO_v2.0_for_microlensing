#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include"global.h"
	
void hocll(int nlens,Ls *Lens,struct Ip ImagePlane,int *ll,int *hoc){

 int i,q1,q2; 
 int nx = ImagePlane.l1nx;
 int ny = ImagePlane.l1ny;
 double xmin = ImagePlane.xmin;
 double ymin = ImagePlane.ymin;
 double dx = ImagePlane.l1xgridscale;
 double dy = ImagePlane.l1ygridscale;

 memset(hoc, -1, sizeof(int) * nx * ny);

 for(i=0;i<nlens;i++){
     q1 = (int)((Lens[i].lensx - xmin) / dx);
     q2 = (int)((Lens[i].lensy - ymin) / dy); 
     
     if(q1<0 || q1>=nx || q2<0 || q2>=ny){
         fprintf(stderr, "[Warning] lens[%d] out of grid: (%.3f, %.3f) → (%d, %d)\n",i, Lens[i].lensx, Lens[i].lensy, q1, q2);
         continue;      
     }
    
     ll[i] = hoc[q1 * ny + q2];
     hoc[q1 * ny + q2] = i;
 }

}
  
  
  


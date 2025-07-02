#include<stdio.h>
#include<stdlib.h>
#include<cuda_runtime.h>
#include<math.h>
#include"calmag.h"
#include"setinterp.h"
#include"global.h"


__global__ void Cal_Coe(
    Double8 *d_far_lenses,int *d_ll,int *d_hoc,Ls *d_Lens,
    double mass_ave, double kstar,
    double Ixmin, double Ixmax, double Iymin, double Iymax,
    int nx1, int ny1, int nx2, int ny2,
    double dx1, double dy1, double dx2, double dy2)
{
  int i=blockDim.x*blockIdx.x+threadIdx.x;
  int j=blockDim.y*blockIdx.y+threadIdx.y;
  
  if (i >= nx2 || j >= ny2) return;

  double pi=acos(-1.0);
  double Ix = Ixmax+1.0e-20;
  double Iy = Iymax+1.0e-20;
 
  double posx = Ixmin + (i + 0.5) * dx2;
  double posy = Iymin + (j + 0.5) * dy2;

  double px[4] = {
         posx + 0.5 * dx2, //px[0]
         posx - 0.5 * dx2, //px[1]
         posx - 0.5 * dx2, //px[2]
         posx + 0.5 * dx2  //px[3]
  };
  double py[4] = {
         posy + 0.5 * dy2,
         posy + 0.5 * dy2,
         posy - 0.5 * dy2, 
         posy - 0.5 * dy2
  };
  
  int cx = (int)((posx - Ixmin) / dx1);
  int cy = (int)((posy - Iymin) / dy1);
 
  double alpha[8] = {0};
  double alphax_corr[4],alphay_corr[4];
  double inv_mass = 1.0/mass_ave;

  for(int ii=0;ii<nx1;ii++){
       for(int jj=0;jj<ny1;jj++){
            if (abs(ii - cx) > 1 || abs(jj - cy) > 1){
                int it = d_hoc[ii*ny1+jj];
                while(it!=-1){                          
                    double far_lensx = d_Lens[it].lensx;
                    double far_lensy = d_Lens[it].lensy;
                    double far_lensm = d_Lens[it].lensm*inv_mass;

                    for (int s = 0; s < 4; ++s) {
                        double dx = px[s] - far_lensx;
                        double dy = py[s] - far_lensy;
                        double r2 = dx * dx + dy * dy;
                        double inv_r2 = 1.0 / r2;
                        alpha[2 * s + 0] += far_lensm * dx * inv_r2;
                        alpha[2 * s + 1] += far_lensm * dy * inv_r2;
                    }
                 
                    it=d_ll[it];
                      
                }
            }
       }
  }
  
  for(int s=0;s<4;s++){
      double p0=px[s];
      double p1=py[s];

      double am_p0 = Ix-p0, ap_p0 = Ix+p0;
      double bm_p1 = Iy-p1, bp_p1 = Iy+p1;

      double ax1 = 0.5*bm_p1*log((am_p0*am_p0+bm_p1*bm_p1)/(ap_p0*ap_p0+bm_p1*bm_p1));
      double ax2 = 0.5*bp_p1*log((am_p0*am_p0+bp_p1*bp_p1)/(ap_p0*ap_p0+bp_p1*bp_p1));
      double ax3 = am_p0*(atan(bm_p1/am_p0)+atan(bp_p1/am_p0));
      double ax4 = ap_p0*(atan(bm_p1/ap_p0)+atan(bp_p1/ap_p0));
     
      double bx1 = 0.5*am_p0*log((am_p0*am_p0+bm_p1*bm_p1)/(am_p0*am_p0+bp_p1*bp_p1));
      double bx2 = 0.5*ap_p0*log((ap_p0*ap_p0+bm_p1*bm_p1)/(ap_p0*ap_p0+bp_p1*bp_p1));
      double bx3 = bm_p1*(atan(am_p0/bm_p1)+atan(ap_p0/bm_p1));
      double bx4 = bp_p1*(atan(am_p0/bp_p1)+atan(ap_p0/bp_p1));

      alphax_corr[s] = (ax1+ax2+ax3-ax4)*(-kstar/pi);
      alphay_corr[s] = (bx1+bx2+bx3-bx4)*(-kstar/pi);
  }

          
      int idx = i*ny2+j;

       d_far_lenses[idx].a0 = alpha[0]-alphax_corr[0];
       d_far_lenses[idx].a1 = alpha[1]-alphay_corr[0]+2*kstar*(py[0]-posy);
  
       d_far_lenses[idx].a2 = alpha[2]-alphax_corr[1];  
       d_far_lenses[idx].a3 = alpha[3]-alphay_corr[1]+2*kstar*(py[1]-posy);  

       d_far_lenses[idx].a4 = alpha[4]-alphax_corr[2];  
       d_far_lenses[idx].a5 = alpha[5]-alphay_corr[2]+2*kstar*(py[2]-posy);
  
       d_far_lenses[idx].a6 = alpha[6]-alphax_corr[3];  
       d_far_lenses[idx].a7 = alpha[7]-alphay_corr[3]+2*kstar*(py[3]-posy);
  
       d_far_lenses[idx].xc = posx;  
       d_far_lenses[idx].yc = posy;
}



void setinterp(
    Double8 *far_lenses,int nlens,double mass_ave,Ls *Lens,
    struct Pa Parameters,struct Ip ImagePlane,
    int *hoc,int *ll,int dev)
{
  cudaSetDevice(dev);
  cudaError_t err;

  int nx1 = ImagePlane.l1nx;
  int ny1 = ImagePlane.l1ny;
  int nx2 = ImagePlane.l2nx;
  int ny2 = ImagePlane.l2ny;

  Ls *d_Lens;
  cudaMalloc(&d_Lens,nlens*sizeof(Ls));

  Double8 *d_far_lenses;
  cudaMalloc(&d_far_lenses,nx2*ny2*sizeof(Double8));

  int *d_hoc;
  cudaMalloc(&d_hoc,nx1*ny1*sizeof(int));

  int *d_ll;
  cudaMalloc(&d_ll,nlens*sizeof(int));


  err=cudaMemcpy(d_Lens,Lens,nlens*sizeof(Ls),cudaMemcpyHostToDevice);
  if(err != cudaSuccess){
      fprintf(stderr, "cudaMemcpy in setinterp Lens to d_Lens: %s\n", cudaGetErrorString(cudaGetLastError()));
      exit(EXIT_FAILURE);
                        }
  
  err=cudaMemcpy(d_hoc,hoc,nx1*ny1*sizeof(int),cudaMemcpyHostToDevice);
  if(err != cudaSuccess){
    fprintf(stderr, "cudaMemcpy in setinterp hoc to d_hoc: %s\n", cudaGetErrorString(cudaGetLastError()));
    exit(EXIT_FAILURE);  
                        }
  
  err=cudaMemcpy(d_ll,ll,nlens*sizeof(int),cudaMemcpyHostToDevice);
  if(err != cudaSuccess){
    fprintf(stderr, "cudaMemcpy in setinterp ll to d_ll: %s\n", cudaGetErrorString(cudaGetLastError()));
    exit(EXIT_FAILURE);  
                        }

  dim3 blocksize(2,32);
  dim3 gridsize((nx2+blocksize.x-1)/blocksize.x,(ny2+blocksize.y-1)/blocksize.y);
  
  Cal_Coe<<<gridsize,blocksize>>>(
      d_far_lenses,d_ll,d_hoc,d_Lens,
      mass_ave,Parameters.kstar,
      ImagePlane.xmin, ImagePlane.xmax,
      ImagePlane.ymin, ImagePlane.ymax,
      ImagePlane.l1nx, ImagePlane.l1ny,
      ImagePlane.l2nx, ImagePlane.l2ny,
      ImagePlane.l1xgridscale, ImagePlane.l1ygridscale,
      ImagePlane.l2xgridscale, ImagePlane.l2ygridscale
  );
  cudaDeviceSynchronize();

  err=cudaMemcpy(far_lenses,d_far_lenses,nx2*ny2*sizeof(Double8),cudaMemcpyDeviceToHost);
  if(err != cudaSuccess){
    fprintf(stderr, "cudaMemcpy in setinterp d_far_lenses to far_lenses: %s\n", cudaGetErrorString(cudaGetLastError()));
    exit(EXIT_FAILURE);  
                        }

  cudaFree(d_Lens);
  cudaFree(d_hoc);
  cudaFree(d_ll);
  cudaFree(d_far_lenses);

  cudaDeviceReset();
  }
  



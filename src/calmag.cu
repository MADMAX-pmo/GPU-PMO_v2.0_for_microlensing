#include<stdio.h>
#include<stdlib.h>
#include<cuda_runtime.h>
#include<math.h>
#include "fits_IO.h"
#include"calmag.h"
#include"setinterp.h"
#include"global.h"

__device__ void detko(double a0,double a1,double a2,double a3,
                      double a4,double a5,double a6,double a7,
                      double *phi,double delx,double dely)
{
  phi[0] = 0.25*(a0+a2+a4+a6);
  phi[1] = 0.25*(a1+a3+a5+a7);
  phi[2] = 0.125*delx*(a0-a1-a2-a3-a4+a5+a6+a7);
  phi[3] = 0.125*dely*(a0+a1+a2-a3-a4-a5-a6+a7);
  phi[4] = 0.25*delx*delx*(-a1+a3-a5+a7);
  phi[5] = 0.25*dely*dely*(a0-a2+a4-a6);
  phi[6] = 0.375*delx*delx*delx*(-a0-a1+a2-a3+a4+a5-a6+a7);
  phi[7] = 0.375*dely*dely*dely*(a0-a1+a2+a3-a4+a5-a6-a7);  
}


__global__ void Map_Mag(
    int *d_num_map,int nlens,double mass_ave,Ls *d_Lens,
    Double8 *d_far_lenses,int *d_hoc,int *d_ll,
    double ktotal,double kstar,double gamma,
    double Ixmin,double Ixmax,double Iymin,double Iymax,
    int nx1, int ny1,int nx2,int ny2,
    long long unsigned int nx3, long long unsigned int ny3,
    double dx1,double dy1,double dx2,double dy2,
    double dx3,double dy3,
    double Sxmin,double Sxmax,double Symin,double Symax,
    int Snpix, double inv_Spixscale)
{
  int i=blockDim.x*blockIdx.x+threadIdx.x;
  int j=blockDim.y*blockIdx.y+threadIdx.y;

  if (i >= nx3 || j >= ny3) return;

  double posx = Ixmin+(i+0.5)*dx3;
  double posy = Iymin+(j+0.5)*dy3; 

  double alphax = 0.0, alphay = 0.0;
  double inv_mass = 1.0/mass_ave;  

  int cx1 = (int)((posx - Ixmin) / dx1);
  int cy1 = (int)((posy - Iymin) / dy1);
   
  for(int ii=cx1-1;ii<=cx1+1;ii++){
      if(ii>=0 && ii<nx1){ 
          for(int jj=cy1-1;jj<=cy1+1;jj++){
              if(jj>=0 && jj<ny1){
                  int it=d_hoc[ii*ny1+jj];
                  while(it!=-1){
                      double dx = posx-d_Lens[it].lensx;
                      double dy = posy-d_Lens[it].lensy;
                      double r2 = (d_Lens[it].lensm*inv_mass)/(dx*dx+dy*dy);  
                      alphax += dx*r2;
                      alphay += dy*r2;
                      it=d_ll[it];
                  }
              }
          }
      }
  }

  int cx2 = (int)((posx-Ixmin)/dx2);
  int cy2 = (int)((posy-Iymin)/dy2);
  int idx = cx2*ny2+cy2;     
  double phi[8];
   
  detko(
      d_far_lenses[idx].a0,d_far_lenses[idx].a1,
      d_far_lenses[idx].a2,d_far_lenses[idx].a3,
      d_far_lenses[idx].a4,d_far_lenses[idx].a5,
      d_far_lenses[idx].a6,d_far_lenses[idx].a7,
      phi,2.0/dx2,2.0/dy2);
       
  double xx = posx-d_far_lenses[idx].xc;
  double yy = posy-d_far_lenses[idx].yc;
  double xx2 = xx*xx;
  double yy2 = yy*yy;
  double z1 = (xx2-yy2)/2.0;
  double z2 = xx*yy;
  double z3 = 0.5*xx*(xx2/3.0-yy2);
  double z4 = 0.5*yy*(xx2-yy2/3.0);
   
  double alphax_far = phi[0]+xx*phi[2]+yy*phi[3]+z1*phi[4]+z2*phi[5]+z3*phi[6]+z4*phi[7];
  double alphay_far = phi[1]-yy*phi[2]+xx*phi[3]-z2*phi[4]+z1*phi[5]-z4*phi[6]+z3*phi[7];
   
  alphax += alphax_far+ktotal*posx+gamma*posx;
  alphay += alphay_far+ktotal*posy-gamma*posy-2.0*kstar*yy;
        
  double betax = posx-alphax;
  double betay = posy-alphay;

  if(betax>=Sxmin && betax<=Sxmax && betay>=Symin && betay<=Symax){
      int sx = (int)((betax-Sxmin)*inv_Spixscale);
      int sy = (int)((betay-Symin)*inv_Spixscale);
      if(sx>=0 && sx<Snpix && sy>=0 && sy<Snpix){
          atomicAdd(&d_num_map[sx*Snpix+sy],1);
      }
  }
} 



void calmag(
    int *num_map, int *hoc, int *ll, Double8 *far_lenses,
    int nlens,Ls *Lens,double mass_ave,
    struct Pa Parameters,struct Sp SourcePlane,struct Ip ImagePlane,int dev)
{
  cudaSetDevice(dev);
  cudaError_t err;

  // allocate device memory
  int Snpix = SourcePlane.npix;
  double inv_Spixscale = 1.0/SourcePlane.pixscale;

  int *d_num_map;
  cudaMalloc(&d_num_map,Snpix*Snpix*sizeof(int));

  int *d_hoc;
  cudaMalloc(&d_hoc,ImagePlane.l1nx*ImagePlane.l1ny*sizeof(int));

  int *d_ll;
  cudaMalloc(&d_ll,nlens*sizeof(int));

  Ls *d_Lens;
  cudaMalloc(&d_Lens,nlens*sizeof(Ls));

  Double8 *d_far_lenses;
  cudaMalloc(&d_far_lenses,ImagePlane.l2nx*ImagePlane.l2ny*sizeof(Double8));  
 
  // copy to device 
  err=cudaMemcpy(d_num_map,num_map, Snpix*Snpix*sizeof(int),cudaMemcpyHostToDevice);
  if(err != cudaSuccess) {
    fprintf(stderr, "cudaMemcpy in calmag num_map to d_num_map: %s\n", cudaGetErrorString(cudaGetLastError()));
    exit(EXIT_FAILURE);
                         }
  
  err=cudaMemcpy(d_hoc,hoc, ImagePlane.l1nx*ImagePlane.l1ny*sizeof(int),cudaMemcpyHostToDevice);
  if(err != cudaSuccess) {
    fprintf(stderr, "cudaMemcpy in calmag hoc to d_hoc: %s\n", cudaGetErrorString(cudaGetLastError()));
    exit(EXIT_FAILURE);
                         }

  err=cudaMemcpy(d_ll,ll,nlens*sizeof(int),cudaMemcpyHostToDevice);
  if(err != cudaSuccess) {
    fprintf(stderr, "cudaMemcpy in calmag ll to d_ll: %s\n", cudaGetErrorString(cudaGetLastError()));
    exit(EXIT_FAILURE);
                         }

  err=cudaMemcpy(d_Lens,Lens,nlens*sizeof(Ls),cudaMemcpyHostToDevice);
  if(err != cudaSuccess) {
     fprintf(stderr, "cudaMemcpy in calmag Lens to d_Lens: %s\n", cudaGetErrorString(cudaGetLastError()));
     exit(EXIT_FAILURE);
                         }

  err=cudaMemcpy(d_far_lenses,far_lenses,ImagePlane.l2nx*ImagePlane.l2ny*sizeof(Double8),cudaMemcpyHostToDevice);
  if(err != cudaSuccess) {
    fprintf(stderr, "cudaMemcpy in calmag far_lenses to d_far_lenses: %s\n", cudaGetErrorString(cudaGetLastError()));
    exit(EXIT_FAILURE);
                         }

  // launch
  dim3 block(2,32);
  dim3 grid((ImagePlane.l3nx+block.x-1)/block.x,
            (ImagePlane.l3ny+block.y-1)/block.y);

  Map_Mag<<<grid,block>>>(
      d_num_map,nlens,mass_ave,d_Lens,d_far_lenses,d_hoc,d_ll,
      Parameters.ktotal, Parameters.kstar, Parameters.gamma,
      ImagePlane.xmin, ImagePlane.xmax, ImagePlane.ymin, ImagePlane.ymax,
      ImagePlane.l1nx, ImagePlane.l1ny, ImagePlane.l2nx, ImagePlane.l2ny,
      ImagePlane.l3nx, ImagePlane.l3ny,
      ImagePlane.l1xgridscale, ImagePlane.l1ygridscale,
      ImagePlane.l2xgridscale, ImagePlane.l2ygridscale,
      ImagePlane.l3xgridscale, ImagePlane.l3ygridscale,
      SourcePlane.xmin, SourcePlane.xmax,
      SourcePlane.ymin, SourcePlane.ymax,
      Snpix, inv_Spixscale);
  
  cudaDeviceSynchronize();

  // copy back to host 
  err=cudaMemcpy(num_map,d_num_map,Snpix*Snpix*sizeof(int),cudaMemcpyDeviceToHost);
  if(err != cudaSuccess) {
    fprintf(stderr, "cudaMemcpy in calmag d_num_map to num_map: %s\n", cudaGetErrorString(cudaGetLastError()));
    exit(EXIT_FAILURE);
                         }
  cudaDeviceReset();
 
  cudaFree(d_num_map);
  cudaFree(d_Lens);
  cudaFree(d_far_lenses);
  cudaFree(d_ll);
  cudaFree(d_hoc);
}




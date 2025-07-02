#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>
#include <math.h>
#include <time.h>
#include "hocll.h"
#include "calmag.h"
//#include "sal_massave.h"
#include "lens_salpetermf.h"
#include "dm.h"
#include "fits_IO.h"
#include "setinterp.h"
#include "global.h"
#include "Parameters.h"


struct Pa init_parameters() {
    struct Pa p;
    p.ktotal = ktotal;
    p.fstar = fstar;
    p.gamma = g;
    p.kstar = p.fstar * p.ktotal;
    p.ksmooth = p.ktotal - p.kstar;
    return p;
}

struct Mf init_massfunction() {
    struct Mf m;
    m.power = salmf_power;
    m.minmass = salmf_minmass;
    m.maxmass = salmf_maxmass;
    return m;
}

struct Sp init_source_plane() {
    struct Sp s;
    s.scale = source_scale;
    s.xmin = -0.5 * s.scale;
    s.xmax = 0.5 * s.scale;
    s.ymin = s.xmin;
    s.ymax = s.xmax;
    s.npix = source_npix;
    s.pixscale = (s.xmax - s.xmin) / s.npix;
    return s;
}

struct Ip init_image_plane(struct Pa p, struct Sp s, double Pro_scale) {
    struct Ip ip;

    double Xpart = 1.0 / fabs(1 - p.ktotal - p.gamma);
    double Ypart = 1.0 / fabs(1 - p.ktotal + p.gamma);

    double xscale0 = (s.scale + 2.0 * Pro_scale) * Xpart;
    double yscale0 = (s.scale + 2.0 * Pro_scale) * Ypart;

    int n_star0 = (int)((xscale0 * yscale0 * p.kstar) / PI);
    double L0 = sqrt((xscale0 * yscale0) / n_star0);
    double min_Image = fmin(xscale0, yscale0);

    int l1nx = (2 * L0 <= yscale0 / 10) ?
               (int)(xscale0 / (2 * L0)) :
               (int)(xscale0 / (min_Image / 10));
    double l1xscale = xscale0 / l1nx;

    ip.l1nx = l1nx;
    ip.l1ny = (g == 0.0) ? l1nx : (int)(yscale0 / l1xscale) + 1;
    ip.l1xgridscale = l1xscale;
    ip.l1ygridscale = l1xscale;

    ip.l2nx = ip.l1nx * 20;
    ip.l2ny = ip.l1ny * 20;
    ip.l2xgridscale = xscale0 / ip.l2nx;
    ip.l2ygridscale = ip.l1ygridscale * ip.l1ny / ip.l2ny;

    ip.xscale = xscale0;
    ip.yscale = ip.l1ygridscale * ip.l1ny;
    ip.xmin = -0.5 * ip.xscale;
    ip.xmax = 0.5 * ip.xscale;
    ip.ymin = -0.5 * ip.yscale;
    ip.ymax = 0.5 * ip.yscale;

    long long Nrays0 = Nave * s.npix * s.npix * pow((s.scale + 2.0 * Pro_scale) / s.scale, 2);
    double XYRatio = ip.xscale / ip.yscale;
    double ny0 = sqrt(Nrays0 / XYRatio);
    ip.l3ny = ((int)ny0) + ((int)ny0 % 2 != 0);
    double nx0 = Nrays0 / ny0;
    ip.l3nx = ((int)nx0) + ((int)nx0 % 2 != 0);
    ip.l3xgridscale = ip.xscale / ip.l3nx;
    ip.l3ygridscale = ip.yscale / ip.l3ny;

    return ip;
}

Ls *init_lenses(int nlens, struct Ip ip, struct Mf m, double *real_massave) {
    Ls *Lens = (Ls *)malloc(nlens * sizeof(Ls));
    double cal_massave;
    lens_salpetermf(nlens,Lens,ip,m,&cal_massave);
    *real_massave = cal_massave;

    FILE *ft = fopen("Output/microlenses/lenspos.dat", "w");
    for (int i = 0; i < nlens; i++)
        fprintf(ft, "%e %e %e\n", Lens[i].lensx, Lens[i].lensy, Lens[i].lensm);
    fclose(ft);

    ft = fopen("Output/microlenses/lenspos.bin", "wb");
    fwrite(Lens, sizeof(Ls), nlens, ft);
    fclose(ft);

    return Lens;
}

/*
void print_lens_scales(double zl, double zs, double mass_ave) {
    const double dl = dm(zl);     // Mpc
    const double ds = dm(zs);     // Mpc
    const double dls = (ds * (1.0 + zs) - dl * (1.0 + zl)) / (1.0 + zs); // Mpc

    // θ_E in radians, then convert to arcsec
    double thetaE_rad = sqrt((4.0 * GRAVITY * M_SUN * mass_ave * dls * MPC) / (LIGHT_SPEED * LIGHT_SPEED * ds * MPC * dl * MPC));
    double thetaE_arcsec = thetaE_rad * (180.0 / PI) * 3600.0;
    
    // Σ_crit in kg/m²
    double sigma_crit = (LIGHT_SPEED * LIGHT_SPEED * ds * MPC) / (4.0 * PI * GRAVITY * dl * MPC * dls * MPC);
    
    printf("\n--- Einstein Radius & Critical Surface Density ---\n");
    printf("theta_E     = %.8e arcsec\n", thetaE_arcsec);
    printf("Σ_critical  = %.8e kg/m²\n", sigma_crit);
    printf("---------------------------------------------------\n\n");
}
*/


void compute_magnification(int *num_map, float **mag, long long *rays_count, 
                           double *Ave_rays, double *Ave_mag,
                           struct Sp s, struct Ip ip) {
    *rays_count = 0;
    double magsum = 0;

    double factor = ip.l3xgridscale * ip.l3ygridscale * (1.0 / (s.pixscale * s.pixscale));
    for(int i=0;i< s.npix;i++){
        for(int j=0;j<s.npix;j++){
            int idx = i * s.npix + j;
            (*rays_count) += num_map[idx];
            mag[i][j] = num_map[idx] * factor;
            (magsum) += mag[i][j];
        }
    }
    *Ave_rays = (double)*rays_count/(s.npix*s.npix);
    *Ave_mag = magsum/(s.npix*s.npix*1.0);

}


float **allocate_mag_map(int npix) {
    float **mag = (float **)malloc(npix * sizeof(float *));
    for (int i = 0; i < npix; i++)
        mag[i] = (float *)malloc(npix * sizeof(float));
    return mag;
}

void write_mag_output(float **mag, int npix) {
    FILE *ft = fopen("Output/Microlensing_MagMap.bin", "wb");
    for (int i = 0; i < npix; i++)
        fwrite(mag[i], sizeof(float), npix, ft);
    fclose(ft);

    int dim[2] = {npix, npix};
    write_fits_2D("Output/Microlensing_MagMap.fits", mag, dim);
}

void print_summary(struct Pa p, struct Sp s, struct Ip ip, double mu, int nlens, 
                   long long Nrays, long long rays_count, double Ave_rays, double Ave_mag, 
                   double real_massave, double mass_ave){
    printf("\n***** PARAMETERS *****\n");
    printf("ktotal = %lf\n", p.ktotal);
    printf("fstar = %lf\n", p.fstar);
    printf("kstar = %lf\n", p.kstar);
    printf("ksmooth = %lf\n", p.ksmooth);
    printf("gamma = %lf\n", p.gamma);
    printf("mu = %lf\n", mu);

    printf("mass_ave = %f,  real_massave = %f\n\n",mass_ave, real_massave);
     
    printf("ImagePlane:\n");
    printf("  xscale = %.6lf, yscale = %.6lf\n", ip.xscale, ip.yscale);
    printf("  l1: %d x %d\n", ip.l1nx, ip.l1ny);
    printf("    l1xgridscale = %.6lf\n", ip.l1xgridscale);
    printf("    l1ygridscale = %.6lf\n", ip.l1ygridscale);
    printf("  l2: %d x %d\n", ip.l2nx, ip.l2ny);
    printf("    l2xgridscale = %.6lf\n", ip.l2xgridscale);
    printf("    l2ygridscale = %.6lf\n", ip.l2ygridscale);
    printf("  l3: %ld x %ld\n", ip.l3nx, ip.l3ny);
    printf("    l3xgridscale = %.6lf\n", ip.l3xgridscale);
    printf("    l3ygridscale = %.6lf\n", ip.l3ygridscale);

    printf("  Nrays = %lld , Nave = %d\n", Nrays, Nave);
    printf("  nlens = %d\n", nlens);

    printf("\n\n");
    printf("*********************************** Results *********************************\n");
    printf("rays_count = %ld = %e\n",rays_count,rays_count*1.0);
    printf("Ave_rays = %lf\n",Ave_rays);
    printf("Ave_mag=%lf\n",Ave_mag);
    printf("*********************************** Results *********************************\n");
    printf("\n\n");

}



int main(int argc, char **argv) {
    clock_t total_start = clock();

    struct Pa p = init_parameters();
    struct Sp s = init_source_plane();
    double mu = 1 / ((1 - p.ktotal) * (1 - p.ktotal) - p.gamma * p.gamma);
    double r0 = sqrt(p.kstar), Pro_scale = 10.0 * r0;
    struct Ip ip = init_image_plane(p, s, Pro_scale);
 
    struct Mf m = init_massfunction();
    double mass_ave = massave(m);
    double real_massave;
    int nlens = (int)((ip.xscale * ip.yscale * p.kstar) / PI);
    Ls *Lens = init_lenses(nlens, ip, m, &real_massave);    

    clock_t t1 = clock();
    int *hoc = (int *)malloc(ip.l1nx * ip.l1ny * sizeof(int));
    int *ll = (int *)malloc(nlens * sizeof(int));
    hocll(nlens, Lens, ip, ll, hoc);
    printf("Level 1: %.2fs\n", (double)(clock() - t1) / CLOCKS_PER_SEC);

    clock_t t2 = clock();
    Double8 *far_lenses = (Double8 *)malloc(ip.l2nx * ip.l2ny * sizeof(Double8));
    setinterp(far_lenses, nlens, mass_ave, Lens, p, ip, hoc, ll, dev);
    printf("Level 2: %.2fs\n", (double)(clock() - t2) / CLOCKS_PER_SEC);

    int *num_map = (int *)malloc(s.npix * s.npix * sizeof(int));
    memset(num_map, 0, s.npix * s.npix * sizeof(int));
    float **mag = allocate_mag_map(s.npix);


    clock_t t3 = clock();
    calmag(num_map, hoc, ll, far_lenses, nlens, Lens, mass_ave, p, s, ip, dev);
    printf("Level 3: %.2fs\n", (double)(clock() - t3) / CLOCKS_PER_SEC);

    long long rays_count;
    double Ave_rays;
    double Ave_mag;
    compute_magnification(num_map, mag, &rays_count, &Ave_rays, &Ave_mag,s, ip);

    write_mag_output(mag, s.npix);

    long long Nrays = (long long)s.npix * s.npix * Nave *
                      pow((s.scale + 2.0 * Pro_scale) / s.scale, 2);
    print_summary(p, s, ip, mu, nlens, Nrays, rays_count,  Ave_rays, Ave_mag, real_massave, mass_ave);

    printf("Total time: %.2fs\n", (double)(clock() - total_start) / CLOCKS_PER_SEC);

    free(hoc); free(ll); free(far_lenses); free(num_map);
    for (int i = 0; i < s.npix; i++) free(mag[i]);
    free(mag);
    free(Lens);

}





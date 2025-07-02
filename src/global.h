#ifndef GLOBAL_H
#define GLOBAL_H

#define H0 7.0e+4
#define LIGHT_SPEED 2.9979e+8
#define OM 0.3
#define OA 0.7

#define GRAVITY 6.67408e-11
#define M_SUN 1.989e+30
#define MPC 3.08568e+22

#define PI 3.14159265358979323846264338327950288


typedef struct
{
  double lensx;
  double lensy;
  double lensm;
}Ls;

struct Mf
{
  double power;
  double minmass;
  double maxmass;
};

struct Pa
{
  double ktotal;
  double ksmooth;
  double kstar;
  double fstar;
  double gamma;
};

struct Sp
{
  double scale;
  double xmin;
  double xmax;
  double ymin;
  double ymax;
  int npix;
  double pixscale;
};

struct Ip
{
  double xscale;
  double yscale;
  double xmin;
  double xmax;
  double ymin;
  double ymax;
  int l1nx;
  int l1ny;
  double l1xgridscale;
  double l1ygridscale;
  int l2nx;
  int l2ny;
  double l2xgridscale;
  double l2ygridscale;
  long long unsigned int l3nx;
  long long unsigned int l3ny;
  double l3xgridscale;
  double l3ygridscale;
};

struct doubleeight
{
 double a0;
 double a1;
 double a2;
 double a3;
 double a4;
 double a5;
 double a6;
 double a7;
 double xc;
 double yc;
};
typedef struct doubleeight Double8;


#endif

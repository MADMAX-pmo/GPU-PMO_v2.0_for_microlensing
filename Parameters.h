#ifndef PARAMETERS_H
#define PARAMETERS_H

int dev=1;  // ID of the GPU device to use

// === Microlens Configuration ===
double salmf_power = -2.35;  // Parameters for Salpeter MF
double salmf_minmass = 0.1;
double salmf_maxmass = 10;

double ktotal   = 0.45; // Total convergence (κ)
double fstar    = 0.5;  // Fraction of κ contributed by microlenses
double g        = 0.45; // External shear

// === Source Plane Setup ===
double source_scale = 20;     // Side length of the source plane (in units of θ_E)
int source_npix     = 2000;  // Number of pixels per axis on the source plane

// === Ray-Tracing Accuracy ===
long int Nave = 1000;   // Target number of light rays per source pixel

#endif

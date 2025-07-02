============================================================================
                  GPU-PMO code,  version 2.0 (2025.6)

           **Authors:** Wenwen Zheng,  Xuechun Chen,  Guoliang Li
           **Language:** C + CUDA
============================================================================

**Description**
This is a high-performance microlensing simulation code based on **inverse ray-shooting**, featuring a **three-level grid scheme** to efficiently optimize microlens deflection calculations. The code is accelerated using **CUDA on NVIDIA GPUs**, achieving significant performance gains in computing magnification maps.

-- For implementation details, refer to:
Zheng, W., Chen, X., Li, G., & Chen, H.-Z. (2022). An improved GPU-based ray-shooting code for gravitational microlensing. The Astrophysical Journal, 931(2), 114
https://ui.adsabs.harvard.edu/abs/2022ApJ...931..114Z/abstract
-- For any help, please contact **wwzheng@pmo.ac.cn**

**Configuration**
-- Build options  → Edit `Makefile` based on your GPU (e.g. CUDA arch, paths)
-- For FITS file output, make sure [CFITSIO](https://heasarc.gsfc.nasa.gov/fitsio/) is installed and paths are correct.
-- Microlensing parameters  → Edit `Parameters.h`
-- Microlenses follow a Salpeter mass distribution by default.

**Output**
Results are saved in the `Output/` directory:
- `Microlensing_MagMap.bin` — raw binary format
- `Microlensing_MagMap.fits` — standard FITS format
- Maps are in units of Einstein radius (θ_E)

** Compile & run**
./run.bash


Don't worry! Be happy!


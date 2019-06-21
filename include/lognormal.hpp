#ifndef _LOGNORMAL_HPP_
#define _LOGNORMAL_HPP_

#include <vector>
#include <string>
#include <fftw3.h>
#include <gsl/gsl_spline.h>
#include "tpods.hpp"
#include "grid.hpp"
#include "galaxy.hpp"

class lognormal : public grid3D{
    double b, f;
    std::vector<double> kx, ky, kz;
    fftw_plan dr2dk, dk2dr;
    std::vector<double> F_i;
    
    pod2<size_t> getComplexIndex(int i, int j, int k);
    
    void fillGridWithPower(gsl_spline *Pk, gsl_interp_accel *acc);
    
    std::vector<double> fftFrequencies(int n, double l);
    
    void getLog();
    
    void getdk_ln();
    
    public:
        
        lognormal(pod3<int> N, pod3<double> L, std::string pkFile, double b = 1.0, double f = 0.0);
        
        ~lognormal();
        
        void init(std::string pkFile);
        
        void sample();
        
        std::vector<galaxy> getGalaxies(cosmology &cosmo, double nbar, pod3<double> r_min);
        
};

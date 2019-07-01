#ifndef _LOGNORMAL_HPP_
#define _LOGNORMAL_HPP_

#include <vector>
#include <string>
#include <random>
#include <fftw3.h>
#include <gsl/gsl_spline.h>
#include "tpods.hpp"
#include "grid.hpp"
#include "galaxy.hpp"

class lognormal : public grid3D{
    double b, f;
    std::vector<double> kx, ky, kz;
    fftw_plan dk2dr;
    std::vector<double> F_i;
    std::mt19937_64 gen;
    std::normal_distribution<double> norm;
    std::uniform_real_distribution<double> U;
    
    size_t getRealIndex(int i, int j, int k);
    
    pod2<size_t> getComplexIndex(int i, int j, int k);
    
    void fillGridWithPower(gsl_spline *Pk, gsl_interp_accel *acc);
    
    std::vector<double> fftFrequencies(int n, double l);
    
    void getLog();
    
    void getdk_ln(fftw_complex *dk);
    
    void updateStats(double val, double &mean, double &var, long &count);
    
    pod3<double> cartToEqua(double x, double y, double z, cosmology &cosmo);
    
    public:
        
        lognormal(pod3<int> N, pod3<double> L, std::string pkFile, double b = 1.0, double f = 0.0);
        
        ~lognormal();
        
        void init(std::string pkFile);
        
        void sample();
        
        std::vector<galaxy> getGalaxies(double nbar, pod3<double> r_min);
        
        std::vector<galaxy> getRandoms(double nbar, pod3<double> r_min, double timesRan);
        
        std::vector<galaxy> getGalaxies(cosmology &cosmo, gsl_spline *NofZ, gsl_interp_accel *acc,
                                        std::vector<int> &map, int nside, pod3<double> r_min,
                                        double z_min, double z_max);
        
        std::vector<galaxy> getRandoms(cosmology &cosmo, gsl_spline *NofZ, gsl_interp_accel *acc,
                                       std::vector<int> &map, int nside, pod3<double> r_min,
                                       double z_min, double z_max, double timeRan);
};

#endif

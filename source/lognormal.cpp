// Standard library includes
#include <vector>
#include <string>
#include <sstream>
#include <cmath>

// 3rd Party library includes
#include <omp.h>
#include <fftw3.h>
#include <gsl/gsl_spline.h>

// Custom library includes
#include "../include/tpods.hpp"
#include "../include/grid.hpp"
#include "../include/galaxy.hpp"
#include "../include/lognormal.hpp"

pod2<size_t> getComplexIndex(int i, int j, int k) {
    pod2<size_t> index;
    index.x = (2*k    ) + 2*(this->N.z/2 + 1)*(j + this->N.y*i);
    index.y = (2*k + 1) + 2*(this->N.z/2 + 1)*(j + this->N.y*i);
    return index;
}

void lognormal::fillGridWithPower(gsl_spline *Pk, gsl_interp_accel *acc) {
    double V = this->L.x*this->L.y*this->L.z;
    for (int i = 0; i < this->N.x; ++i) {
        double kxsq = this->kx[i]*this->kx[i];
        for (int j = 0; j < this->N.y; ++j) {
            double kysq = this->ky[j]*this->ky[j];
            for (int k = 0; k <= this->N.z/2; ++k) {
                double kzsq = this->kz[k]*this->kz[k];
                double k_mag = std::sqrt(kxsq + kysq + kzsq);
                
                pod2<size_t> index = lognormal(getComplexIndex(i, j, k));
                
                if (k_mag > 0) {
                    double mu = kx[i]/k_mag;
                    double mono = (this->b + mu*mu*this->f)*(this->b + mu*mu*this->f);
                    this->F_i[index.x] = mono*gsl_spline_eval(Pk, k_mag, acc);
                    this->F_i[index.y] = 0.0;
                } else {
                    this->F_i[index.x] = 0.0;
                    this->F_i[index.y] = 0.0;
                }
            }
        }
    }
}

std::vector<double> lognormal::fftFrequencies(int n, double l) {
    std::vector<double> k(n);
    double dk = 2.0*M_PI/l;
    for (int i = 0; i <= n/2; ++i)
        k[i] = i*dk;
    for (int i = n/2 + 1; i < n; ++i)
        k[i] = (i - n)*dk;
    return k;
}

void lognormal::getLog() {
#pragma omp parallel for
    for (int i = 0; i < this->F_i.size(); ++i)
        this->F_i[i] = std::log(1.0 + F_i[i]);
}

void lognormal::getdk_ln(fftw_complex *dk) {
    int N_tot = this->N.x*this->N.y*(this->N.z/2 + 1);
#pragma omp parallel for
    for (int i = 0; i < N_tot; ++i) {
        if (dk[i][0] > 0) {
            dk[i][0] /= N_tot;
            dk[i][1] = 0.0;
        } else {
            dk[i][0] = 0.0;
            dk[i][1] = 0.0;
        }
    }
}

void lognormal::updateStats(double val, double &mean, double &var, long &count) {
    count += 1L;
    double delta = val - mean;
    mean += delta/count;
    double delta2 = val - mean;
    double var += delta*delta2;
}

lognormal::lognormal(pod3<int> N, pod3<double> L, std::string pkFile, double b, double f) : 
gen((std::random_device())()), norm(0.0, 1.0), U(0.0, 1.0) {
    this->b = b;
    this->f = f;
    this->N = N;
    this->L = L;
    this->F.resize(N.x*N.y*2*(N.z/2 + 1));
    this->F_i.resize(N.x*N.y*2*(N.z/2 + 1));
    
    this->kx = lognormal::fftFrequencies(N.x, L.x);
    this->ky = lognormal::fftFrequencies(N.y, L.y);
    this->kz = lognormal::fftFrequencies(N.z, L.z);
    
    fftw_init_threads();
    fftw_import_wisdom_from_filename("fftwWisdom.dat");
    fftw_plan_with_nthreads(omp_get_max_threads());
    this->dr2dk = fftw_plan_dft_r2c_3d(N.x, N.y, N.z, this->F.data(), (fftw_complex *)this->F.data(),
                                       FFTW_MEASURE);
    this->dk2dr = fftw_plan_dft_c2r_3d(N.x, N.y, N.z, (fftw_complex *)this->F.data(), this->F.data(),
                                       FFTW_MEASURE);
    fftw_export_wisfrom_from_filename("fftwWisdom.dat"):
    
    lognormal::init(pkFile);
}

lognormal::~lognormal{
    fftw_destroy_plan(this->dr2dk);
    fftw_destroy_plan(this->dk2dr);
    fftw_cleanup_threads();
}

void lognormal::init(std::string pkFile) {
    if (std::ifstream(pkFile)) {
        std::vector<double> k, P;
        
        fftw_plan dr2dk_i = fftw_plan_dft_r2c_3d(this->N.x, this->N.y, this->N.z, this->F_i.data(),
                                                 (fftw_complex *)this->F_i.data(), FFTW_MEASURE);
        fftw_plan dk2dr_i = fftw_plan_dft_c2r_3d(this->N.x, this->N.y, this->N.z, 
                                                 (fftw_complex *)this->F_i.data(), this->F_i.data(),
                                                 FFTW_MEASURE);
        
        std::ifstream fin(pkFile);
        while (!fin.eof()) {
            double kt, Pt;
            fin >> kt >> Pt;
            if (!fin.eof()) {
                k.push_back(kt);
                P.push_back(Pt);
            }
        }
        fin.close();
        
        gsl_spline *Pk = gsl_spline_alloc(gsl_interp_cspline, P.size());
        gsl_interp_accel *acc = gsl_interp_accel_alloc();
        
        gsl_spline_init(Pk, k.data(), P.data(), P.size());
        
        lognormal::fillGridWithPower(Pk, acc);
        fftw_execute(dk2dr_i);
        lognormal::getLog();
        fftw_execute(dr2dk_i);
        lognormal::getdk_ln((fftw_complex *)this->F_i.data());
        
        gsl_spline_free(Pk);
        gsl_interp_accel_free(acc);
        fftw_destroy_plan(dr2dk_i);
        fftw_destroy_plan(dk2dr_i);
    } else {
        std::stringstream errMsg;
        errMsg << "Could not open " << pkFile << "\n";
        throw std::runtime_error(errMsg.str());
    }
}

void lognormal::sample() {
    for (int i = 0; i < this->N.x; ++i) {
        int i2 - (2*this->N.x - i) % this->N.x;
        for (int j = 0; j < this->N.y; ++j) {
            int j2 = (2*this->N.y - j) % this->N.y;
            for (int k = 0; k <= this->N.z/2; ++k) {
                pod2<size_t> index1 = getComplexIndex(i, j, k);
                
                if ((i == 0 or i == this->N.x/2) and (j == 0 or j == this->N.y/2) and
                    (k == 0 or k == this->N.z/2)) {
                    this->F[index1.x] = this->norm(this->gen)*std::sqrt(this->F_i[index1.x]);
                    this->F[index1.y] = 0.0;
                } else if (k == 0 or k == this->N.z/2) {
                    pod2<size_t> index2 = getComplexIndex(i2, j2, k);
                    this->F[index1.x] = this->norm(this->gen)*std::sqrt(this->F_i[index1.x]/2);
                    this->F[index1.y] = this->norm(this->gen)*std::sqrt(this->F_i[index1.x]/2);
                    
                    this->F[index2.x] =  this->F[index1.x];
                    this->F[index2.y] = -this->F[index1.y];
                } else {
                    this->F[index1.x] = this->norm(this->gen)*std::sqrt(this->F_i[index1.x]/2);
                    this->F[index1.y] = this->norm(this->gen)*std::sqrt(this->F_i[index1.x]/2);
                }
            }
        }
    }
    
    fftw_execute(dk2dr);
}

std::vector<galaxy> lognormal::getGalaxies(double nbar, pod3<double> r_min) {
    std::vector<galaxy> gals;
    double n = nbar*this->Delta_r.x*this->Delta_r.y*this->Delta_r.z;
    long count = 0L;
    double mean = 0.0, var = 0.0;
    for (int i = 0; i < this->N.x; ++i) {
        for (int j = 0; j < this->N.y; ++j) {
            for (int k = 0; k < this->N.z; ++k) {
                size_t index = getRealIndex(i, j, k);
                
                lognormal::updateStats(this->F[index], mean, var, count);
            }
        }
    }
    var /= (count - 1L);
    
    for (int i = 0; i < this->N.x; ++i) {
        double r_x = i*this->Delta_r.x + r_min.x;
        for (int j = 0; j < this->N.y; ++j) {
            double r_y = j*this->Delta_r.y + r_min.y;
            for (int k = 0; k < this->N.z; ++k) {
                double r_z = k*this->Delta_r.z + r_min.z;
                size_t index = getRealIndex(i, j, k);
                this->F[index] -= mean;
                
                double density = n*exp(this->F[index] - var/2.0);
                std::poisson_distribution<int> p_dist(density);
                int num_gals = p_dist(this->gen);
                
                
                for (int gal = 0; gal < num_gals; ++gal) {
                    galaxy gal(r_x + this->U(this->gen)*this->Delta_r.x, 
                               r_y + this->U(this->gen)*this->Delta_r.y,
                               r_z + this->U(this->gen)*this->Delta_r.z,
                               nbar, 1.0, this->b);
                    gals.push_back(gal);
                }
            }
        }
    }
    
    return gals;
}

std::vector<galaxy> lognormal::getGalaxies(cosmology &cosmo, gsl_spline *NofZ, gsl_interp_accel *acc,
                                           std::vector<int> &map, int nside, pod3<double> r_min,
                                           double z_min, double z_max) {
    std::vector<galaxy> gals;
    long count = 0L;
    double mean = 0.0, var = 0.0;
    for (int i = 0; i < this->N.x; ++i) {
        for (int j = 0; j < this->N.y; ++j) {
            for (int k = 0; k < this->N.z; ++k) {
                size_t index = getRealIndex(i, j, k);
                
                lognormal::updateStats(this->F[index], mean, var, count);
            }
        }
    }
    var /= (count - 1L);
    
    for (int i = 0; i < this->N.x; ++i) {
        double r_x = i*this->Delta_r.x + r_min.x;
        for (int j = 0; j < this->N.y; ++j) {
            double r_y = j*this->Delta_r.y + r_min.y;
            for (int k = 0; k < this->N.z; ++k) {
                double r_z = k*this->Delta_r.z + r_min.z;
                size_t index = getRealIndex(i, j, k);
                this->F[index] -= mean;
                
                double r = std::sqrt(r_x*r_x + r_y*r_y + r_z*r_z);
                double z = get_redshift_from_comoving_distance(r);
                double n = gsl_spline_eval(NofZ, z, acc)*this->Delta_r.x*this->Delta_r.y*this->Delta_r.z;
                double dec = std::asin(r_z/r);
                
            }
        }
    }

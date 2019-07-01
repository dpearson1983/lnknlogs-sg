#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <gsl/gsl_spline.h>
#include <chealpix.h>
#include "../include/cosmology.hpp"
#include "../include/file_io.hpp"
#include "../include/galaxy.hpp"
#include "../include/grid.hpp"
#include "../include/harppi.hpp"
#include "../include/lognormal.hpp"
#include "../include/tpods.hpp"

void getNofZSpline(std::string file, std::vector<double> &n, std::vector<double> &z) {
    if (std::ifstream(file)) {
        std::ifstream fin(file);
        while (!fin.eof()) {
            double zt, nt;
            fin >> zt >> nt;
            if (!fin.eof()) {
                z.push_back(zt);
                n.push_back(nt);
            }
        }
        fin.close();
    } else {
        std::stringstream errMsg;
        errMsg << "Could not open " << file << "\n";
        throw std::runtime_error(errMsg.str());
    }
}

std::vector<int> getMap(std::string file, int nside) {
    long n_pix = nside2npix(nside);
    std::vector<int> map(n_pix);
    if (std::ifstream(file)) {
        std::ifstream fin(file);
        for (size_t i = 0; i < map.size(); ++i) {
            int pix;
            double ra, dec;
            fin >> pix >> ra >> dec >> map[i];
        }
        fin.close();
    } else {
        std::stringstream errMsg;
        errMsg << "Could not open " << file << "\n";
        throw std::runtime_error(errMsg.str());
    }
    return map;
}

int main(int argc, char *argv[]) {
    parameters p(argv[1]);
    p.print();
    
    if (p.getb("surveyMode")) {
        std::cout << "Get grid properties..." << std::endl;
        pod3<int> N = {p.geti("Nx"), p.geti("Ny"), p.geti("Nz")};
        pod3<double> L = {p.getd("Lx"), p.getd("Ly"), p.getd("Lz")};
        pod3<double> r_min = {p.getd("x_min"), p.getd("y_min"), p.getd("z_min")};
        
        std::cout << "Read in the map..." << std::endl;
        std::vector<int> map = getMap(p.gets("mapFile"), p.geti("nside"));
        std::cout << "    Number of pixels in map: " << map.size() << std::endl;
        
        std::cout << "Initialize the n(z) spline..." << std::endl;
        std::vector<double> n, z;
        getNofZSpline(p.gets("NofZFile"), n, z);
        gsl_spline *NofZ = gsl_spline_alloc(gsl_interp_cspline, n.size());
        gsl_interp_accel *acc = gsl_interp_accel_alloc();
        gsl_spline_init(NofZ, z.data(), n.data(), n.size());
        std::cout << "    Spline test: " << gsl_spline_eval(NofZ, 0.5, acc) << std::endl;
        
        std::cout << "Initialize cosmology..." << std::endl;
        cosmology cosmo(p.getd("H_0"), p.getd("Omega_M"), p.getd("Omega_L"));
        std::cout << "    Test redshift: " << cosmo.get_redshift_from_comoving_distance(1875.0) << std::endl;
        
        std::cout << "Initialize file writer..." << std::endl;
        fileIO writer(p.gets("outBase"), p.gets("outExt"), p.geti("digits"), p.geti("startNum"));
        
        std::cout << "Initialize lognormal object..." << std::endl;
        lognormal lnknlog(N, L, p.gets("pkFile"), p.getd("b"), p.getd("f"));
        
        if (p.getb("generateRandoms")) {
            std::cout << "Generating randoms..." << std::endl;
            std::vector<galaxy> rans = lnknlog.getRandoms(cosmo, NofZ, acc, map, p.geti("nside"), r_min, 
                                                          p.getd("red_min"), p.getd("red_max"), 
                                                          p.getd("timesRan"));
            std::cout << "    Total number of randoms: " << rans.size() << std::endl;
            writer.write(rans, p.gets("ransFile"));
        }
        
        std::cout << "Generating mocks..." << std::endl;
        for (int mock = p.geti("startNum"); mock < p.geti("startNum") + p.geti("numMocks"); ++mock) {
            std::cout << "    Mock #" << mock << "..." << std::endl;
            lnknlog.sample();
            
            std::vector<galaxy> gals = lnknlog.getGalaxies(cosmo, NofZ, acc, map, p.geti("nside"), r_min,
                                                           p.getd("red_min"), p.getd("red_max"));
            std::cout << "        Number of galaxies: " << gals.size() << std::endl;
            writer.write(gals);
        }
        
        gsl_spline_free(NofZ);
        gsl_interp_accel_free(acc);
    } else {
        std::cout << "Get grid properties..." << std::endl;
        pod3<int> N = {p.geti("Nx"), p.geti("Ny"), p.geti("Nz")};
        pod3<double> L = {p.getd("Lx"), p.getd("Ly"), p.getd("Lz")};
        pod3<double> r_min = {p.getd("x_min"), p.getd("y_min"), p.getd("z_min")};
        
        std::cout << "Initialize file writer..." << std::endl;
        fileIO writer(p.gets("outBase"), p.gets("outExt"), p.geti("digits"), p.geti("startNum"));
        
        std::cout << "Initialize cosmology..." << std::endl;
        cosmology cosmo(p.getd("H_0"), p.getd("Omega_M"), p.getd("Omega_L"));
        std::cout << "    Test redshift: " << cosmo.get_redshift_from_comoving_distance(1875.0) << std::endl;
        
        std::cout << "Initialize lognormal object..." << std::endl;
        lognormal lnknlog(N, L, p.gets("pkFile"), p.getd("b"), p.getd("f"));
        
        if (p.getb("generateRandoms")) {
            std::cout << "Generating randoms..." << std::endl;
            std::vector<galaxy> rans = lnknlog.getRandoms(cosmo, p.getd("nbar"), r_min, p.getd("timesRan"));
            std::cout << "    Total number of randoms: " << rans.size() << std::endl;
            writer.write(rans, p.gets("ransFile"));
        }
        
        std::cout << "Generating mocks..." << std::endl;
        for (int mock = p.geti("startNum"); mock < p.geti("startNum") + p.geti("numMocks"); ++mock) {
            std::cout << "    Mock #" << mock << "..." << std::endl;
            lnknlog.sample();
            
            std::vector<galaxy> gals = lnknlog.getGalaxies(cosmo, p.getd("nbar"), r_min);
            std::cout << "        Number of galaxies: " << gals.size() << std::endl;
            writer.write(gals);
        }
    }
    
    return 0;
}

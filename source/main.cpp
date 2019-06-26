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

void getNofZSpline(std::string file, gsl_spline *NofZ) {
    std::vector<double> n, z;
    
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
    
    NofZ = gsl_spline_alloc(gsl_interp_cspline, n.size());
    gsl_spline_init(NofZ, z.data(), n.data(), n.size());
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
    
    std::cout << "Get grid properties..." << std::endl;
    pod3<int> N = {p.geti("Nx"), p.geti("Ny"), p.geti("Nz")};
    pod3<double> L = {p.getd("Lx"), p.getd("Ly"), p.getd("Lz")};
    pod3<double> r_min = {p.getd("x_min"), p.getd("y_min"), p.getd("z_min")};
    
    std::cout << "Read in the map..." << std::endl;
    std::vector<int> map = getMap(p.gets("mapFile"), p.geti("nside"));
    
    std::cout << "Initialize the n(z) spline..." << std::endl;
    gsl_spline *NofZ;
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    getNofZSpline(p.gets("NofZFile"), NofZ);
    
    std::cout << "Initialize cosmology..." << std::endl;
    cosmology cosmo(p.getd("H_0"), p.getd("Omega_M"), p.getd("Omega_L"));
    
    std::cout << "Initialize file writer..." << std::endl;
    fileIO writer(p.gets("outBase"), p.gets("outExt"), p.geti("digits"), p.geti("startNum"));
    
    std::cout << "Initialize lognormal object..." << std::endl;
    lognormal lnknlog(N, L, p.gets("pkFile"), p.getd("b"), p.getd("f"));
    
    for (int mock = p.geti("startNum"); mock < p.geti("startNum") + p.geti("numMocks"); ++mock) {
        lnknlog.sample();
        
        std::vector<galaxy> gals = lnknlog.getGalaxies(cosmo, NofZ, acc, map, p.geti("nside"), r_min,
                                                       p.getd("red_min"), p.getd("red_max"));
        
        writer.write(gals);
    }
    
    gsl_spline_free(NofZ);
    gsl_interp_accel_free(acc);
    
    return 0;
}

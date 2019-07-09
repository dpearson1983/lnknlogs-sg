#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <vector>
#include <limits>
#include "../include/galaxy.hpp"
#include "../include/file_io.hpp"

std::string fileIO::filename() {
    std::stringstream file;
    file << this->base << "_" << std::setw(this->digits) << std::setfill('0') << this->N << "." << this->ext;
    N++;
    return file.str();
}

fileIO::fileIO(std::string base, std::string ext, int digits, int startNum) {
    this->base = base;
    this->ext = ext;
    this->digits = digits;
    this->N = startNum;
}

void fileIO::writeTXT(std::vector<galaxy> &gals) {
    std::string file = fileIO::filename();
    std::ofstream fout(file);
    fout.precision(std::numeric_limits<double>::digits10);
    for (size_t i = 0; i < gals.size(); ++i) {
        fout << gals[i].ra << " " << gals[i].dec << " " << gals[i].red << " " << gals[i].x << " ";
        fout << gals[i].y << " " << gals[i].z << " " << gals[i].nbar << " " << gals[i].w << " ";
        fout << gals[i].bias << "\n";
    }
    fout.close();
}

void fileIO::writeTXT(std::vector<galaxy> &gals, std::string file) {
    std::ofstream fout(file);
    fout.precision(std::numeric_limits<double>::digits10);
    for (size_t i = 0; i < gals.size(); ++i) {
        fout << gals[i].ra << " " << gals[i].dec << " " << gals[i].red << " " << gals[i].x << " ";
        fout << gals[i].y << " " << gals[i].z << " " << gals[i].nbar << " " << gals[i].w << " ";
        fout << gals[i].bias << "\n";
    }
    fout.close();
}

void fileIO::writeBIN(std::vector<galaxy> &gals) {
    long N_gals = gals.size();
    std::string file = fileIO::filename();
    std::ofstream fout(file, std::ios::out|std::ios::binary);
    fout.write((char *)&N_gals, sizeof(long));
    fout.write((char *)gals.data(), gals.size()*sizeof(galaxy));
    fout.close();
}

void fileIO::writeBIN(std::vector<galaxy> &gals, std::string file) {
    long N_gals = gals.size();
    std::ofstream fout(file, std::ios::out|std::ios::binary);
    fout.write((char *)&N_gals, sizeof(long));
    fout.write((char *)gals.data(), gals.size()*sizeof(galaxy));
    fout.close();
}

std::vector<galaxy> fileIO::readTXT() {
    std::vector<galaxy> gals;
    std::string file = fileIO::filename();
    if (std::ifstream(file)) {
        std::ifstream fin(file);
        while (!fin.eof()) {
            galaxy gal;
            fin >> gal.ra >> gal.dec >> gal.red >> gal.x >> gal.y >> gal.z >> gal.nbar >> gal.w >> gal.bias;
            if (!fin.eof()) {
                gals.push_back(gal);
            }
        }
        fin.close();
    } else {
        std::stringstream errMsg;
        errMsg << "Could not open " << file << "\n";
        throw std::runtime_error(errMsg.str());
    }
    return gals;
}

std::vector<galaxy> fileIO::readTXT(std::string file) {
    std::vector<galaxy> gals;
    if (std::ifstream(file)) {
        std::ifstream fin(file);
        while (!fin.eof()) {
            galaxy gal;
            fin >> gal.ra >> gal.dec >> gal.red >> gal.x >> gal.y >> gal.z >> gal.nbar >> gal.w >> gal.bias;
            if (!fin.eof()) {
                gals.push_back(gal);
            }
        }
        fin.close();
    } else {
        std::stringstream errMsg;
        errMsg << "Could not open " << file << "\n";
        throw std::runtime_error(errMsg.str());
    }
    return gals;
}

std::vector<galaxy> fileIO::readBIN() {
    std::string file = fileIO::filename();
    long N_gals;
    std::ifstream fin(file, std::ios::in|std::ios::binary);
    fin.read((char *)&N_gals, sizeof(long));
    std::vector<galaxy> gals(N_gals);
    fin.read((char *)gals.data(), N_gals*sizeof(galaxy));
    fin.close();
    return gals;
}

std::vector<galaxy> fileIO::readBIN(std::string file) {
    long N_gals;
    std::ifstream fin(file, std::ios::in|std::ios::binary);
    fin.read((char *)&N_gals, sizeof(long));
    std::vector<galaxy> gals(N_gals);
    fin.read((char *)gals.data(), N_gals*sizeof(galaxy));
    fin.close();
    return gals;
}

void fileIO::write(std::vector<galaxy> &gals) {
    if (this->ext == "bin") {
        fileIO::writeBIN(gals);
    } else {
        fileIO::writeTXT(gals);
    }
}

void fileIO::write(std::vector<galaxy> &gals, std::string file) {
    if (this->ext == "bin") {
        fileIO::writeBIN(gals, file);
    } else {
        fileIO::writeTXT(gals, file);
    }
}

std::vector<galaxy> fileIO::read() {
    std::vector<galaxy> gals;
    if (this->ext == "bin") {
        gals = fileIO::readBIN();
    } else {
        gals = fileIO::readTXT();
    }
    return gals;
}

std::vector<galaxy> fileIO::read(std::string file) {
    std::vector<galaxy> gals;
    if (this->ext == "bin") {
        gals = fileIO::readBIN(file);
    } else {
        gals = fileIO::readTXT(file);
    }
    return gals;
}

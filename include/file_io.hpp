#ifndef _FILE_IO_HPP_
#define _FILE_IO_HPP_

#include <string>
#include <vector>
#include "galaxy.hpp"

class fileIO{
    int N, digits;
    std::string base, ext;
    
    std::string filename();
    
    public:
        fileIO(std::string base, std::string ext, int digits, int startNum = 1);
        
        void write(std::vector<galaxy> &gals);
        
        void write(std::vector<galaxy> &gals, std::string file);
        
        std::vector<galaxy> read();
        
};

#endif

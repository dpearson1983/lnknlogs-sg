#ifndef _FILE_IO_HPP_
#define _FILE_IO_HPP_

#include <string>
#include <vector>
#include "galaxy.hpp"

class fileIO{
    int N, digits;
    std::string base, ext;
    
    public:
        fileIO(std::string base, std::string ext, int digits);
        
        void write(std::vector<galaxy> &gals);
        
        std::vector<galaxy> read();
        
};

#endif

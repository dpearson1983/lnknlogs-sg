#ifndef _FILE_IO_HPP_
#define _FILE_IO_HPP_

#include <string>
#include <vector>
#include "galaxy.hpp"

class fileIO{
    int N, digits;
    std::string base, ext;
    
    std::string filename();
    
    void writeTXT(std::vector<galaxy> &gals);
    
    void writeTXT(std::vector<galaxy> &gals, std::string file);
    
    void writeBIN(std::vector<galaxy> &gals);
    
    void writeBIN(std::vector<galaxy> &gals, std::string file);
    
    std::vector<galaxy> readTXT();
    
    std::vector<galaxy> readTXT(std::string file);
    
    std::vector<galaxy> readBIN();
    
    std::vector<galaxy> readBIN(std::string file);
    
    public:
        fileIO(std::string base, std::string ext, int digits, int startNum = 1);
        
        void write(std::vector<galaxy> &gals);
        
        void write(std::vector<galaxy> &gals, std::string file);
        
        std::vector<galaxy> read();
        
        std::vector<galaxy> read(std::string file);
        
};

#endif

#ifndef _GRID_HPP_
#define _GRID_HPP_

#include <vector>
#include "tpods.hpp"

class grid3D{
    protected:
        pod3<double> L, Delta_r;
        pod3<int> N;
    
        std::vector<double> F;
        
        size_t getIndex(int i, int j, int k);
        
    public:
        grid3D(pod3<int> N, pod3<double> L);
        
        void set(int i, int j, int k, double val);
        
        double get(int i, int j, int k);
        
};

#endif
    

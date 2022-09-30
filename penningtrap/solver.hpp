//
//  solver.hpp
//  penningtrap
//
//  Created by Felix Aarekol Forseth on 30/09/2022.
//

#ifndef solver_hpp
#define solver_hpp

#include <stdio.h>
#include <iostream>
#include <vector>

class Solver{
public:
    double time;           // Length of simulation.
    double N;              // Steps of simulation.
    std::string filename;  // Name of files.
    
    // Declaration function.
    Solver(double time, double N, std::string filename);
    
    // Solve with forward-Euler.
    // Needs to be solved for each particle.
    void forwardEuler(vector<Particle> particles);
    
    // Function to print to datafiles. 
    
};

#endif /* solver_hpp */

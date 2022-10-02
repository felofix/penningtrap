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
#include "particle.hpp"
#include "penningTrap.hpp"

class Solver{
public:
    double time;           // Length of simulation.
    double N;              // Steps of simulation.
    std::string filename;  // Name of files.
    double dt;    // timestep.
    
    // Declaration function.
    Solver(double time, double N, std::string filename);
    
    // Solve with forward-Euler.
    // Needs to be solved for each particle.
    void SolveforwardEuler(PenningTrap pt);
    
    // Stepping forward like Euler would.
    arma::mat forwardEulerStep(PenningTrap pt, Particle p);
    
    // Function to print to datafiles. 
    void writetofilefloat(arma::vec timer, arma::mat xyz, std::string direc);
    
};

#endif /* solver_hpp */

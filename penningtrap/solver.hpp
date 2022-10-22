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
#include <complex> 

class Solver{
public:
    double time;           // Length of simulation.
    double N;              // Steps of simulation.
    std::string filename;  // Name of files.
    double dt;    // timestep.
    bool etimedep; // time dependant electric field.
    double f; // amplitude
    double omegaV; // angular frequency.
    bool writetofile; // If you want to write to file.
    bool errorplot; //
    
    // Declaration function.
    Solver(double time, double N, std::string filename, bool etimedep, double f, double omegaV, bool wrfile, bool errorplot);
    
    // Solve with forward-Euler.
    
    // Needs to be solved for each particle.
    void SolveforwardEuler(PenningTrap pt);
    
    // Rung-kutta solver. 
    void RungeKuttaW(PenningTrap &pt);
    
    // Runge-kutta step.
    arma::mat forwardRKStepW(Particle p, double ddt, arma::vec vel, arma::vec pos, arma::vec a);
    
    void updateposition(PenningTrap &pt, arma::mat velpos);
    
    // Stepping forward like Euler would.
    arma::mat forwardEulerStep(PenningTrap pt, Particle p);
    
    // Function to print to datafiles. 
    void writetofilefloat(arma::vec timer, arma::mat xyz, std::string direc);
    
    // Analytical soulution.
    void solve_analytically(PenningTrap pt, arma::mat numeric, std::string euorru);
    
    //
    arma::mat create_analytical(PenningTrap pt);
};

#endif /* solver_hpp */

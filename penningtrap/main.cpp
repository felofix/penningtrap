//
//  main.cpp
//  penningtrap
//
//  Created by Felix Aarekol Forseth on 30/09/2022.
//

#include <iostream>
#include "particle.hpp"
#include "penningTrap.hpp"
#include "solver.hpp"

int main(int argc, const char * argv[]) {
    double c = 5.0;
    double m = 4.0; 
    arma::vec posi1 = {0, 0, 1};
    arma::vec posi2 = {1, 0, 0};
    arma::vec vel = {0, 0, 0};
    Particle calcium1 = Particle(c, m, posi1, vel);
    Particle calcium2 = Particle(c, m, posi2, vel);
    PenningTrap pp = PenningTrap(4.0, 4.0, 4.0);
    pp.add_particle(calcium1);
    pp.add_particle(calcium2);
    Solver solve = Solver(1, 10, "felix.txt");
    
}


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
    arma::vec posi1 = {0, 0, 5};
    arma::vec vel = {0, 0, 0};
    Particle calcium1 = Particle(c, m, posi1, vel);
    PenningTrap pp = PenningTrap(4.0, 4.0, 4.0);
    pp.add_particle(calcium1);
    Solver solve = Solver(10, 1000, "felix.txt");
    std::cout << pp.total_acceleration(calcium1) << std::endl;
    arma::mat velpos = solve.forwardEulerStep(pp, calcium1);
    std::cout << velpos << std::endl;
}


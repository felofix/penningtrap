//
//  particle.hpp
//  penningtrap
//
//  Created by Felix Aarekol Forseth on 30/09/2022.
//

#ifndef particle_hpp
#define particle_hpp

#include <stdio.h>
#include <armadillo>

class Particle{
public:
    double charge;    // Charge of the particle.
    double mass;      // Mass of particle.
    arma::vec position;   // Poisition of particle.
    arma::vec velocity;   // Velocity of particle.
    
    // Declaration function.
    Particle(double c, double m, arma::vec pos, arma::vec vel);
    
    // A nice function to print out all the values.
    void printToScreen(); 
    
};

#endif /* penningTrap_hpp */

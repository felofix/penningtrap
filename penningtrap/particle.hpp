//
//  particle.hpp
//  penningtrap
//
//  Created by Felix Aarekol Forseth on 30/09/2022.
//

#define particle_hpp

#include <stdio.h>
#include <armadillo>

class Particle{
public:
    double charge;    // Charge of the particle.
    double mass;      // Mass of particle.
    arma::vec position;   // Poisition of particle.
    arma::vec velocity;   // Velocity of particle.
    arma::vec force; // Force experienced by the particle.
    arma::vec acceleration; // Acelleration experienced by the particle.
    
    // Declaration function.
    Particle(double c, double m, arma::vec pos, arma::vec vel, arma::vec forc, arma::vec accel);
    
    // A nice function to print out all the values.
    void printToScreen(); 
    
};


//
//  particle.cpp
//  penningtrap
//
//  Created by Felix Aarekol Forseth on 30/09/2022.
//

#include "particle.hpp"

// Declaration function.
Particle::Particle(double c, double m, arma::vec pos, arma::vec vel, arma::vec forc, arma::vec accel){
    charge = c;   // Charge of the particle,
    mass = m;     // Mass of particle.
    position = pos;  // Position of particle.
    velocity = vel;  // Velocity of particle. 
    force = forc; // Force experienced by the particle.
    acceleration  = accel;// Acelleration experienced by the particle.
}
    // A nice function to print out all the values.
void Particle::printToScreen(){
    std::cout << "Charge: " << charge << "." << std::endl;
    std::cout << "Mass: " << mass << "." << std::endl;
    std::cout << "Position: " << std::endl;
    std::cout << position;
    std::cout << "Velocity: " << std::endl;
    std::cout << velocity;
    std::cout << "Force: " << std::endl;
    std::cout << force;
    std::cout << "Acceleration: " << std::endl;
    std::cout << acceleration;
}


//
//  particle.cpp
//  penningtrap
//
//  Created by Felix Aarekol Forseth on 30/09/2022.
//

#include "particle.hpp"

// Declaration function.
Particle::Particle(double c, double m, arma::vec pos, arma::vec vel){
    charge = c;   // Charge of the particle,
    mass = m;     // Mass of particle.
    position = pos;  // Position of particle.
    velocity = vel;  // Velocity of particle. 
}
    // A nice function to print out all the values.
void Particle::printToScreen(){
    std::cout << "Charge: " << charge << "." << std::endl;
    std::cout << "Mass: " << mass << "." << std::endl;
    std::cout << "Position: " << std::endl;
    std::cout << position;
    std::cout << "Velocity: " << std::endl;
    std::cout << velocity;
}


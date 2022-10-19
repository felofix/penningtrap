//
//  penningTrap.hpp
//  penningtrap
//
//  Created by Felix Aarekol Forseth on 30/09/2022.
//

#ifndef penningTrap_hpp
#define penningTrap_hpp

#include <stdio.h>
#include <armadillo>
#include <vector>
#include "particle.hpp"    

class PenningTrap{
public:
    double B0;     // Applied magnetic field.
    double V0;     // Applied electric field.
    double d;      // characteristic dimension.
    std::vector<Particle> particles; // Particles in penning trap.
    double gamma;  // Difference.
    double ke = 1.38935333e5;
    bool coulumb;
    
    // Declaration function.
    PenningTrap(double b, double v, double d, bool coul);
    
    // A nice function to print out all the values.
    void printToScreen();
    
    // Calculate the external electric field at certain position, where a specific potential is used(outlined in the README file).
    arma::vec find_efield(Particle P);
    
    // Calculate the external magnetic field, where a specific field is used(outlined in the README file).
    arma::vec find_mfield(Particle P);
    
    // Adding particle.
    void add_particle(Particle add);
    
    // Calculate force due to intneractions between particles.
    arma::vec force_from_particles(Particle pcheck);
    
    // Total force of particle.
    arma::vec total_force(Particle pcheck);
    
    // Acceleration of particle.
    arma::vec total_acceleration(Particle pcheck);
    
    // counting number of particles in trap.
    int count_particles();
    
    // filling the penningtrap with random particles.
    void fill_with_particles(int number, double charge, double mass);
};

#endif /* penningTrap_hpp */

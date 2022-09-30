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
#include "particle.hpp"     // Spørre hvordan jeg kan inkludere en folder med header filer.

class PenningTrap{
public:
    double B0;     // Applied magnetic field.
    double V0;     // Applied electric field.
    double d;      // characteristic dimension.
    std::vector<Particle> particles; // Particles in penning trap.
    double gamma;  // Difference. 
    
    // Declaration function.
    PenningTrap(double b, double v, double d);
    
    // A nice function to print out all the values.
    void printToScreen();
    
    // Calculate the external electric field at certain position, where a specific potential is used(outlined in the README file).
    arma::vec find_efield(arma::vec position);
    
    // Calculate the external magnetic field, where a specific field is used(outlined in the README file).
    arma::vec find_mfield(arma::vec position);
    
    // Calculate force due to intneractions between particles.
    
};

#endif /* penningTrap_hpp */

//
//  penningTrap.cpp
//  penningtrap
//
//  Created by Felix Aarekol Forseth on 30/09/2022.
//

#include "penningTrap.hpp"

PenningTrap::PenningTrap(double b, double v, double d2){
    B0 = b;
    V0 = v;
    d = d2;
    gamma = V0/pow(d, 2);
}

void PenningTrap::printToScreen(){
    std::cout << "B0: " << B0 << "." << std::endl;
    std::cout << "V0: " << V0 << "." << std::endl;
    std::cout << "d: " << d << "." << std::endl;
}

arma::vec PenningTrap::find_efield(Particle P){
    arma::vec e_field = {0, 0, 0};
    e_field(0) = P.position(0)*gamma;  // x
    e_field(1) = P.position(1)*gamma;  // y
    e_field(2) = P.position(2)*gamma*(-2); // z
    return e_field;
}

arma::vec PenningTrap::find_mfield(Particle P){
    arma::vec m_field = {0, 0, 0};
    m_field(2) = P.position(2)*B0; // z
    return m_field;
}

void PenningTrap::add_particle(Particle add){
    particles.push_back(add);
}

arma::vec PenningTrap::force_from_particles(Particle pcheck){
    arma::vec force_from_others = {0, 0, 0};
    for (int i = 0; i < particles.size(); i++){
        // Finding distance between particles.
        arma::vec posdiff = pcheck.position - particles[i].position; // Difference in position vector.
        double distance2 = arma::dot(posdiff, posdiff);   // Distance raised to second power

        if (distance2 != 0){
            force_from_others += pcheck.charge*ke*particles[i].charge*posdiff/(pow(distance2, 2));
        }
        else{
            continue;
        }
    }
    return (force_from_others);
}

arma::vec PenningTrap::total_force(Particle pcheck){
    arma::vec total = {0, 0, 0};
    total += force_from_particles(pcheck);   // Force from other particles.
    total += pcheck.charge*find_efield(pcheck); // Force from external electric field.
    total += arma::cross(pcheck.velocity, find_mfield(pcheck)); // Force from magnetic field.
    return total;
}

arma::vec PenningTrap::total_acceleration(Particle pcheck){
    arma::vec force = total_force(pcheck);
    arma::vec accel = force/(pcheck.mass);
    return accel;
}

//
//  main.cpp
//  penningtrap
//
//  Created by Felix Aarekol Forseth on 30/09/2022.
//

#include <iostream>
#include "particle.hpp"
#include "penningTrap.hpp"

int main(int argc, const char * argv[]) {
    double cha = 5.0;
    double mas = 6.0;
    arma::vec posi = arma::linspace(1, 2, 3);
    arma::vec velo = arma::linspace(2, 3, 3);
    arma::vec forc = arma::linspace(2, 3, 3);
    arma::vec accel = arma::linspace(2, 3, 3);

    Particle calcium = Particle(cha, mas, posi, velo, forc, accel);
    PenningTrap pp = PenningTrap(0.2, 0.3, 0.4);
}

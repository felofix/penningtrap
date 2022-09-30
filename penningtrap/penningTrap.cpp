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

arma::vec PenningTrap::find_efield(arma::vec position){
    position(0) = position(0)*gamma;  // x
    position(1) = position(1)*gamma;  // y
    position(2) = position(2)*gamma*(-2); // z
    return position;
}

arma::vec PenningTrap::find_mfield(arma::vec position);
    position(0) = 0;  // x
    position(1) = 0;  // y
    position(2) = position(2)*B0; // z
    return position;
}

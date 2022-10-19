#include <iostream>
#include <string>
#include <array>
#include <armadillo>
#include <math.h>

arma::vec analytical(arma::vec t, double w_z, double z0){
    arma::vec analytical = z0*arma::cos(t*w_z);
    return analytical;
}

int main(){
    arma::vec t = arma::linspace(0,1,100);
    arma::vec test = analytical(t, 1,1);

    std::cout << test << std::endl;

};

// Hmmm... this
// HALLO
// Didn't mean to make your cry 
// MOMMMAAAA
// UUuuUUUUUu
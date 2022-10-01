//
//  solver.cpp
//  penningtrap
//
//  Created by Felix Aarekol Forseth on 30/09/2022.
//

#include "solver.hpp"

Solver::Solver(double t, double NN, std::string f){
    time = t;           // Length of simulation.
    N = NN;             // Steps of simulation.
    filename = f;  // Outputting data to datafiles folder.
    dt = time/N;    // timestep.
}


void Solver::SolveforwardEuler(PenningTrap pt){
    arma::vec timer = arma::linspace(0, time, N);        // timer.
    arma::mat positions(N, pt.particles.size()*3); // Create a position matrix for all particles.
    
    for (int i = 0; i < N; i++){ // For all timesteps loop.
        
        
        
    }
    
    
    // update positions int he very end, so that  we get the correct force between particles. 
}

// Stepping forward like Euler would.
arma::mat Solver::forwardEulerStep(PenningTrap pt, Particle p){
    // Number refers to which particle we are looking at.
    // Returns both velocity and position so that we can use it for Runge-Kutta aswell.
    arma::mat velpos(3, 2);
    velpos.col(0) = p.velocity + pt.total_acceleration(p)*dt;
    velpos.col(1) = p.position + p.velocity*dt;
    return velpos;
}

void Solver::writetofilefloat(arma::vec timer, arma::mat xyz){
    // Writing to file with float values.
    std::ofstream fw(filename, std::ofstream::out);  // Setting the stream to output.
    if (fw.is_open())
    {
      for (int i = 0; i < xyz.n_cols; i++) {
          fw << timer(i) << " " << xyz.row(i) << "\n";
      }
      fw.close();
    }
    else std::cout << "The file couldnt be opened. " << std::endl;
}



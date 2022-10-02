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
    arma::mat velocities(N, pt.particles.size()*3); // Create a position matrix for all particles.
    
    for (int i = 0; i < N; i++){ // For all timesteps.
        // Put it here so we get zeros if something isnt filled.
        arma::mat temp(3, pt.particles.size()*2);      // temporary particle positions and veloicities.
        for (int j = 0; j < pt.particles.size(); j++){ // For all particles.
            arma::mat velpos = forwardEulerStep(pt, pt.particles[j]); // Forward-Euler step for paarticle j.
            temp.col(j*2) = velpos.col(0);
            temp.col(j*2+1) = velpos.col(1);
        }
        
        // Actually changing velocities and positions.
        for (int p = 0; p < pt.particles.size(); p++){
            velocities(i, arma::span(3*p, 3*p + 2)) = pt.particles[p].velocity.as_row(); // Changing values in positions matrix.
            positions(i, arma::span(3*p, 3*p + 2)) = pt.particles[p].position.as_row(); // Changing values in positions matrix.
            pt.particles[p].velocity = temp.col(p*2);
            pt.particles[p].position = temp.col(p*2+1);
        }
    }
    
    // Printing data to datafiles.
    for (int p = 0; p < pt.particles.size(); p++){
        std::string postring = "datafiles/pos_part_" + std::to_string(p) + filename;
        std::string velstring = "datafiles/vel_part_" + std::to_string(p) + filename;
        writetofilefloat(timer, positions(arma::span(0, N-1), arma::span(3*p, 3*p + 2)), postring);
        writetofilefloat(timer, velocities(arma::span(0, N-1), arma::span(3*p, 3*p + 2)), velstring);
    }
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


void Solver::writetofilefloat(arma::vec timer, arma::mat xyz, std::string direc){
    // Writing to file with float values.
    std::ofstream fw(direc, std::ofstream::out);  // Setting the stream to output.
    if (fw.is_open())
    {
      for (int i = 0; i < xyz.n_cols; i++) {
          fw << timer(i) << " " << xyz.row(i) << "\n";
      }
      fw.close();
    }
    else std::cout << "The file couldnt be opened. " << std::endl;
}



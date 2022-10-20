//
//  solver.cpp
//  penningtrap
//
//  Created by Felix Aarekol Forseth on 30/09/2022.
//

#include "solver.hpp"

Solver::Solver(double t, double NN, std::string fn, bool td, double amp, double oV, bool wrfile){
    time = t;           // Length of simulation.
    N = NN;             // Steps of simulation.
    filename = fn;  // Outputting data to datafiles folder.
    dt = time/N;    // timestep.
    etimedep = td;
    f = amp;
    omegaV = oV;
    writetofile = wrfile;
}


void Solver::RungeKuttaW(PenningTrap pt){
    arma::vec timer = arma::linspace(0, time, N);
    arma::mat positions(N, pt.particles.size()*3);
    arma::mat velocities(N, pt.particles.size()*3);
    double initV0 = pt.V0; // For alternativ efield.
    
    for (int i = 0; i < N; i++){
        arma::mat clone(3, pt.particles.size()*2); // Original position and velocity.
        arma::mat k1(3, pt.particles.size()*2);
        arma::mat k2(3, pt.particles.size()*2);
        arma::mat k3(3, pt.particles.size()*2);
        arma::mat k4(3, pt.particles.size()*2);
        arma::mat temp(3, pt.particles.size()*2);
        
        if (etimedep == true){
            pt.V0 = initV0*(1 + f*cos(omegaV*(timer(i))));
        }
        
        for (int p = 0; p < pt.particles.size(); p++){ //r_i v_i
            clone.col(p*2) = pt.particles[p].velocity;     // velocity
            clone.col(p*2 + 1) = pt.particles[p].position; // position
        }
        
        for (int p = 0; p < pt.particles.size(); p++){
            k1.col(p*2) = clone.col(p*2);                           // Velocity.
            k1.col(p*2+1) = pt.total_acceleration(pt.particles[p]); // acceleration.
            arma::mat velpos = forwardRKStepW(pt.particles[p], dt/2, clone.col(p*2), clone.col(p*2+1), k1.col(p*2+1)); // first step.
            temp.col(p*2) = velpos.col(0);   // Velocity.
            temp.col(p*2+1) = velpos.col(1); // Position.
        }
        
        updateposition(pt, temp);
        
        for (int p = 0; p < pt.particles.size(); p++){
            k2.col(p*2) = pt.particles[p].velocity; // updating k2 velocity.
            k2.col(p*2+1) = pt.total_acceleration(pt.particles[p]); // updating k2 accel.
            arma::mat velpos = forwardRKStepW(pt.particles[p], dt/2, clone.col(p*2), clone.col(p*2+1), k2.col(p*2+1)); // second step
            temp.col(p*2) = velpos.col(0);   // Velocity
            temp.col(p*2+1) = velpos.col(1); // Position.
        }
        
        updateposition(pt, temp);
        
        for (int p = 0; p < pt.particles.size(); p++){
            k3.col(p*2) = pt.particles[p].velocity; // updating k3 velocity.
            k3.col(p*2+1) = pt.total_acceleration(pt.particles[p]);
            arma::mat velpos = forwardRKStepW(pt.particles[p], dt, clone.col(p*2), clone.col(p*2+1), k3.col(p*2+1)); // third step
            temp.col(p*2) = velpos.col(0); // Velocity
            temp.col(p*2+1) = velpos.col(1); // Position.
        }
        
        updateposition(pt, temp);
        
        for (int p = 0; p < pt.particles.size(); p++){
            k4.col(p*2) = pt.particles[p].velocity;
            k4.col(p*2+1) = pt.total_acceleration(pt.particles[p]);
        }
        
        arma::mat godupdate = (k1 + 2*k2 + 2*k3 + k4)/6;
        
        for (int p = 0; p < pt.particles.size(); p++){   // Chaning the values of the particles.
            pt.particles[p].position = clone.col(p*2+1) + godupdate.col(0)*dt;
            pt.particles[p].velocity = clone.col(p*2) + godupdate.col(1)*dt;
            velocities(i, arma::span(3*p, 3*p + 2)) = pt.particles[p].velocity.as_row(); // Changing values in positions matrix.
            positions(i, arma::span(3*p, 3*p + 2)) = pt.particles[p].position.as_row(); // Changing values in positions matrix.
        }
    }
    
    if (writetofile == true){
        // Printing data to datafiles.
        for (int p = 0; p < pt.particles.size(); p++){
            std::string postring = "datafiles/pos_runge_part_" + std::to_string(p) + filename;
            std::string velstring = "datafiles/vel_runge_part_" + std::to_string(p) + filename;
            writetofilefloat(timer, positions(arma::span(0, N-1), arma::span(3*p, 3*p + 2)), postring);
            writetofilefloat(timer, velocities(arma::span(0, N-1), arma::span(3*p, 3*p + 2)), velstring);
        }
    }
    
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
            pt.particles[p].velocity = temp.col(p*2);
            pt.particles[p].position = temp.col(p*2+1);
            velocities(i, arma::span(3*p, 3*p + 2)) = pt.particles[p].velocity.as_row(); // Changing values in positions matrix.
            positions(i, arma::span(3*p, 3*p + 2)) = pt.particles[p].position.as_row(); // Changing values in positions matrix.
        }
    }
    // Printing data to datafiles.
    
    for (int p = 0; p < pt.particles.size(); p++){
        std::string postring = "datafiles/pos_euler_part_" + std::to_string(p) + filename;
        std::string velstring = "datafiles/vel_euler_part_" + std::to_string(p) + filename;
        writetofilefloat(timer, positions(arma::span(0, N-1), arma::span(3*p, 3*p + 2)), postring);
        writetofilefloat(timer, velocities(arma::span(0, N-1), arma::span(3*p, 3*p + 2)), velstring);
    }
}

void Solver::updateposition(PenningTrap &pt, arma::mat velpos){
    for (int p = 0; p < pt.particles.size(); p++){ // Changing the values of the particles.
        pt.particles[p].velocity = velpos.col(p*2);
        pt.particles[p].position = velpos.col(p*2+1);
    }
}

arma::mat Solver::forwardRKStepW(Particle p, double ddt, arma::vec vel, arma::vec pos, arma::vec a){
    arma::mat velpos(3, 2);
    velpos.col(1) = pos + p.velocity*ddt; // Velocity.
    velpos.col(0) = vel + a*ddt; // Position.
    return velpos;
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
      for (int i = 0; i < xyz.n_rows; i++) {
          fw << timer(i) << " " << xyz.row(i) << "\n";
      }
      fw.close();
    }
    else std::cout << "The file couldnt be opened. " << std::endl;
}



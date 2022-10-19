//
//  solver.cpp
//  penningtrap
//
//  Created by Felix Aarekol Forseth on 30/09/2022.
//

#include "solver.hpp"

Solver::Solver(double t, double NN, std::string fn, bool td, double amp, double oV){
    time = t;           // Length of simulation.
    N = NN;             // Steps of simulation.
    filename = fn;  // Outputting data to datafiles folder.
    dt = time/N;    // timestep.
    etimedep = td;
    f = amp;
    omegaV = oV;
}


void Solver::RungeKutta(PenningTrap pt){
    arma::vec timer = arma::linspace(0, time, N);        // timer.
    arma::mat positions(N, pt.particles.size()*3); // Create a position matrix for all particles.
    arma::mat velocities(N, pt.particles.size()*3); // Create a position matrix for all particles.
    double initV0 = pt.V0; // For alternativ efield.
    
    for (int i = 0; i < N; i++){ // For all timesteps.
        
        arma::mat clonevelpos(3, pt.particles.size()*2);      // first temp step of runge.
        
        for (int p = 0; p < pt.particles.size(); p++){   // Chaning the values of the particles.
            clonevelpos.col(p*2) = pt.particles[p].velocity; // velocity
            clonevelpos.col(p*2 + 1) = pt.particles[p].position; // position
        }
        arma::mat velacc1(3, pt.particles.size()*2);      // first temp step of runge.
        arma::mat velacc2(3, pt.particles.size()*2);      // second temp step of runge.
        arma::mat velacc3(3, pt.particles.size()*2);      // third temp step of runge.
        arma::mat velacc4(3, pt.particles.size()*2);      // fourth temp step of runge.
        
        // Put it here so we get zeros if something isnt filled.
        arma::mat temp2(3, pt.particles.size()*2);      // second temp step of runge.
        arma::mat temp3(3, pt.particles.size()*2);      // third temp step of runge.
        arma::mat temp4(3, pt.particles.size()*2);      // fourth temp step of runge.
        
        if (etimedep == true){
            pt.V0 = initV0*(1 + f*cos(omegaV*(timer(i))));
        }
        
        // Runge kutta.
        for (int r = 0; r < 3; r++){
            switch (r) {
                case 0: // calculating values at ti.
                    for (int j = 0; j < pt.particles.size(); j++){ // For all particles.
                        velacc1.col(j*2) = clonevelpos.col(0);                        // first vel.
                        velacc1.col(j*2+1) = pt.total_acceleration(pt.particles[j]);  // first accel.
                        
                        // First runge half-step for ti + dt/2.
                        arma::mat velpos = forwardRKStep(pt, pt.particles[j], dt/2, clonevelpos.col(j*2), clonevelpos.col(j*2 + 1));
                        temp2.col(j*2) = velpos.col(0);                              // tempvel
                        temp2.col(j*2+1) = velpos.col(1);                            // tempvelzq

                        velacc2.col(j*2) = velpos.col(0);                            // seocond vel.
                        velacc2.col(j*2+1) = pt.total_acceleration(pt.particles[j]); // seocond accel.
                    }
                    
                    for (int p = 0; p < pt.particles.size(); p++){ // Chaning the values of the particles.
                        pt.particles[p].velocity = temp2.col(p*2);
                        pt.particles[p].position = temp2.col(p*2+1);
                    }
                    break;

                case 1: // second time step, using runge ti + dt/2.
                    for (int j = 0; j < pt.particles.size(); j++){ // For all particles.
                        arma::mat velpos = forwardRKStep(pt, pt.particles[j], dt/2, clonevelpos.col(j*2), clonevelpos.col(j*2 + 1));      // Runge-kutta step for particle j.
                        temp3.col(j*2) = velpos.col(0);
                        temp3.col(j*2+1) = velpos.col(1);
                        velacc3.col(j*2) = velpos.col(0);                           // third vel.
                        velacc3.col(j*2+1) = pt.total_acceleration(pt.particles[j]);// third accel.
                    }
                    
                    for (int p = 0; p < pt.particles.size(); p++){// Chaning the values of the particles.
                        pt.particles[p].velocity = temp3.col(p*2);
                        pt.particles[p].position = temp3.col(p*2+1);
                    }
                    break;
                    
                case 2: // Fourth time step, using regular euler for ti + dt. 
                    for (int j = 0; j < pt.particles.size(); j++){ // For all particles.
                        arma::mat velpos = forwardRKStep(pt, pt.particles[j], dt, clonevelpos.col(j*2), clonevelpos.col(j*2 + 1)); // Forward-Euler step for particle j.
                        temp4.col(j*2) = velpos.col(0);
                        temp4.col(j*2+1) = velpos.col(1);
                        velacc4.col(j*2) = velpos.col(0);                           // fourth vel.
                        velacc4.col(j*2+1) = pt.total_acceleration(pt.particles[j]);// fourth accel.
                    }
                    break;
            }
        }
        
        arma::mat rungetemp = (velacc1 + 2*velacc2 + 2*velacc3 + velacc4)/6; // Runge-step.
        
        for (int p = 0; p < pt.particles.size(); p++){   // Chaning the values of the particles.
            pt.particles[p].velocity = clonevelpos.col(p*2) + rungetemp.col(1)*dt;
            pt.particles[p].position = clonevelpos.col(p*2 + 1) + rungetemp.col(0)*dt;
            velocities(i, arma::span(3*p, 3*p + 2)) = pt.particles[p].velocity.as_row(); // Changing values in positions matrix.
            positions(i, arma::span(3*p, 3*p + 2)) = pt.particles[p].position.as_row(); // Changing values in positions matrix.
        }
    }
    
    // Printing data to datafiles.
    for (int p = 0; p < pt.particles.size(); p++){
        std::string postring = "datafiles/pos_runge_part_" + std::to_string(p) + filename;
        std::string velstring = "datafiles/vel_runge_part_" + std::to_string(p) + filename;
        writetofilefloat(timer, positions(arma::span(0, N-1), arma::span(3*p, 3*p + 2)), postring);
        writetofilefloat(timer, velocities(arma::span(0, N-1), arma::span(3*p, 3*p + 2)), velstring);
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

arma::mat Solver::forwardRKStep(PenningTrap pt, Particle p, double ddt, arma::vec vel, arma::vec pos){
    // Returns both velocity and position so that we can use it for Runge-Kutta aswell.
    arma::mat velpos(3, 2);
    velpos.col(1) = pos + p.velocity*ddt;
    velpos.col(0) = vel + pt.total_acceleration(p)*ddt;
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



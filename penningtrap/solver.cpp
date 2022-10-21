//
//  solver.cpp
//  penningtrap
//
//  Created by Felix Aarekol Forseth on 30/09/2022.
//

#include "solver.hpp"

Solver::Solver(double t, double NN, std::string fn, bool td, double amp, double oV, bool wrfile, bool err){
    time = t;           // Length of simulation.
    N = NN;             // Steps of simulation.
    filename = fn;  // Outputting data to datafiles folder.
    dt = time/N;    // timestep.
    etimedep = td;
    f = amp;
    omegaV = oV;
    writetofile = wrfile;
    errorplot = err;
}

void Solver::RungeKuttaW(PenningTrap &pt){
    arma::vec timer = arma::linspace(0, time, N);
    arma::mat positions(N, pt.particles.size()*3);
    arma::mat velocities(N, pt.particles.size()*3);
    double initV0 = pt.V0; // For alternativ efield.
    for (int p = 0; p < pt.particles.size(); p++){ // initial values.
        velocities(0, arma::span(3*p, 3*p + 2)) = pt.particles[p].velocity.as_row(); // Changing values in positions matrix.
        positions(0, arma::span(3*p, 3*p + 2)) = pt.particles[p].position.as_row(); // Changing values in positions matrix.
    }
    
    for (int i = 1; i < N; i++){
        arma::mat clone(3, pt.particles.size()*2); // Original position and velocity.
        arma::mat k1(3, pt.particles.size()*2);
        arma::mat k2(3, pt.particles.size()*2);
        arma::mat k3(3, pt.particles.size()*2);
        arma::mat k4(3, pt.particles.size()*2);
        arma::mat temp(3, pt.particles.size()*2);
        
        if (etimedep == true){
            pt.gamma = initV0*(1 + f*cos(omegaV*(timer(i))))/pow(pt.d, 2);
        }
    
        for (int p = 0; p < pt.particles.size(); p++){
            clone.col(p*2) = pt.particles[p].velocity;     // velocity
            clone.col(p*2 + 1) = pt.particles[p].position; // position
        }
        
        for (int p = 0; p < pt.particles.size(); p++){
            // Updating k1 for all particles.
            k1.col(p*2) = clone.col(p*2);                           // Velocity.
            k1.col(p*2+1) = pt.total_acceleration(pt.particles[p]); // acceleration.
            // Forward half-step.
            arma::mat velpos = forwardRKStepW(pt.particles[p], dt/2, clone.col(p*2), clone.col(p*2+1), k1.col(p*2+1)); // first step.
            // Changing temp values.
            temp.col(p*2) = velpos.col(0);       // Velocity.
            temp.col(p*2+1) = velpos.col(1); // Position.
        }
        
        // Updating position of the particles.
        updateposition(pt, temp);
        
        for (int p = 0; p < pt.particles.size(); p++){
            // Updating k2 for all particles.
            k2.col(p*2) = pt.particles[p].velocity; // updating k2 velocity.
            k2.col(p*2+1) = pt.total_acceleration(pt.particles[p]); // updating k2 accel.
            // Forward half-step.
            arma::mat velpos = forwardRKStepW(pt.particles[p], dt/2, clone.col(p*2), clone.col(p*2+1), k2.col(p*2+1)); // second step
            // Changing temp values.
            temp.col(p*2) = velpos.col(0);   // Velocity
            temp.col(p*2+1) = velpos.col(1); // Position.
        }
        
        // Updating position of the particles.
        updateposition(pt, temp);
        
        for (int p = 0; p < pt.particles.size(); p++){
            // Updating k3 for all particles.
            k3.col(p*2) = pt.particles[p].velocity; // updating k3 velocity.
            k3.col(p*2+1) = pt.total_acceleration(pt.particles[p]);
            // Forward half-step.
            arma::mat velpos = forwardRKStepW(pt.particles[p], dt, clone.col(p*2), clone.col(p*2+1), k3.col(p*2+1)); // third step
            // Changing temp values.
            temp.col(p*2) = velpos.col(0); // Velocity
            temp.col(p*2+1) = velpos.col(1); // Position.
        }
        
        // Updating position of the particles.
        updateposition(pt, temp);
        
        for (int p = 0; p < pt.particles.size(); p++){
            // Updating k4 for all particles.
            k4.col(p*2) = pt.particles[p].velocity;
            k4.col(p*2+1) = pt.total_acceleration(pt.particles[p]);
        }
        
        arma::mat godupdate = (k1 + 2*k2 + 2*k3 + k4)/6;
        
        for (int p = 0; p < pt.particles.size(); p++){   // Chaning the values of the particles.
            pt.particles[p].position = clone.col(p*2+1) + godupdate.col(p*2)*dt;
            pt.particles[p].velocity = clone.col(p*2) + godupdate.col(p*2+1)*dt;
            velocities(i, arma::span(3*p, 3*p + 2)) = pt.particles[p].velocity.as_row(); // Changing values in positions matrix.
            positions(i, arma::span(3*p, 3*p + 2)) = pt.particles[p].position.as_row(); // Changing values in positions matrix.
        }
    }
    
    if (errorplot == true){
        solve_analytically(pt, positions, "runge");
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
    
    // initial values.
    for (int p = 0; p < pt.particles.size(); p++){
        velocities(0, arma::span(3*p, 3*p + 2)) = pt.particles[p].velocity.as_row(); // Changing values in positions matrix.
        positions(0, arma::span(3*p, 3*p + 2)) = pt.particles[p].position.as_row(); // Changing values in positions matrix.
    }
    
    for (int i = 1; i < N; i++){ // For all timesteps.
        // Put it here so we get zeros if something isnt filled.
        arma::mat temp(3, pt.particles.size()*2);      // temporary particle positions and veloicities.
        for (int j = 0; j < pt.particles.size(); j++){ // For all particles.
            arma::mat velpos = forwardEulerStep(pt, pt.particles[j]); // Forward-Euler step for particle j.
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
    
    if (errorplot == true){
        solve_analytically(pt, positions, "euler");
    }
    
    if (writetofile == true){
        for (int p = 0; p < pt.particles.size(); p++){
            std::string postring = "datafiles/pos_euler_part_" + std::to_string(p) + filename;
            std::string velstring = "datafiles/vel_euler_part_" + std::to_string(p) + filename;
            writetofilefloat(timer, positions(arma::span(0, N-1), arma::span(3*p, 3*p + 2)), postring);
            writetofilefloat(timer, velocities(arma::span(0, N-1), arma::span(3*p, 3*p + 2)), velstring);
        }
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

void Solver::solve_analytically(PenningTrap pt, arma::mat numeric, std::string euorru){
    arma::vec timearray = arma::linspace(0, time, N);
    arma::mat analytical = create_analytical(pt);
    std::string abserrorstring = "datafiles/abserrorN" + euorru + std::to_string((int)N) + ".txt";
    std::string relerrorstring = "datafiles/relerrorN" + euorru + std::to_string((int)N) + ".txt";
    arma::vec abserror(N);
    arma::vec relerror(N);
    
    for (int i = 0; i < N; i++){
        abserror[i] = arma::norm(analytical.row(i) - numeric.row(i), 2);
        relerror[i] = abserror[i]/(arma::norm(analytical.row(i), 2));
    }
    
    writetofilefloat(timearray, abserror, abserrorstring);
    writetofilefloat(timearray, relerror, relerrorstring);
}

arma::mat Solver::create_analytical(PenningTrap pt){
    // constants
    double omega2z = (2*pt.particles[0].charge*pt.V0)/(pt.particles[0].mass*pow((pt.d),2));
    double omega0 = (pt.particles[0].charge*pt.B0)/pt.particles[0].mass;
    double v0 = 25;
    double x0 = 20;
    double z0 = 20;
    double omegaplus = (omega0 + sqrt(pow(omega0,2) - 2*omega2z))/2;
    double omegaminus = (omega0 - sqrt(pow(omega0,2) - 2*omega2z))/2;
    double Ap = (v0 + omegaminus*x0)/(omegaminus - omegaplus);
    double Am = -(v0 + omegaplus*x0)/(omegaminus - omegaplus);
    arma::vec timearray = arma::linspace(0, time, N);
    arma::mat posvel(N, 3);
    std::complex<double> i(0.0, 1.0);
    
    for (int j = 0; j < N; j++){
        posvel.row(j)(2) = z0*cos(sqrt(omega2z)*timearray[j]);
        std::complex<double> xy = Ap*exp(-i*omegaplus*timearray[j]) + Am*exp(-i*omegaminus*timearray[j]);
        posvel.row(j)(0) = xy.real();
        posvel.row(j)(1) = xy.imag();
    }
    return posvel;
}

void Solver::writetofilefloat(arma::vec timer, arma::mat xyz, std::string direc){
    // Writing to file with float values.
    std::ofstream fw(direc, std::ofstream::out);  // Setting the stream to output.
    if (fw.is_open())
    {
      for (int i = 0; i < xyz.n_rows; i++) {
      fw << timer(i) << xyz.row(i) << "\n";
      }
      fw.close();
    }
    else std::cout << "The file couldnt be opened. " << std::endl;
}



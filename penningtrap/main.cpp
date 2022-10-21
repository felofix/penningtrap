//
//  main.cpp
//  penningtrap
//
//  Created by Felix Aarekol Forseth on 30/09/2022.
//

#include <iostream>
#include "particle.hpp"
#include "penningTrap.hpp"
#include "solver.hpp"

int main(int argc, const char * argv[]) {
    // Constants.
    double c = 1.0;   // charge of calcium.
    double m = 40.078; // mass of calcium.
    double T = 9.64852558e1; // Tesla
    double V = 9.64852558e7;
    double B0 = 1*T;
    double V0 = 25e-3*V;
    double d = 500;
    double p1x0 = 20;
    double p1z0 = 20;
    double p1ydot = 25;
    double p2x0 = 25;
    double p2y0 = 25;
    double p2ydot = 40;
    double p2zdot = 5;
    double time = 50; // microseconds.
    
    // number of steps
    double n = 15000;
    double n1 = 4000;
    double n2 = 8000;
    double n3 = 16000;
    double n4 = 32000;
    
    // Classes.
    arma::vec posi1 = {p1x0, 0, p1z0};
    arma::vec posi2 = {p2x0, p2y0, 0};
    arma::vec vel = {0, p1ydot, 0};
    arma::vec vel2 = {0, p2ydot, p2zdot};
    Particle calcium1 = Particle(c, m, posi1, vel);
    Particle calcium2 = Particle(c, m, posi2, vel2);
    PenningTrap pp = PenningTrap(B0, V0, d, true);  // For just the first particle.
    PenningTrap pp2 = PenningTrap(B0, V0, d, true); // For just the second particle.
    PenningTrap pp3 = PenningTrap(B0, V0, d, true); // For both particles.
    PenningTrap pperror1 = PenningTrap(B0, V0, d, true); // for error.
    PenningTrap pperror2 = PenningTrap(B0, V0, d, true); // for error.
    PenningTrap pperror3 = PenningTrap(B0, V0, d, true); // for error.
    PenningTrap pperror4 = PenningTrap(B0, V0, d, true); // for error.
    
    // Adding particles.
    pp.add_particle(calcium1);
    pp2.add_particle(calcium2);
    pp3.add_particle(calcium1);
    pp3.add_particle(calcium2);
    pperror1.add_particle(calcium1);
    pperror2.add_particle(calcium1);
    pperror3.add_particle(calcium1);
    pperror4.add_particle(calcium1);
    
    // Solving for particles.
    Solver solve = Solver(time, n, "Single1.txt", false, 0, 0, true, false);
    Solver solve2 = Solver(time, n, "Single2.txt", false, 0, 0, true, false);
    Solver solve3 = Solver(time, n, "Both.txt", false, 0, 0, true, false);
    
    // Different stepsizes to compare with analytical soulution.
    Solver solveerror1 = Solver(time, n1, "n1.txt", false, 0, 0, false, true);
    Solver solveerror2 = Solver(time, n2, "n2.txt", false, 0, 0, false, true);
    Solver solveerror3 = Solver(time, n3, "n3.txt", false, 0, 0, false, true);
    Solver solveerror4 = Solver(time, n4, "n4.txt", false, 0, 0, false, true);
    
    // Actually solving for one and two particle.
    /*
    solve.SolveforwardEuler(pp);
    solve.RungeKuttaW(pp);
    solve2.SolveforwardEuler(pp2);
    solve2.RungeKuttaW(pp2);
    solve3.SolveforwardEuler(pp3);
    solve3.RungeKuttaW(pp3);
    */
    
    // Actually solving for the single particle to compare with analytical solution.
    solveerror1.SolveforwardEuler(pperror1);
    solveerror1.RungeKuttaW(pperror1);
    solveerror2.SolveforwardEuler(pperror2);
    solveerror2.RungeKuttaW(pperror2);
    solveerror3.SolveforwardEuler(pperror3);
    solveerror3.RungeKuttaW(pperror3);
    solveerror4.SolveforwardEuler(pperror4);
    solveerror4.RungeKuttaW(pperror4);
    
    /*
    // Simulation with 100 particles.
    int newtime = 500;
    int steps = 10;
    int particles = 1;
    std::vector<double> f = {0.1, 0.4, 0.7};
    for (int j = 0; j < f.size(); j++){
        arma::vec wz = arma::linspace(0.2, 2.5, steps);
        arma::vec pleft(steps, arma::fill::zeros);
        PenningTrap pt100 =  PenningTrap(B0, V0, d, false);
        pt100.fill_with_particles(particles, c, m);

        for (int i = 0; i < steps; i++){ // For f = 0.1.
            std::cout<< i << std::endl;
            Solver s100 = Solver(newtime, 20000, "null", true, f[j], wz[i], false);
            s100.RungeKuttaW(pt100);
            pleft[i] = pt100.count_particles();
        }
        
        std::ofstream fw("datafiles/f" + std::to_string(f[j] + ".txt", std::ofstream::out);  // Setting the stream to output.
        if (fw.is_open())
        {
          for (int g = 0; g < wz.n_rows; g++) {
              fw << wz(g) << " " << pleft(g) << "\n";
          }
          fw.close();
        }
        else std::cout << "The file couldnt be opened. " << std::endl;
    }
    */
}


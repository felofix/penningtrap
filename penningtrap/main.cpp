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
    
    // Adding particles.
    pp.add_particle(calcium1);
    pp2.add_particle(calcium2);
    pp3.add_particle(calcium1);
    pp3.add_particle(calcium2);
    
    // Solving for particles.
    Solver solve = Solver(time, n, "Single1.txt", false, 0, 0);
    Solver solve2 = Solver(time, n, "Single2.txt", false, 0, 0);
    Solver solve3 = Solver(time, n, "Both.txt", false, 0, 0);
    
    // Different stepsizes to compare with analytical soulution.
    Solver solveerror1 = Solver(time, n1, "n1.txt", false, 0, 0);
    Solver solveerror2 = Solver(time, n2, "n2.txt", false, 0, 0);
    Solver solveerror3 = Solver(time, n3, "n3.txt", false, 0, 0);
    Solver solveerror4 = Solver(time, n4, "n4.txt", false, 0, 0);
    
    // Actually solving for one and two particle.
    solve.SolveforwardEuler(pp);
    solve.RungeKutta(pp);
    solve2.SolveforwardEuler(pp2);
    solve2.RungeKutta(pp2);
    solve3.SolveforwardEuler(pp3);
    solve3.RungeKutta(pp3);
    
    // Actually solving for the single particle to compare with analytical solution.
    solveerror1.SolveforwardEuler(pp);
    solveerror1.RungeKutta(pp);
    solveerror2.SolveforwardEuler(pp);
    solveerror2.RungeKutta(pp);
    solveerror3.SolveforwardEuler(pp);
    solveerror3.RungeKutta(pp);
    solveerror4.SolveforwardEuler(pp);
    solveerror4.RungeKutta(pp);
}


/*
    FLIP FLuid Simulation: PBS Project
    Created by Tom Lausberg & Safira Piasko
    19.12.2018
*/

#include <iostream>
#include <eigen3/Eigen/Dense>
#include <eigen3/unsupported/Eigen/CXX11/Tensor>

#include "simulation.hpp"
#include "initialiser.hpp"
#include "util.hpp"
#include "grid.hpp"
#include "writer.hpp"

int main() {
    // timestep
    double dt = 0.1;
    // length of a cell
    double dx = 1;
    // density
    double rho = 1;

    Simulation sim(dt, dx, rho);

    // initialise particle to start the simulation with (5x5x5 box at positon (10,10,10) )
    std::vector<Particle> parts = initialise_NPerCell(5, 5, 5, 2, dx, Vec3(10*dx,10*dx,10*dx));

    // add particles to simulation
    sim.add_particles(parts);

    //run simulation for 200 timesteps and output data in every timestep
    sim.run(200,1);

    return 0;
}

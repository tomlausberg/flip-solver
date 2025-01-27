/*
    FLIP FLuid Simulation: PBS Project
    Created by Tom Lausberg & Safira Piasko
    19.12.2018
*/

#ifndef INITIALISER_H
#define INITIALISER_H
#include <eigen3/Eigen/Dense>
#include "util.hpp"
#include "simulation.hpp"
#include <vector>

std::vector<Particle> initialise_rectangle(int x, int y, int z, double spacing, Vec3 origin);

// void initialise_cylinder(double height, double radius, Coord origin, int nop);

//starts grid with n by n particles per cell ( n is usually 2,3)
std::vector<Particle> initialise_NPerCell(int x, int y, int z, int n, double dx, Vec3 origin);

#endif

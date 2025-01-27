/*
    FLIP FLuid Simulation: PBS Project
    Created by Tom Lausberg & Safira Piasko
    19.12.2018
*/

#include "util.hpp"
#include "simulation.hpp"
#include <cassert>
//Vec3 class constructors
Vec3::Vec3(void) : x(0), y(0), z(0) {}

Vec3::Vec3(double n) : x(n), y(n), z(n) {}

Vec3::Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

//Particle constructors
Particle::Particle(Vec3 pos_) : pos(pos_) {}

Particle::Particle(Vec3 pos_, Vec3 vel_) : pos(pos_), velocity(vel_) {}

// --------------------------------- helper functions --------------------------------

int clamp(int a, int lower, int upper) {
     if(a<lower) return lower;
     else if(a>upper) return upper;
     else return a;
}

double clamp(double a, double lower, double upper) {
     if(a<lower) return lower;
     else if(a>upper) return upper;
     else return a;
}

// find grid cell in which the particle is and find the difference to the grid points
void find_cell(double dx, double x, double& du, int& i) {
    double sx = x/dx;
    i = (int)sx;
    du = sx - floor(sx);
}

// same as above but doesn't return du
void find_cell(double dx, double x, int& i) {
    i = int(x/dx);
}

// for MACgrid
void find_cell_center(double dx, double x, double& du, int& i){
    double sx = (x/dx) - 0.5;
    i = (int)sx;
    // boundary conditions
    if(i<0){
        i=0;
        du=0;
    } else if(i>X_DIM-2) {
        i=X_DIM-2;
        du=1.;
    } else {
        du = sx-floor(sx);
    }
}

double cubicSplineKernel(double distance){
    double weight;
    if(distance < 0.5) {
        //falls 0 < dist < 0.5 --> weight = 1 - 6*d^2 - 6*d^3
        weight = 6 * distance*distance;
        weight -= weight * distance;
        return 1 - weight;
    }
    if(distance < 1) {
        //falls 0.5 < dist < 1 --> weight = 2*(1-d)^3
        weight = 1 - distance;
        weight = weight*weight*weight;
        return 2*weight;
    }
    // falls 1 <= dist --> weight = 0
    return 0;
}

/*
    FLIP FLuid Simulation: PBS Project
    Created by Tom Lausberg & Safira Piasko
    19.12.2018
*/

#ifndef SIMULATION_H
#define SIMULATION_H

#include "util.hpp"
#include "grid.hpp"
#include "timer.hpp"

#include <vector>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/SparseCore>
#include <eigen3/Eigen/IterativeLinearSolvers>
#include <algorithm>

enum cell_type {
  AIR,
  FLUID,
  SURFACE,
  SOLID
};

// Simulation: holds all nessesary data and members for interfacin
class Simulation {
private:
    /* data */
    Grid<double> pressure; // stores pressures in cells
    Grid<double> velx; // stores velocity in x direction in cells
    Grid<double> vely; // stores velocity in x direction in cells
    Grid<double> velz; // stores velocity in x direction in cells

    Grid<cell_type> cell;     // stores celltype in cells
    Grid<int> fluid_idx;     // stores fluid cell index
    Grid<int> particle_count; //stores number of particles in each cell;

    Grid<double> sum;
    Grid<double> sumx;
    Grid<double> sumy;
    Grid<double> sumz;
    Grid<double> x_dist_front; // stores distance to fluid cells
    Grid<double> x_dist_back;
    Grid<double> y_dist_front; // stores distance to fluid cells
    Grid<double> y_dist_back;
    Grid<double> z_dist_front; // stores distance to fluid cells
    Grid<double> z_dist_back;

    std::vector<Particle> particles;

    Grid<double> pre_velx;
    Grid<double> pre_vely;
    Grid<double> pre_velz;

    Grid<double> diff_velx;
    Grid<double> diff_vely;
    Grid<double> diff_velz;

    Grid<double> distance;
    Grid<int> index;

    // solves of poisson pressure equation
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double> > cg_solver;
    // laplacian matrix
    Eigen::SparseMatrix<double> A;

    int max_step;
    double dt;
    double dx;
    double density;

    Timer tTotal;
    Timer tPressure;
    Timer textendVelocity;
    Timer tAdvection;
    Timer tStep;

    int write_frequency;

public:
    //constructor/destructor = Simulation
    int current_step;
    Simulation (double dt_, double dx_, double density_);

    void source_particles(int x, int y, int z, int n, Vec3 origin);
    void source_cylinder(int r, Vec3 origin);
    void add_particles(std::vector<Particle> parts);
    std::vector<Particle> get_particles();

    // run simulation
    void run(int steps, int steps_per_frame);
    // advance one step in simulation
    void step();

    //vtk writer functions;
    void write_frame();
    void write_pressure();
    void write_cells();
    void write_particle_count();
    void write_scalars();
    void write_velocity();
    // functions to calculate new velocity and pressure
    void particleToGrid();
    void cubicSplineToGrid();
    void addGravity();
    void gridVelocityToParticles();
    void updateParticlePositions();
    void DirichletBoundary();
    void extendVelocity();
    void pressure_solver();
    void simpleAdvect();
    void levelset();
    void extend_velocity();
    void applyBC();

    void marchingCubes();

    // include helper functions
    void extend_u(int x_start, int x_end, int y_start, int y_end, int z_start, int z_end);
    void extend_v(int x_start, int x_end, int y_start, int y_end, int z_start, int z_end);
    void extend_w(int x_start, int x_end, int y_start, int y_end, int z_start, int z_end);

    bool IsInGrid(Particle&);
    void checkParticles();

    void swapVelocity();
    void swapDistance();

    void getValidNeighborValues(double, int, int, int, double&, int&);
    template <class T>
    void velAtEdges(Grid<T>&, double, int, int, int, double, double, double);
};
void solve_distance(double, double, double, double&, double);

#endif

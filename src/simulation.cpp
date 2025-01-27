/*
    FLIP FLuid Simulation: PBS Project
    Created by Tom Lausberg & Safira Piasko
    19.12.2018
*/

#include <vector>
#include <iostream>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <cstdio>

#include <eigen3/Eigen/Dense>
// #include <eigen3/unsupported/Eigen/CXX11/Tensor>
#include <eigen3/Eigen/SparseCore>
#include <eigen3/Eigen/IterativeLinearSolvers>

#include "initialiser.hpp"
#include "simulation.hpp"
#include "util.hpp"
#include "grid.hpp"
#include "writer.hpp"
#include "timer.hpp"
#include "makePLY.hpp"

// constructor for simulation
Simulation::Simulation(double dt_, double dx_, double density_) : dt(dt_), dx(dx_), density(density_) {
    pressure = Grid<double>();
    velx = Grid<double>();
    vely = Grid<double>();
    velz = Grid<double>();

    cell = Grid<cell_type>();
    fluid_idx = Grid<int>();
    particle_count = Grid<int>();

    sum = Grid<double>();
    x_dist_front = Grid<double>();
    x_dist_back =Grid<double>();
    y_dist_front = Grid<double>();
    y_dist_back = Grid<double>();
    z_dist_front = Grid<double>();
    z_dist_back = Grid<double>();

    diff_velx = Grid<double>();
    diff_vely = Grid<double>();
    diff_velz = Grid<double>();

    current_step = 0;
    write_frequency = 1;
    max_step = 300;
}

// add & get particles
void Simulation::add_particles(std::vector<Particle> parts){
    for(auto part : parts){
        if(IsInGrid(part)){
            particles.push_back(part);
        }
        else {
            std::cout << "particle not in grid" << std::endl;
        }
    }
  std::cout << "added particles: " << parts.size() << '\n';
}

std::vector<Particle> Simulation::get_particles(){
    return particles;
}

// advection step
void Simulation::simpleAdvect() {
  for(auto &part : particles) {
    part.pos.x = part.pos.x + part.velocity.x * dt;
    part.pos.y = part.pos.y + part.velocity.y * dt;
    part.pos.z = part.pos.z + part.velocity.z * dt;

  }
}

// the bounaries are solid and have no velocities
void Simulation::applyBC(){
    for(int k = 0; k < Z_DIM; ++k){
      for(int j = 0; j < Y_DIM; ++j) {
        for(int i = 0; i < X_DIM; ++i) {
            if(i <= 1 || i >= X_DIM - 1){
              if(i == 0 || i == X_DIM-1)
                cell(i,j,k) = SOLID;
              velx(i,j,k) = 0;
            }

            if(j <= 1 || j >= Y_DIM - 1){
              if(j == 0 || j == Y_DIM-1)
                cell(i,j,k) = SOLID;
              vely(i,j,k) = 0;
            }

            if(k <= 1 || k >= Z_DIM - 1){
              if(k == 0 || k == Z_DIM-1)
                cell(i,j,k) = SOLID;
              velz(i,j,k) = 0;
            }

          }
       }
    }
    // move particles back into the valid space, when they end up in the boundary
    for(auto& part : particles){
        if(part.pos.x < dx) {
            part.pos.x = dx*1.001;
        }
        else if(part.pos.x > dx*(X_DIM-1)) {
            part.pos.x = dx*(X_DIM-1.001);
        }

        if(part.pos.y < dx) {
            part.pos.y = dx*1.001;
        }
        else if(part.pos.y > dx*(Y_DIM-1)) {
            part.pos.y = dx*(Y_DIM-1.001);
        }

        if(part.pos.z < dx) {
            part.pos.z = dx*1.001;
        }
        else if(part.pos.z > dx*(Z_DIM-1)) {
            part.pos.z = dx*(Z_DIM-1.001);
        }
    }
}

// delete particles if they are out of the grid (for debugging)
void Simulation::checkParticles(){
    for( auto part = particles.begin(); part != particles.end();){
        if(part->pos.x > dx*X_DIM || part->pos.x < 0
        || part->pos.y > dx*Y_DIM || part->pos.y < 0
        || part->pos.z > dx*Z_DIM || part->pos.z < 0
        || std::isnan(part->pos.x)|| std::isnan(part->pos.x)|| std::isnan(part->pos.x)){
            std::cout << "particle out of grid" << '\n';
             part = particles.erase(part);
        }
        else {
            ++part;
        }
    }
}

// calculate velocities of grid cells with known velocities of particles
void Simulation::particleToGrid() {
  // du: Abstand von x Koordinate zu Grid koorinate von u (Am Rande vom Cell)
  // dpx: Abstand von x Koordinate zu Grid Koordinate von p (Mittelpunkt)
  double du, dv, dw;
  int i = 0;
  int j = 0;
  int k = 0;

  velx.zero();
  sum.zero();

  // velocity in x direction
  for (auto part : particles){
    find_cell(dx, part.pos.x, du, i);
    find_cell_center(dx, part.pos.y, dv, j);
    find_cell_center(dx, part.pos.z, dw, k);
    // find weights
    velAtEdges(velx, part.velocity.x, i, j, k, du, dv, dw);
  }

  for(k=1; k<Z_DIM-1; ++k) for(j=1; j<Y_DIM-1; ++j) for(i=1; i<X_DIM-1; ++i){
     if(sum(i,j,k)!=0)
        velx(i,j,k) /= sum(i,j,k);
  }

  vely.zero();
  sum.zero();

  // velocity in y direction
  for (auto part : particles){
    find_cell_center(dx, part.pos.x, du, i);
    find_cell(dx, part.pos.y, dv, j);
    find_cell_center(dx, part.pos.z, dw, k);
    // find weights
    velAtEdges(vely, part.velocity.y, i, j, k, du, dv, dw);
  }

  for(k=1; k<Z_DIM-1; ++k) for(j=1; j<Y_DIM-1; ++j) for(i=1; i<X_DIM-1; ++i){
     if(sum(i,j,k)!=0)
        vely(i,j,k) =  vely(i,j,k) / sum(i,j,k);
  }

  velz.zero();
  sum.zero();

  // velocity in z direction
  for (auto part : particles){
    find_cell_center(dx, part.pos.x, du, i);
    find_cell_center(dx, part.pos.y, dv, j);
    find_cell(dx, part.pos.z, dw, k);
    // find weights
    velAtEdges(velz, part.velocity.z, i, j, k, du, dv, dw);
  }

  for(k=1; k<Z_DIM-1; ++k) for(j=1; j<Y_DIM-1; ++j) for(i=1; i<X_DIM-1; ++i){
     if(sum(i,j,k)!=0) velz(i,j,k) = velz(i,j,k)/sum(i,j,k);
  }

  cell.zero();
  // find celltype grid
  for(auto part : particles){
    find_cell(dx, part.pos.x, i);
    find_cell(dx, part.pos.y, j);
    find_cell(dx, part.pos.z, k);

    cell(i,j,k) = FLUID;
  }
}

void Simulation::addGravity(){
    // accelerates the grid according to gravity.
    double g = 1;
    double gdt = dt*g;
    for(int i = 0; i<X_DIM; ++i){
        for(int j = 0; j<Y_DIM; ++j){
            for(int k = 0; k<Z_DIM; ++k){
                // if(cell(i,j,k) == FLUID)
                    velz(i,j,k) = velz(i,j,k)-gdt;
            }
        }
    }
}

void Simulation::levelset(){
  // initialize distance grid with maximum number
  int max_d = (X_DIM + Y_DIM + Z_DIM) * dx;
  distance.set_constant(max_d);

  // distance from point to cell face
  int i, j, k;
  double xdist, ydist, zdist;

  for(int it = 0; it < particles.size(); ++it){
    find_cell(dx, particles[it].pos.x, xdist, i);
    find_cell(dx, particles[it].pos.y, ydist, j);
    find_cell(dx, particles[it].pos.z, zdist, k);

    double dist = std::sqrt(xdist*xdist + ydist*ydist + zdist*zdist);

    if(dist < distance(i,j,k)){
        distance(i,j,k) = -dist;
    }
}
    // all fluid cells get distance -0.5 negative, because its inside the fluid
    // for(int i = 1; i < X_DIM-1; i++){
    //   for (int j  = 1; j  < Y_DIM-1; j++) {
    //     for (size_t k = 1; k < Z_DIM-1; k++) {
    //       if(cell(i,j,k) == FLUID)
    //         // indicator inside fluid
    //         distance(i,j,k) = -0.5;
    //     }
    //   }
    // }

    // propagate closest point and distance etimate to the rest of the grid points
    for(int l = 0; l < 2; ++l){
      // 1. sweep order
      for(int i = 1; i < X_DIM; ++i){
        for(int j = 1; j < Y_DIM; ++j){
          for (int k = 1; k < Z_DIM; k++) {
            if(cell(i,j,k) != FLUID)
              solve_distance(distance(i-1,j,k), distance(i,j-1,k), distance(i,j,k-1), distance(i,j,k),dx);
          }
        }
      }

      // 2. sweep order
      for(int i = 1; i < X_DIM; ++i){
        for(int j = 1; j < Y_DIM; ++j){
          for (int k = Z_DIM-2; k >= 0; k--) {
            if(cell(i,j,k) != FLUID)
              solve_distance(distance(i-1,j,k), distance(i,j-1,k), distance(i,j,k+1), distance(i,j,k),dx);
          }
        }
      }

      // 3. sweep order
      for(int i = 1; i < X_DIM; ++i){
        for(int j = Y_DIM-2; j >= 0; --j){
          for (int k = 1; k < Z_DIM; k++) {
            if(cell(i,j,k) != FLUID)
              solve_distance(distance(i-1,j,k), distance(i,j+1,k), distance(i,j,k-1), distance(i,j,k),dx);
          }
        }
      }

      // 4. sweep order
      for(int i = X_DIM-2; i >= 0; --i){
        for(int j = 1; j < Y_DIM; ++j){
          for (int k = 1; k < Z_DIM; k++) {
            if(cell(i,j,k) != FLUID)
              solve_distance(distance(i+1,j,k), distance(i,j-1,k), distance(i,j,k-1), distance(i,j,k),dx);
          }
        }
      }

      // 5. sweep order
      for(int i = 1; i < X_DIM; ++i){
        for(int j = Y_DIM-2; j >= 0; --j){
          for (int k = Z_DIM-2; k >= 0; k--) {
            if(cell(i,j,k) != FLUID)
              solve_distance(distance(i-1,j,k), distance(i,j+1,k), distance(i,j,k+1), distance(i,j,k),dx);
          }
        }
      }

      // 6. sweep order
      for(int i = X_DIM-2; i >= 0; --i){
        for(int j = 1; j < Y_DIM; ++j){
          for (int k = Z_DIM-2; k >= 0; k--) {
            if(cell(i,j,k) != FLUID)
              solve_distance(distance(i+1,j,k), distance(i,j-1,k), distance(i,j,k+1), distance(i,j,k),dx);
          }
        }
      }

      // 7. sweep order
      for(int i = X_DIM-2; i >= 0; --i){
        for(int j = Y_DIM-2; j >= 0; --j){
          for (int k = 1; k < Z_DIM; k++) {
            if(cell(i,j,k) != FLUID)
              solve_distance(distance(i+1,j,k), distance(i,j+1,k), distance(i,j,k-1), distance(i,j,k),dx);
          }
        }
      }

      // 8. sweep order
      for(int i = X_DIM-2; i >= 0; --i){
        for(int j = Y_DIM-2; j >= 0; --j){
          for (int k = Z_DIM-2; k >= 0; --k) {
            if(cell(i,j,k) != FLUID)
              solve_distance(distance(i+1,j,k), distance(i,j+1,k), distance(i,j,k+1), distance(i,j,k),dx);
          }
        }
      }
    }
}

void refined_levelset(Grid<double>& distance, std::vector<Particle>& particles, double dx, int ref){
  // initialize distance grid with maximum number
  int X_D = X_DIM*ref;
  int Y_D = Y_DIM*ref;
  int Z_D = Z_DIM*ref;
  dx /= ref;
  Grid<int> cell = Grid<int>(ref);

  int t,r,s;
  for(auto part : particles){
    find_cell(dx, part.pos.x, t);
    find_cell(dx, part.pos.y, r);
    find_cell(dx, part.pos.z, s);

    cell(t,r,s) = FLUID;
  }

  int max_d = (X_D + Y_D + Z_D) * dx;
  distance.set_constant(max_d);

  // distance from point to cell face
  int i, j, k;
  double xdist, ydist, zdist;

  for(int it = 0; it < particles.size(); ++it){
    find_cell(dx, particles[it].pos.x, xdist, i);
    find_cell(dx, particles[it].pos.y, ydist, j);
    find_cell(dx, particles[it].pos.z, zdist, k);

    double dist = std::sqrt(xdist*xdist + ydist*ydist + zdist*zdist);

    if(dist < distance(i,j,k)){
        distance(i,j,k) = -dist;
    }
}

    // all fluid cells get distance -0.5 negative, because its inside the fluid
    // for(int i = 1; i < X_DIM-1; i++){
    //   for (int j  = 1; j  < Y_DIM-1; j++) {
    //     for (size_t k = 1; k < Z_DIM-1; k++) {
    //       if(cell(i,j,k) == FLUID)
    //         // indicator inside fluid
    //         distance(i,j,k) = -0.5;
    //     }
    //   }
    // }

    // propagate closest point and distance etimate to the rest of the grid points
    for(int l = 0; l < 2; ++l){
      // 1. sweep order
      for(int i = 1; i < X_D; ++i){
        for(int j = 1; j < Y_D; ++j){
          for (int k = 1; k < Z_D; k++) {
            if(cell(i,j,k) != FLUID)
              solve_distance(distance(i-1,j,k), distance(i,j-1,k), distance(i,j,k-1), distance(i,j,k),dx);
          }
        }
      }

      // 2. sweep order
      for(int i = 1; i < X_D; ++i){
        for(int j = 1; j < Y_D; ++j){
          for (int k = Z_D-2; k >= 0; k--) {
            if(cell(i,j,k) != FLUID)
              solve_distance(distance(i-1,j,k), distance(i,j-1,k), distance(i,j,k+1), distance(i,j,k),dx);
          }
        }
      }

      // 3. sweep order
      for(int i = 1; i < X_D; ++i){
        for(int j = Y_D-2; j >= 0; --j){
          for (int k = 1; k < Z_D; k++) {
            if(cell(i,j,k) != FLUID)
              solve_distance(distance(i-1,j,k), distance(i,j+1,k), distance(i,j,k-1), distance(i,j,k),dx);
          }
        }
      }

      // 4. sweep order
      for(int i = X_D-2; i >= 0; --i){
        for(int j = 1; j < Y_D; ++j){
          for (int k = 1; k < Z_D; k++) {
            if(cell(i,j,k) != FLUID)
              solve_distance(distance(i+1,j,k), distance(i,j-1,k), distance(i,j,k-1), distance(i,j,k),dx);
          }
        }
      }

      // 5. sweep order
      for(int i = 1; i < X_D; ++i){
        for(int j = Y_D-2; j >= 0; --j){
          for (int k = Z_D-2; k >= 0; k--) {
            if(cell(i,j,k) != FLUID)
              solve_distance(distance(i-1,j,k), distance(i,j+1,k), distance(i,j,k+1), distance(i,j,k),dx);
          }
        }
      }

      // 6. sweep order
      for(int i = X_D-2; i >= 0; --i){
        for(int j = 1; j < Y_D; ++j){
          for (int k = Z_D-2; k >= 0; k--) {
            if(cell(i,j,k) != FLUID)
              solve_distance(distance(i+1,j,k), distance(i,j-1,k), distance(i,j,k+1), distance(i,j,k),dx);
          }
        }
      }

      // 7. sweep order
      for(int i = X_D-2; i >= 0; --i){
        for(int j = Y_D-2; j >= 0; --j){
          for (int k = 1; k < Z_D; k++) {
            if(cell(i,j,k) != FLUID)
              solve_distance(distance(i+1,j,k), distance(i,j+1,k), distance(i,j,k-1), distance(i,j,k),dx);
          }
        }
      }

      // 8. sweep order
      for(int i = X_D-2; i >= 0; --i){
        for(int j = Y_D-2; j >= 0; --j){
          for (int k = Z_D-2; k >= 0; --k) {
            if(cell(i,j,k) != FLUID)
              solve_distance(distance(i+1,j,k), distance(i,j+1,k), distance(i,j,k+1), distance(i,j,k),dx);
          }
        }
      }
    }
}

void Simulation::extend_velocity(){
  // expand to max 4 cells
  for(int it = 0; it < 4; ++it){

    extend_u(1,X_DIM-1,1,Y_DIM-1,1,Z_DIM-1);
    extend_u(X_DIM-2,0,1,Y_DIM-1,1,Z_DIM-1);
    extend_u(1,X_DIM-1,Y_DIM-2,0,1,Z_DIM-1);
    extend_u(1,X_DIM-1,1,Y_DIM-1,Z_DIM-2,0);
    extend_u(X_DIM-2,0,Y_DIM-2,0,1,Z_DIM-1);
    extend_u(X_DIM-2,0,1,Y_DIM-1,Z_DIM-2,0);
    extend_u(1,X_DIM-1,Y_DIM-2,0,Z_DIM-2,0);
    extend_u(X_DIM-2,0,Y_DIM-2,0,Z_DIM-2,0);

    // boundary conditions
    velx.clampBoundary();

    extend_v(1,X_DIM-1,1,Y_DIM-1,1,Z_DIM-1);
    extend_v(X_DIM-2,0,1,Y_DIM-1,1,Z_DIM-1);
    extend_v(1,X_DIM-1,Y_DIM-2,0,1,Z_DIM-1);
    extend_v(1,X_DIM-1,1,Y_DIM-1,Z_DIM-2,0);
    extend_v(X_DIM-2,0,Y_DIM-2,0,1,Z_DIM-1);
    extend_v(X_DIM-2,0,1,Y_DIM-1,Z_DIM-2,0);
    extend_v(1,X_DIM-1,Y_DIM-2,0,Z_DIM-2,0);
    extend_v(X_DIM-2,0,Y_DIM-2,0,Z_DIM-2,0);

    // boundary conditions
    vely.clampBoundary();

    extend_w(1,X_DIM-1,1,Y_DIM-1,1,Z_DIM-1);
    extend_w(X_DIM-2,0,1,Y_DIM-1,1,Z_DIM-1);
    extend_w(1,X_DIM-1,Y_DIM-2,0,1,Z_DIM-1);
    extend_w(1,X_DIM-1,1,Y_DIM-1,Z_DIM-2,0);
    extend_w(X_DIM-2,0,Y_DIM-2,0,1,Z_DIM-1);
    extend_w(X_DIM-2,0,1,Y_DIM-1,Z_DIM-2,0);
    extend_w(1,X_DIM-1,Y_DIM-2,0,Z_DIM-2,0);
    extend_w(X_DIM-2,0,Y_DIM-2,0,Z_DIM-2,0);

    // boundary conditions
    velz.clampBoundary();
  }
}

void Simulation::gridVelocityToParticles(){
    // get velocities from grid cells and interpolates to all particles
    //cell in which the particle lies
    int i,j,k;
    int ui,vj,wk;
    // distance from cell to the velocity which will be interpolated
    double fx, fy, fz;
    double ufx, vfy, wfz;

    for (auto &part : particles){

        find_cell(dx, part.pos.x, ufx, ui);
        find_cell_center(dx, part.pos.x, fx, i);
        find_cell(dx, part.pos.y, vfy, vj);
        find_cell_center(dx, part.pos.y, fy, j);
        find_cell(dx, part.pos.z, wfz, wk);
        find_cell_center(dx, part.pos.z, fz, k);

        // part.velocity.x += diff_velx.interpolate(ui,j,k,ufx,fy,fz,dx);
        // part.velocity.y += diff_vely.interpolate(i,vj,k,fx,vfy,fz,dx);
        // part.velocity.z += diff_velz.interpolate(i,j,wk,fx,fy,wfz,dx);


        //PIC to FLIP ration
        double alpha = 0.1;

        part.velocity.x = (part.velocity.x + diff_velx.interpolate(ui,j,k,ufx,fy,fz,dx))*(1-alpha) + alpha * velx.interpolate(ui,j,k,ufx,fy,fz,dx);
        part.velocity.y = (part.velocity.y + diff_vely.interpolate(i,vj,k,fx,vfy,fz,dx))*(1-alpha) + alpha * vely.interpolate(i,vj,k,fx,vfy,fz,dx);
        part.velocity.z = (part.velocity.z + diff_velz.interpolate(i,j,wk,fx,fy,wfz,dx))*(1-alpha) + alpha * velz.interpolate(i,j,wk,fx,fy,wfz,dx);
  }
}

void Simulation::run(int steps, int steps_per_frame){
    max_step = steps;
    write_frequency = steps_per_frame;
    std::cout << "Starting Simulation" << std::endl;
    std::cout << "Grid Size: " << X_DIM << "x" << Y_DIM << "x" << Z_DIM << std::endl;
    std::cout << "Number of Cells: " << X_DIM*Y_DIM*Z_DIM<< std::endl;
    std::cout << "Steps: " << max_step<< std::endl;
    char filename[40];

    tTotal.start();

    for(;current_step < max_step; ++current_step) {
        std::cout << "=========================================================================" << std::endl;
        std::cout << "Step: " << current_step << "of" << max_step << std::endl;
        std::cout << "#Particles:" << particles.size() << std::endl;

        step();
        tTotal.lap();

        std::cout << "step time: " << tStep.duration() << std::endl;
        std::cout << "advection time: " << tAdvection.duration() << std::endl;
        std::cout << "velocity extension time: " << textendVelocity.duration() << std::endl;
        std::cout << "pressure time: " << tPressure.duration() << std::endl;

        sprintf(filename,"data/test%03d.csv",current_step);
        if(current_step%write_frequency==0){
            write_to_csv(get_particles(),filename);
            write_frame();
        }
        std::cout << "=========================================================================" << std::endl;

    }

    std::cout << "Total simulation time: " << tTotal.duration() << std::endl;
    std::cout << "Average step time: " << tTotal.mean() << std::endl;
    std::cout << "Average advection time: " << tAdvection.duration() << std::endl;
    std::cout << "Average velocity extension time: " << textendVelocity.duration() << std::endl;
    std::cout << "Average pressure time: " << tPressure.duration() << std::endl;

}

void Simulation::step(){
    tStep.start();
    //// FLIP Algorithm
    //update particle positions with RK2 or simple advection
    tAdvection.start();
    simpleAdvect();
    applyBC();
    checkParticles();
    tAdvection.stop();

    // Splat particles velocities to grid to find intermediate velocities v*
    particleToGrid();
    // boundary conditions
    applyBC();

    // save velocities in pre grid

    pre_velx = velx;
    pre_vely = vely;
    pre_velz = velz;

    // resolve forces on grid
    addGravity();
    textendVelocity.start();
    levelset();
    marchingCubes();

    //creates refined level set grid
    Grid<double> ref_dist = Grid<double>(2);
    refined_levelset(ref_dist,particles,dx,2);
    refined_marchingCubes(ref_dist,dx,2,current_step);

    extend_velocity();
    textendVelocity.stop();
    applyBC();

    // algorithm to find the surface of the fluid
    marchingCubes();

    // solve pressures (Nabla)² p = Nabla v*
    tPressure.start();
    pressure_solver();
    tPressure.stop();

    extend_velocity();
    // use diff for flip
    diff_velx = velx - pre_velx;
    diff_vely = vely - pre_vely;
    diff_velz = velz - pre_velz;

    // update particles velocities
    gridVelocityToParticles();
    tStep.stop();

    // initialiser for new particles
    // source_cylinder(5,Vec3(X_DIM/2,Y_DIM/2,Z_DIM/2));
}

// another splat option which uses an SPH cubic spline kernel
void Simulation::cubicSplineToGrid(){
    int i,j,k;
    int ci,cj,ck;
    double du,dv,dw;
    double cdu,cdv,cdw;
    double xWeight,yWeight,zWeight;
    double xDist,yDist,zDist;

    sumx.zero();
    velx.zero();
    sumy.zero();
    vely.zero();
    sumz.zero();
    velz.zero();

    //get velocity from particles with a cublic spline kernel and add it to the grid.
    for(auto& part : particles) {
        if(IsInGrid(part)) {
            //find which cell the partice is in
            find_cell(dx, part.pos.x, du, i);
            find_cell_center(dx, part.pos.x, cdu, ci);
            find_cell(dx, part.pos.y, dv, j);
            find_cell_center(dx, part.pos.y, cdv, cj);
            find_cell(dx, part.pos.z, dw, k);
            find_cell_center(dx, part.pos.z, cdw, ck);


            //loop over 3x3x3 grid around particle
            for(int x = -1; x < 2; ++x) for(int y = -1; y < 2; ++y) for(int z = -1; z < 2; ++z) {

                xDist = std::sqrt(std::abs(x-du)*std::abs(x-du)+std::abs(y-cdv)*std::abs(y-cdv)+std::abs(z-cdw)*std::abs(z-cdw));
                yDist = std::sqrt(std::abs(x-cdu)*std::abs(x-cdu)+std::abs(y-dv)*std::abs(y-dv)+std::abs(z-cdw)*std::abs(z-cdw));
                zDist = std::sqrt(std::abs(x-cdu)*std::abs(x-cdu)+std::abs(y-cdv)*std::abs(y-cdv)+std::abs(z-dw)*std::abs(z-dw));

                xWeight = cubicSplineKernel(xDist);
                yWeight = cubicSplineKernel(yDist);
                zWeight = cubicSplineKernel(zDist);

                sumx(i+x,cj+y,ck+z) += xWeight;
                sumy(ci+x,j+y,ck+z) += yWeight;
                sumz(ci+x,cj+y,k+z) += zWeight;
                velx(i+x,cj+y,ck+z) += xWeight * part.velocity.x;
                vely(ci+x,j+y,ck+z) += yWeight * part.velocity.y;
                velz(ci+x,cj+y,k+z) += zWeight * part.velocity.z;
            }
        }
    }
    for(int i = 0; i < X_DIM; ++i){
      for(int j = 0; j < Y_DIM; ++j){
        for(int k = 0; k < Z_DIM; ++k){
          if(sumx(i,j,k) != 0){
            velx(i,j,k) /= sumx(i,j,k);
            vely(i,j,k) /= sumy(i,j,k);
            velz(i,j,k) /= sumz(i,j,k);
          }
        }
      }
    }
}

// another advection option
void Simulation::updateParticlePositions(){
    //Uses RK2 time integration to update the particle posiotions.

    //loop over all particles in the set
    for(auto &part : particles){
        double x_temp, y_temp, z_temp;
        double u_temp, v_temp, w_temp;
        double u, v, w;
        double x_min = 1.001 *dx;
        double x_max = (X_DIM-1.001) *dx;
        double y_min = 1.001 *dx;
        double y_max = (Y_DIM-1.001) *dx;
        double z_min = 1.001 *dx;
        double z_max = (Z_DIM-1.001) *dx;

        // u*
        //gets estimated interpolated velocities at (x,y,z)
        u_temp = velx.interpolate_vel(dx, part.pos.x, part.pos.y, part.pos.z,0);
        v_temp = vely.interpolate_vel(dx, part.pos.x, part.pos.y, part.pos.z,1);
        w_temp = velz.interpolate_vel(dx, part.pos.x, part.pos.y, part.pos.z,2);

        // x*
        x_temp = clamp(part.pos.x + (dt * u_temp) / 2.,x_min,x_max);
        y_temp = clamp(part.pos.y + (dt * v_temp) / 2.,y_min,y_max);
        z_temp = clamp(part.pos.z + (dt * w_temp) / 2.,z_min,z_max);

        //u
        //gets interpolated velocities at (x*,y*,z*)
        u = velx.interpolate_vel(dx, x_temp,y_temp, z_temp,0);
        v = vely.interpolate_vel(dx, x_temp,y_temp, z_temp,1);
        w = velz.interpolate_vel(dx, x_temp,y_temp, z_temp,2);

        //x
        part.pos.x += dt * u;
        part.pos.y += dt * v;
        part.pos.z += dt * w;
    }

}

void Simulation::DirichletBoundary() {
  for(unsigned int i = 0; i < X_DIM; i++){
    for(unsigned int j = 0; j < Y_DIM; j++){
      for(unsigned int k = 0; k < Z_DIM; k++){

        //make top and bottom of simulation solid
        if(k == 0 || k == Z_DIM-1){
            cell(i,j,k) = SOLID;
        }
        // X velocity
        if(	(cell.cval(i-1,j,k) == SOLID && velx(i,j,k) < 0) ||
          (cell(i, j, k) == SOLID && velx(i,j,k) > 0) )
          velx(i, j, k) = 0;

        // Y velocity
        if(	(cell.cval(i,j-1,k) == SOLID && vely(i,j,k) < 0) ||
          (cell(i, j, k) == SOLID && vely(i,j,k) > 0) )
          vely(i, j, k) = 0;

        // Z velocity
        if(	(cell.cval(i,j,k-1) == SOLID && velz(i,j,k) < 0) ||
          (cell(i, j, k) == SOLID && velz(i,j,k) > 0) )
          velz(i, j, k)  = 0;
      }
    }
  }
}

// -------------------------helper functions-------------------------

// velocity extention in x direction
void Simulation::extend_u(int x_start, int x_end, int y_start, int y_end, int z_start, int z_end){
    // if x,y,z is decreasing in loop --> -1, for increasing +1
    int di = -1;
    int dj = -1;
    int dk = -1;

    if(x_start<x_end)
        di = 1;
    if(y_start<y_end)
        dj = 1;
    if(z_start<z_end)
        dk = 1;

  double da, db, dc, alpha, beta;
  // itaretes through all possible directions
  for(int i = x_start; i != x_end; i += di){
    for(int j = y_start; j!= y_end; j += dj){
      for(int k = z_start; k!= z_end; k += dk){
        // calculate new extended velocities for air cells
        if(cell(i,j,k) == AIR && cell.cval(i-1,j,k) == AIR){
          // compute distance difference in x direction
          da = di * (distance(i,j,k) - distance(i-1,j,k));
          if(da < 0)
            continue; // this direction is unnessecairy
          // compute distance difference in y direction
          db = 0.5 * (distance(i,j,k) + distance(i-1,j,k)) - 0.5 * (distance(i,j-dj,k) + distance(i-1,j-dj,k));
          if(db < 0)
            continue;

          // compute distance difference in y direction
          dc = 0.5 * (distance(i,j,k) + distance(i-1,j,k)) - 0.5 * (distance(i,j,k-dk) + distance(i-1,j,k-dk));
          if(dc < 0)
            continue;
          // so dividing by 0 wont  be a problem
          if(da+db+dc == 0){
            alpha = 1./3.;
            beta = 1./3.;
          }
          else {
            alpha = da / (da+db+dc);
            beta = db / (da+db+dc);
          }

          velx(i,j,k) = alpha * velx(i-di,j,k) + beta * velx(i,j-dj,k) + (1-alpha-beta) * velx(i,j,k-dk);
        }
      }
    }
  }
}

// velocity extention in y direction
void Simulation::extend_v(int x_start, int x_end, int y_start, int y_end, int z_start, int z_end){
    int di = -1;
    int dj = -1;
    int dk = -1;

    if(x_start<x_end)
        di = 1;
    if(y_start<y_end)
        dj = 1;
    if(z_start<z_end)
        dk = 1;

  double da, db, dc, alpha, beta;

for(int k = z_start; k!= z_end; k += dk){
  for(int j = y_start; j!= y_end; j += dj){
     for(int i = x_start; i != x_end; i += di){
        // calculate new extended velocities for air cells
        if(cell(i,j,k) == AIR && cell.cval(i,j-1,k) == AIR){
          // compute distance difference in y direction
          da = dj * (distance(i,j,k) - distance(i,j-1,k));
          if(da < 0)
            continue; // this direction is unnessecairy
          // compute distance difference in x direction
          db = 0.5 * (distance(i,j,k) + distance(i,j-1,k)) - 0.5*(distance(i-di,j,k) + distance(i-di,j-1,k));
          if(db < 0)
            continue;

          // compute distance difference in z direction
          dc = 0.5 * (distance(i,j,k) + distance(i,j-1,k)) - 0.5*(distance(i,j,k-dk) + distance(i,j-1,k-dk));
          if(dc < 0)
            continue;

          if(da+db+dc == 0){
            alpha = 1./3.;
            beta = 1./3.;
          }
          else {
            alpha = da / (da+db+dc);
            beta = db / (da+db+dc);
          }
          // TODO: eventuell brauchts hier noch eine if bedingung
          vely(i,j,k) = alpha * vely(i,j-dj,k) + beta * vely(i-di,j,k) + (1-alpha-beta) * vely(i,j,k-dk);
        }
      }
    }
  }
}

// velocity extention in z direction
void Simulation::extend_w(int x_start, int x_end, int y_start, int y_end, int z_start, int z_end){
    int di = -1;
    int dj = -1;
    int dk = -1;

    if(x_start<x_end)
        di = 1;
    if(y_start<y_end)
        dj = 1;
    if(z_start<z_end)
        dk = 1;

    double da, db, dc, alpha, beta;

    for(int i = z_start; i != x_end; i += di){
      for(int j = y_start; j!= y_end; j += dj){
          for(int k = x_start; k!= z_end; k += dk){
          // calculate new extended velocities for air cells
          if(cell(i,j,k) == AIR && cell.cval(i,j,k-1) == AIR){
            // compute distance difference in z direction
            da = dk * (distance(i,j,k) - distance(i,j,k-1));
            if(da < 0)
              continue; // this direction is unnessecairy
            // compute distance difference in x direction
            db = 0.5 * (distance(i,j,k) + distance(i,j,k-1)) - 0.5*(distance(i-di,j,k) + distance(i-di,j,k-1));
            if(db < 0)
              continue;

            // compute distance difference in y direction
            dc = 0.5 * (distance(i,j,k) + distance(i,j,k-1)) - 0.5*(distance(i,j-dj,k) + distance(i,j-dj,k-1));
            if(dc < 0)
              continue;

            if(da+db+dc == 0){
              alpha = 1./3.;
              beta = 1./3.;
            }
            else {
              alpha = da / (da+db+dc);
              beta = db / (da+db+dc);
            }
            // TODO: eventuell brauchts hier noch eine if bedingung
            velz(i,j,k) = alpha * velz(i,j,k-dk) + beta * velz(i-di,j,k) + (1-alpha-beta) * velz(i,j-dj,k);
          }
        }
      }
    }
  }

void solve_distance(double fx, double fy, double fz, double& current, double dx){
  double phi_min;
  double phi_max;
  double phi_mid;

  phi_min = std::min({fx, fy, fz});
  phi_max = std::max({fx, fy, fz});
  phi_mid = std::max(std::min(fx,fy),std::min(std::max(fx,fy),fz));

  // try just the closest neighbor
  double d = phi_min + dx;

  if(d > phi_mid){
    // try the two closest neighbors
    d = 0.5 * (phi_min + phi_mid + std::sqrt(2*dx*dx - (phi_mid - phi_min)*(phi_mid - phi_min)));

    if(d > phi_max)
      d = (phi_min + phi_mid + phi_max + std::sqrt(std::max(0., std::pow(phi_min + phi_mid + phi_max, 2) - 3*(phi_min*phi_min + phi_mid*phi_mid + phi_max*phi_max - dx*dx)))) / 3.;
  }

  if (d < current)
    current = d;
}

// calculate the the velocity of each edge of one cell
template<class T>
void Simulation::velAtEdges(Grid<T> &v_grid, double v, int i, int j, int k, double du, double dv, double dw){
  float weight;

  du = du / dx;
  dv = dv / dx;
  dw = dw / dx;

  weight=(1-du)*(1-dv)*(1-dw);
  v_grid(i,j,k)+=weight*v;
  sum(i,j,k)+=weight;

  weight=du*(1-dv)*(1-dw);
  v_grid(i+1,j,k)+=weight*v;
  sum(i+1,j,k)+=weight;

  weight=(1-du)*dv*(1-dw);
  v_grid(i,j+1,k)+=weight*v;
  sum(i,j+1,k)+=weight;

  weight=(1-du)*(1-dv)*dw;
  v_grid(i,j,k+1)+=weight*v;
  sum(i,j,k+1)+=weight;

  weight=(1-du)*dv*dw;
  v_grid(i,j+1,k+1)+=weight*v;
  sum(i,j+1,k+1)+=weight;

  weight=du*(1-dv)*dw;
  v_grid(i+1,j,k+1)+=weight*v;
  sum(i+1,j,k+1)+=weight;

  weight=du*dv*(1-dw);
  v_grid(i+1,j+1,k)+=weight*v;
  sum(i+1,j+1,k)+=weight;

  weight=du*dv*dw;
  v_grid(i+1,j+1,k+1)+=weight*v;
  sum(i+1,j+1,k+1)+=weight;
}

//checks that the particle is in the grid
bool Simulation::IsInGrid(Particle& part) {
    if(part.pos.x < 0 || part.pos.y < 0 || part.pos.z < 0) {
        return false;
    }
    if(part.pos.x > X_DIM*dx || part.pos.y > Y_DIM * dx || part.pos.z > Z_DIM * dx) {
        return false;
    }
    if(std::isnan(part.pos.x) || std::isnan(part.pos.x) || std::isnan(part.pos.x)){
        return false;
    }
    return true;
}

/*
    FLIP FLuid Simulation: PBS Project
    Created by Tom Lausberg & Safira Piasko
    19.12.2018
*/

#include "simulation.hpp"
#include "util.hpp"
#include "grid.hpp"
#include <iostream>
#include <cstdio>

// index function
int idx(int i, int j, int k){
  return (k-1)*(Z_DIM-2)*(Y_DIM-2) + (j-1)*(Y_DIM-2) + i - 1;
}

void Simulation::pressure_solver(){
    // laplacian matrix
    A.resize((X_DIM-2)*(Y_DIM-2)*(Z_DIM-2),(X_DIM-2)*(Y_DIM-2)*(Z_DIM-2));
    A.setZero();
    A.reserve(particles.size()); //reserve space for nnz

    // divergence vector
    Eigen::VectorXd rhs((X_DIM-2)*(Y_DIM-2)*(Z_DIM-2));

    // calculate discrete laplacian
    for(int k = 1; k < Z_DIM - 1; ++k) {
      for(int j = 1; j < Y_DIM - 1; ++j) {
        for(int i = 1; i < X_DIM - 1; ++i){
            // Create function to get the poisson matrix cell index
            // Add entries on the pressure matrix by quering xyz indices to the poisson row calculation

            if(cell(i,j,k) == FLUID){
                int non_solid_neighbor_count = 0;

                // X direction
                if(cell(i-1,j,k) != SOLID) {
                    if(cell(i-1,j,k) == FLUID){
                        A.insert(idx(i,j,k),idx(i-1,j,k)) = -1;
                    }
                    non_solid_neighbor_count++;
                }
                if(cell(i+1,j,k) != SOLID) {
                    if(cell(i+1,j,k) == FLUID){
                        A.insert(idx(i,j,k),idx(i+1,j,k)) = -1;
                    }
                    non_solid_neighbor_count++;
                }

                // Y direction
                if(cell(i,j-1,k) != SOLID) {
                    if(cell(i,j-1,k) == FLUID){
                        A.insert(idx(i,j,k),idx(i,j-1,k)) = -1;
                    }
                    non_solid_neighbor_count++;
                }

                if(cell(i,j+1,k) != SOLID) {
                    if(cell(i,j+1,k) == FLUID){
                        A.insert(idx(i,j,k),idx(i,j+1,k)) = -1;
                    }
                    non_solid_neighbor_count++;
                }

                // Z direction
                if(cell(i,j,k-1) != SOLID) {
                    if(cell(i,j,k-1) == FLUID){
                        A.insert(idx(i,j,k),idx(i,j,k-1)) = -1;
                    }
                    non_solid_neighbor_count++;
                }
                if(cell(i,j,k+1) != SOLID) {
                    if(cell(i,j,k+1) == FLUID){
                        A.insert(idx(i,j,k),idx(i,j,k+1)) = -1;
                    }
                    non_solid_neighbor_count++;
                }

                //set diagonal
                A.insert(idx(i,j,k),idx(i,j,k)) = non_solid_neighbor_count;
            }

            else if(cell(i,j,k) == AIR){
                A.insert(idx(i,j,k),idx(i,j,k)) = 1;
            }

            if(cell(i,j,k)==FLUID){
                // store div in rhs
                rhs(idx(i,j,k)) = -(velx(i + 1,j,k) - velx(i,j,k)
                                  + vely(i,j + 1,k) - vely(i,j,k)
                                  + velz(i,j,k + 1) - velz(i,j,k))*(dx)/(dt);
            }
            else {
                // store div in rhs
                rhs(idx(i,j,k)) = 0;
            }
        }
      }
    }

    //vector for solution of p
    Eigen::initParallel();
    int n = Eigen::nbThreads();
    std::cout << "Number of threads = " << n << std::endl;
    Eigen::VectorXd p_sol((X_DIM-2)*(Y_DIM-2)*(Z_DIM-2));
    cg_solver.compute(A);
    p_sol = cg_solver.solve(rhs);

    // Loop through all cells to set new velocities form pressures
    for(int k = 1; k < Z_DIM; ++k){
      for(int j = 1; j < Y_DIM; ++j) {
        for(int i = 1; i < X_DIM; ++i) {

          if(cell(i,j,k)==FLUID){
              pressure(i,j,k) = p_sol(idx(i,j,k)); // used for writer vars
          }

          if(cell(i,j,k) == SOLID || cell(i-1,j,k) == SOLID){
            velx(i,j,k) = 0;
          }
          else if(cell(i,j,k) == FLUID || cell(i-1,j,k) == FLUID || cell(i,j-1,k) == FLUID || cell(i,j,k-1) == FLUID){
            double p_diff_x = p_sol(idx(i,j,k)) - p_sol(idx(i-1,j,k));
            velx(i,j,k) -= dt * p_diff_x / dx;
          }

          if(cell(i,j,k) == SOLID || cell(i,j-1,k) == SOLID){
            vely(i,j,k) = 0;
          }
          else if(cell(i,j,k) == FLUID || cell(i-1,j,k) == FLUID || cell(i,j-1,k) == FLUID || cell(i,j,k-1) == FLUID){
            double p_diff_y = p_sol(idx(i,j,k)) - p_sol(idx(i,j-1,k));
            vely(i,j,k) -= dt * p_diff_y / dx;
          }

          if(cell(i,j,k) == SOLID || cell(i,j,k-1) == SOLID){
            velz(i,j,k) = 0;
          }
          else if(cell(i,j,k) == FLUID || cell(i-1,j,k) == FLUID || cell(i,j-1,k) == FLUID || cell(i,j,k-1) == FLUID){
            double p_diff_z = p_sol(idx(i,j,k)) - p_sol(idx(i,j,k-1));
            velz(i,j,k) -= dt * p_diff_z / dx;
          }
        }
      }
    }
}

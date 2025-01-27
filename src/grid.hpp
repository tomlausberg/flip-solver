/*
    FLIP FLuid Simulation: PBS Project
    Created by Tom Lausberg & Safira Piasko
    19.12.2018
*/

#ifndef GRID_H
#define GRID_H

#include <eigen3/Eigen/Dense>
#include <eigen3/unsupported/Eigen/CXX11/Tensor>
#include "util.hpp"
#include <iostream>
#include <cassert>


template<class T>
class Grid {
private:
    typedef T value_type;
    typedef T& ref_type;
    typedef const T& const_ref_type;
    Eigen::Tensor<T, 3> data;

public:
    // constructors
    Grid(void){
        data = Eigen::Tensor<T,3>(X_DIM,Y_DIM,Z_DIM);
        data.setZero();
    }
    Grid(int i){
        data = Eigen::Tensor<T,3>(i*X_DIM,i*Y_DIM,i*Z_DIM);
        data.setZero();
    }

    // copy constructors
    Grid(const Grid<T> &copy) : data(copy.data) {};
    Grid(const Eigen::Tensor<T, 3> &data_copy) : data(data_copy) {};

    // set all entries in grid to value
    void set_constant(const T& value){
        data.setConstant(value);
    }

    // returns rank of dimension d
    int dim(int d){
        return data.dimension(d);
    }

    // Boundary get the same values as their inner neighbors
    void clampBoundary(){

      for(int i = 0; i < X_DIM; i++){
        for (int j = 0; j < Y_DIM; j++) {
          for(int k = 0; k < Z_DIM; k++){
            if(i==0)
              data(i,j,k) = data(i+1,j,k);

            if(i==X_DIM-1)
              data(i,j,k) = data(i-1,j,k);

            if(j==0)
              data(i,j,k) = data(i,j+1,k);

            if(j==Y_DIM-1)
              data(i,j,k) = data(i,j-1,k);

            if(k==0)
              data(i,j,k) = data(i,j,k+1);

            if(k==Z_DIM-1)
              data(i,j,k) = data(i,j,k-1);
          }
        }
      }
      // clamping edges of the grid
      data(0,0,0) = data(1,1,1);
      data(X_DIM-1,Y_DIM-1,Z_DIM-1) = data(X_DIM-2,X_DIM-2,X_DIM-2);
      data(0,0,Z_DIM-1) = data(1,1,X_DIM-2);
      data(X_DIM-1,0,0) = data(X_DIM-2,1,1);
      data(0,Y_DIM-1,0) = data(1,X_DIM-2,1);
      data(0,Y_DIM-1,Z_DIM-1) = data(1,X_DIM-2,X_DIM-2);
      data(X_DIM-1,0,Z_DIM-1) = data(X_DIM-2,1,X_DIM-2);
      data(X_DIM-1,Y_DIM-1,0) = data(X_DIM-2,X_DIM-2,1);
    }

    ref_type operator()(int i, int j, int k) {return data(i,j,k);}
    const_ref_type operator()(int i, int j, int k) const {return data(i,j,k);}

    // clamping when using an entry of a grid
    T cval(int i,int j,int k){
        return data(clamp(i,0,X_DIM),clamp(j,0,Y_DIM),clamp(k,0,Z_DIM));
    }

    void tensor_to_array(float array[]){
        for(int k = 0; k < Z_DIM; ++k){
            for(int j = 0; j < Y_DIM; ++j) {
                for(int i = 0;i < X_DIM; ++i) {
                    array[i+j*Y_DIM+k*Z_DIM*Z_DIM] = float(data(i,j,k));
                }
            }
        }
    }

    void tensor_to_vec_array(float array[], int n){
        for(int k = 0; k < Z_DIM; ++k){
            for(int j = 0; j < Y_DIM; ++j) {
                for(int i = 0;i < 3 * X_DIM; i+=3) {
                    array[n +i*3+ j*Y_DIM*3+ k*Z_DIM*Y_DIM*3] = float(data(i,j,k));
                }
            }
        }
    }

    // interpolation from grid to particle
    double interpolate(int i, int j, int k, double gdx, double gdy, double gdz, double delta_x){

        double dx = gdx/delta_x;
        double dy = gdy/delta_x;
        double dz = gdz/delta_x;

        //Trilinearinterpolation
        double z1 = (1-dz) * cval(i,j,k) + dz * cval(i,j,k+1);             //0,0,0, + 0,0,1
        double z2 = (1-dz) * cval(i,j+1,k) + dz * cval(i,j+1,k+1);         //0,1,0  + 0,1,1

        double z3 = (1-dz) * cval(i+1,j,k) + dz * cval(i+1,j,k+1);         //1,0,0, + 1,0,1
        double z4 = (1-dz) * cval(i+1,j+1,k) + dz * cval(i+1,j+1,k+1);     //1,1,0, + 1,1,1

        double y1 = (1-dy) * z1 + dy * z2;
        double y2 = (1-dy) * z3 + dy * z4;

        return (1-dx) * y1 + dx * y2;
    }

    // interpolation used for Runge Kutta 2
    // finds interpolated velocity at x,y,z for the current grid.
    // only works for velocity grids with coorect dim! (for u: dim = 0, v: dim=1, w: dim=2)
    double interpolate_vel(double dx, double x, double y, double z, int dim){
        //cell in which the current coordinates lie
        int i,j,k;
        //distance in cell from velocity value
        double du,dv,dw;

        if(dim==0){
            find_cell(dx,x,du,i);
            find_cell_center(dx,y,dv,j);
            find_cell_center(dx,z,dw,k);
        }
        else if(dim==1){
            find_cell_center(dx,x,du,i);
            find_cell(dx,y,dv,j);
            find_cell_center(dx,z,dw,k);
        }
        else if(dim==2){
            find_cell_center(dx,x,du,i);
            find_cell_center(dx,y,dv,j);
            find_cell(dx,z,dw,k);
        }
        double vel = interpolate(i,j,k,du,dv,dw,dx); //du,dv or dx,dy
        return vel;
    }

    // sets the grid to zero
    void zero(){
        data.setZero();
    }

    // minus operator for grid
    Grid<T> operator-(const Grid<T>& b){
        assert(this->data.dimension(0) == b.data.dimension(0) &&
               this->data.dimension(1) == b.data.dimension(1) &&
               this->data.dimension(2) == b.data.dimension(2) &&
               "grid operator-: dimensions do not match");
        Grid<T> val(this->data - b.data);
        return val;
        //TODO optimise
    }

    // dividing operator for grid
    void divideBy(const Grid<T>& b){
        assert(this->data.dimension(0) == b.data.dimension(0) &&
               this->data.dimension(1) == b.data.dimension(1) &&
               this->data.dimension(2) == b.data.dimension(2) &&
               "grid operator/=: dimensions do not match");
        this->data = this->data / b.data;
    }

    // output operator for grid (for debugging)
    std::ostream& operator<<(std::ostream& os){
        for (int i = 0; i < X_DIM; i++) {
            for (size_t j = 0; j < Y_DIM; j++) {
                for (size_t k = 0; k < Z_DIM; k++) {
                    os << data(i,j,k) << " ";
                }
                os << "\n";
            }
            os << "\n";
        }
        return os;
    }
};


#endif

/*
    FLIP FLuid Simulation: PBS Project
    Created by Tom Lausberg & Safira Piasko
    19.12.2018
*/

#include <vector>
#include <iostream>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <string>
#include <typeinfo>
#include <cstdio>

#include "initialiser.hpp"
// #include "simulation.hpp"
#include "util.hpp"
#include "grid.hpp"
#include "writer.hpp"
#include "timer.hpp"
#include "makePLY.hpp"
#include "MCtable.hpp"

// find point with value zero between two points
Vec3 linInterpolation(Vec3 place0, Vec3 place1, double d0, double d1){
  Vec3 point;

  if(d0 != d1){
    point.x = place0.x + (place1.x -place0.x)/(d1 - d0)*(-d0);
    point.y = place0.y + (place1.y -place0.y)/(d1 - d0)*(-d0);
    point.z = place0.z + (place1.z -place0.z)/(d1 - d0)*(-d0);

    return point;
  }
  else
    return place0;
}


// give all points an index number
void get_verticies(std::vector<Vec3>& points, std::vector<Triangle>& triangles){
    for(auto& tri : triangles){
      int i = 0;
      for(auto point : points){
        // check if one of the points in the current triangle already exists in the point vector
        if(tri.p[0] == point){
          tri.v1 = i;
        }
        else if(tri.p[1] == point){
          tri.v2 = i;
        }
        else if(tri.p[2] == point){
          tri.v3 = i;
        }
        ++i;
      }

      // new indices for new points
      if(tri.v1==-1){
        points.push_back(tri.p[0]);
        tri.v1 = i;
        ++i;
      }
      if(tri.v2==-1){
        points.push_back(tri.p[1]);
        tri.v2 = i;
        ++i;
      }
      if(tri.v3==-1){
        points.push_back(tri.p[2]);
        tri.v3 = i;
        ++i;
    }
  }

}

// write koordinates of the points and triangles int ply file
void write_ply(int frame_no, std::vector<Vec3>& points, std::vector<Triangle>& triangles, int ref){
    char name[40];
    if(ref){
        sprintf(name, "ply_ref/surface_%03d.ply",frame_no);
    } else {
        sprintf(name, "ply/surface_%03d.ply",frame_no);
    }
    std::string fileName(name);
    std::ofstream plyFile(fileName.c_str());

    //write header;
    plyFile << "ply" << std::endl;
    plyFile << "format ascii 1.0" << std::endl;
    plyFile << "element vertex " << points.size() << std::endl;
    plyFile << "property float x" << std::endl;
    plyFile << "property float y" << std::endl;
    plyFile << "property float z" << std::endl;
    plyFile << "element face " << triangles.size() << std::endl;
    plyFile << "property list uchar int vertex_index" << std::endl;
    plyFile << "end_header" << std::endl;

    //write points
    for ( auto point : points ) {
        plyFile << point.x << " "
                << point.y << " "
                << point.z << std::endl;
    }

    //write triangles
    for (auto triangle : triangles){
        plyFile << "3 "
                << triangle.v1 << " "
                << triangle.v2 << " "
                << triangle.v3 << std::endl;
    }

    return;
}

// includes the algorithm marching cubes and formats the triangles into a ply format
void Simulation::marchingCubes(){
  std::vector<Triangle> triangles;
  std::vector<Vec3> glob_zero_points;

  // go through all the points and initialize the vertices
  for (size_t k = 0; k < Z_DIM-1; k++) {
    for (size_t j = 0; j < Y_DIM-1; j++) {
      for (size_t i = 0; i < X_DIM-1; i++) {
        std::vector<double> vert; // saves the distance value
        std::vector<Vec3> place; // saves where the nodes are
        vert.push_back(distance(i,j,k));
        vert.push_back(distance(i,j,k+1));
        vert.push_back(distance(i+1,j,k+1));
        vert.push_back(distance(i+1,j,k));
        vert.push_back(distance(i,j+1,k));
        vert.push_back(distance(i,j+1,k+1));
        vert.push_back(distance(i+1,j+1,k+1));
        vert.push_back(distance(i+1,j+1,k));

        place.push_back(Vec3(i*dx,j*dx,k*dx));
        place.push_back(Vec3(i*dx,j*dx,(k+1)*dx));
        place.push_back(Vec3((i+1)*dx,j*dx,(k+1)*dx));
        place.push_back(Vec3((i+1)*dx,j*dx,k*dx));
        place.push_back(Vec3(i*dx,(j+1)*dx,k*dx));
        place.push_back(Vec3(i*dx,(j+1)*dx,(k+1)*dx));
        place.push_back(Vec3((i+1)*dx,(j+1)*dx,(k+1)*dx));
        place.push_back(Vec3((i+1)*dx,(j+1)*dx,k*dx));

        // determine weather vertex is inside or outside
        // save in a number of bits 1 for inside 0 for outside
        std::vector<int> inside;
        for (size_t n = 0; n < 8; n++) {
          if(vert[n] <= 0)
            inside.push_back(1);
          else
            inside.push_back(0);
        }

        // check if cube is completly inside or outside with summing up the values of the vector
        int sum = 0;
        for (size_t i = 0; i < inside.size(); i++) {
            sum += inside[i];
        }
        if(sum == 0 || sum == 8){
            continue;
        }

        // find zero vertices  at edges through interpolation and save them in z_points vector
        std::vector<Vec3> z_points(12);
        // for edges, 0:(0,1),1:(1,2),2:(2,3),4:(4,5),5:(5,6),6:(6,7)
        for(int n = 0; n < 7; n++){
          if(inside[n] + inside[n+1] == 1)
            z_points[n] = linInterpolation(place[n], place[n+1], vert[n], vert[n+1]);
          else
            z_points[n] = Vec3(0,0,0);
          if(n == 2)
            n++;
        }

        // edges 3:(3,0) and 7:(7,4)
        if(inside[3]+ inside[0] == 1)
          z_points[3] = linInterpolation(place[3], place[0], vert[3], vert[0]);
        else z_points[3] = Vec3(0,0,0);
        if(inside[7]+ inside[4] == 1)
          z_points[7] = linInterpolation(place[7], place[4], vert[7], vert[4]);
        else z_points[7] = Vec3(0,0,0);

        // edges 8:(0,4), 9:(1,5), 10:(2,6), 11:(3,7)
        for (size_t n = 0; n < 4; n++) {
          if(inside[n] + inside[n+4] == 1)
            z_points[n+8] = linInterpolation(place[n], place[n+4], vert[n], vert[n+4]);
          else
            z_points[n+8] = Vec3(0,0,0);
        }

        // make "binary" vector inside to dezimal number
        int dezimal = 0;
        for (size_t n = 0; n < 8; n++) {
          if(inside[n] == 1)
            dezimal |= (1 << n);
        }

        // build the triangles using triTable
        for(int n = 0; triTable[dezimal][n]!= -1; n+=3){
          Triangle triangle;
          triangle.p[0] = z_points[triTable[dezimal][n+2]];
          triangle.p[1] = z_points[triTable[dezimal][n+1]];
          triangle.p[2] = z_points[triTable[dezimal][n]];

          triangles.push_back(triangle);
        }
      }
    }
  }
  //formating for ply
  get_verticies(glob_zero_points, triangles);
  write_ply(current_step, glob_zero_points, triangles, 0);
}
  void refined_marchingCubes(Grid<double>& distance, double dx, int ref, int current_step){
    std::vector<Triangle> triangles;
    std::vector<Vec3> glob_zero_points;

    int X_D = X_DIM*ref;
    int Y_D = Y_DIM*ref;
    int Z_D = Z_DIM*ref;
    dx /= ref;



    // go through all the points and initialize the vertices
    for (size_t k = 0; k < Z_D-1; k++) {
      for (size_t j = 0; j < Y_D-1; j++) {
        for (size_t i = 0; i < X_D-1; i++) {
          std::vector<double> vert; // saves the distance value
          std::vector<Vec3> place; // saves where the nodes are
          vert.push_back(distance(i,j,k));
          vert.push_back(distance(i,j,k+1));
          vert.push_back(distance(i+1,j,k+1));
          vert.push_back(distance(i+1,j,k));
          vert.push_back(distance(i,j+1,k));
          vert.push_back(distance(i,j+1,k+1));
          vert.push_back(distance(i+1,j+1,k+1));
          vert.push_back(distance(i+1,j+1,k));

          place.push_back(Vec3(i*dx,j*dx,k*dx));
          place.push_back(Vec3(i*dx,j*dx,(k+1)*dx));
          place.push_back(Vec3((i+1)*dx,j*dx,(k+1)*dx));
          place.push_back(Vec3((i+1)*dx,j*dx,k*dx));
          place.push_back(Vec3(i*dx,(j+1)*dx,k*dx));
          place.push_back(Vec3(i*dx,(j+1)*dx,(k+1)*dx));
          place.push_back(Vec3((i+1)*dx,(j+1)*dx,(k+1)*dx));
          place.push_back(Vec3((i+1)*dx,(j+1)*dx,k*dx));

          // determine weather vertex is inside or outside
          // save in a number of bits 1 for inside 0 for outside
          std::vector<int> inside;
          for (size_t n = 0; n < 8; n++) {
            if(vert[n] <= 0)
              inside.push_back(1);
            else
              inside.push_back(0);
          }

          // check if cube is completly inside or outside with summing up the values of the vector
          int sum = 0;
          for (size_t i = 0; i < inside.size(); i++) {
              sum += inside[i];
          }
          if(sum == 0 || sum == 8){
              continue;
          }

          // find zero vertices  at edges through interpolation and save them in z_points vector
          std::vector<Vec3> z_points(12);
          // for edges, 0:(0,1),1:(1,2),2:(2,3),4:(4,5),5:(5,6),6:(6,7)
          for(int n = 0; n < 7; n++){
            if(inside[n] + inside[n+1] == 1)
              z_points[n] = linInterpolation(place[n], place[n+1], vert[n], vert[n+1]);
            else
              z_points[n] = Vec3(0,0,0);
            if(n == 2)
              n++;
          }

          // edges 3:(3,0) and 7:(7,4)
          if(inside[3]+ inside[0] == 1)
            z_points[3] = linInterpolation(place[3], place[0], vert[3], vert[0]);
          else z_points[3] = Vec3(0,0,0);
          if(inside[7]+ inside[4] == 1)
            z_points[7] = linInterpolation(place[7], place[4], vert[7], vert[4]);
          else z_points[7] = Vec3(0,0,0);

          // edges 8:(0,4), 9:(1,5), 10:(2,6), 11:(3,7)
          for (size_t n = 0; n < 4; n++) {
            if(inside[n] + inside[n+4] == 1)
              z_points[n+8] = linInterpolation(place[n], place[n+4], vert[n], vert[n+4]);
            else
              z_points[n+8] = Vec3(0,0,0);
          }

          // make "binary" vector inside to dezimal number
          int dezimal = 0;
          for (size_t n = 0; n < 8; n++) {
            if(inside[n] == 1)
              dezimal |= (1 << n);
          }

          // build the triangles using triTable
          for(int n = 0; triTable[dezimal][n]!= -1; n+=3){
            Triangle triangle;
            triangle.p[0] = z_points[triTable[dezimal][n+2]];
            triangle.p[1] = z_points[triTable[dezimal][n+1]];
            triangle.p[2] = z_points[triTable[dezimal][n]];

            triangles.push_back(triangle);
          }
        }
      }
    }

  //formating for ply
  get_verticies(glob_zero_points, triangles);
  write_ply(current_step, glob_zero_points, triangles, ref);

}

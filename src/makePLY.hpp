/*
    FLIP FLuid Simulation: PBS Project
    Created by Tom Lausberg & Safira Piasko
    19.12.2018
*/

#ifndef MAKEPLY_H_
#define MAKEPLY_H_

#include "grid.hpp"
#include "simulation.hpp"

struct Triangle{
  // coordinates
  Vec3 p[3];
  int v1;
  int v2;
  int v3;

  Triangle(){
      v1 = -1; //initialise with invalid vertex in point vector
      v2 = -1;
      v3 = -1;
  }
};

void write_ply(int frame_no, std::vector<Vec3>& points, std::vector<Triangle>& triangles, int ref);
void get_verticies(std::vector<Vec3>& points, std::vector<Triangle>& triangles);
void refined_marchingCubes(Grid<double>& distance, double dx, int ref, int current_step);

#endif

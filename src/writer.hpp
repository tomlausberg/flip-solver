/*
    FLIP FLuid Simulation: PBS Project
    Created by Tom Lausberg & Safira Piasko
    19.12.2018
*/

#ifndef VTK_H
#define VTK_H

#include <iostream>
#include <vector>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<cassert>

#include "util.hpp"
#include "simulation.hpp"

void write_to_csv( std::vector<Particle> particles, const char *filename );

FILE *vtk_open(char* filename, int no);

//Create vtk header
void vtk_header(FILE *data);

// Inserts point data into Vtk file
void vtk_grid_points(FILE *data, int dim,int x_length, int y_length, int z_length, double dx, Vec3 origin);

void vtk_poly_vortex_and_celltype(FILE *data, int x_length, int y_length, int z_length);

//Inserts point data into Vtk file
void vtk_pointdata(FILE *data,char *val_name, int x_length, int y_length, int z_length);

void write_xvel_vtk(int frame);



#endif

/*
    FLIP FLuid Simulation: PBS Project
    Created by Tom Lausberg & Safira Piasko
    19.12.2018
*/


#include <fstream>
#include <iostream>
#include <vector>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include <string>
#include<cassert>

#include "writer.hpp"
#include "visit_writer.hpp"
#include "simulation.hpp"
#include "util.hpp"

void write_to_csv( std::vector<Particle> particles, const char *filename ){
    std::ofstream pointCloudCsv(filename);
    pointCloudCsv << "x coord, y coord, z coord, scalar\n";

    for(auto part : particles){
         pointCloudCsv << part.pos.x << ", " << part.pos.y << ", " << part.pos.z << ", 0\n";
    }

    pointCloudCsv.close();
}


void Simulation::write_frame(){
    // writes different variables into seperate vtk files
    write_pressure();
    write_cells();
    return;
}

void Simulation::write_pressure(){
    float* pressureData = new float[X_DIM*Y_DIM*Z_DIM];
    pressure.tensor_to_array(pressureData);

    char pressureFile[40];
    sprintf(pressureFile,"data/pressure.vtk.%03d",current_step);
    int dims[] = {X_DIM,Y_DIM,Z_DIM};

    float x[X_DIM];
    float y[Y_DIM];
    float z[Z_DIM];
    int  i;
    for (i = 0 ; i < dims[0]; ++i) x[i] = i*dx;
    for (i = 0 ; i < dims[1]; ++i) y[i] = i*dx;
    for (i = 0 ; i < dims[2]; ++i) z[i] = i*dx;

    int nvar = 1;
    int vardim[]={1};
    int centering[]={1};
    const char *pressureVar[] = {"pressure"};
    float *vars[] = {(float*)pressureData};

    write_rectilinear_mesh(pressureFile,0,dims,x,y,z,nvar,vardim,centering,pressureVar,vars);
    delete(pressureData);
}

void Simulation::write_cells(){
    //create buffer for data to be written
    float* data = new float[X_DIM*Y_DIM*Z_DIM];
    //write data from Eigen::Tensor grid to buffer
    cell.tensor_to_array(data);

    char file[40];
    sprintf(file,"data/cells.vtk.%03d",current_step);
    int dims[] = {X_DIM,Y_DIM,Z_DIM};

    float x[X_DIM];
    float y[Y_DIM];
    float z[Z_DIM];
    int  i;
    for (i = 0 ; i < dims[0]; ++i) x[i] = i*dx;
    for (i = 0 ; i < dims[1]; ++i) y[i] = i*dx;
    for (i = 0 ; i < dims[2]; ++i) z[i] = i*dx;

    int nvar = 1; //number of varibles
    int vardim[]={1}; //dimension of varibles
    int centering[]={1}; //centering?
    const char *var[] = {"cell_type"}; //varible name
    float *vars[] = {(float*)data}; //create pointer to data

    write_rectilinear_mesh(file,0,dims,x,y,z,nvar,vardim,centering,var,vars); //write to vtk file
    delete(data); //free space
}

void Simulation::write_particle_count(){
    //create buffer for data to be written
    float* data = new float[X_DIM*Y_DIM*Z_DIM];
    //write data from Eigen::Tensor grid to buffer
    particle_count.tensor_to_array(data);

    char file[40];
    sprintf(file,"data/count.vtk.%03d",current_step);
    int dims[] = {X_DIM,Y_DIM,Z_DIM};

    float x[X_DIM];
    float y[Y_DIM];
    float z[Z_DIM];
    int  i;
    for (i = 0 ; i < dims[0]; ++i) x[i] = i*dx;
    for (i = 0 ; i < dims[1]; ++i) y[i] = i*dx;
    for (i = 0 ; i < dims[2]; ++i) z[i] = i*dx;

    int nvar = 1; //number of varibles
    int vardim[]={1}; //dimension of varibles
    int centering[]={0}; //centering?
    const char *var[] = {"particle count"}; //varible name
    float *vars[] = {(float*)data}; //create pointer to data

    write_rectilinear_mesh(file,0,dims,x,y,z,nvar,vardim,centering,var,vars); //write to vtk file
    delete(data); //free space
}

void Simulation::write_scalars(){
    //create buffer for data to be written
    float* dataPressure = new float[X_DIM*Y_DIM*Z_DIM];
    float* dataCell = new float[X_DIM*Y_DIM*Z_DIM];
    float* dataCount = new float[X_DIM*Y_DIM*Z_DIM];

    //write data from Eigen::Tensor grid to buffer
    pressure.tensor_to_array(dataPressure);
    cell.tensor_to_array(dataCell);
    particle_count.tensor_to_array(dataCount);

    char file[40];
    sprintf(file,"data/scalar.vtk.%03d",current_step);
    int dims[] = {X_DIM,Y_DIM,Z_DIM};

    float x[X_DIM];
    float y[Y_DIM];
    float z[Z_DIM];

    //create spacings for grid
    int  i;
    for (i = 0 ; i < dims[0]; ++i) x[i] = i*dx;
    for (i = 0 ; i < dims[1]; ++i) y[i] = i*dx;
    for (i = 0 ; i < dims[2]; ++i) z[i] = i*dx;

    int nvar = 3; //number of varibles
    int vardim[]={1,1,1}; //dimension of varibles
    int centering[]={1,1,1}; //1: value stored at cell center, 0: value stored at node
    const char *var[] = {"pressure", "cell_type", "particle count"}; //varible name
    float *vars[] = {(float*)dataPressure,(float*)dataCell,(float*)dataCount}; //create pointer to data

    //write to vtk file
    write_rectilinear_mesh(file,0,dims,x,y,z,nvar,vardim,centering,var,vars);
    delete(dataPressure); //free space
    delete(dataCell); //free space
    delete(dataCount); //free space
}

void Simulation::write_velocity(){
    //create buffer for data to be written
    float* data = new float[3*X_DIM*Y_DIM*Z_DIM];


    //write data from Eigen::Tensor grid to buffer
    velx.tensor_to_vec_array(data,0);
    vely.tensor_to_vec_array(data,1);
    velz.tensor_to_vec_array(data,2);

    char file[40];
    sprintf(file,"data/velocity.vtk.%03d",current_step);
    int dims[] = {X_DIM,Y_DIM,Z_DIM};

    float x[X_DIM];
    float y[Y_DIM];
    float z[Z_DIM];

    //create spacings for grid
    int  i;
    for (i = 0 ; i < dims[0]; ++i) x[i] = i*dx;
    for (i = 0 ; i < dims[1]; ++i) y[i] = i*dx;
    for (i = 0 ; i < dims[2]; ++i) z[i] = i*dx;

    int nvar = 1; //number of varibles
    int vardim[]={3}; //dimension of varibles
    int centering[]={1}; //1: value stored at cell center, 0: value stored at node
    const char *var[] = {"velocity"}; //varible name
    float *vars[] = {(float*)data}; //create pointer to data
    //TODO fix vars definition

    //write to vtk file
    write_rectilinear_mesh(file,0,dims,x,y,z,nvar,vardim,centering,var,vars);
    delete(data); //free space
}

//Opens data stream to vtk timestamp
FILE *vtk_open(const char* filename, int no) {
	//char name1[]= var.filename; //creates name
	char name_t[40]; //INIT total name
	FILE *data;
	snprintf(name_t,sizeof(name_t),"%s%i",filename,no); //combines name and no.1 Bsp. "2d_sim.vtk.13"
	data=fopen(name_t,"w"); //opens/creates writeable file with name
	return data;
}

//Create vtk header
void vtk_header(FILE *data) {
    fprintf(data,"# vtk DataFile Version 3.1\n FLIP Simulation\nASCII\nDATASET UNSTRUCTURED_GRID\n");
    return;
}

// Inserts point data into Vtk file
void vtk_grid_points(FILE *data, int dim,int x_length, int y_length, int z_length, double dx, Vec3 origin) {
    //dim 0: cell center(0.5, 0.5, 0.5)
    //    1: xvel point (0  , 0.5, 0.5)
    //    2: yvel point (0.5, 0  , 0.5)
    //    3: zvel point (0.5, 0.5, 0  )

    fprintf(data,"POINTS %i DOUBLE\n", x_length*y_length*z_length); //Point Header, number of points, type
    int i,j,k;
    double x_shift,y_shift,z_shift;
    if(dim == 0){
        x_shift = 0.5;
        y_shift = 0.5;
        z_shift = 0.5;
    } else if(dim == 1){
        x_shift = 0;
        y_shift = 0.5;
        z_shift = 0.5;
    } else if(dim == 2){
        x_shift = 0.5;
        y_shift = 0;
        z_shift = 0.5;
    } else if(dim == 3){
        x_shift = 0.5;
        y_shift = 0.5;
        z_shift = 0;
    } else {
        assert("wrong dimension centering");
    }


    for(i=0;i<x_length;i++) {
        for(j=0;j<y_length;j++) {
            for(k=0;k<z_length;k++){
                fprintf(data,"%f %f %f \n", origin.x + (i+x_shift) * dx,
                                            origin.y + (j+y_shift) * dx,
                                            origin.z + (k+z_shift) * dx); // print to grid
            }
        }
    }
    fprintf(data,"\n");
    return;
}

void vtk_poly_vortex_and_celltype(FILE *data, int x_length, int y_length, int z_length) {

    fprintf(data,"CELLS %i %i\n", 1, 1+x_length*y_length*z_length); //Cell Header, number of Cells, total data( no. of rows of triangle strips*no. of points per strip + number of cells)
    int i;
    fprintf(data,"%i",x_length*y_length*z_length); //number of data points in row
    for(i=0;i<x_length*y_length*z_length;i++) {
        fprintf(data," %i",i); // adds all point to poly_vortex
    }
    fprintf(data,"\n\nCELL_TYPES 1\n2\n"); //def CELL_TYPES as 1 poly_vortex (=2)
    return;
}

//Inserts point data into Vtk file
void vtk_pointdata(FILE *data,char* val_name, int x_length, int y_length, int z_length){

    fprintf(data,"POINT_DATA %i\nSCALARS %s FLOAT 1\nLOOKUP_TABLE Info\n",x_length*y_length*z_length,val_name);
    int i;
    for(i=0;i<x_length*y_length*z_length;i++) {
            fprintf(data,"%i\n",1); //TODO: add data
    }
    fprintf(data,"\n");
}


void write_xvel_vtk(int frame){
    FILE *data;
    char varname[] = {'v','e','l','x'};
    const std::string& filename("velx.vtk.");

    data = vtk_open(filename.c_str(),frame);

    vtk_header(data);
    vtk_grid_points(data,1,X_DIM,Y_DIM,Z_DIM,1,Vec3(0,0,0));
    vtk_poly_vortex_and_celltype(data,X_DIM,Y_DIM,Z_DIM);
    vtk_pointdata(data,varname,X_DIM,Y_DIM,Z_DIM);
    fclose(data);
}

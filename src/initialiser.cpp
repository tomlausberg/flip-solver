/*
    FLIP FLuid Simulation: PBS Project
    Created by Tom Lausberg & Safira Piasko
    19.12.2018
*/

#include <vector>
#include "initialiser.hpp"
#include "simulation.hpp"

// initialiser a single rectangle,that drops on the floor
// origin = starting point, spacing = distance between the particles
std::vector<Particle> initialise_rectangle(int x, int y, int z, double spacing, Vec3 origin){
    std::vector<Particle> particles;
    Vec3 pos = origin;
    for(int i = 0; i < x; i++) {
        pos.y =origin.y;
        for(int j = 0; j < y; j++) {
            pos.z = origin.z;
            for(int k = 0; k < z; k++) {
                particles.push_back(Particle(pos));
                pos.z += spacing;
            }
            pos.y += spacing;
        }
        pos.x += spacing;
    }
    return particles;
}

std::vector<Particle> initialise_NPerCell(int x, int y, int z, int n, double dx, Vec3 origin){
    std::vector<Particle> particles;
    double dist = dx/n;

    Vec3 start_pos(origin.x*dx,origin.y*dx,origin.z*dx);
    Vec3 grid_pos = start_pos;
    Vec3 current_pos;
    for(int i = 0; i < x; i++) {
        for(int j = 0; j < y; j++) {
            for(int k = 0; k < z; k++) {
                grid_pos.x = start_pos.x + dx*i;
                grid_pos.y = start_pos.y + dx*j;
                grid_pos.z = start_pos.z + dx*k;
                current_pos = grid_pos;
                for(double a = dist/2; a < dx; a+=dist){
                    for(double b = dist/2; b < dx; b+=dist){
                        for(double c = dist/2; c < dx; c+=dist){
                            current_pos.x = grid_pos.x + a;
                            current_pos.y = grid_pos.y + b;
                            current_pos.z = grid_pos.z + c;
                            particles.push_back(Particle(current_pos));
                        }
                    }
                }

            }
        }
    }
    return particles;
}

void Simulation::source_particles(int x, int y, int z, int n, Vec3 origin){
    std::vector<Particle> particles;
    double dist = dx/n;

    Vec3 start_pos(origin.x*dx,origin.y*dx,origin.z*dx);
    Vec3 grid_pos = start_pos;
    Vec3 current_pos;
    int k = int(origin.z);
    for(int i = 0; i < x; i++) {
        for(int j = 0; j < y; j++) {
            if(cell(int(i+origin.x*dx),int(j+origin.y*dx),int(k+origin.z*dx))==AIR){
                grid_pos.x = start_pos.x + dx*i + dist/2;
                grid_pos.y = start_pos.y + dx*j + dist/2;
                grid_pos.z = start_pos.z + dx*k + dist/2;
                current_pos = grid_pos;
                for(int a = 0; a < n; a++){
                    current_pos.y = grid_pos.y;
                    for(int b = 0; b < n; b++){
                        current_pos.z = grid_pos.z;
                        for(int c = 0; c < n; c++){
                            particles.push_back(Particle(current_pos,Vec3(a/n,b/n,c/n)));
                            current_pos.x += dist;
                        }
                        current_pos.y += dist;
                    }
                    current_pos.z += dist;
                }
            }
        }
    }
    add_particles(particles);
}

void Simulation::source_cylinder(int r, Vec3 origin){
    std::vector<Particle> particles;
    int i,j,k;
    double z = origin.z;
    for(int x = -r; x <= r; ++x){
        for(int y = -r; y <= r; ++y){
            if(x*x+y*y<r*r){
                find_cell(dx,origin.x+x*dx*0.5,i);
                find_cell(dx,origin.y+y*dx*0.5,j);
                find_cell(dx,origin.z,k);
                if(cell(i,j,k) != FLUID) {
                    particles.push_back(Particle(Vec3(origin.x+x*dx*0.5,origin.y+y*dx*0.5,z)));
                    particles.push_back(Particle(Vec3(origin.x+x*dx*0.5,origin.y+y*dx*0.5,z+0.5*dx)));
                }
            }
        }
    }
    add_particles(particles);
}

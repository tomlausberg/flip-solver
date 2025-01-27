#ifndef UTIL_H
#define UTIL_H

#include <eigen3/Eigen/Dense>
#define X_DIM 20
#define Y_DIM 20
#define Z_DIM 20

struct Vec3 {
    double x;
    double y;
    double z;

    Vec3();
    Vec3(double);
    Vec3(double x_, double y_, double z_);

    Vec3 operator+(const double& h) {
        return {x + h, y + h, z + h};
    }

    Vec3 operator+(const Vec3& v) {
        return {x + v.x, y + v.y, z + v.z};
    }

    Vec3 operator*(const double& s) {
        return {x * s, y * s, z * s};
    }

    Vec3 operator/(const double& s) {
        return {x / s, y / s, z / s};
    }

    bool operator==(const Vec3& v){
        return (x==v.x && y==v.y && z==v.z);
    }

};

// Particle Class
class Particle {
private:

public:
    // Coord pos;
    Vec3 pos;
    Vec3 velocity;

    // constructors
    Particle (Vec3);
    Particle (Vec3, Vec3);
};

void find_cell(double, double, double&, int&);
void find_cell(double, double, int&);
void find_cell_center(double, double, double&, int&);

int clamp(int a, int lower, int upper);

double clamp(double a, double lower, double upper);

double cubicSplineKernel(double distance);

#endif

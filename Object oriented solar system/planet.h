#ifndef PLANET_H
#define PLANET_H
#include <cmath>

#include <armadillo>
#include <cstdlib>
#include <math.h>
#include "time.h"
using namespace arma;
using namespace std;

class planet
{
public:
    double mass;
    double position[2];
    double velocity[2];
    double fourpi2 = 4*M_PI*M_PI;
    double c = 63239.7263; // Speed of light in AU per year
    //Initializers
    planet();
    planet(double m, double x, double y, double velo_x, double velo_y);

    double distance(planet other_planet);
    double kinetic_energy();
    double potential_energy(planet other_planet);
    double angular_momentum();
    double forcebutnoxory(planet other_planet);
    double perihelion(planet other_planet);
    void find_perihelion(planet other_planet, int i, double step, std::ofstream &output);
};

#endif

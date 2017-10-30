#ifndef SOLVER_H
#define SOLVER_H
#include "planet.h"

#include <armadillo>
#include <cstdlib>
#include <math.h>
using namespace arma;
using namespace std;

class solver
{
public:
    friend class planet;

    int total_planets;
    vector <planet> all_planets;
    //Initializers
    solver();

    void add(planet newplanet);
    void print_positions(std::ofstream &output);
    void verlet(int n, double step);
    void verlet_perihelion(int n, double step);
    void euler(int n, double step);
    void center_of_mass();

};


#endif // SOLVER_H

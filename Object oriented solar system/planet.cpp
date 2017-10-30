#include "planet.h"

planet::planet()
{
    mass = 1;
    position[0] = 0;
    position[1] = 1;
    velocity[0] = -2*M_PI;
    velocity[1] = 0;


}

planet::planet(double m, double x, double y, double velo_x, double velo_y)
{

    mass = m;
    position[0] = x;
    position[1] = y;
    velocity[0] = velo_x;
    velocity[1] = velo_y;


}

double planet::distance(planet other_planet)
{
    double x1 = this->position[0];
    double y1 = this->position[1];

    double x2 = other_planet.position[0];
    double y2 = other_planet.position[1];


    return pow((pow(x1-x2,2) + pow(y1-y2,2)), 0.5);
}

double planet::kinetic_energy(){
    double vx = this->velocity[0];
    double vy = this->velocity[1];
    double vsquared = vx*vx + vy*vy;
    return 0.5 * this->mass * vsquared;
}

double planet::potential_energy(planet other_planet){
    double radius = this->distance(other_planet);
    return -fourpi2 *this->mass* other_planet.mass / radius;
}

double planet::angular_momentum(){
    double x = this->position[0];
    double y = this->position[1];
    double vx = this->velocity[0];
    double vy = this->velocity[1];
    return this->mass * (x*vy - y*vx);;
}



double planet::perihelion(planet other_planet){
    double r = this->distance(other_planet);
    double vx = this->velocity[0];
    double vy = this->velocity[1];
    double x = this->position[0];
    double y = this->position[1];
    double l = (x*vy - y*vx);
    return 3*l*l / (r*r*c*c);
}

void planet::find_perihelion(planet other_planet, int i, double step, std::ofstream &output){
    double tol = pow(10,-12);
    if (this->distance(other_planet) < 0.3075 + tol && this->distance(other_planet) > 0.3075 - tol){
        output << atan(this->position[1]/this->position[0]) << "  " << i*step << endl;
    }
}


double planet::forcebutnoxory(planet other_planet){
    double r = this->distance(other_planet);
    return fourpi2 * other_planet.mass / (r*r*r);
}




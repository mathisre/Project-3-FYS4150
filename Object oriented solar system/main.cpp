#include "planet.h"
#include "solver.h"
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <time.h>
using namespace std;


int main()
{
    int n = 10000;
    double final_time = 50;
    double step = final_time / n;
    cout << "Step size = " << step <<" years" << endl;
    //Planets mass in units of solar masses
    double earth_mass= 3*pow(10,-6); double jupiter_mass = 9.5*pow(10,-4);
    double merc_mass = 1.65 *pow(10,-7); double venus_mass = 2.45 * pow(10,-6);
    double mars_mass = 3.3*pow(10,-7); double saturn_mass = 2.75*pow(10,-4);
    double uranus_mass = 4.4*pow(10,-5); double neptune_mass = 5.15*pow(10,-5);
    double pluto_mass = 6.55*pow(10,-9);

    //Initial positions and velocities of planets + Pluto
    //double earth_x = 0; double earth_y = 1;  double earth_vx = -2*M_PI;  double earth_vy = 0;
    double earth_x = 0.988608; double earth_y = 0.225297;  double earth_vx = -1.483794;  double earth_vy = 6.10647;
    double jupiter_x = -4.62783; double jupiter_y = -2.8547; double jupiter_vx = 1.4149; double jupiter_vy = -2.2149;
    double mercx = -0.38964; double mercy = -0.0281635; double merc_vx = -1.25572; double merc_xy = -9.79417;
    double venus_x = -0.512269; double venus_y = 0.505719; double venus_vx = -5.17741; double venus_vy = -5.33447;
    double mars_x = -1.51055; double mars_y = 0.701684; double mars_vx = -1.94638; double mars_vy = -4.20472;
    double saturn_x = -0.410687; double saturn_y = -10.0466; double saturn_vx = 1.9243; double saturn_vy = 0.089789;
    double uranus_x = 17.87905; double uranus_y = 8.77302; double uranus_vx = -0.643322; double uranus_vy = 1.22269;
    double neptune_x = 28.6038; double neptune_y = -8.85453; double neptune_vx = 0.33138; double neptune_vy = 1.10195;
    double pluto_x = 10.5135; double pluto_y = -31.7162; double pluto_vx = 1.11339; double pluto_vy = 0.123435;
    double sun_x = 2.28064*pow(10,-3); double sun_y = 5.66903*pow(10,-3); double sun_vx = -0.001871801; double sun_vy = 0.002027733;


    //Create planet objects
    planet Earth = planet(earth_mass,earth_x, earth_y, earth_vx, earth_vy);
    planet Jupiter= planet(jupiter_mass, jupiter_x, jupiter_y, jupiter_vx, jupiter_vy);
    planet Mercury = planet(merc_mass, mercx, mercy, merc_vx, merc_xy);
    planet Venus = planet(venus_mass, venus_x, venus_y, venus_vx, venus_vy);
    planet Mars = planet(mars_mass, mars_x, mars_y, mars_vx, mars_vy);
    planet Saturn = planet(saturn_mass, saturn_x, saturn_y, saturn_vx, saturn_vy);
    planet Uranus = planet(uranus_mass, uranus_x, uranus_y, uranus_vx, uranus_vy);
    planet Neptune = planet(neptune_mass, neptune_x, neptune_y, neptune_vx, neptune_vy);
    planet Pluto = planet(pluto_mass, pluto_x, pluto_y, pluto_vx, pluto_vy);    
    planet Sun = planet(1, sun_x, sun_y, sun_vx, sun_vy);

    //Creating solver object and adding planets to solver
    solver verlet_solver = solver();
    verlet_solver.add(Sun);
    verlet_solver.add(Mercury);
    verlet_solver.add(Venus);
    verlet_solver.add(Earth);
    verlet_solver.add(Mars);
    verlet_solver.add(Jupiter);
    verlet_solver.add(Saturn);
    verlet_solver.add(Uranus);
    verlet_solver.add(Neptune);
    verlet_solver.add(Pluto);

    verlet_solver.center_of_mass();
    verlet_solver.verlet(n, step);
    verlet_solver.center_of_mass();


    //Coordinates for Mercury and Sun using Sun body center coordinates.
    //These are used for perihelion precision calculations
    planet Sun2 = planet(1, 0, 0, 0, 0);
    planet Mercury2 = planet(merc_mass, 0.3075, 0, 0, 12.44);

    solver perihelion_precision = solver();
    perihelion_precision.add(Sun2);
    perihelion_precision.add(Mercury2);

    //perihelion_precision.verlet_perihelion(n,step);

    //Forward Euler solver with just the stationary sun and Earth
    planet Earth2 = planet(earth_mass, 0, 1, -2*M_PI, 0);
    solver euler_solver = solver();
    euler_solver.add(Sun2);
    euler_solver.add(Earth2);
    //euler_solver.euler(n,step);

    return 0;
}

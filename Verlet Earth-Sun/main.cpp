#include <iostream>
#include <armadillo>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <string>
#include <time.h>
using namespace arma;
using namespace std;
ofstream ofile;

int main(int argc, char* argv[])
{
    if (argc <= 2){
        cout << "Bad usage, should be file then n in command line" << endl;
        exit(1);
    }
    string filename;
    filename = argv[1]; //importing filename and n from command line
    ofile.open(filename);
    int n = atoi(argv[2]);
    ofile << n << endl;
    double final_time = 10;
    double step = final_time / n;

    //Earth initial conditions
    vec x = zeros<vec>(n); vec y = zeros<vec>(n); vec velo_x = zeros<vec>(n); vec velo_y = zeros<vec>(n);
    x(0) = 0.0; y(0) = 1.0; velo_x(0) = -2*M_PI; velo_y(0) = 0.0;
    vec acc_x = zeros<vec>(n); vec acc_y = zeros<vec>(n);

    double fourpi2 = 4*M_PI*M_PI; //Precalculation factor
    acc_x(0) = -fourpi2*x(0) ; //Initial accelerations
    acc_y(0) = -fourpi2*y(0) ;


    for (int i = 0; i<n-1; i++){
        double radius = pow(x(i)*x(i) + y(i)*y(i), 0.5);
        x(i+1) = x(i) + velo_x(i) * step + 0.5*step*step*acc_x(i);
        y(i+1) = y(i) + velo_y(i) * step + 0.5*step*step*acc_y(i);

        acc_x(i+1) = -fourpi2*x(i+1) * (pow(radius,-3));
        acc_y(i+1) = -fourpi2*y(i+1) * (pow(radius,-3));

        velo_x(i+1) = velo_x(i) + 0.5*step*(acc_x(i) + acc_x(i+1));
        velo_y(i+1) = velo_y(i) + 0.5*step*(acc_y(i) + acc_y(i+1));
    }
    for (int i = 0; i<n-1; i++){
        ofile << x(i) <<"  " <<  y(i) << "  " << velo_x(i) << "  " <<  velo_y(i) << endl;
    }

    ofile.close();
    return 0;
}

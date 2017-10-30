#include "solver.h"
#include <vector>

solver::solver(){
    total_planets = 0;
}
void solver::add(planet newplanet){
    total_planets += 1;
    all_planets.push_back(newplanet);
}


void solver::print_positions(std::ofstream &output){
    for (int i=0; i< total_planets; i++){
        output << all_planets[i].position[0] << "  ";
        output << all_planets[i].position[1] << "  ";

        }
    output << endl;
}

void solver::center_of_mass(){
    double x_center = 0; double y_center = 0; double m_total = 0;
    for (int k=0; k<total_planets; k++){
        m_total += all_planets[k].mass;
        x_center+= all_planets[k].mass * all_planets[k].position[0];
        y_center+= all_planets[k].mass * all_planets[k].position[1];
    }
    x_center = x_center / m_total;
    y_center = y_center / m_total;
    cout << "x center of mass =  " <<  x_center << endl;
    cout << "y center of mass =  " <<  y_center << endl;
    for (int k=0; k<total_planets; k++){
        all_planets[k].position[0] =   all_planets[k].position[0] - x_center;
        all_planets[k].position[1] =   all_planets[k].position[1] - y_center;
    }

}

void solver::verlet(int n, double step){

    char *filename = new char[1000];
    sprintf(filename, "verleto.txt");
    std::ofstream output_file(filename);
    output_file << n << endl;
    print_positions(output_file);

    double tol = pow(10,-3);

    vector<double> acc_curr_x(total_planets,0);
    vector<double> acc_curr_y(total_planets,0);
    vector<double> acc_new_x(total_planets,0);
    vector<double> acc_new_y(total_planets,0);

    vector<double> total_energy(total_planets,0);
    vector<double> potential_energy(total_planets,0);
    vector<double> initial_total_energy(total_planets,0);
    vector<double> initial_potential_energy(total_planets,0);


    vector<double> angular_momentum(total_planets,0);
    vector<double> initial_angular_momentum(total_planets,0);



    //Initial values for acceleration and energy
    for (int m=0; m< total_planets; m++){        

        for (int k = 0; k < total_planets; k++){
            if (k != m){
                acc_curr_x[m] += -all_planets[m].forcebutnoxory(all_planets[k]) * (all_planets[m].position[0] - all_planets[k].position[0]);
                acc_curr_y[m] += -all_planets[m].forcebutnoxory(all_planets[k]) * (all_planets[m].position[1] - all_planets[k].position[1]);
                initial_potential_energy[m] += all_planets[m].potential_energy(all_planets[k]);
            }
            initial_total_energy[m] = all_planets[m].kinetic_energy() + initial_potential_energy[m];
            initial_angular_momentum[m] = all_planets[m].angular_momentum();
        }
    }
    clock_t start, finish;
    start = clock();
    for (int i =0; i<n-1; i++){
        for (int k = 0; k< total_planets; k++){
            all_planets[k].position[0] += step*all_planets[k].velocity[0] + 0.5*step*step*acc_curr_x[k];
            all_planets[k].position[1] += step*all_planets[k].velocity[1] + 0.5*step*step*acc_curr_y[k];
        }
        for (int m=0; m< total_planets; m++){
            //Put to zero because of the += instead of just =
            acc_new_x[m] = 0;
            acc_new_y[m] = 0;
            for (int k = 0; k< total_planets; k++){
                if (k != m){
                    acc_new_x[m] += -all_planets[m].forcebutnoxory(all_planets[k]) * (all_planets[m].position[0] - all_planets[k].position[0]);
                    acc_new_y[m] += -all_planets[m].forcebutnoxory(all_planets[k]) * (all_planets[m].position[1] - all_planets[k].position[1]);
            }}
        }
        for (int k = 0; k < total_planets; k++){
            all_planets[k].velocity[0] += 0.5*step*(acc_curr_x[k] + acc_new_x[k]);
            all_planets[k].velocity[1] += 0.5*step*(acc_curr_y[k] + acc_new_y[k]);
            acc_curr_x[k] = acc_new_x[k];
            acc_curr_y[k] = acc_new_y[k];

            //Calculate energy and angular momentum with updated position and velocity
            potential_energy[k] = 0;
            for (int l = 0; l < total_planets; l++){
                if (k != l){
                    potential_energy[k] += all_planets[k].potential_energy(all_planets[l]);
                }
            }
            total_energy[k] = all_planets[k].kinetic_energy() + potential_energy[k];
            angular_momentum[k] = all_planets[k].angular_momentum();
        }

        // Unit tests for energy and angular momentum
        for (int k = 0; k < total_planets; k++){
            if (total_energy[k] > initial_total_energy[k] + tol || total_energy[k] <  initial_total_energy[k] - tol ){
                cout << "Energy not conserved after " << (i+1)*step << " years in planet " << k << endl;
                cout << total_energy[k]  << "  and  " << initial_total_energy[k] ;
                exit(EXIT_FAILURE);
            }

            if (angular_momentum[k] > initial_angular_momentum[k] + tol || angular_momentum[k] < initial_angular_momentum[k] - tol){
                cout << "Angular momentum not conserved after " << (i+1)*step << " years in planet " << k << endl;
                cout << angular_momentum[k]  << "  and  " << initial_angular_momentum[k] ;
                exit(EXIT_FAILURE);
            }
        }
        print_positions(output_file);
    }


    finish = clock();
    double ftime = double (finish-start)/CLOCKS_PER_SEC;
    cout<< endl << "Verlet time = " << ftime<< "s"<< endl << endl;
}

void solver::verlet_perihelion(int n, double step){

    char *filename = new char[1000];
    sprintf(filename, "perihelion.txt");
    std::ofstream output_file(filename);
    output_file << n << endl;
    all_planets[1].find_perihelion(all_planets[0],0, step, output_file);

    double tol = pow(10,-3);

    vector<double> acc_curr_x(total_planets,0);
    vector<double> acc_curr_y(total_planets,0);
    vector<double> acc_new_x(total_planets,0);
    vector<double> acc_new_y(total_planets,0);

    vector<double> total_energy(total_planets,0);
    vector<double> potential_energy(total_planets,0);
    vector<double> initial_total_energy(total_planets,0);
    vector<double> initial_potential_energy(total_planets,0);


    vector<double> angular_momentum(total_planets,0);
    vector<double> initial_angular_momentum(total_planets,0);



    //Initial values for acceleration and energy
    for (int m=1; m< total_planets; m++){

        for (int k = 0; k < total_planets; k++){
            if (k != m){
                acc_curr_x[m] += -all_planets[m].forcebutnoxory(all_planets[k]) * (all_planets[m].position[0] - all_planets[k].position[0])*(1+all_planets[m].perihelion(all_planets[k]));
                acc_curr_y[m] += -all_planets[m].forcebutnoxory(all_planets[k]) * (all_planets[m].position[1] - all_planets[k].position[1])*(1+all_planets[m].perihelion(all_planets[k]));
                initial_potential_energy[m] += all_planets[m].potential_energy(all_planets[k]);
            }
            initial_total_energy[m] = all_planets[m].kinetic_energy() + initial_potential_energy[m];
            initial_angular_momentum[m] = all_planets[m].angular_momentum();
        }
    }
    clock_t start, finish;
    start = clock();
    for (int i =0; i<n-1; i++){
        for (int k = 1; k< total_planets; k++){
            all_planets[k].position[0] += step*all_planets[k].velocity[0] + 0.5*step*step*acc_curr_x[k];
            all_planets[k].position[1] += step*all_planets[k].velocity[1] + 0.5*step*step*acc_curr_y[k];
        }
        for (int m=0; m< total_planets; m++){
            //Put to zero because of the += instead of just =
            acc_new_x[m] = 0;
            acc_new_y[m] = 0;
            for (int k = 0; k< total_planets; k++){

                if (k != m){
                    acc_new_x[m] += -all_planets[m].forcebutnoxory(all_planets[k]) * (all_planets[m].position[0] - all_planets[k].position[0])*(1+ all_planets[m].perihelion(all_planets[k]));
                    acc_new_y[m] += -all_planets[m].forcebutnoxory(all_planets[k]) * (all_planets[m].position[1] - all_planets[k].position[1])*(1+ all_planets[m].perihelion(all_planets[k]));
            }}
        }
        for (int k = 0; k < total_planets; k++){
            all_planets[k].velocity[0] += 0.5*step*(acc_curr_x[k] + acc_new_x[k]);
            all_planets[k].velocity[1] += 0.5*step*(acc_curr_y[k] + acc_new_y[k]);
            acc_curr_x[k] = acc_new_x[k];
            acc_curr_y[k] = acc_new_y[k];

            //Calculate energy and angular momentum with updated position and velocity
            potential_energy[k] = 0;
            for (int l = 0; l < total_planets; l++){
                if (k != l){
                    potential_energy[k] += all_planets[k].potential_energy(all_planets[l]);
                }
            }
            total_energy[k] = all_planets[k].kinetic_energy() + potential_energy[k];
            angular_momentum[k] = all_planets[k].angular_momentum();
        }

        // Unit tests for energy and angular momentum
        for (int k = 0; k < total_planets; k++){
            if (total_energy[k] > initial_total_energy[k] + tol || total_energy[k] <  initial_total_energy[k] - tol ){
                cout << "Energy not conserved after " << (i+1)*step << " years in planet " << k << endl;
                cout << total_energy[k]  << "  and  " << initial_total_energy[k] ;
                exit(EXIT_FAILURE);
            }

            if (angular_momentum[k] > initial_angular_momentum[k] + tol || angular_momentum[k] < initial_angular_momentum[k] - tol){
                cout << "Angular momentum not conserved after " << (i+1)*step << " years in planet " << k << endl;
                cout << angular_momentum[k]  << "  and  " << initial_angular_momentum[k] ;
                exit(EXIT_FAILURE);
            }
        }
        all_planets[1].find_perihelion(all_planets[0],i, step, output_file);
        }

    finish = clock();
    double ftime = double (finish-start)/CLOCKS_PER_SEC;
    cout<< endl << "Verlet time = " << ftime<< "s"<< endl << endl;
}



void solver::euler(int n, double step){

    char *filename = new char[1000];
    sprintf(filename, "euleroo.txt");
    std::ofstream output_file(filename);
    output_file << n << endl;
    print_positions(output_file);

    double tol = pow(10,-1);

    vector<double> acc_curr_x(total_planets,0);
    vector<double> acc_curr_y(total_planets,0);

    clock_t start, finish;
    start = clock();
    for (int i =0; i<n-1; i++){
        for (int m=0; m< total_planets; m++){
            //Put to zero because of the += instead of just =
            acc_curr_x[m] = 0;
            acc_curr_y[m] = 0;
            for (int k = 0; k< total_planets; k++){
                if (k != m){
                    acc_curr_x[m] += -all_planets[m].forcebutnoxory(all_planets[k]) * (all_planets[m].position[0] - all_planets[k].position[0]);
                    acc_curr_y[m] += -all_planets[m].forcebutnoxory(all_planets[k]) * (all_planets[m].position[1] - all_planets[k].position[1]);
            }}
        }
        for (int k = 0; k< total_planets; k++){
            all_planets[k].position[0] += step*all_planets[k].velocity[0];
            all_planets[k].position[1] += step*all_planets[k].velocity[1];
            all_planets[k].velocity[0] += step*acc_curr_x[k];
            all_planets[k].velocity[1] += step*acc_curr_y[k];
        }
        print_positions(output_file);
        }
    finish = clock();
    double ftime = double (finish-start)/CLOCKS_PER_SEC;
    cout<< endl << "Euler time = " << ftime<< "s"<< endl << endl;
}


/********************************************************************************
 * Author:   Gerrit Wellecke
 * Project:  Chimera states in populations of Kuramoto oscillators
 * Date:     March 2020
 *
 * Main function for the project on chimera states. As this project grew
 * gradually I only wrote one main for both systems (explicit and reduced). In
 * order to run the program properly only one of the two integration blocks
 * should be active at once, see below.
 * The integration functions were left here to simplify adaptations of code. All
 * functions that are "set in stone" are outsourced to other files.
 *
 * If no parallelisation is wanted the #pragma's in main can simply be commented
 * out.
 * Further the makefile might need to be adapted, depending on whether or not
 * parallelisation is wanted and what system architecture is used.
 *
 * OUTPUT is written to binary files. Each timestep of pop*osci oscillators is
 * written in sequence using 8-byte double precision floats.
 *******************************************************************************/
#include "functions.hpp"        // functions used in integration() and main()
#include "kuramoto.hpp"         // kuramoto model class
#include "ottantonsen.hpp"      // Ott-Antonsen reduction

#include <boost/numeric/odeint.hpp>
#include <omp.h>

#include <chrono>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace boost::numeric::odeint;

using state_type = std::vector<double>;

/* *************************************************************************** */

// integration of explicit Kuramoto system
void integration(const double omega, const double alpha, const int pop, const int osci,
        std::vector< std::vector<double> > &K, const uint maxIteration, std::string simName,
        const std::vector<double> phase, const std::vector<double> variance, const double dt,
        const Distribution dist, const long seedGiven
        ) {
    // state vector with initial phases
    state_type x;
    x.resize(pop * osci);

    // set initial conditions and return seed of RNG
    long seed {setInitCond(x, pop, osci, variance, phase, dist, seedGiven)};


    // filename w/ given parameters
    std::string fileName {"data/" + simName};

    // fileIO
    std::ofstream output(fileName, std::ios::binary);
    std::ofstream ordParamOut("data/" + simName + ".order", std::ios::binary);

    // initialize stepper
    runge_kutta4<state_type> stepper;

    // initialize system
    Kuramoto kuramoto_func{omega, K, alpha, pop, osci};

    // write log-file
    kuramoto_func.writeLog(simName, maxIteration, dt, seed);

    // write initial time step
    output.write(reinterpret_cast<char*>(x.data()), x.size() * sizeof(double));

    // write initial order parameter
    std::vector<double> ordP {orderParam(x, pop, osci)};
    ordParamOut.write(reinterpret_cast<char*>(ordP.data()), ordP.size() * sizeof(double));

    // integration loop
    uint iteration {0};

    while (iteration < maxIteration) {
        // as t isn't explicitly considered it is just set 0 here
        stepper.do_step(kuramoto_func, x, 0., dt);

        // renormalize phases in state vector
        renorm(x);

        // write time step to file
        output.write(reinterpret_cast<char*>(x.data()), x.size() * sizeof(double));

        // compute order parameter -- write to separate file
        ordP = orderParam(x, pop, osci);
        ordParamOut.write(reinterpret_cast<char*>(ordP.data()), ordP.size() * sizeof(double));

        // iterate integration step
        iteration++;
    }

    output.close();
    ordParamOut.close();
}

// integration of reduced OA system (analogous to integration of Kuramoto sys.)
void integrationOA(const double omega, const double alpha, const int pop, const int osci,
        std::vector< std::vector<double> > &K, const uint maxIteration, std::string simName,
        const double dt, /* const std::string Kname, */ const long seed
        ) {
    // state vector with init. cond.
    state_type x;
    x.resize(3);

    // set initial conditions and get seed
    // uncomment line below if wanting to reproduce initial state of Kuramoto sys.

    // long seed {setInitCondOAfromKuramoto(x, pop, osci, Kname)};

    // if using above line, comment the line below
    setInitCondOA(x, seed);

    // filename w/ given parameters
    std::string fileName {"data/" + simName + ".order"};

    // fileIO
    std::ofstream output(fileName, std::ios::binary);

    // initialize stepper
    runge_kutta4<state_type> stepper;

    // initialize system
    OttAntonsen OA_func{omega, K, alpha, pop, osci};

    // write log-file
    OA_func.writeLog(simName, maxIteration, dt, seed);

    // write initial time step
    output.write(reinterpret_cast<char*>(x.data()), x.size() * sizeof(double));

    // integration loop
    uint iteration {0};

    while (iteration < maxIteration) {
        // as t isn't explicitly considered it is just set 0 here
        stepper.do_step(OA_func, x, 0., dt);

        // renormalize mean phase difference
        x[2] = fmod(x[2], 2*M_PI);
        if (x[2] < 0)
            x[2] += 2*M_PI;

        // write time step
        output.write(reinterpret_cast<char*>(x.data()), x.size() * sizeof(double));

        // iterate integration step
        iteration++;
    }

    output.close();
}

/* *************************************************************************** */

int main(void) {
    // system parameters as needed by Kuramoto and OttAntonsen
    double omega {0.};
    double beta {.01}; // for beta -> 0, P(Chimera) -> 1;
    double alpha {.5 * M_PI - beta};
    int pop {2};
    int osci {128};

    // parameters for initial conditions
    std::vector<double> phase{0., 0.}; // (0, 0) -- in phase
    std::vector<double> variance{.2, M_PI};
    Distribution distribution {Distribution::NORMAL}; // only for explicit case

    // parameters for integration
    uint maxIteration {500000};
    double dt{.1};

    // vectors to iterate over
    std::vector<double> Avec {.05, .15, .20, .40};
    std::vector<long> seedVec {4578, 130547, 1337, 4578};

    std::vector<std::string> Kname {"data/seed_4578_A_0.050000",
        "data/seed_130547_A_0.150000", "data/seed_1337_A_0.200000",
        "data/seed_4578_A_0.400000"};

    // Numerical solution of explicit Kuramoto system
    /*
#pragma omp parallel for
    for (int par = 0; par < static_cast<int>(Avec.size()); par++) {
        // set seed
        long seedGiven {seedVec[par]};

        // coupling matrix
        double A {Avec[par]};
        double mu {.5 + A/2.};
        double nu {.5 - A/2.};

        std::vector< std::vector<double> > K {{mu, nu}, {nu, mu}};

        // name for given integration
        // std::string simName {"seed_" + std::to_string(seedGiven) + "_A_" + std::to_string(A)};
        std::string simName{"OttAntonsen_1"};

        // just benchmarking fun
        auto start {std::chrono::high_resolution_clock::now()};

        // integration of Kuramoto system w/ above setup
        integration(omega, alpha, pop, osci, K, maxIteration, simName,
                    phase, variance, dt, distribution, seedGiven);

        // continue benchmarking fun
        auto stop {std::chrono::high_resolution_clock::now()};

        auto benchmark {std::chrono::duration_cast<std::chrono::microseconds>(stop - start)};

        std::cout.precision(6);
        // std::cout << "Thread " << omp_get_thread_num() << " has finished calculation "
                  // << par + 1 << '\n';
        std::cout << "Execution time: " << std::fixed << benchmark.count() * 1e-6 << "s\n";
    }
    */

    // Numerical solution of Ott-Antonsen reduction
    // set seed
    long seedGiven {std::time(nullptr)};

#pragma omp parallel for
    //for (int par = 0; par < static_cast<int>(Avec.size()); par++) {
    for (int par = 0; par <= 50 ; par++) {
        //seedGiven = seedVec[par];

        // coupling matrix
        //double A {Avec[par]};
        double A {.01 * par};
        double mu {.5 + A/2.};
        double nu {.5 - A/2.};

        std::vector< std::vector<double> > K {{mu, nu}, {nu, mu}};

        // name for given integration
        std::string simName{"OttAntonsen_A" + std::to_string(A)};

        // just benchmarking fun
        auto start {std::chrono::high_resolution_clock::now()};

        integrationOA(omega, alpha, pop, osci, K, maxIteration, simName,
                        dt, /* Kname[par], */ seedGiven);

        // continue benchmarking fun
        auto stop {std::chrono::high_resolution_clock::now()};

        auto benchmark {std::chrono::duration_cast<std::chrono::microseconds>(stop - start)};

        std::cout.precision(6);
        std::cout << "Execution time: " << std::fixed << benchmark.count() * 1e-6 << "s\n";
    }

    return 0;
}

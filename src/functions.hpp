/********************************************************************************
 * Author:   Gerrit Wellecke
 * Project:  Chimera states in populations of Kuramoto oscillators
 * Date:     March 2020
 * 
 * Functions used by main.cpp in separate file to enhance readability of main.
 *******************************************************************************/
#ifndef FUNCTIONS_HPP
#define FUNCTIONS_HPP

#include <cmath>
#include <fstream>
#include <vector>
#include <random>

using state_type = std::vector<double>;


enum class Distribution {
    UNIFORM,
    NORMAL,
    CAUCHY,
};


// adaptive coupling. actually DEPRECATED now
void adaptCoupling(std::vector< std::vector<double> > &K, const std::vector<double> &ordP,
                   const int &pop);

// renormalization of the state vector onto [0, 2pi)
void renorm(state_type &x);

// compute order parameter
std::vector<double> orderParam(const state_type &x, const int numPop, const int numOsci);

// set initial state for numerical system
double initCond(std::mt19937_64 &mersenne, Distribution distribution, const double phase,
                const double variance);

// set initial state for Ott-Antonsen reduced system
double initCondOA(std::mt19937_64 &mersenne);

// set initial phases of Kuramoto-system
long setInitCond(state_type &x, const int pop, const int osci,
                 const std::vector<double> variance, const std::vector<double> phase,
                 const Distribution distribution, const long seed);

// set initial order parameters and phase difference of reduced system
long setInitCondOA(state_type &x, const long seed);

// copy inital order parameters and phase difference from Kuramoto-system for the
// reduced case of the Ott-Antonsen (OA) system
long setInitCondOAfromKuramoto(state_type &x, const int pop, const int osci,
        const std::string fileTS);
#endif

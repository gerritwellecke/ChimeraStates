/********************************************************************************
 * Author:   Gerrit Wellecke
 * Project:  Chimera states in populations of Kuramoto oscillators
 * Date:     March 2020
 * 
 * Kuramoto model with arbitrary number of oscillators in point-like clusters.
 * Each coupled to the others as given by a coupling matrix m_K.
 * 
 * OOP approach to boost::odeint by overloading the operator(). 
 * Method to write log-file of current simulation to ensure reproducibility.
 *******************************************************************************/
#ifndef KURAMOTO_HPP
#define KURAMOTO_HPP

#include <cmath>
#include <vector>
#include <fstream>
#include <iomanip>

using state_type = std::vector<double>;

// Put this class into its own files (*.hpp, *.cpp)
// maybe another thought: all parameters could just be const
class Kuramoto {
protected:
    // parameters
    double m_omega;                          // natural frequency
    std::vector< std::vector<double> > m_K;  // coupling within/between populations (this is a matrix)
    double m_alpha;                          // phase lag
    int m_numPop;                            // number of populations
    int m_numOsci;                           // number of oscillators per population

public:
    // default constructor
    Kuramoto(double omega, std::vector< std::vector<double> > K, double alpha, int pop, int osci);

    // ()-operator for RHS of ODE
    void operator() (const state_type &x, state_type &dxdt, const double /* t */);

    // generate logfile for given simulation
    void writeLog(const std::string &simName, int maxIteration, double dt, int seed);
};

#endif

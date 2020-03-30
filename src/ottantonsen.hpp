/********************************************************************************
 * Author:   Gerrit Wellecke
 * Project:  Chimera states in populations of Kuramoto oscillators
 * Date:     March 2020
 * 
 * Inherited class from Kuramoto, to save code repetition (this way the logs are
 * also identical in format).
 * Different ODE through overloading operator().
 *******************************************************************************/
#ifndef OTTANTONSEN_HPP
#define OTTANTONSEN_HPP

#include "kuramoto.hpp"

#include <cmath>
#include <vector>

using state_type = std::vector<double>;

class OttAntonsen : public Kuramoto {
public:
    // RHS of ODE
    OttAntonsen(double omega, std::vector< std::vector<double> > K, double alpha,
                int pop, int osci)      : Kuramoto(omega, K, alpha, pop, osci) {};

    void operator() (const state_type &x, state_type &dxdt, const double /* t */);
};

#endif

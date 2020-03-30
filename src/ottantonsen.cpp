#include "ottantonsen.hpp"

// RHD of ODE in reduced system
void OttAntonsen::operator() (const state_type &x, state_type &dxdt,
                              const double /* t */) {
    // in literature mostly $\beta$
    double b {M_PI / 2. - m_alpha};

    // straightforward ODE
    dxdt[0] = (1 - x[0]*x[0]) / 2. * (m_K[0][0] * x[0] * sin(b) + m_K[0][1] * x[1]
              * sin(b - x[2]));
    dxdt[1] = (1 - x[1]*x[1]) / 2. * (m_K[1][1] * x[1] * sin(b) + m_K[1][0] * x[0]
              * sin(b + x[2]));
    dxdt[2] = (1 + x[1]*x[1]) / (2. * x[1])
              * (m_K[1][1] * x[1] * cos(b) + m_K[1][0] * x[0] * cos(b + x[2])
              - (1 + x[0]*x[0]) / (2. * x[0])
              * (m_K[0][0] * x[0] * cos(b) + m_K[0][1] * x[1] * cos(b - x[2])));
}

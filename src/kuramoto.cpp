#include "kuramoto.hpp"


// default constructor
Kuramoto::Kuramoto(double omega, std::vector< std::vector<double> > K, double alpha, int pop, int osci)
    : m_omega{omega}, m_K{K}, m_alpha{alpha}, m_numPop{pop}, m_numOsci{osci} {};


// ()-operator w/ 3 args for RHS of ODE (as required by boost::odeint)
void Kuramoto::operator() (const state_type &x, state_type &dxdt, const double /* t */) {
    // dimension of state vector
    const double N_ss {static_cast<double>(m_numOsci)};

    // compute RHS
    // iterate over all populations
    for (int s=0; s<m_numPop; s++) {

        // iterate over oscillators
        for (int i=0; i<m_numOsci; i++) {
            dxdt[s * m_numOsci + i] = m_omega;

            // iterate over all populations for coupling
            for (int ss=0; ss<m_numPop; ss++) {
                // iterate over all oscillators for coupling
                for (int j=0; j<m_numOsci; j++) {
                    double sinarg = x[ss*m_numOsci + j] - x[s * m_numOsci + i];

                    dxdt[s * m_numOsci + i] += m_K[s][ss] / N_ss * sin(sinarg - m_alpha);
                }
            }
        }
    }
}

// function to write logfile with all relevant parameters
void Kuramoto::writeLog(const std::string &simName, int maxIteration, double dt, int seed) {
    int w1{60}, w2{20};

    std::ofstream log("data/" + simName + ".log");

    log << std::setw(w1) << std::setfill(' ') << std::left << "Simulation:"
        << std::right << std::setw(w2) << simName << "\n\n";
    log << std::setw(w1) << std::setfill(' ') << std::left << "Number of populations:"
        << std::right << std::setw(w2) << m_numPop << '\n';
    log << std::setw(w1) << std::setfill(' ') << std::left << "Number of oscillators per population:"
        << std::right << std::setw(w2) << m_numOsci << '\n';
    log << std::setw(w1) << std::setfill(' ') << std::left << "Natural frequency:"
        << std::right << std::setw(w2) << m_omega << '\n';
    log << std::setw(w1) << std::setfill(' ') << std::left << "Phase lag:"
        << std::right << std::setw(w2) << m_alpha << "\n\n";
    log << std::setw(w1) << std::setfill(' ') << std::left << "Integration timestep:"
        << std::right << std::setw(w2) << dt << "\n";
    log << std::setw(w1) << std::setfill(' ') << std::left << "Total timesteps:"
        << std::right << std::setw(w2) << maxIteration + 1 << "\n\n";
    log << std::setw(w1) << std::setfill(' ') << std::left << "Random seed:"
        << std::right << std::setw(w2) << seed << "\n\n";

    log << std::left << "Coupling matrix:\n";

    for (auto i : m_K) {
        for (auto j : i)
            log << j << ' ';
        log << '\n';
    }

    log.close();
}

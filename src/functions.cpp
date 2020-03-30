#include "functions.hpp"

// apply adaptive coupling
void adaptCoupling(std::vector< std::vector<double> > &K, const std::vector<double> &ordP,
                   const int &pop) {
    for (int i=0; i<pop; i++) {
        for (auto &coupl : K[i]) {
            coupl *= ordP[i];
        }
    }
}


// renormalise all phases to interval [0, 2\pi)
void renorm(state_type &x) {
    for (auto &val : x) {
        val = fmod(val, 2*M_PI);
        if (val < 0)
            val += 2*M_PI;
    }
}


// compute order parameter
std::vector<double> orderParam(const state_type &x, const int numPop, const int numOsci) {
    double N_s {static_cast<double>(numOsci)}, cosSum{};

    // order parameter defined as: r(t) = |<exp(i theta(t))>|
    std::vector<double> ordParam;
    ordParam.resize(numPop);

    // fill ordParam w/ zeros
    for (auto p : ordParam)
        p = 0;

    // compute order parameter for each population
    for (int s=0; s<numPop; s++) {
        cosSum = 0;
        for (int i=0; i<numOsci; i++) {
            for (int j=0; j<i; j++) {
                cosSum += cos(x[s * numOsci + i] - x[s * numOsci + j]);
            }
        }

        ordParam[s] += sqrt(N_s + 2 * cosSum) / (N_s);
    }

    return ordParam;
}

// use given distribution to set initial phase
double initCond(std::mt19937_64 &mersenne, Distribution distribution, const double phase,
                const double variance) {
    switch (distribution) {
        case Distribution::CAUCHY: {
            std::cauchy_distribution<double> initCondC{M_PI + phase, variance};
            return initCondC(mersenne);
        }
        case Distribution::NORMAL: {
            std::normal_distribution<double> initCondN{M_PI + phase, variance};
            return initCondN(mersenne);
        }
        // use uniform distribution as default
        case Distribution::UNIFORM:
        default: {
            std::uniform_real_distribution<double> initCondU{0, 2*M_PI};
            return initCondU(mersenne);
        }
    }
}

// returns random value [0, 1] for initial value of order parameter
double initCondOA(std::mt19937_64 &mersenne) {
    std::uniform_real_distribution<double> initCondU{0, 1};
    return initCondU(mersenne);
}

// set initial phases of Kuramoto-system
long setInitCond(state_type &x, const int pop, const int osci,
                 const std::vector<double> variance, const std::vector<double> phase,
                 const Distribution distribution, const long seed) {
    // initialise RNG
    std::mt19937_64 mersenne{static_cast<std::mt19937_64::result_type>(seed)};

    // set initial conditions
    for (int s=0; s<pop; s++) {
        for (int i=0; i<osci; i++) {
            do {
                x[s * osci + i] = initCond(mersenne, distribution, phase[s], variance[s]);
            } while (x[s * osci + i] > 2 * M_PI || x[s * osci + i] < 0);
        }
    }

    return seed;
}

// set initial order parameters and phase difference of reduced system
long setInitCondOA(state_type &x, const long seed) {
    std::mt19937_64 mersenne{static_cast<std::mt19937_64::result_type>(seed)};
    x[0] = initCondOA(mersenne);
    x[1] = initCondOA(mersenne);
    x[2] = initCond(mersenne, Distribution::UNIFORM, 0., 0.);

    return seed;
}

// copy inital order parameters and phase difference from Kuramoto-system for the
// reduced case of the Ott-Antonsen (OA) system
long setInitCondOAfromKuramoto(state_type &x, const int pop, const int osci,
        const std::string fileTS) {
    // open fstreams to files
    std::ifstream timeSeries(fileTS, std::ios::binary);
    std::ifstream orderParameter(fileTS + ".order", std::ios::binary);

    // read init state to vector
    std::vector<double> initState(pop * osci);
    timeSeries.read(reinterpret_cast<char*>(initState.data()),
                    initState.size() * sizeof(double));


    // read r values from *.order file
    std::vector<double> initR(pop);
    orderParameter.read(reinterpret_cast<char*>(initR.data()),
                        initR.size() * sizeof(double));

    // close fstreams
    timeSeries.close();
    orderParameter.close();

    // compute psi from init state of time series
    std::vector<double> phi {0., 0.};

    for (int s=0; s<pop; s++) {
        for (int i=0; i<osci; i++){
            phi[s] += initState[s * osci + i];
        }
        phi[s] /= osci;
    }

    // compute mean phase difference between clusters/populations
    double psi{phi[0] - phi[1]};
    if (psi < 0)
        psi += 2 * M_PI;

    // set initial state
    x[0] = initR[0];
    x[1] = initR[1];
    x[2] = psi;

    return 0;
}

#include <Eigen/Core>
#include <vector>
#include <iostream>
#include "writer.hpp"
#include <cmath>

/// Uses the Heun's method to compute u from time 0 to time T 
/// for the ODE $u'=e^{-2t}-2u$
///
/// @param[out] u at the end of the call (solution at all time steps)
/// @param[out] time contains the time levels
/// @param[in] u0 the initial data
/// @param[in] dt the step size
/// @param[in] T the final time up to which to compute the solution.
///

// function f=u'=e^(-2t)-2u to avoid writting too many times the same function
double f(double t, double u) {

    return std::exp(-2. * t) - 2. * u;
}

void Heun(std::vector<double> & u, std::vector<double> & time,
          const double & u0, double dt, double T) {

    const unsigned int nsteps = std::round(T/dt);

// (write your solution here)

    /// Start of my solution ///

    u.resize(nsteps+1);
    time.resize(nsteps+1);

    // initialize the first steps
    double y1;
    double y2;

    // initialize the initial values
    u[0] = u0;
    time[0] = 0.;

    for (int i = 0; i < nsteps; i++) {

        // I declare the function f above to avoid writting too many times the same function.
        // compute the first steps
        y1 = u[i];
        y2 = u[i] + dt * f(time[i], y1);

        // compute U_n+1
        u[i+1] = u[i] + dt * ( 0.5 * f(time[i], y1) + 0.5 * f(time[i]+dt, y2));

        // compute t_n+1
        time[i+1] = time[i] + dt;
    }

    /// End of my solution ///
}

int main(int argc, char** argv) {

    double T = 10.0;

    double dt = 0.2;


    // To make some plotting easier, we take the dt parameter in as an optional
    // parameter.
    if (argc == 2) {
        dt = atof(argv[1]);
    }

    const double u0 = 0.;
    std::vector<double> time;
    std::vector<double> u;
    Heun(u,time,u0,dt,T);

    writeToFile("solution.txt", u);
    writeToFile("time.txt",time);

    return 0;
}

// (write your solution here)

/// Start of my solution ///

#include <Eigen/Core>
#include <vector>
#include <iostream>
#include "writer.hpp"
#include <cmath>

// function f=u'=e^{-2t}-2u to avoid writting too many times the same function
double f(double t, double u) {

    return std::exp(-2. * t) - 2. * u;
}

// exact solution of the ODE u'=e^{-2t}-2u
double exact(double t) {

    return t * std::exp(-2.*t);
}

// Heun's method (same as in heun.cpp)
void Heun(std::vector<double> & u, std::vector<double> & time,
          const double & u0, double dt, double T) {

    const unsigned int nsteps = std::round(T/dt);

    u.resize(nsteps+1);
    time.resize(nsteps+1);

    // initialize the first steps
    double y1;
    double y2;

    // initialize the initial values
    u[0] = u0;
    time[0] = 0.;

    for (int i = 0; i < nsteps; i++) {

        // compute the first steps
        y1 = u[i];
        y2 = u[i] + dt * f(time[i], y1);

        // compute U_n+1
        u[i+1] = u[i] + dt * ( 0.5 * f(time[i], y1) + 0.5 * f(time[i]+dt, y2));

        // compute t_n+1
        time[i+1] = time[i] + dt;
    }
}

int main() {

    double T = 10.0;

    // initialize the vector dt and error
    std::vector<double> dt(8);
    std::vector<double> error(8);

    const double u0 = 0.;
    std::vector<double> time;
    std::vector<double> u;

    // store the value in dt, calculate the error and store it
    for (int i = 0; i < 8; i++) {

        dt[i] = std::pow(2, -i);

        Heun(u,time,u0,dt[i],T);

        error[i] = std::abs( u.back() - exact(T) );

    }

    writeToFile("dt.txt", dt);
    writeToFile("error.txt",error);

    return 0;
}

/// End of my solution ///


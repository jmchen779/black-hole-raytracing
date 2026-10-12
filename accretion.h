#include <cmath>
#include <vector>
#include <iostream>
#include <print>
#include "rk45_dromandprince.h"
#include "metric.h"
#include "newton_raphson.h"

class accretion {
public:
    double r1, r2, r3, r4, r5;
 const double d1 = -12715105075.0 / 11282082432.0;
  const double d2 = 0.0;
  const double d3 = 87487479700.0 / 32700410799.0;
  const double d4 = -10690763975.0 / 1880347072.0;
  const double d5 = 701980252875.0 / 199316789632.0;
  const double d6 = -1453857185.0 / 822651844.0;
  const double d7 = 69997945.0 / 29380423.0;
    void compute_dense_rs(const double y,
                        const double y_prev, double h, double k1, double k2, double k3, double k4, double k5, double k6, double k7) {
     r1 = y_prev;
     r2 = y - y_prev;
     r3 = y_prev + h*k1 - y;
     r4 = 2*( y - y_prev) - h*(k1 + k7);
     r5 = h*(d1*k1 + d3*k3 + d4*k4 + d5*k5 + d6*k6 + d7*k7);
  
  }
  //a method for checking if a photon hits the accretion disk, and calculating what the resulting inensity is.
  //The parameters it takes are a and M, the spin and massof the black hole respectively, r_in and r_out which are
  //the inner and outer radii of the accretion disk. Finally, y0 is a std::vector of the initial conditions
   double hit_accretion(double a, double M, double r_in, double r_out, const std::vector<double> &y0) {
    boyer_lindquist_metric metric(a, M);
    //our integrand, set of six coupled differential equations we are solving
    auto dydx = [&metric](double x, const std::vector<double> &y) {
        metric.compute_metric(y[0], y[1]); 
        double u_lt = ((std::sqrt(metric.gamma11 * y[3] * y[3] + metric.gamma22 * y[4] * y[4] + metric.gamma33 * y[5] * y[5])) / metric.alpha);
        return std::vector<double>{
            metric.gamma11 * y[3] / u_lt,
            metric.gamma22 * y[4] / u_lt,
            metric.gamma33 * y[5] / u_lt - metric.beta3,
            -metric.alpha * u_lt * metric.d_alpha_dr + y[5] * metric.d_beta3_dr - (1.0 / (2.0 * u_lt)) * (y[3] * y[3] * metric.d_gamma11_dr + y[4] * y[4] * metric.d_gamma22_dr + y[5] * y[5] * metric.d_gamma33_dr),
            -metric.alpha * u_lt * metric.d_alpha_dth + y[5] * metric.d_beta3_dth - (1.0 / (2.0 * u_lt)) * (y[3] * y[3] * metric.d_gamma11_dth + y[4] * y[4] * metric.d_gamma22_dth + y[5] * y[5] * metric.d_gamma33_dth),
            0
        };
    };
    double x0 = 0.0;
    double x_end = 10000.0;
    //event horizon
    double r_horizon = M + std::sqrt(M*M - a*a);
    
    bool hit_disk_during_integration = false;
    //custom stopping condition: checks if our photon has crossed the equator (means it hit the accretions disk)
    //or if it has gone too close to the black hole or excaped to infinity
    //If any of these three conditions are satisfied then the integrator stops, and we know if our photon hit the disk or not
    auto stop = [r_in, r_out, M, a, r_horizon, &hit_disk_during_integration](double x, const std::vector<double>& y, const std::vector<double>& y0) {
        bool crossed_equator = (y[1] - M_PI/2.0) * (y0[1] - M_PI/2.0) <= 0;
        bool in_disk_radius = y[0] >= r_in && y[0] <= r_out;
        
        if (crossed_equator && in_disk_radius) {
            hit_disk_during_integration = true;
        }
        
        bool hit_horizon = y[0] <= 1.01 * r_horizon;
        bool escaped = y[0] > 500 * M;
        
        return hit_disk_during_integration || hit_horizon || escaped;
    };
    rk45_dormand_prince solver(1e-10, 1e-10);
    solver.integrate(dydx, stop, x0, x_end, y0, false);

    //if we didnt hit the disk, then our intensity is 0
    if (!hit_disk_during_integration){return 0.0;}
    //if we did hit the disk, we need to use dense output and newton raphson iterations to determine the exact time 
    //we hit the disk, as well as the values of our positions and momenta at that time
    if (hit_disk_during_integration){
        //computing dense output coefficients in the theta direction
    compute_dense_rs(
    solver.result[solver.xs.size() - 1][1],
    solver.result[solver.xs.size() - 2][1],
    solver.hs.back(),
    solver.k1[1], solver.k2[1], solver.k3[1],
    solver.k4[1], solver.k5[1], solver.k6[1], solver.k7[1]
     );
     double phi_target =  M_PI / 2.0;
    double phi_prev = solver.result[solver.xs.size() - 2][1];
    double phi_last = solver.result[solver.xs.size() - 1][1];

    //std::print("phi_prev = {}, phi_last = {}, target = {}\n", phi_prev, phi_last, phi_target);
    double x_value = 0.0;

    if (phi_prev < phi_target && phi_last > phi_target) {
        //difference between our interpolated phi value and our desired phi value
    auto F = [&](double x) {
        double phi_interp =
            r1 + x * (r2 + (1 - x) * (r3 + x * (r4 + (1 - x) * r5)));
        return phi_interp - phi_target;
    };
    //derivative of the interpolation polynomial
    auto Fprime = [&](double x) {
        return r2 + r3 * (1 - 2 * x)
            + r4 * (3 * x * x - 4 * x + 1)
            + r5 * (3 * x * x - 2 * x);
    };

    std::vector<double> guesses = {0.2, 0.5, 0.8};
    bool converged = false;
    //go through all guesses and see if newton raphson iterations converge and give us the time our photon crossed disk
    for (double guess : guesses) {
        auto result = newton_raphson(F, Fprime, guess, 1e-10);
        if (std::get<0>(result)) {
            x_value = std::get<1>(result);
            converged = true;
            //std::print("Converged with initial guess {}: theta = {}\n", guess, x_value);
            break;
        } else {
           // std::print("Guess {} failed to converge\n", guess);
        }
    }

    if (converged) {
        double x_interp = solver.xs[solver.xs.size() - 2] + solver.hs.back() * x_value;
       // std::print("Interpolated x = {}, phi ≈ {}\n", x_interp, phi_target);
    } else {
        std::print("Newton–Raphson failed for all guesses.\n");
    }
    } 
    //Now that we have our crossing time, we pplug this crossing time into our interpolation polynomails
    //for each of our positions and momentums
    compute_dense_rs(solver.result.back()[0], solver.result[solver.xs.size()-2][0],
                 solver.hs.back(), solver.k1[0], solver.k2[0], solver.k3[0],
                 solver.k4[0], solver.k5[0], solver.k6[0], solver.k7[0]);
    double r_final = r1 + x_value * (r2 + (1 - x_value) * (r3 + x_value * (r4 + (1 - x_value) * r5)));
    compute_dense_rs(solver.result.back()[2], solver.result[solver.xs.size()-2][2],
                 solver.hs.back(), solver.k1[2], solver.k2[2], solver.k3[2],
                 solver.k4[2], solver.k5[2], solver.k6[2], solver.k7[2]);
    double phi_final = r1 + x_value * (r2 + (1 - x_value) * (r3 + x_value * (r4 + (1 - x_value) * r5)));
    compute_dense_rs(solver.result.back()[3], solver.result[solver.xs.size()-2][3],
                 solver.hs.back(), solver.k1[3], solver.k2[3], solver.k3[3],
                 solver.k4[3], solver.k5[3], solver.k6[3], solver.k7[3]);
    double ur_final = r1 + x_value * (r2 + (1 - x_value) * (r3 + x_value * (r4 + (1 - x_value) * r5)));
    compute_dense_rs(solver.result.back()[4], solver.result[solver.xs.size()-2][4],
                 solver.hs.back(), solver.k1[4], solver.k2[4], solver.k3[4],
                 solver.k4[4], solver.k5[4], solver.k6[4], solver.k7[4]);
    double u_theta_final = r1 + x_value * (r2 + (1 - x_value) * (r3 + x_value * (r4 + (1 - x_value) * r5)));
    compute_dense_rs(solver.result.back()[5], solver.result[solver.xs.size()-2][5],
                 solver.hs.back(), solver.k1[5], solver.k2[5], solver.k3[5],
                 solver.k4[5], solver.k5[5], solver.k6[5], solver.k7[5]);
    double u_phi_final = r1 + x_value * (r2 + (1 - x_value) * (r3 + x_value * (r4 + (1 - x_value) * r5)));

    //  ur_final = -ur_final;
    //  u_theta_final = -u_theta_final;
    // u_phi_final = -u_phi_final;
    //calculate gravitational redshift and doppler boosting effects here
    metric.compute_metric(r_final, M_PI / 2.0);
    double omega = 1.0 / (a + (std::pow(r_final, 3.0 / 2.0) / std::sqrt(M) ));
        double u_lt = ((std::sqrt(metric.gamma11 * ur_final * ur_final + metric.gamma22 * u_theta_final * u_theta_final + metric.gamma33 * u_phi_final *u_phi_final ) ) / metric.alpha);
    double z = (1.0 + omega * u_phi_final / ( (-metric.alpha * metric.alpha  * u_lt + u_phi_final*metric.beta3) )) / std::sqrt( -metric.g_00 - omega*omega*metric.g_33-2*omega*metric.g_03);
    
    return 1.0 / (z*z*z);}
    //return phi_final;
//}
   // else{
   // return 0.0;}

    }
    

};
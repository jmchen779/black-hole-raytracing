#include <cmath>
#include <vector>
#include <print>
#include <fstream>
#include <tuple>
#include <iostream>
#include "metric.h"
#include "rk45_dormandprince2.h"
#include "newton_raphson.h"

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


int main(){
double a = 0.0;
double M = 1.0;
boyer_lindquist_metric metric(a, M);
auto dydx = [&metric] (double x, const std::vector<double> &y) {
metric.compute_metric(y[0], y[1]); 
double u_lt = ((std::sqrt(metric.gamma11 * y[3] * y[3] + metric.gamma22 * y[4] * y[4] + metric.gamma33 * y[5] *y[5] ) ) / metric.alpha);
return std::vector<double>{
    metric.gamma11 * y[3] / u_lt,
    metric.gamma22 * y[4] / u_lt,
    metric.gamma33 * y[5] / u_lt - metric.beta3,
    -metric.alpha * u_lt * metric.d_alpha_dr + y[5] * metric.d_beta3_dr - ( 1.0 / (2.0 * u_lt)) * (y[3] * y[3] *  metric.d_gamma11_dr +  y[4] * y[4] *  metric.d_gamma22_dr + y[5] * y[5] *  metric.d_gamma33_dr),
     -metric.alpha * u_lt * metric.d_alpha_dth + y[5] * metric.d_beta3_dth - ( 1.0 / (2.0 * u_lt)) * (y[3] * y[3] *  metric.d_gamma11_dth +  y[4] * y[4] *  metric.d_gamma22_dth + y[5] * y[5] *  metric.d_gamma33_dth),
     0};
};
metric.compute_metric(10.0*M, M_PI/2.0);
double xi = 0.48857874;
double x0 = 0.0;
double x_end = 100.0;
double phi_end = 2*M_PI;
auto stop = [phi_end, M](double x, const std::vector<double>& y) {
    return (y[2] >= phi_end) || (x > 1e4)|| (y[0] <= 2.01*M);
};
std::vector<double> y0 = {10.0*M, M_PI / 2.0, 0, -std::sqrt(1.0 / metric.gamma11) *std::cos(xi), 0,  std::sqrt(1.0 / metric.gamma33) *std::sin(xi)};
rk45_dormand_prince solver(1e-4, 1e-4);
std::cerr << "About to integrate...\n";
solver.integrate(dydx, stop, x0, x_end, y0, false);
std::cerr << "Finished integration call!\n";
    std::ofstream output_file("e_ring.csv");
for (int i = 0; i < solver.xs.size(); i++) {
     output_file <<  solver.result[i][0] << "," << solver.result[i][2] << std::endl;
    //std::print(output_file, "{},{}\n", rk4.xs[i], rk4.result[i][0]);
  }

  output_file.close();

compute_dense_rs(
    solver.result[solver.xs.size() - 1][2],
    solver.result[solver.xs.size() - 2][2],
    solver.hs.back(),
    solver.k1[2], solver.k2[2], solver.k3[2],
    solver.k4[2], solver.k5[2], solver.k6[2], solver.k7[2]
);

double phi_target = 2 * M_PI;
double phi_prev = solver.result[solver.xs.size() - 2][2];
double phi_last = solver.result[solver.xs.size() - 1][2];

std::print("phi_prev = {}, phi_last = {}, target = {}\n", phi_prev, phi_last, phi_target);
    double x_value = 0.0;

if (phi_prev < phi_target && phi_last > phi_target) {
    auto F = [&](double x) {
        double phi_interp =
            r1 + x * (r2 + (1 - x) * (r3 + x * (r4 + (1 - x) * r5)));
        return phi_interp - phi_target;
    };

    auto Fprime = [&](double x) {
        return r2 + r3 * (1 - 2 * x)
            + r4 * (3 * x * x - 4 * x + 1)
            + r5 * (3 * x * x - 2 * x);
    };

    std::vector<double> guesses = {0.2, 0.5, 0.8};
    bool converged = false;

    for (double guess : guesses) {
        auto result = newton_raphson(F, Fprime, guess, 1e-10);
        if (std::get<0>(result)) {
            x_value = std::get<1>(result);
            converged = true;
            std::print("Converged with initial guess {}: theta = {}\n", guess, x_value);
            break;
        } else {
            std::print("Guess {} failed to converge\n", guess);
        }
    }

    if (converged) {
        double x_interp = solver.xs[solver.xs.size() - 2] + solver.hs.back() * x_value;
        std::print("Interpolated x = {}, phi ≈ {}\n", x_interp, phi_target);
    } else {
        std::print("Newton–Raphson failed for all guesses.\n");
    }
} 
compute_dense_rs(solver.result.back()[0], solver.result[solver.xs.size()-2][0],
                 solver.hs.back(), solver.k1[0], solver.k2[0], solver.k3[0],
                 solver.k4[0], solver.k5[0], solver.k6[0], solver.k7[0]);
double r_obs = r1 + x_value * (r2 + (1 - x_value) * (r3 + x_value * (r4 + (1 - x_value) * r5)));

std::print("Error: {}\n", std::abs(r_obs - 10.0*M));
return 0;
}


//Error: 1.6412281066635614e-05


#pragma once
#include <fstream>
#include <vector>
#include <cmath>
#include <print>

class rk45_dormand_prince {
public:
  rk45_dormand_prince(double tolerance_abs, double tolerance_rel) {
    // constructor, initialize the data members
    atol = tolerance_abs;
    rtol = tolerance_rel;
  }

  template <typename F, typename U>
  std::vector<double> integrate(const F &f, const U &stop_condition,  double x0, double x_end,
                                const std::vector<double> &y0,
                                bool dense_output = false,
                                const std::vector<double> &x_out = {}) {
    // resize the temporary arrays according to the given y0
    n_eq = y0.size();
    y_temp.resize(n_eq);
    y_err.resize(n_eq);

    k1.resize(n_eq);
    k2.resize(n_eq);
    k3.resize(n_eq);
    k4.resize(n_eq);
    k5.resize(n_eq);
    k6.resize(n_eq);
    k7.resize(n_eq);

    r1.resize(n_eq);
    r2.resize(n_eq);
    r3.resize(n_eq);
    r4.resize(n_eq);
    r5.resize(n_eq);

    // clear the arrays so that we start fresh every time we call this function
    hs.clear();
    xs.clear();
    result.clear();

    // initial conditions
    double h = 0.1; // arbitrary initial step size
    double x = x0;
    std::vector<double> y = y0;
    double err = 1.0;
    // err_prev is for you to implement the more advanced step size control
    double err_prev = 1.0;
    
    // Only push the initial condition into the results if it is not a dense output
    if (!dense_output) {  
      xs.push_back(x);
      result.push_back(y);
    }

    auto y_prev = y; // This is used for computing dense output coefficients
    std::vector<double> y_out(n_eq); // This is used to store the value for dense output
    int xout_idx = 0; // This is the current index for the x_out array

    while (!stop_condition(x,y)) {
      std::vector<double> y_trial = step(f, h, x, y);
      // Compute the rescaled error
      err = error(y_trial);
      // If err is fortuitously too small, set it to some lower bound
      err = std::max(err, 1.0e-12);

      // If the error is below 1, we accept the step by pushing the new x, y,
      // and h to the arrays, and advancing x by the current value of h
      if (err < 1.0) {
        double x_prev = x;
         x += h;
         y = y_trial;
        if (dense_output) {
          compute_dense_rs(y, y_prev, h);
          for (int i = xout_idx; i < x_out.size(); i++) {
    if (x_out[i] > x_prev && x_out[i] < x) {
        xs.push_back(x_out[i]);
        theta.push_back((x_out[i] - x_prev) / h);
        xout_idx = i;

        std::vector<double> row; 
        double th = theta.back();
        for (int j = 0; j < n_eq; j++) {
            row.push_back(r1[j] + th * (r2[j] + (1 - th) * (r3[j] + th * (r4[j] + (1 - th) * r5[j]))));
        }
        result.push_back(std::move(row));
    }
}
         // TODO: Carry out dense output here, in two steps. First, compute the
          // coefficients r1, r2, r3, r4, r5. You can use the already-defined
          // function compute_dense_rs to do this, and implement the "TODO" part
          // in the function. Then, go through the x_out array, find the x_out
          // values that are between x and x + h, interpolate the corresponding
          // y values using the coefficients. Push the x values to the xs array
          // and y values to the result array.

          // Note that you probably don't want to loop through the whole x_out
          // array every time you compute a step. You can keep track of the
          // index of the last x_out value you used as xout_idx, and start from
          // there next time. You can assume that the x_out array is already
          // sorted in ascending order.
        } else {
          // This is what we do when not using dense output
          hs.push_back(h);
          y = y_trial;
          xs.push_back(x);
          result.push_back(y);
        }
       
        y_prev = y;
      }

      // Adjust h regardless of whether the step was accepted or rejected
      const double safety = 0.9;
      h = h * safety * std::pow(err, -0.7/5.0) * std::pow(err_prev, 0.4/5.0);

      if ( h > hmax){
        h = hmax;
      }
      else if ( h < hmin){
        h = hmin;
      }
      if (err < 1.0){
      err_prev = err;}
      // TODO: Implement this part where we adjust h according to the error.
      // Note that using the formula given in the lecture, it will increase h if
      // err is small and decrease it if err is large. Do not let h go below
      // hmin or above hmax. You can either use the simple step size control or
      // the more advanced one that uses err_prev. If you do use err_prev, make
      // sure to update it accordingly.

      // Uncomment the following line to see the step size and error at each step
      // std::println("x = {}, h = {}, err = {}", x, h, err);
    }
    return y;
  }

  template <typename F>
  std::vector<double> step(const F& f, double h, double x, const std::vector<double> &y) {
    // Compute the next step in y, given x and y of the current step
    std::vector<double> y_next(n_eq);
     k1 = f(x, y);
    for (int i = 0; i < n_eq; i++){
      y_temp[i] = y[i] + h * k1[i]*a21;
    }
    k2 = f(x + h*c2, y_temp);
    for (int i = 0; i < n_eq; i++){
      y_temp[i] = y[i] + h * k1[i]*a31 + h*k2[i]* a32;
    }
    k3 = f(x + h*c3, y_temp);
    for (int i = 0; i < n_eq; i++){
      y_temp[i] = y[i] + h * k1[i]*a41 + h*k2[i]* a42 + h*k3[i]* a43;
    }
    k4 = f(x + h*c4, y_temp);
    for (int i = 0; i < n_eq; i++){
      y_temp[i] = y[i] + h * k1[i]*a51 + h*k2[i]* a52 + h*k3[i]* a53 + h*k4[i]*a54;
    }
    k5 = f(x + h*c5, y_temp);
    for (int i = 0; i < n_eq; i++){
      y_temp[i] = y[i] + h * k1[i]*a61 + h*k2[i]* a62 + h*k3[i]* a63 + h*k4[i] * a64 + h*k5[i]*a65;
    }
    
    k6 = f(x + h, y_temp);
    for (int i = 0; i < n_eq; i++){
     y_temp[i] = y[i] + h * k1[i]*a71 + h*k2[i]* a72 + h*k3[i]* a73 + h*k4[i] * a74 + h*k5[i]*a75 + h*k6[i] * a76;
      //y_temp[i] = y[i] + h * (k1[i]*5179.0/57600.0 + 
                                // k3[i]*7571.0/16695.0 + 
                                // k4[i]*393.0/640.0 + 
                                // k5[i]*(-92097.0/339200.0) + 
                                // k6[i]*187.0/2100.0);
     }
    k7 = f(x + h, y_temp);

    for (int i = 0; i < n_eq; i++){
      y_next[i] = y[i] + h*(k1[i]*a71 + k2[i]*a72 + k3[i]*a73 + k4[i]*a74+k5[i]*a75+k6[i]*a76);
      y_err[i] = std::abs((h*(k1[i]*e1 + k2[i]*e2 + k3[i]*e3 + k4[i]*e4+k5[i]*e5+k6[i]*e6 + k7[i] * e7) ));
    }
    // TODO: Compute 5th order y_next using the Dormand-Prince method. Then,
    // compute a y_err vector as the difference between the 5th order and 4th
    // order estimates of DP45. Note that the coefficients e1, e2, ..., e7 are
    // already defined for you, which are the differences between b_i
    // (coefficients of the 5th order solution) and b_i^* (coefficients of the
    // 4th order solution).

    return y_next;
  }

  double error(const std::vector<double> &y) {
    // TODO: Compute the rescaled scalar error from the y_err vector computed in
    // the step function, and return it.
    double sum = 0;
    for (int i = 0; i < n_eq; i++){
      sum += (y_err[i] / (atol + rtol * std::abs(y[i]))) * (y_err[i] / (atol + rtol * std::abs(y[i])));
    }
    return std::sqrt((1.0 / n_eq) * sum);
  }

  void compute_dense_rs(const std::vector<double> &y,
                        const std::vector<double> &y_prev, double h) {
    // TODO: Compute the coefficients r1, r2, r3, r4, r5 here
    for (int i = 0; i < n_eq; i++){
    r1[i] = y_prev[i];
    r2[i] = y[i] - y_prev[i];
    r3[i] = y_prev[i] + h*k1[i] - y[i];
    r4[i] = 2*( y[i] - y_prev[i]) - h*(k1[i] + k7[i]);
    r5[i] = h*(d1*k1[i] + d3*k3[i] + d4*k4[i] + d5*k5[i] + d6*k6[i] + d7*k7[i]);
  }
  }

  int n_eq; // number of equations
  double atol, rtol; // absolute and relative tolerances

  // We impose a minimum and maximum step size
  const double hmin = 1.0e-10;
  const double hmax = 1.0;

  // These are temporary variables used to store the coefficients. You are allowed to
  // define additional temporary variables if you need them.
  std::vector<double> k1, k2, k3, k4, k5, k6, k7, y_temp, y_err, theta, row;

  // vectors that store the results
  std::vector<double> hs;
  std::vector<double> xs;
  std::vector<std::vector<double>> result;

  // These are temporary variables used to store the coefficients
  // used in dense output
  std::vector<double> r1, r2, r3, r4, r5;

  // c1 is zero, c6 and c7 are 1.0
  const double c2 = 1.0 / 5.0;
  const double c3 = 3.0 / 10.0;
  const double c4 = 4.0 / 5.0;
  const double c5 = 8.0 / 9.0;

  const double a21 = 1.0 / 5.0;
  const double a31 = 3.0 / 40.0;
  const double a32 = 9.0 / 40.0;
  const double a41 = 44.0 / 45.0;
  const double a42 = -56.0 / 15.0;
  const double a43 = 32.0 / 9.0;
  const double a51 = 19372.0 / 6561.0;
  const double a52 = -25360.0 / 2187.0;
  const double a53 = 64448.0 / 6561.0;
  const double a54 = -212.0 / 729.0;
  const double a61 = 9017.0 / 3168.0;
  const double a62 = -355.0 / 33.0;
  const double a63 = 46732.0 / 5247.0;
  const double a64 = 49.0 / 176.0;
  const double a65 = -5103.0 / 18656.0;

  // Note that a71, a72, a73, a74, a75, a76 are essentially the coefficients of the 5th order solution
  const double a71 = 35.0 / 384.0;
  const double a72 = 0.0;
  const double a73 = 500.0 / 1113.0;
  const double a74 = 125.0 / 192.0;
  const double a75 = -2187.0 / 6784.0;
  const double a76 = 11.0 / 84.0;

  // These coefficients are used to estimate the error in the solution. They are essentially
  // the coefficients of the 5th order solution minus the 4th order solution, i.e. b_i - b_i^*.
  const double e1 = 71.0 / 57600.0;
  const double e2 = 0.0;
  const double e3 = -71.0 / 16695.0;
  const double e4 = 71.0 / 1920.0;
  const double e5 = -17253.0 / 339200.0;
  const double e6 = 22.0 / 525.0;
  const double e7 = -1.0 / 40.0;

  // These are the coefficients for dense output
  const double d1 = -12715105075.0 / 11282082432.0;
  const double d2 = 0.0;
  const double d3 = 87487479700.0 / 32700410799.0;
  const double d4 = -10690763975.0 / 1880347072.0;
  const double d5 = 701980252875.0 / 199316789632.0;
  const double d6 = -1453857185.0 / 822651844.0;
  const double d7 = 69997945.0 / 29380423.0;
};

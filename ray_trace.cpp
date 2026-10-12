#include <cmath>
#include <vector>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <thread>
#include "accretion.h"
#include "metric.h"


int main(){
    int n = 100;
    int m = 200;
    std::vector<std::vector<double>> picture(n, std::vector<double>(m));
    double a = 0.99;
    double M = 1;
    double D = 500.0;
    double theta0 = 85.0 * M_PI / 180.0;  
    //inner and outer radii of accretion disk
    double r_in = 5.0;
    double r_out = 20.0;
    double expected_disk_size = 2.0 * r_out;
    double screen_width = expected_disk_size * 2.0;
    double screen_height = screen_width * n / m;

        //loop through each pixel and compute intensity
        for (int i = 0; i < n; ++i){
            for (int j = 0; j < m; ++j){
                double alpha_screen = -(screen_width * (j - m/2.0) / m) / M;
                double beta_screen = -(screen_height * (i - n/2.0) / n) / M;

                double u_t = -1.0; 
                double u_phi = -alpha_screen * sin(theta0) * u_t;
                double u_theta = -beta_screen * u_t;
                
                if (std::abs(beta_screen) < 1e-10) {
                    u_theta = -0.001;
                }
                if (std::abs(alpha_screen) < 1e-10) {
                    u_phi = -0.001;
                }
                
                boyer_lindquist_metric metric(a, M);
                metric.compute_metric(D, theta0);
                
        
                double gu_tt = -1.0 / (metric.alpha * metric.alpha);
                double gu_t3 = metric.beta3 / (metric.alpha * metric.alpha);
                double gu_rr = metric.gamma11;
                double gu_22 = metric.gamma22;
                double gu_33 = metric.gamma33 - (metric.beta3 / metric.alpha) * (metric.beta3 / metric.alpha);
                
                double expression = gu_tt * u_t * u_t + 2.0 * gu_t3 * u_t * u_phi + 
                                   gu_22 * u_theta * u_theta + gu_33 * u_phi * u_phi;

                if (expression > 0) {
                    picture[i][j] = 0.0;
                    continue;
                }
                
                double u_r = -sqrt(-expression / gu_rr);
                
                std::vector<double> y0 = {D, theta0, 0.0, u_r, u_theta, u_phi};
                
                accretion checker;
                double result = checker.hit_accretion(a, M, r_in, r_out, y0);
                
                if (result > 0) {
                    picture[i][j] = result;
                } else {
                    
                    picture[i][j] = 0.0;
                }
            }
        }
  


  
    double min_val = picture[0][0], max_val = picture[0][0];
    for (auto &row : picture) {
        for (double val : row) {
            if (val < min_val) min_val = val;
            if (val > max_val) max_val = val;
        }
    }

    std::ofstream out("bh_image.pgm");
    out << "P2\n" << m << " " << n << "\n255\n";  

    for (auto &row : picture) {
        for (double val : row) {
            int pixel = static_cast<int>(255 * (val - min_val) / (max_val - min_val) + 0.5);
            out << pixel << " ";
        }
        out << "\n";
    }

    return 0;
}
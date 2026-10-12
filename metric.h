#pragma once
#include <cmath>

// a simple class for computing elements of the metric tensor
class boyer_lindquist_metric {
public:
    boyer_lindquist_metric(double a0, double M0) {
        a = a0;
        M = M0;
    }

    void compute_metric(double r, double th) {
    double sin_th = std::sin(th);
    double cos_th = std::cos(th);
    double sin2_th = sin_th * sin_th;
    double cos2_th = cos_th * cos_th;

    double rho2 = r*r + a*a*cos2_th;
    double delta = r*r + a*a - 2*M*r;
    double sigma = (r*r + a*a)*(r*r + a*a) - a*a * delta * sin2_th;

    alpha = std::sqrt(rho2 * delta / sigma);
    beta3 = - (2.0 * M * a * r) / sigma; 
    
    gamma11 = delta / rho2;      
    gamma22 = 1.0 / rho2;         
    gamma33 = rho2 / (sigma * sin2_th); 

    g_00 = (2.0 * M * r / rho2) - 1.0;                   
    g_03 = -2.0 * M * a * r * sin2_th / rho2;             
    g_11 = 1.0 / gamma11; 
    g_22 = 1.0 / gamma22; 
    g_33 = 1.0 / gamma33; 

    double d_rho2_dr = 2.0 * r;
    double d_rho2_dth = -2.0 * a*a * sin_th * cos_th;
    
    double d_delta_dr = 2.0 * (r - M);
    double d_delta_dth = 0.0;
    
    double d_sigma_dr = 4.0 * r * (r*r + a*a) - a*a * sin2_th * d_delta_dr;
    double d_sigma_dth = -2.0 * a*a * delta * sin_th * cos_th;

    d_alpha_dr = 0.5 * alpha * (d_rho2_dr/rho2 + d_delta_dr/delta - d_sigma_dr/sigma);
    d_alpha_dth = 0.5 * alpha * (d_rho2_dth/rho2 - d_sigma_dth/sigma);

    d_gamma11_dr = (d_delta_dr * rho2 - delta * d_rho2_dr) / (rho2 * rho2);
    d_gamma11_dth = -delta * d_rho2_dth / (rho2 * rho2);
    
    d_gamma22_dr = -d_rho2_dr / (rho2 * rho2);
    d_gamma22_dth = -d_rho2_dth / (rho2 * rho2);
    
    double term_r = d_rho2_dr * sigma - rho2 * d_sigma_dr;
    d_gamma33_dr = term_r / (sigma * sigma * sin2_th);
    
    double term_th = d_rho2_dth * sigma - rho2 * d_sigma_dth;
    d_gamma33_dth = (term_th * sin2_th - 2.0 * rho2 * sigma * sin_th * cos_th) / (sigma * sigma * sin2_th * sin2_th);

    d_beta3_dr = -2.0 * M * a * (sigma - r * d_sigma_dr) / (sigma * sigma);
    d_beta3_dth = 2.0 * M * a * r * d_sigma_dth / (sigma * sigma);
}


    double a, M;
    double alpha, beta3;
    double gamma11, gamma22, gamma33; // components of upper gamma^ij
    double g_00, g_03, g_11, g_22, g_33; // components of lower g_\mu\nu
    double d_alpha_dr, d_beta3_dr, d_gamma11_dr, d_gamma22_dr, d_gamma33_dr;
    double d_alpha_dth, d_beta3_dth, d_gamma11_dth, d_gamma22_dth, d_gamma33_dth;
};
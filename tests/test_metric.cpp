#include <iostream>
#include "metric.h"

int main() {
    double a = 0.99;
    double M = 4.0;
    double D = 500.0;
    double theta0 = 85.0 * M_PI / 180.0;

    boyer_lindquist_metric metric(a, M);
    metric.compute_metric(D, theta0);

    std::cout << "Metric at observer (r=" << D << ", theta=" << theta0 << "):\n";
    std::cout << "alpha=" << metric.alpha << "\n";
    std::cout << "beta3=" << metric.beta3 << "\n";
    std::cout << "gamma11=" << metric.gamma11 << " gamma22=" << metric.gamma22 
              << " gamma33=" << metric.gamma33 << "\n";

    double gu_tt = -1.0 / (metric.alpha * metric.alpha);
    double gu_t3 = metric.beta3 / (metric.alpha * metric.alpha);
    double gu_rr = metric.gamma11;
    double gu_22 = metric.gamma22;
    double gu_33 = metric.gamma33 - (metric.beta3 / metric.alpha)*(metric.beta3 / metric.alpha);

    std::cout << "gu_tt=" << gu_tt << " gu_t3=" << gu_t3 << " gu_rr=" << gu_rr 
              << " gu_22=" << gu_22 << " gu_33=" << gu_33 << "\n";
}


#include <print>
#include <cmath>
#include <tuple>

template <typename F>
std::tuple<bool, double> fixed_point_iteration(const F& g, double x0, double tolerance){
    auto OG = [&g](double x) {return g(x);};
    auto G = [&g](double x) {return g(x) + x;};
    auto Ge = [&g](double x) { return -(g(x) - x); };
    int counter = 0;
    int counter2 = 0;
    double xNext = x0;
    const double initial = x0;
    while (counter < 100){
       //std::print("{}\n", x0);
        xNext = G(x0);
        std::print("{}\n", xNext);
        counter++;
        if (std::abs(OG(xNext)) < tolerance){
            //std::print("{}\n", counter);
            return std::make_tuple(true,xNext);
        }
        x0 = xNext;
    }
    x0 = initial;
     while (counter2 < 100){
       //std::print("{}\n", x0);
        xNext = Ge(x0);
       // std::print("{}\n", xNext);
        counter2++;
        if (std::abs(OG(xNext)) < tolerance){
            //std::print("{}\n", counter);
            return std::make_tuple(true,xNext);
        }
        x0 = xNext;
    }
    return std::make_tuple(false, xNext);
}
template <typename F>
std::tuple<bool, double> bisection(const F& f, double a, double b, double tolerance){
    auto Ff = [&f](double x){ return f(x);};
    double midpoint = 0;
    int counter = 0;
    while (Ff(a)*Ff(b) < 0 && counter < 100){
        midpoint = (a+b)/2;
        if (std::abs(Ff(midpoint)) < tolerance){
            return std::make_tuple(true, midpoint);
        }
        counter++;
        if ((Ff(a) > 0 && Ff(midpoint) > 0) || (Ff(a) < 0 && Ff(midpoint) < 0)){
            a = midpoint;
        }
        else{
            b = midpoint;
        }
    }
    return std::make_tuple(false, midpoint);
}
template <typename F, typename Fprime>
std::tuple<bool, double> newton_raphson(const F& f, const Fprime& fprime, double x0, double tolerance){
    auto Ff = [&f](double x){ return f(x);};
    auto Fpprime = [&fprime](double x){ return fprime(x);};
    double xNext = x0;
    int counter = 0;
    const double initial = x0;
    while (counter < 1000){
        if (Fpprime(x0) == 0){
            return std::make_tuple(false, xNext);
        }
        xNext = x0 - Ff(x0)/ Fpprime(x0);
       // std::print("{}\n", xNext);
        //if (std::abs(xNext - initial  ) > 50){
           // return std::make_tuple(false, x0);
       // }
        if (std::abs(Ff(xNext)) < tolerance){
            return std::make_tuple(true, xNext);
        }
        counter++;
        x0 = xNext;
    }
    return std::tuple(false, xNext);
}

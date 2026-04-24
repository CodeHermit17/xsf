#pragma once

#include <cmath>
#include <limits>

#include "gamma.h"

namespace xsf {
namespace detail {
    
    // Returns the shape parameter a of a gamma distribution with rate b
    // such that the CDF at x equals p.
    // ADD: additional input check for x and p
    inline double gdtria(double p, double b, double x) { 
        if ((x <= 0) || (p < 0.0) || (p > 1.0)) {
            return std::numeric_limits<double>::quiet_NaN();
        }
        if ((b == 0) && (p == 0)) {
            if (std::isinf(x) && (x > 0)) {
                return std::numeric_limits<double>::quiet_NaN();
            }
            return 0.0;
        }
        return gammaincinv(b, p) / x;
    }

    // Returns the x value of a gamma distribution with shape b and rate a
    // such that the CDF at x equals p.
    // ADD: additional input check for p
    inline double gdtrix(double a, double b, double p) { 
        if ((p < 0.0) || (p > 1.0)) {
            return std::numeric_limits<double>::quiet_NaN();
        }        
        if ((a == 0) && (b == 0)) {
            return std::numeric_limits<double>::quiet_NaN();
        }
        // if a or b is positive infinite, return NaN
        if ((std::isinf(a) || std::isinf(b)) && (a >= 0 && b >= 0)) {
            return std::numeric_limits<double>::quiet_NaN();
        }
        return gammaincinv(b, p) / a;
    }

} // namespace detail

inline double gdtria(double p, double b, double x) { return detail::gdtria(p, b, x); }
inline double gdtrix(double a, double b, double p) { return detail::gdtrix(a, b, p);}

inline float gdtria(float p, float b, float x) {
    return static_cast<float>(detail::gdtria(static_cast<double>(p), static_cast<double>(b), static_cast<double>(x))
    );
}
inline float gdtrix(float a, float b, float p){
    return static_cast<float>(detail::gdtrix(static_cast<double>(a), static_cast<double>(b), static_cast<double>(p))
    );
}   

} //namespace xsf
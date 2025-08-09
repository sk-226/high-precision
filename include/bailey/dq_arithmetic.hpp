#pragma once

#include <iostream>
#include <string>
#include <cstring>
#include <cmath>
#include <Eigen/Sparse>

// DQ (Quad-Double) precision arithmetic using Bailey's DQFUN library
extern "C" {
    void dqadd_(const long double* a, const long double* b, long double* c);      // c = a + b
    void dqsub_(const long double* a, const long double* b, long double* c);      // c = a - b
    void dqmul_(const long double* a, const long double* b, long double* c);      // c = a * b
    void dqdiv_(const long double* a, const long double* b, long double* c);      // c = a / b
    void dqdqd_(const double* d, long double* a);                                  // a = (double)d
    void dqsqrt_(const long double* a, long double* b);                            // b = sqrt(a)
    void dqtoqd_(const long double* a, int* n, char* c, int cl);
}

namespace bailey {

struct DQNumber {
    long double dq[2] = {0.0L, 0.0L};
    
    DQNumber() = default;
    DQNumber(double val) { dqdqd_(&val, dq); }
};

// Basic Arithmetic Operators
inline DQNumber operator+(const DQNumber& a, const DQNumber& b) { 
    DQNumber r; dqadd_(a.dq, b.dq, r.dq); return r; 
}

inline DQNumber operator-(const DQNumber& a, const DQNumber& b) { 
    DQNumber r; dqsub_(a.dq, b.dq, r.dq); return r; 
}

inline DQNumber operator*(const DQNumber& a, const DQNumber& b) { 
    DQNumber r; dqmul_(a.dq, b.dq, r.dq); return r; 
}

inline DQNumber operator/(const DQNumber& a, const DQNumber& b) { 
    DQNumber r; dqdiv_(a.dq, b.dq, r.dq); return r; 
}

// Assignment Operators
inline DQNumber& operator+=(DQNumber& a, const DQNumber& b) { 
    a = a + b; return a; 
}

inline DQNumber& operator-=(DQNumber& a, const DQNumber& b) { 
    a = a - b; return a; 
}

inline DQNumber& operator*=(DQNumber& a, const DQNumber& b) { 
    a = a * b; return a; 
}

inline DQNumber& operator/=(DQNumber& a, const DQNumber& b) { 
    a = a / b; return a; 
}

// Mathematical Functions
inline DQNumber sqrt(const DQNumber& a) { 
    DQNumber r; dqsqrt_(a.dq, r.dq); return r; 
}

// Type Conversion (avoid narrowing to double for precision-sensitive output)
inline std::string to_string(const DQNumber& a, int digits=64) {
    char s[128] = {0};
    dqtoqd_(a.dq, &digits, s, sizeof(s));
    return std::string(s, strnlen(s, sizeof(s)));
}

inline double to_double(const DQNumber& a) {
    char s[128] = {0};
    int d = 64;
    dqtoqd_(a.dq, &d, s, sizeof(s));
    try {
        return std::stod(s);
    } catch (...) {
        return a.dq[0];
    }
}

// Stream Output
inline std::ostream& operator<<(std::ostream& os, const DQNumber& dq) {
    char s[128] = {0}; 
    int digits = 64;
    dqtoqd_(dq.dq, &digits, s, sizeof(s));
    os << s;
    return os;
}

} // namespace bailey

// Eigen Integration
namespace Eigen {
    template<> struct NumTraits<bailey::DQNumber> : GenericNumTraits<bailey::DQNumber> {
        typedef bailey::DQNumber Real; 
        typedef bailey::DQNumber NonInteger; 
        typedef bailey::DQNumber Nested;
        enum { 
            IsComplex = 0, 
            IsInteger = 0, 
            IsSigned = 1, 
            RequireInitialization = 1, 
            ReadCost = 4, 
            AddCost = 32, 
            MulCost = 64 
        };
    };
}
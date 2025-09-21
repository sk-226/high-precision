#pragma once

#include <iostream>
#include <string>
#include <string_view>
#include <cstring>
#include <cmath>
#include <Eigen/Sparse>

// DD (Double-Double) precision arithmetic using Bailey's DDFUN library
extern "C" {
    void ddadd_(const double* a, const double* b, double* c);      // c = a + b
    void ddsub_(const double* a, const double* b, double* c);      // c = a - b
    void ddmul_(const double* a, const double* b, double* c);      // c = a * b
    void dddiv_(const double* a, const double* b, double* c);      // c = a / b
    void dddqd_(const double* d, double* a);                       // a = (double)d
    void ddsqrt_(const double* a, double* b);                      // b = sqrt(a)
    void ddabs_(const double* a, double* b);                       // b = abs(a)
    void dd_to_string(const double* a, int* n, char* c, int cl);
    void ddfromstr_(const char* s, int len, double* a);            // a = parse(s)
    void ddcpr_(const double* a, const double* b, int* ic);        // compare(a,b)
    void dd_read_line(const char* s, int len, double* a);          // parse via ddread
}

namespace bailey {

struct DDNumber {
    double dd[2] = {0.0, 0.0};
    
    DDNumber() = default;
    DDNumber(double val) { dddqd_(&val, dd); }
};

// Basic Arithmetic Operators
inline DDNumber operator+(const DDNumber& a, const DDNumber& b) { 
    DDNumber r; ddadd_(a.dd, b.dd, r.dd); return r; 
}

inline DDNumber operator-(const DDNumber& a, const DDNumber& b) { 
    DDNumber r; ddsub_(a.dd, b.dd, r.dd); return r; 
}

inline DDNumber operator*(const DDNumber& a, const DDNumber& b) { 
    DDNumber r; ddmul_(a.dd, b.dd, r.dd); return r; 
}

inline DDNumber operator/(const DDNumber& a, const DDNumber& b) { 
    DDNumber r; dddiv_(a.dd, b.dd, r.dd); return r; 
}

// Assignment Operators
inline DDNumber& operator+=(DDNumber& a, const DDNumber& b) { 
    a = a + b; return a; 
}

inline DDNumber& operator-=(DDNumber& a, const DDNumber& b) { 
    a = a - b; return a; 
}

inline DDNumber& operator*=(DDNumber& a, const DDNumber& b) { 
    a = a * b; return a; 
}

inline DDNumber& operator/=(DDNumber& a, const DDNumber& b) { 
    a = a / b; return a; 
}

// Mathematical Functions
inline DDNumber sqrt(const DDNumber& a) { 
    DDNumber r; ddsqrt_(a.dd, r.dd); return r; 
}

inline DDNumber abs(const DDNumber& a) {
    DDNumber r; ddabs_(a.dd, r.dd); return r;
}

// Comparison helpers (exact ordering using Bailey routines)
inline int compare(const DDNumber& a, const DDNumber& b) {
    int ic = 0;
    ddcpr_(a.dd, b.dd, &ic);
    return ic;
}

inline bool operator==(const DDNumber& a, const DDNumber& b) {
    return compare(a, b) == 0;
}

inline bool operator!=(const DDNumber& a, const DDNumber& b) {
    return compare(a, b) != 0;
}

inline bool operator<(const DDNumber& a, const DDNumber& b) {
    return compare(a, b) < 0;
}

inline bool operator<=(const DDNumber& a, const DDNumber& b) {
    return compare(a, b) <= 0;
}

inline bool operator>(const DDNumber& a, const DDNumber& b) {
    return compare(a, b) > 0;
}

inline bool operator>=(const DDNumber& a, const DDNumber& b) {
    return compare(a, b) >= 0;
}

// Type Conversion 
inline std::string to_string(const DDNumber& a, int digits=32) {
    char s[80] = {0};
    dd_to_string(a.dd, &digits, s, sizeof(s));
    return std::string(s, strnlen(s, sizeof(s)));
}

inline double to_double(const DDNumber& a) {
    char s[80] = {0};
    int d = 32;
    dd_to_string(a.dd, &d, s, sizeof(s));
    try {
        return std::stod(s);
    } catch (...) {
        return a.dd[0];
    }
}

// Stream Output
inline std::ostream& operator<<(std::ostream& os, const DDNumber& d) {
    char s[80] = {0};
    int digits = 32;
    dd_to_string(d.dd, &digits, s, sizeof(s));
    os << s;
    return os;
}

} // namespace bailey

// Eigen Integration
namespace Eigen {
    template<> struct NumTraits<bailey::DDNumber> : GenericNumTraits<bailey::DDNumber> {
        using Real = bailey::DDNumber; 
        using NonInteger = bailey::DDNumber; 
        using Nested = bailey::DDNumber;
        // Notes: the DD mantissa holds ~106 bits (~31.9 decimal digits).
        // Bailey's docs describe "about 30 digits" to stay conservative;
        // we report 31 here to match Eigen's digits10 semantics.

        static inline int digits10() { return 31; }
        static inline int digits() { return 106; }
        enum { 
            IsComplex = 0, 
            IsInteger = 0, 
            IsSigned = 1, 
            RequireInitialization = 1, 
            ReadCost = 2, 
            AddCost = 16, 
            MulCost = 32 
        };
    };
} // namespace Eigen

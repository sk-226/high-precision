#pragma once

#include <iostream>
#include <string>
#include <cmath>
#include <iomanip>
#include <cstring>
#include <algorithm>
#include <Eigen/Sparse>
#include <Eigen/Core>

// ==============================================================================
//  Bailey QX高精度算術ライブラリとの連携のためのQXNumber型定義
//  公式仕様: QXFUN uses single real(qxknd) values (~33 digit precision)
// ==============================================================================

extern "C" {
    void qxadd_(const long double* a, const long double* b, long double* c);      // c = a + b
    void qxsub_(const long double* a, const long double* b, long double* c);      // c = a - b  
    void qxmul_(const long double* a, const long double* b, long double* c);      // c = a * b
    void qxdiv_(const long double* a, const long double* b, long double* c);      // c = a / b
    void qxdqd_(const double* d, long double* a);                                 // a = (double)d
    void qxsqrt_(const long double* a, long double* b);                           // b = sqrt(a)
    void qxtoqd_(const long double* a, int* n, char* c, int cl);
}

namespace bailey {

/// QX (Extended Quad) precision number (~33 decimal digits)
/// Based on Bailey's QXFUN library using extended precision arithmetic
struct QXNumber {
    long double qx = 0.0L;  // Use long double for better precision than double
    
    QXNumber() = default;
    QXNumber(double val) : qx(static_cast<long double>(val)) {}  // Remove explicit to allow implicit conversion
    QXNumber(long double val) : qx(val) {}
    QXNumber(int val) : qx(static_cast<long double>(val)) {}     // Add int constructor for Eigen
    
    // Copy constructor and assignment
    QXNumber(const QXNumber& other) : qx(other.qx) {}
    QXNumber& operator=(const QXNumber& other) {
        if (this != &other) {
            qx = other.qx;
        }
        return *this;
    }
    
    // Direct access to long double for Bailey Fortran interface
    const long double* get_qx_ptr() const {
        return &qx;
    }
    
    long double* get_qx_ptr() {
        return &qx;
    }
};

// --- Basic Arithmetic Operators ---
inline QXNumber operator+(const QXNumber& a, const QXNumber& b) { 
    QXNumber result;
    qxadd_(a.get_qx_ptr(), b.get_qx_ptr(), result.get_qx_ptr());
    return result; 
}

inline QXNumber operator-(const QXNumber& a, const QXNumber& b) { 
    QXNumber result;
    qxsub_(a.get_qx_ptr(), b.get_qx_ptr(), result.get_qx_ptr());
    return result; 
}

inline QXNumber operator*(const QXNumber& a, const QXNumber& b) { 
    QXNumber result;
    qxmul_(a.get_qx_ptr(), b.get_qx_ptr(), result.get_qx_ptr());
    return result; 
}

inline QXNumber operator/(const QXNumber& a, const QXNumber& b) { 
    QXNumber result;
    qxdiv_(a.get_qx_ptr(), b.get_qx_ptr(), result.get_qx_ptr());
    return result; 
}

// --- Assignment Operators ---
inline QXNumber& operator+=(QXNumber& a, const QXNumber& b) { 
    a = a + b; 
    return a; 
}

inline QXNumber& operator-=(QXNumber& a, const QXNumber& b) { 
    a = a - b; 
    return a; 
}

inline QXNumber& operator*=(QXNumber& a, const QXNumber& b) { 
    a = a * b; 
    return a; 
}

inline QXNumber& operator/=(QXNumber& a, const QXNumber& b) { 
    a = a / b; 
    return a; 
}

// --- Comparison Operators (Essential for CG algorithm) ---
inline bool operator==(const QXNumber& a, const QXNumber& b) {
    // Use high-precision epsilon for QX (~33 digits)
    constexpr long double QX_EPSILON = 1e-31L;
    long double diff = std::abs(a.qx - b.qx);
    long double max_val = std::max(std::abs(a.qx), std::abs(b.qx));
    
    if (max_val < 1e-15L) {
        return diff < 1e-15L; // Absolute comparison for near-zero
    }
    return diff / max_val < QX_EPSILON;
}

inline bool operator!=(const QXNumber& a, const QXNumber& b) {
    return !(a == b);
}

inline bool operator<(const QXNumber& a, const QXNumber& b) {
    constexpr long double QX_EPSILON = 1e-31L;
    return (b.qx - a.qx) > QX_EPSILON * std::max(std::abs(a.qx), std::abs(b.qx));
}

inline bool operator<=(const QXNumber& a, const QXNumber& b) {
    return (a < b) || (a == b);
}

inline bool operator>(const QXNumber& a, const QXNumber& b) {
    return b < a;
}

inline bool operator>=(const QXNumber& a, const QXNumber& b) {
    return (a > b) || (a == b);
}

// --- Mathematical Functions ---
inline QXNumber sqrt(const QXNumber& a) { 
    QXNumber result;
    qxsqrt_(a.get_qx_ptr(), result.get_qx_ptr());
    return result; 
}

inline QXNumber abs(const QXNumber& a) {
    return QXNumber(std::abs(a.qx));
}

// --- Type Conversion ---
inline std::string to_string(const QXNumber& a, int digits=33) {
    char s[128] = {0};
    qxtoqd_(&a.qx, &digits, s, sizeof(s));
    return std::string(s, strnlen(s, sizeof(s)));
}

inline double to_double(const QXNumber& a) {
    char s[128] = {0};
    int d = 33;
    qxtoqd_(&a.qx, &d, s, sizeof(s));
    try {
        return std::stod(s);
    } catch (...) {
        return static_cast<double>(a.qx);
    }
}

// --- Stream Output ---
inline std::ostream& operator<<(std::ostream& os, const QXNumber& q) {
    char s[128] = {0};
    int digits = 33;
    qxtoqd_(&q.qx, &digits, s, sizeof(s));
    os << s;
    return os;
}

} // namespace bailey

// --- Eigen Integration ---
namespace Eigen {
    template<> struct NumTraits<bailey::QXNumber> : GenericNumTraits<bailey::QXNumber> {
        typedef bailey::QXNumber Real; 
        typedef bailey::QXNumber NonInteger; 
        typedef bailey::QXNumber Nested;
        enum { 
            IsComplex = 0, 
            IsInteger = 0, 
            IsSigned = 1, 
            RequireInitialization = 1, 
            ReadCost = 1,      // Single scalar value
            AddCost = 8,       // Estimated cost for QX operations
            MulCost = 16       // Estimated cost for QX operations
        };
    };
}

// --- Type Aliases ---
using SpMat_QX = Eigen::SparseMatrix<bailey::QXNumber>;
using Vec_QX = Eigen::Vector<bailey::QXNumber, Eigen::Dynamic>;
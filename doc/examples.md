# Examples and Use Cases

This document provides practical examples for using the Bailey high-precision numerical analysis library.

## Basic Arithmetic Examples

### Simple Calculations

```cpp
#include "high_precision.h"  // Your header file

int main() {
    // Basic arithmetic
    QuadDouble a(1.0);
    QuadDouble b(3.0);
    
    QuadDouble sum = a + b;
    QuadDouble product = a * b;
    QuadDouble quotient = b / a;
    
    std::cout << "1 + 3 = " << sum << std::endl;
    std::cout << "1 * 3 = " << product << std::endl;
    std::cout << "3 / 1 = " << quotient << std::endl;
    
    return 0;
}
```

### High-Precision Constants

```cpp
// Computing π with high precision
QuadDouble pi_approx(3.1415926535897932384626433832795);
QuadDouble radius(1.0);
QuadDouble area = pi_approx * radius * radius;

std::cout << "π ≈ " << pi_approx << std::endl;
std::cout << "Area of unit circle: " << area << std::endl;
```

### Mathematical Functions

```cpp
// Square root calculations
QuadDouble x(2.0);
QuadDouble sqrt_2 = sqrt(x);

std::cout << "√2 = " << sqrt_2 << std::endl;

// Pythagorean theorem with high precision
QuadDouble a(3.0), b(4.0);
QuadDouble c_squared = a*a + b*b;
QuadDouble c = sqrt(c_squared);

std::cout << "Hypotenuse: " << c << std::endl;
```

## Matrix Operations

### Creating Sparse Matrices

```cpp
#include <Eigen/Sparse>
#include <vector>

// Create a tridiagonal matrix
SpMat_QD createTridiagonalMatrix(int n, QuadDouble diag, QuadDouble off_diag) {
    SpMat_QD A(n, n);
    std::vector<Eigen::Triplet<QuadDouble>> triplets;
    
    for (int i = 0; i < n; ++i) {
        // Diagonal elements
        triplets.push_back(Eigen::Triplet<QuadDouble>(i, i, diag));
        
        // Off-diagonal elements
        if (i > 0) {
            triplets.push_back(Eigen::Triplet<QuadDouble>(i, i-1, off_diag));
        }
        if (i < n-1) {
            triplets.push_back(Eigen::Triplet<QuadDouble>(i, i+1, off_diag));
        }
    }
    
    A.setFromTriplets(triplets.begin(), triplets.end());
    return A;
}

// Usage
int main() {
    int n = 100;
    QuadDouble diagonal(4.0);
    QuadDouble off_diagonal(-1.0);
    
    SpMat_QD A = createTridiagonalMatrix(n, diagonal, off_diagonal);
    std::cout << "Created " << n << "×" << n << " tridiagonal matrix" << std::endl;
    std::cout << "Non-zeros: " << A.nonZeros() << std::endl;
    
    return 0;
}
```

### Vector Operations

```cpp
// Create and manipulate high-precision vectors
Vec_QD createTestVector(int n, QuadDouble start_value, QuadDouble increment) {
    Vec_QD v(n);
    
    for (int i = 0; i < n; ++i) {
        v(i) = start_value + QuadDouble(i) * increment;
    }
    
    return v;
}

// Vector operations example
int main() {
    int n = 10;
    
    // Create vectors
    Vec_QD x = createTestVector(n, QuadDouble(1.0), QuadDouble(0.1));
    Vec_QD y = createTestVector(n, QuadDouble(2.0), QuadDouble(0.2));
    
    // Dot product
    QuadDouble dot_product = x.dot(y);
    std::cout << "Dot product: " << dot_product << std::endl;
    
    // Vector norm
    QuadDouble norm_x = sqrt(x.dot(x));
    std::cout << "||x|| = " << norm_x << std::endl;
    
    // Element-wise operations
    Vec_QD sum = x + y;
    Vec_QD scaled = QuadDouble(2.0) * x;
    
    std::cout << "First element of x+y: " << sum(0) << std::endl;
    std::cout << "First element of 2*x: " << scaled(0) << std::endl;
    
    return 0;
}
```

## Linear System Solving

### Simple CG Example

```cpp
#include "conjugate_gradient.h"  // Your CG implementation

int main() {
    // Create a simple 2D Laplacian matrix (5-point stencil)
    int n = 100;  // 10×10 grid
    SpMat_QD A(n, n);
    std::vector<Eigen::Triplet<QuadDouble>> triplets;
    
    for (int i = 0; i < n; ++i) {
        // Diagonal: 4
        triplets.push_back(Eigen::Triplet<QuadDouble>(i, i, QuadDouble(4.0)));
        
        // Off-diagonals: -1 (simplified, assumes grid structure)
        if (i % 10 != 0) {  // Not left boundary
            triplets.push_back(Eigen::Triplet<QuadDouble>(i, i-1, QuadDouble(-1.0)));
        }
        if (i % 10 != 9) {  // Not right boundary  
            triplets.push_back(Eigen::Triplet<QuadDouble>(i, i+1, QuadDouble(-1.0)));
        }
        if (i >= 10) {      // Not top boundary
            triplets.push_back(Eigen::Triplet<QuadDouble>(i, i-10, QuadDouble(-1.0)));
        }
        if (i < 90) {       // Not bottom boundary
            triplets.push_back(Eigen::Triplet<QuadDouble>(i, i+10, QuadDouble(-1.0)));
        }
    }
    
    A.setFromTriplets(triplets.begin(), triplets.end());
    
    // Create right-hand side (all ones)
    Vec_QD b = Vec_QD::Ones(n);
    
    // Initial guess (zero)
    Vec_QD x = Vec_QD::Zero(n);
    
    // Solve
    int max_iterations = 1000;
    double tolerance = 1e-12;
    
    std::cout << "Solving " << n << "×" << n << " system..." << std::endl;
    int iterations = conjugateGradient(A, b, x, max_iterations, tolerance);
    
    if (iterations < max_iterations) {
        std::cout << "Converged in " << iterations << " iterations" << std::endl;
        
        // Check solution quality
        Vec_QD residual = b - A * x;
        QuadDouble residual_norm = sqrt(residual.dot(residual));
        std::cout << "Residual norm: " << residual_norm << std::endl;
    } else {
        std::cout << "Failed to converge" << std::endl;
    }
    
    return 0;
}
```

### Iterative Refinement Example

```cpp
// High-precision iterative refinement
Vec_QD iterativeRefinement(const SpMat_QD& A, const Vec_QD& b, 
                          Vec_QD x, int max_refinements) {
    
    for (int refinement = 0; refinement < max_refinements; ++refinement) {
        // Compute residual: r = b - A*x
        Vec_QD residual = b - A * x;
        
        // Check if refinement needed
        QuadDouble residual_norm = sqrt(residual.dot(residual));
        std::cout << "Refinement " << refinement 
                  << ", residual norm: " << residual_norm << std::endl;
        
        if (to_double(residual_norm) < 1e-15) {
            std::cout << "Refinement converged" << std::endl;
            break;
        }
        
        // Solve correction: A * delta_x = residual
        Vec_QD delta_x = Vec_QD::Zero(A.cols());
        conjugateGradient(A, residual, delta_x, 100, 1e-12);
        
        // Update solution
        x += delta_x;
    }
    
    return x;
}
```

## Numerical Analysis Applications

### Computing π Using Series

```cpp
// Compute π using Machin's formula: π/4 = 4*arctan(1/5) - arctan(1/239)
QuadDouble arctan_series(QuadDouble x, int terms) {
    QuadDouble result(0.0);
    QuadDouble x_power = x;
    QuadDouble x_squared = x * x;
    
    for (int n = 0; n < terms; ++n) {
        QuadDouble term = x_power / QuadDouble(2*n + 1);
        
        if (n % 2 == 0) {
            result += term;
        } else {
            result -= term;
        }
        
        x_power *= x_squared;
    }
    
    return result;
}

QuadDouble compute_pi_machin(int terms) {
    QuadDouble one_fifth(0.2);        // 1/5
    QuadDouble one_239(1.0/239.0);    // 1/239
    
    QuadDouble term1 = QuadDouble(4.0) * arctan_series(one_fifth, terms);
    QuadDouble term2 = arctan_series(one_239, terms);
    
    QuadDouble pi_quarter = term1 - term2;
    return QuadDouble(4.0) * pi_quarter;
}

int main() {
    int terms = 100;
    QuadDouble pi_computed = compute_pi_machin(terms);
    
    std::cout << "π computed with " << terms << " terms:" << std::endl;
    std::cout << pi_computed << std::endl;
    
    return 0;
}
```

### High-Precision Root Finding

```cpp
// Newton's method with high precision
QuadDouble newton_sqrt(QuadDouble x, int max_iterations) {
    QuadDouble guess = x;  // Initial guess
    
    for (int i = 0; i < max_iterations; ++i) {
        QuadDouble f = guess * guess - x;           // f(guess) = guess² - x
        QuadDouble fprime = QuadDouble(2.0) * guess; // f'(guess) = 2*guess
        
        QuadDouble new_guess = guess - f / fprime;
        
        // Check convergence
        QuadDouble diff = new_guess - guess;
        if (to_double(diff * diff) < 1e-30) {
            std::cout << "Newton's method converged in " << i+1 << " iterations" << std::endl;
            return new_guess;
        }
        
        guess = new_guess;
    }
    
    std::cout << "Newton's method did not converge" << std::endl;
    return guess;
}

int main() {
    QuadDouble x(2.0);
    QuadDouble sqrt_2_newton = newton_sqrt(x, 20);
    QuadDouble sqrt_2_builtin = sqrt(x);
    
    std::cout << "√2 (Newton's method): " << sqrt_2_newton << std::endl;
    std::cout << "√2 (built-in):       " << sqrt_2_builtin << std::endl;
    
    // Compare accuracy
    QuadDouble difference = sqrt_2_newton - sqrt_2_builtin;
    std::cout << "Difference: " << difference << std::endl;
    
    return 0;
}
```

## Performance Optimization Examples

### Minimizing Conversions

```cpp
// Inefficient: multiple conversions
double compute_sum_inefficient(const std::vector<double>& data) {
    QuadDouble sum(0.0);
    
    for (double value : data) {
        QuadDouble qd_value(value);     // Conversion for each element
        sum += qd_value;
        double temp = to_double(sum);   // Unnecessary conversion
        // ... some logic with temp
    }
    
    return to_double(sum);
}

// Efficient: minimize conversions
double compute_sum_efficient(const std::vector<double>& data) {
    QuadDouble sum(0.0);
    
    for (double value : data) {
        sum += QuadDouble(value);       // Single conversion per element
    }
    
    return to_double(sum);              // Single final conversion
}
```

### Batch Operations

```cpp
// Process multiple similar problems efficiently
void solve_multiple_systems(const std::vector<SpMat_QD>& matrices,
                           const std::vector<Vec_QD>& rhs_vectors,
                           std::vector<Vec_QD>& solutions) {
    
    int max_iter = 1000;
    double tolerance = 1e-12;
    
    for (size_t i = 0; i < matrices.size(); ++i) {
        std::cout << "Solving system " << i+1 << "/" << matrices.size() << std::endl;
        
        Vec_QD x = Vec_QD::Zero(matrices[i].cols());
        int iterations = conjugateGradient(matrices[i], rhs_vectors[i], x, max_iter, tolerance);
        
        if (iterations < max_iter) {
            solutions[i] = x;
            std::cout << "  Converged in " << iterations << " iterations" << std::endl;
        } else {
            std::cout << "  Failed to converge" << std::endl;
        }
    }
}
```

## Error Analysis Examples

### Condition Number Estimation

```cpp
// Simple condition number estimation using power method
QuadDouble estimate_spectral_radius(const SpMat_QD& A, int iterations = 50) {
    int n = A.rows();
    Vec_QD v = Vec_QD::Random(n);  // Random initial vector
    
    // Normalize
    QuadDouble norm = sqrt(v.dot(v));
    v = v / norm;
    
    QuadDouble eigenvalue(0.0);
    
    for (int i = 0; i < iterations; ++i) {
        Vec_QD Av = A * v;
        eigenvalue = v.dot(Av);
        
        QuadDouble norm_Av = sqrt(Av.dot(Av));
        v = Av / norm_Av;
    }
    
    return eigenvalue;
}
```

### Accuracy Verification

```cpp
// Test accuracy against known solution
void verify_solution_accuracy() {
    // Create system with known solution
    int n = 10;
    SpMat_QD A = SpMat_QD::Identity(n, n);  // Identity matrix
    Vec_QD x_true = Vec_QD::Ones(n);        // True solution: all ones
    Vec_QD b = A * x_true;                  // RHS: should also be all ones
    
    // Solve
    Vec_QD x_computed = Vec_QD::Zero(n);
    int iterations = conjugateGradient(A, b, x_computed, 100, 1e-15);
    
    // Analyze error
    Vec_QD error = x_computed - x_true;
    QuadDouble error_norm = sqrt(error.dot(error));
    QuadDouble relative_error = error_norm / sqrt(x_true.dot(x_true));
    
    std::cout << "Solution accuracy test:" << std::endl;
    std::cout << "Absolute error norm: " << error_norm << std::endl;
    std::cout << "Relative error: " << relative_error << std::endl;
    std::cout << "Iterations: " << iterations << std::endl;
    
    // Component-wise analysis
    for (int i = 0; i < n; ++i) {
        QuadDouble component_error = x_computed(i) - x_true(i);
        std::cout << "Error in component " << i << ": " << component_error << std::endl;
    }
}
```

## Integration with Existing Code

### Wrapper for Double-Precision Code

```cpp
// Wrapper to use high-precision solver with double-precision interface
std::vector<double> solve_with_high_precision(
    const std::vector<std::vector<double>>& matrix_dense,
    const std::vector<double>& rhs) {
    
    int n = matrix_dense.size();
    
    // Convert to high-precision sparse matrix
    SpMat_QD A(n, n);
    std::vector<Eigen::Triplet<QuadDouble>> triplets;
    
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (std::abs(matrix_dense[i][j]) > 1e-15) {
                triplets.push_back(Eigen::Triplet<QuadDouble>(
                    i, j, QuadDouble(matrix_dense[i][j])));
            }
        }
    }
    A.setFromTriplets(triplets.begin(), triplets.end());
    
    // Convert RHS
    Vec_QD b(n);
    for (int i = 0; i < n; ++i) {
        b(i) = QuadDouble(rhs[i]);
    }
    
    // Solve
    Vec_QD x = Vec_QD::Zero(n);
    conjugateGradient(A, b, x, 1000, 1e-15);
    
    // Convert back to double
    std::vector<double> result(n);
    for (int i = 0; i < n; ++i) {
        result[i] = to_double(x(i));
    }
    
    return result;
}
```

These examples demonstrate the practical usage patterns and capabilities of the Bailey high-precision library integration. Each example can be adapted for specific numerical analysis problems requiring extended precision.
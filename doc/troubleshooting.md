# Troubleshooting

This document covers common issues and their solutions when working with the Bailey high-precision libraries.

## Build Issues

### CMake Configuration Errors

#### Error: `Cannot find Bailey library files`

```
CMake Error: Could not find library /opt/qxfun/fortran/libqxfun.a
```

**Cause**: Bailey libraries not built or environment variables not set.

**Solution**:
```bash
# Check environment variables
echo $QXFUN_DIR
echo $DQFUN_DIR  
echo $DDFUN_DIR

# Verify library files exist
ls -la $QXFUN_DIR/lib*.a

# If missing, rebuild Docker image
docker build -t bailey-hp . --no-cache
```

#### Error: `Eigen3 not found`

```
CMake Error at CMakeLists.txt:5 (find_package):
  Could not find a package configuration file provided by "Eigen3"
```

**Solution**:
```bash
# Rocky Linux/RHEL
dnf install eigen3-devel

# Ubuntu/Debian  
apt install libeigen3-dev

# Verify installation
pkg-config --modversion eigen3
```

### Compilation Errors

#### Error: `Fortran compilation failed`

```
Error: Cannot convert REAL(8) to TYPE(dq_real)
```

**Cause**: Type mismatch in Fortran wrapper code.

**Solution**: Check wrapper file for correct type conversions:
```fortran
! Correct: Manual field assignment
da%dqr(1) = real(d, dqknd)
da%dqr(2) = 0.0_dqknd

! Incorrect: Direct assignment
da = d  ! This will fail
```

#### Error: `C++ template instantiation failed`

```
error: no match for 'operator+=' (operand types are 'QuadDouble' and 'QuadDouble')
```

**Cause**: Missing operator overloads in QuadDouble struct.

**Solution**: Ensure all required operators are defined:
```cpp
// Required operators for Eigen3 integration
QuadDouble& operator+=(QuadDouble& a, const QuadDouble& b);
QuadDouble& operator-=(QuadDouble& a, const QuadDouble& b);
QuadDouble& operator*=(QuadDouble& a, const QuadDouble& b);
QuadDouble& operator/=(QuadDouble& a, const QuadDouble& b);
```

### Linking Errors

#### Error: `undefined reference to Bailey functions`

```
undefined reference to `qxadd_'
```

**Cause**: Fortran wrapper not linked or function name mangling issues.

**Solution**:
```bash
# Check that wrapper libraries exist
ls -la $QXFUN_DIR/libqxwrap.a

# Verify function names in library
nm $QXFUN_DIR/libqxwrap.a | grep qxadd

# Check CMake linking order
target_link_libraries(sample_qx PRIVATE qxwrap qxfun gfortran)
```

#### Error: `cannot find -lquadmath`

```
/usr/bin/ld: cannot find -lquadmath
```

**Solution**: Remove quadmath dependency or install development package:
```bash
# Option 1: Remove from CMakeLists.txt (recommended)
target_link_libraries(sample_qx PRIVATE qxwrap qxfun gfortran)

# Option 2: Install libquadmath-devel (if available)
dnf install libquadmath-devel
```

## Runtime Issues

### Execution Errors

#### Error: `Fortran runtime error: Positive exponent width required`

```
Fortran runtime error: Positive exponent width required in format string
(ES0.15E)
```

**Cause**: Invalid format specification in string conversion.

**Solution**: Fix format string in Fortran wrapper:
```fortran
! Incorrect
write(tmp, '(ES0.' // trim(fmt_str) // 'E)') a

! Correct  
write(tmp, '(ES0.' // trim(fmt_str) // ')') a
```

#### Error: `Segmentation fault in matrix operations`

**Cause**: Usually memory access issues or uninitialized QuadDouble values.

**Solution**:
```cpp
// Ensure proper initialization
QuadDouble q;                    // May contain garbage
QuadDouble q_safe(0.0);         // Properly initialized

// Check matrix dimensions
assert(A.rows() == b.size());
assert(A.cols() == x.size());

// Initialize solution vector
Vec_QD x = Vec_QD::Zero(n);     // Not just Vec_QD x(n);
```

### Numerical Issues

#### Issue: `Poor convergence in CG solver`

**Symptoms**: 
- Takes many iterations to converge
- Fails to converge within iteration limit
- Large residual norms

**Diagnosis**:
```cpp
// Check matrix condition
SpMat_QD A = createMatrix();
// Compute condition number (requires additional implementation)

// Check residual manually
Vec_QD residual = b - A * x;
QuadDouble residual_norm = sqrt(residual.dot(residual));
std::cout << "Residual norm: " << residual_norm << std::endl;
```

**Solutions**:
1. **Improve conditioning**: Use better matrix construction
2. **Preconditioning**: Implement preconditioned CG
3. **Alternative algorithms**: Consider other solvers for ill-conditioned systems

#### Issue: `Loss of precision in mixed arithmetic`

**Problem**: Mixing double and QuadDouble operations loses precision.

```cpp
// Problematic: precision loss
QuadDouble a(very_precise_value);
double b = 2.0;
QuadDouble result = a + b;  // b converted to double precision first
```

**Solution**: Use consistent precision:
```cpp
// Better: maintain full precision
QuadDouble a(very_precise_value);
QuadDouble b(2.0);          // Explicit QuadDouble conversion
QuadDouble result = a + b;  // Full precision maintained
```

## Performance Issues

### Slow Execution

#### Issue: `Extremely slow matrix operations`

**Expected**: High-precision operations are 50-200x slower than double precision.

**Optimization strategies**:
```cpp
// 1. Minimize precision conversions
QuadDouble sum(0.0);
for (int i = 0; i < n; ++i) {
    sum += vector[i];           // Good: stay in QuadDouble
}
double final = to_double(sum);  // Convert only once at end

// 2. Use sparse matrix benefits
SpMat_QD A_sparse;              // Only store non-zeros
// vs
Matrix<QuadDouble> A_dense;     // Stores all n² elements

// 3. Algorithm selection
// Sometimes better algorithm > higher precision
```

### Memory Issues

#### Issue: `Excessive memory usage`

**Cause**: QuadDouble uses 32 bytes vs 8 bytes for double.

**Memory usage comparison**:
```
1M × 1M sparse matrix:
- Double: ~80 MB (1% density)
- QuadDouble: ~320 MB (1% density)
```

**Solutions**:
1. **Increase sparsity**: Use sparser matrix representations
2. **Problem decomposition**: Split large problems into smaller ones
3. **Streaming algorithms**: Process data in chunks

## Docker Issues

### Container Problems

#### Issue: `Docker build fails on package installation`

```
Error: Problem with package dependencies
```

**Solution**: Use cache-busting and update package lists:
```bash
# Force rebuild without cache
docker build -t bailey-hp . --no-cache

# Check base image updates
docker pull rockylinux:9
```

#### Issue: `Container runs but executable not found`

```
exec: "./sample_qx": no such file or directory
```

**Solution**: Use correct path and verify build:
```bash
# Check executable location
docker run bailey-hp ls -la ./build/

# Use correct path
docker run bailey-hp ./build/sample_qx

# Or set working directory
docker run -w /work/build bailey-hp ./sample_qx
```

## Development Issues

### IDE and Editor Problems

#### Issue: `C++ intellisense errors with QuadDouble`

**Cause**: IDE doesn't understand custom type definitions.

**Solution**: 
1. Ensure proper include paths for Eigen3
2. Add compile commands database:
```bash
cmake -S . -B build -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
ln -s build/compile_commands.json .
```

#### Issue: `Fortran syntax highlighting missing`

**Solution**: Install Fortran language support in your editor.

### Debugging

#### Issue: `GDB doesn't show QuadDouble values correctly`

**Solution**: Add GDB pretty printer or use manual inspection:
```gdb
# Manual inspection
print q.qd[0]
print q.qd[1] 
print q.qd[2]
print q.qd[3]

# Or convert to double for approximation
print to_double(q)
```

## Environment Issues

### Missing Dependencies

#### Issue: `Required tools not available`

```bash
# Check required tools
which cmake
which ninja  
which gfortran
which g++

# Install missing tools
dnf install cmake ninja-build gcc-gfortran gcc-c++
```

### Library Path Issues

#### Issue: `Library not found at runtime`

```
error while loading shared libraries: cannot open shared object file
```

**Solution**: Set library paths (though we use static linking):
```bash
export LD_LIBRARY_PATH=$QXFUN_DIR:$DQFUN_DIR:$DDFUN_DIR:$LD_LIBRARY_PATH
```

## Getting Help

### Diagnostic Information

When reporting issues, include:

1. **System Information**:
```bash
uname -a
gcc --version
gfortran --version
cmake --version
```

2. **Build Information**:
```bash
echo $QXFUN_DIR
ls -la $QXFUN_DIR/
```

3. **Error Messages**: Full error output, not truncated

4. **Minimal Example**: Smallest code that reproduces the problem

### Common Diagnostic Commands

```bash
# Check Bailey library compilation
nm $QXFUN_DIR/libqxfun.a | grep -i add

# Check wrapper compilation  
nm $QXFUN_DIR/libqxwrap.a | grep qxadd

# Verify CMake configuration
cmake -S . -B debug_build -DCMAKE_VERBOSE_MAKEFILE=ON

# Check Eigen3 installation
find /usr -name "Eigen" -type d 2>/dev/null

# Test basic Fortran compilation
echo 'program test; end program' | gfortran -x f95 -
```

### Resources

- **Bailey Libraries**: https://www.davidhbailey.com/dhbsoftware/
- **Eigen3 Documentation**: https://eigen.tuxfamily.org/
- **CMake Documentation**: https://cmake.org/documentation/
- **Project Repository**: Check CLAUDE.md for specific project guidance
# Dependencies

This document describes all external dependencies and their installation requirements.

## Bailey High-Precision Libraries

### Overview

The project depends on three libraries from David H. Bailey:

| Library | Version | Precision | Decimal Digits | Description |
|---------|---------|-----------|----------------|-------------|
| DDFUN   | v03     | Double-double | ~32 | Uses 2 double-precision numbers |
| DQFUN   | v03     | Quad-double | ~64 | Uses 2 quad-precision numbers |
| QXFUN   | v01     | Extended quad | ~128+ | Single quad-precision with extensions |

### Download Sources

All libraries are automatically downloaded from:
- https://www.davidhbailey.com/dhbsoftware/ddfun-v03.tar.gz
- https://www.davidhbailey.com/dhbsoftware/dqfun-v03.tar.gz
- https://www.davidhbailey.com/dhbsoftware/qxfun-v01.tar.gz

### Environment Variables

The build system requires these environment variables:

```bash
export DDFUN_DIR=/opt/ddfun/fortran
export DQFUN_DIR=/opt/dqfun/fortran  
export QXFUN_DIR=/opt/qxfun/fortran
```

## System Dependencies

### Required Packages (Rocky Linux 9)

```bash
# Core development tools
dnf -y groupinstall "Development Tools"

# Specific packages
dnf -y install gcc-gfortran gcc-c++ cmake ninja-build wget tar git eigen3-devel
```

### Compiler Requirements

- **Fortran**: gfortran with Fortran-2008 support
- **C++**: C++17 compatible compiler (GCC 11.5.0+)
- **Build System**: CMake 3.21+, Ninja

### Mathematical Libraries

- **Eigen3**: Version 3.4+ for sparse matrix operations
- **BLAS/LAPACK**: Provided by system packages (optional)

## Library Structure

### Bailey Library Components

Each Bailey library produces:

```
/opt/{dd,dq,qx}fun/fortran/
├── lib{dd,dq,qx}fun.a    # Main library
├── lib{dd,dq,qx}wrap.a   # C interface wrapper
├── *.o                   # Object files
└── *.mod                 # Fortran modules
```

### Build Process

1. **Download**: Libraries are fetched from Bailey's website
2. **Compile**: Using library-specific build scripts (`gnu-complib-*.scr`)
3. **Wrap**: Custom Fortran-to-C interface layers are compiled
4. **Link**: Static libraries are created and installed

## Version Compatibility

### Tested Combinations

| Component | Version | Status |
|-----------|---------|--------|
| Rocky Linux | 9 | ✅ Tested |
| GCC | 11.5.0 | ✅ Tested |
| gfortran | 11.5.0 | ✅ Tested |
| CMake | 3.20+ | ✅ Tested |
| Eigen3 | 3.4+ | ✅ Tested |

### Known Issues

- **Eigen SparseExtra**: Not available in all Eigen3 distributions
  - Solution: Matrix Market loading is disabled, using direct matrix construction
- **libquadmath**: May require additional setup on some systems
  - Solution: Removed from CMake dependencies (included with gfortran)

## Build Verification

### Check Dependencies

```bash
# Verify compilers
gcc --version
gfortran --version
cmake --version

# Check Eigen3
pkg-config --modversion eigen3

# Verify library installation
ls -la /opt/*/fortran/lib*.a
```

### Environment Setup

```bash
# Set library paths
export LD_LIBRARY_PATH=$DDFUN_DIR:$DQFUN_DIR:$QXFUN_DIR:$LD_LIBRARY_PATH

# Verify environment
echo $QXFUN_DIR
echo $DQFUN_DIR  
echo $DDFUN_DIR
```
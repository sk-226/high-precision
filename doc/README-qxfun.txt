
QXFUN: A quad precision package with special functions

Revision date: 18 Mar 2023

AUTHOR:
  David H. Bailey
  Lawrence Berkeley National Lab (retired) and University of California, Davis
  Email: dhbailey@lbl.gov
 
COPYRIGHT AND DISCLAIMER

All software in this package (c) 2023 David H. Bailey. By downloading or using this software you agree to the copyright, disclaimer and license agreement in the accompanying file DISCLAIMER.txt.

INDEX OF THIS README FILE:

I. PURPOSE OF PACKAGE
II. INSTALLING CODING ENVIRONMENT (FOR MAC OS X SYSTEMS, IF NEEDED)
III. INSTALLING FORTRAN COMPILER (IF NEEDED)
IV. DOWNLOADING AND COMPILING QXFUN 
V. BRIEF SUMMARY OF CODING INSTRUCTIONS AND USAGE
VI. SAMPLE APPLICATION PROGRAMS AND TESTS
VII. RECENT UPDATES
VIII. APPENDIX: LIST OF SPECIAL FUNCTIONS


+++++
I.  PURPOSE OF PACKAGE:
This package enhances an IEEE quad precision floating-point facility (approx. 33 digit accuracy) to include a library of numerous special functions, all by making only very minor changes to existing Fortran programs.  The package should run correctly on any Unix-based system supporting a Fortran-2008 compiler and IEEE 128-bit floating-point arithmetic, in hardware or software (for example, the GNU gfortran compiler and the Intel ifort compiler). Note however that results are NOT guaranteed to the last bit.

In addition to fast execution times, one key feature of this package is a 100% THREAD-SAFE design, which means that user-level applications can be easily converted for parallel execution, say using a threaded parallel environment such as OpenMP.

Three related software packages by the same author are DDFUN (double-double real and complex), DQFUN (double-quad real and complex) and MPFUN2020 (arbitrary precision real and complex). They are available from the same site as this suite:

  http://www.davidhbailey.com/dhbsoftware/


II. INSTALLING CODING ENVIRONMENT (FOR MAC OS X SYSTEMS, IF NEEDED)

For Apple Mac OS X systems (highly recommended for QXFUN), first install the latest supported version of Xcode, which is available for free from the Apple Developer website:

  https://developer.apple.com/

Here click on Account, then enter your Apple ID and password, then go to 

  https://developer.apple.com/download/more/?=xcode

From this list, select Xcode. As of the above date, the latest version of Xcode is 13.3. Download this package on your system, double click to decompress, and place the resulting Xcode app somewhere on your system (typically in the Applications folder). Double-click on the app to run Xcode, allowing it install additional components, then quit Xcode. Open a terminal window, using the Terminal application in the Utilities folder, and type

  xcode-select --install

which installs various command-line tools. The entire process of downloading Xcode and installing command-line tools takes roughly 30 minutes. When this is completed, you should be ready to continue with the installation.


III. INSTALLING FORTRAN COMPILER (IF NEEDED)

Running QXFUN is relatively straightforward, provided that one has a Unix-based system, such as Linux or Apple OS X, and a Fortran-2008 compliant compiler that supports IEEE 128-bit arithmetic. These requirements are met by the GNU gfortran compiler, the Intel ifort compiler, and others.

The gfortran compiler (highly recommended for QXFUN) is available free for a variety of systems at this website:

  https://gcc.gnu.org/wiki/GFortranBinaries

For Apple Mac OS X systems, download the installer file here:

  https://github.com/fxcoudert/gfortran-for-macOS/releases

The gfortran compiler is normally placed in /usr/local/lib and /usr/local/bin. Thus before one uses gfortran, one must insert a line in one's shell initialization file (if the Z shell is used, as on most Apple OS X systems, the shell initialization file is ~/.zshrc). The line to be included is:

  PATH=/usr/local/lib:/usr/local/bin:$PATH

The following line is also recommended for gfortran compiler users:

  GFORTRAN_UNBUFFERED_ALL=yes; export GFORTRAN_UNBUFFERED_ALL

The following line is recommended for inclusion in the shell initialization file, no matter what compiler is used (it prevents stack overflow system errors):

  ulimit -s unlimited

On most Unix systems (including Apple Mac OS X systems), the shell initialization file must be manually executed upon initiating a terminal shell, typically by typing "source .zshrc".


IV. DOWNLOADING AND COMPILING QXFUN

From the website https://www.davidhbailey.com/dhbsoftware, download the file "qxfun-vnn.tar.gz" (replace "vnn" by whatever is the current version on the website, such as "v03"). If the file is not decompressed by your browser, use gunzip at the shell level to do this. Some browsers (such as the Apple Safari browser) do not drop the ".gz" suffix after decompression; if so, remove this suffix manually at the shell level. Then type

  tar xfv qxfun-vnn.tar 

(where again "vnn" is replaced by the downloaded version). This should create the directory and unpack all files. Compiling the library (in directory "fortran") is straightforward -- just compile the three library code files: qxfune.f90, qxmodule.f90 and second.f90. Compile/link scripts are available for the GNU gfortran compiler. To compile the library with this script, type

  ./gnu-complib-qx.scr

Then to compile and link the application program tpslqm1qx.f90, using the GNU gfortran compiler, producing the executable file tpslqm1qx, type

  ./gnu-complink-qx.scr tpslqm1qx

To execute the program, with output to tpslqm1qx.txt, type

  ./tpslqm1qx > tpslqm1qx.txt

These scripts assume that the user program is in the same directory as the library files; this can easily be changed by editing the script files.

Several sample test programs, together with reference output files, are included in fortran directory -- see Section VI below.


V. BRIEF SUMMARY OF CODING INSTRUCTIONS AND USAGE

First place the following line in every subprogram of the user's application code that contains a quad variable or array, at the beginning of the declaration section, before any implicit or type statements:

  use qxmodule

To designate a variable or array as quad real (QXR) in your application code, use the Fortran-90 real type statement with the kind "qxknd", as in this example:

  real (qxknd) a, b(m), c(m,n)

Here "qxknd" is the kind parameter for IEEE quad precision, which is set automatically in qxfune.f90. Similarly, to designate a variable or array as quad complex (QXC), use a complex type statement with kind "qxknd", as in

  complex (qxknd) x, y(m), z(m,n)

Most common mixed-mode combinations (arithmetic operations, comparisons and assignments) involving QXR, QXC, double precision real (DP) and double complex (DC) variables should be supported by the compiler. Explicit conversion can be done using the real and cmplx functions with the kind parameter qxknd, as in

  a = real (1.d0, qxknd)
  x = cmplx (1.d0, 2.d0, qxknd)

However, there are some hazards when combining say double precision and quad precision entities. For example, the code r1 = 3.14159d0, where r1 is QXR, does NOT produce the true quad equivalent of 3.14159, unless the numerical value is a modest-sized whole number or exact binary fraction. To obtain the full quad converted value, write this as r1 = 3.14159q0 instead. Similarly, the code r2 = r1 + 3.d0 * sqrt (2.d0), where r1 and r2 are QXR, does NOT produce the true quad value, since the expression 3.d0 * sqrt (2.d0) will be performed in double precision (according to standard Fortran-90 precedence rules). To obtain the fully accurate QXR result, write this as r2 = r1 + 3.q0 * sqrt (2.q0).

As noted above, quad precision constants should be written with a q exponent, such as 3.25q0. Technically speaking, the q exponent, while almost universally supported by compilers that support quad precision, is nonstandard; some fully standard-conforming options are 3.25_qxknd or 3.25e0_qxknd or real(3.25d0,qxknd). Here "qxknd" is the kind parameter for IEEE quad precision, which is set automatically in qxfune.f90.

Most Fortran-2008 intrinsic functions should be supported with QXR and QXC arguments by the compiler, as appropriate. In addition, as mentioned above, the QXFUN package supports numerous special functions. A complete list of special functions and subroutines supported by the package is summarized in the Appendix below (section VIII).


VI. SAMPLE APPLICATION PROGRAMS AND TESTS

Four application programs are included, together with corresponding reference output files (e.g., tpphixqx.ref.txt) for comparison with user results:

testqxfun.f90  Tests most special functions.

tpslqm1qx.f90  Performs the 1-level multipair PSLQ integer relation algorithm.

tpphixqx.f90  Performs a Poisson polynomial application, using 1-level multipair PSLQ.

tquadqx.f90   Evaluates definite integrals, using tanh-sinh, exp-sinh and sinh-sinh algorithms.

In addition, the package includes a test script that compiles the library and run each of the above sample programs above. To run this script, type:

  ./gnu-qxfun-tests.scr

For each test program, the script outputs either TEST PASSED or TEST FAILED. If all tests pass, then one can be fairly confident that the QXFUN software is working properly. 


VII. RECENT UPDATES

16 Mar 2023: Initial release.


VIII. APPENDIX: SPECIAL FUNCTIONS

In these tables, "F" denotes function, "S" denotes subroutine, "QXR" denotes quad precision real. The variable names r1,r2 are QXR; n,nb,np,nq are integers; rb is an QXR vector of length nb; ss is an QXR vector of length n; aa is an QXR vector of length np; bb is an QXR vector of length nq.

Special functions:

Type   Name                 Description 
hline
F:QXR  agm(r1,r2)           Arithmetic-geometric mean of r1 and r2
S      qxberne(nb,rb)       Initialize array rb of length nb with first nb even 
                              Bernoulli numbers [3] 
F:QXR  bessel_i(r1,r2)      BesselI function, order r1, of r2  
F:QXR  bessel_in(n,r1)      BesselI function, order n, of r1 
F:QXR  bessel_j(r1,r2)      BesselJ function, order r1, of r2  
F:QXR  bessel_jn(n,r1)      BesselJ function, order n, of r1 
F:QXR  bessel_i(r1,r2)      BesselK function, order r1, of r2  
F:QXR  bessel_kn(n,r1)      BesselK function, order n, of r1 
F:QXR  bessel_y(r1,r2)      BesselY function, order r1, of r2  
F:QXR  bessel_yn(n,r1)      BesselY function, order n, of r1 
F:QXR  digamma_be(nb,rb,r1) Digamma function of r1, using nb even Bernoulli 
                              numbers in rb [3] 
F:QXR  erf(r1)              Error function 
F:QXR  erfc(r1)             Complementary error function 
F:QXR  expint(r1)           Exponential integral function 
F:QXR  gamma(r1)            Gamma function 
F:QXR  hurwitz_zetan(k,r1)  Hurwitz zeta function, order n >= 2, of r1 [4] 
F:QXR  hurwitz_zetan_be (nb,rb,n,r1)  Hurwitz zeta function, order n >= 2, of r1, 
                              using nb even Bernoulli numbers in rb [3] 
F:QXR  hypergeom_pfq(np,nq,aa,bb,r1)  Hypergeometric pFq function of aa, bb and r1; 
                              dimensions are aa(np) and bb(nq) [5] 
F:QXR  incgamma(r1,r2)      Incomplete gamma function [6] 
F:QXR  polygamma(k,r1)      Polygamma function, order k >= 1, of r1 [4] 
F:QXR  polygamma(nb,rb,k,r1)  Polygamma function, order k >= 1, of r1, using 
                              nb even Bernoulli numbers in rb [3] 
S      polylog_ini(n,ss)    Initialize array ss, of size |n|, for computing
                              polylogarithms of order n when n < 0 [6]
F:QXR  polylog_neg(n,ss,r1) Polylogarithm function of r1, for n < 0, using 
                              precomputed data in ss [7]
F:QXR  polylog_pos(n,r1)    Polylogarithm function, order n >= 0, of r1 [8]
F:QXR  struve_hn(n,r1)      StruveH function, order n >= 0, of r1 [9]
F:QXR  zeta(r1)             Zeta function of r1
F:QXR  zeta_be(nb,rb,r1)    Zeta function of r1, using nb even Bernoulli
                              numbers in rb [3]
F:QXR  zeta_int(n)          Zeta function of integer argument n

Notes:
[3]: For most applications, set nb >= 60; see qxberne above.
[4]: For hurwitz_zetan and polygamma, the argument r1 is limited to the range (0, 1).
[5]: For hypergeom_pfq, the integers np and nq must not exceed 10.
[6]: For incgamma, r1 must not be zero, and must not be negative unless r2 = 0.
[7]: For polylog_ini and polylog_neg, the integer n is limited to the range [-1000, -1]. 
[8]: For polylog_pos, the argument r1 is limited to the range (-1,1). 
[9]: For struve_hn, the argument r1 is limited to the range [-1000, 1000].


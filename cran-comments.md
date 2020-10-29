## Test environments
* ubuntu 14.04 (on travis-ci), R 4.0.2, 3.6.1
* ubuntu 14.04 (on circle-ci), R 4.0.2, 3.6.1
* win-builder (devel and release) Windows Server 2008 R2 SP1, R-devel, 32/64 bit
* rhub (release) Ubuntu Linux 16.04 LTS, R-release, GCC
* Fedora (devel) Linux, R-devel, clang, gfortran
* MacOS 10.14.6 R 3.6.1 
* MacOS 10.15.7 R 4.0.2

## R CMD check results

0 errors ✔ | 0 warnings ✔ | 0 notes ✔

## Updates

Minor update to implement functions in python version.

Improvements to install of python libraries upon loading if python
or conda environments are available (checks suppressed on systems
where they are not).

## Python integration

Vignettes are disabled when python is not available. This works on Linux test environments (with python available) and windows test environments (without python or pip packages).

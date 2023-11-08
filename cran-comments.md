## Test environments
* ubuntu 14.04 (on travis-ci), R 4.1.0, 3.6.1
* ubuntu 14.04 (on circle-ci), R 4.1.0, 3.6.1
* win-builder (devel and release) Windows Server 2008 R2 SP1, R-devel, 32/64 bit
* rhub (release) Ubuntu Linux 16.04 LTS, R-release, GCC
* Fedora (devel) Linux, R-devel, clang, gfortran
* MacOS 10.14.6 R 3.6.1 
* Red Hat Enterprise Linux 8.5 R 4.1.2
* MacOS 10.15.7 R 4.2.0
* CentoOS 7 R 4.2.1
* RockyOS 8.6 R 4.2.2

## R CMD check results

0 errors ✔ | 0 warnings ✔ | 0 notes ✔

## Updates

Avoids archiving as a reverse dependency of knitr on upcoming release. Updates documentation to ensure matching chunk delimiters.

## Python integration

Python is a soft dependency which is still required for some functions but is not essential for core functionality any longer. It is retained for backwards compatibility.

Vignettes are disabled when python is not available. This works on Linux test environments (with python available) and windows test environments (without python or pip packages).

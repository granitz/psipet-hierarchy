# External MATLAB dependencies

Two functions used by the fitting code are not redistributed here, because they
were written by others and are not covered by this repository's MIT license.
Place them in this directory (or anywhere on the MATLAB path) before running
`run_hopf_all_psipet.sh`.

## `HopfModel_LNA.m`

Computes the stationary functional connectivity, lagged covariance, power
spectral density and cross-spectrum of a network of Stuart-Landau (Hopf)
oscillators using the linear noise approximation, by solving the continuous-time
Lyapunov equation. Header attribution: **Adrián Ponce-Alvarez, 06-07-2022**.

Obtain it from the author, or from the code release accompanying the paper that
introduces the method. It carries no license in the copy we used, so we cannot
redistribute it.

## `gaussfilt.m`

Gaussian smoothing of an unevenly sampled time series. Header attribution:
**James Conder, 22 August 2013**. Available from the MATLAB Central File
Exchange. Check the license on the File Exchange page; if it is BSD, the file
may be vendored here together with its license text.

## Checking

The driver scripts assume both are on the path:

```matlab
assert(exist('HopfModel_LNA','file')==2, 'HopfModel_LNA.m not found -- see matlab/external/README.md')
assert(exist('gaussfilt','file')==2,     'gaussfilt.m not found -- see matlab/external/README.md')
```

# PHOENIXgrid.jl

![Lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)<!--
![Lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-stable-green.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-retired-orange.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-archived-red.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-dormant-blue.svg) -->
[![Build Status](https://travis-ci.com/jbrinchmann/PHOENIXgrid.jl.svg?branch=master)](https://travis-ci.com/jbrinchmann/PHOENIXgrid.jl)
[![codecov.io](http://codecov.io/github/jbrinchmann/PHOENIXgrid.jl/coverage.svg?branch=master)](http://codecov.io/github/jbrinchmann/PHOENIXgrid.jl?branch=master)
<!--
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://jbrinchmann.github.io/PHOENIXgrid.jl/stable)
[![Documentation](https://img.shields.io/badge/docs-master-blue.svg)](https://jbrinchmann.github.io/PHOENIXgrid.jl/dev)
-->


A very basic routine to fit PHOENIX grid spectra to MUSE spectra of stars.
The aim is to have this as a complement to the `spexxy` code for santiy
checks but the code is Bayesian and uses a coarse grid which might be
suboptimal in most cases.

# A Scalable Approach for Analysing Multi-Agent Systems with Heterogeneous Stochastic Packet Loss

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7787750.svg)](https://doi.org/10.5281/zenodo.7787750)

## General

This repository contains an implementation of the algorithms and simulations described in the paper

> C. Hespe and H. Werner, "A Scalable Approach for Analysing Multi-Agent Systems with Heterogeneous Stochastic Packet Loss," submitted to *62nd Conference on Decision and Control*. IEEE, 2023.

It may be used to recreate and validate the figures from the paper.
To do so, run either of the two main entry points in the repository, the scripts `scalability.m`, and `uncertainty_sweep.m`.
Be advised that the scripts have a runtime of a few hours.
The raw data used in the figures is available in the subdirectory `figures`.

## Prerequisites

To run the scripts in this repository, you will need a working copy of [*Yalmip*](https://yalmip.github.io/) together with a suitable SDP solver in your *Matlab* path.

The code in this repository was tested in the following environment:

* *Windows 10* Version 21H2
* *Matlab* 2021a
* *Yalmip* 31-March-2021

The *Matlab* [`parfor`](https://de.mathworks.com/help/parallel-computing/parfor.html) feature from the *Parallel Computing Toolbox* is used to speed up the calculations.
*Matlab* should automatically detect if that toolbox is not available and run the iterations sequentially in that case.
However, this will drastically prolong the runtime of the scripts.
You may want to reduce the number of sampling points for the figures or run the calculations for smaller networks.

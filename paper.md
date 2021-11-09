---
title: 'WAVI.jl: An ice sheet model, written in Julia'
tags:
  - Julia
  - glaciology
  - ice sheet modelling
  - ice shelves
  - Antarctica
  - Greenland
  - climate
authors:
  - name: Alexander T. Bradley, Robert Arthern
    orcid: 0000-0001-8381-5317
    affiliation: 1
affiliations:
 - name: British Antarctic Survey, Cambridge, UK
   index: 1
date: 8 October 2021
bibliography: paper.bib
---

# Summary

Why do we need ice sheet models?
The Antarctic and Greenland ice sheets represent two of the largest contributors to future sea level rise [ref], and most uncertainty in future projections arises from these sources. Predictions of future sea level rise come primarily from physically-based ice sheet models. (More of why we need ice sheet models?)

What is WAVI and what does it do? Part one: solution of governing equations
On relevant length scales, ice shelves behave approximately as a viscous fluid with a shear-thinning rheology. WAVI.jl (Wavelet-based Adaptive-grid Vertically-integrated Ice-sheet-model) is a Julia package for the numerical solution of a first-order approximation to equations describing conservation of mass and momentum in such a fluid (the Stokes equations). This approximation, which is appropriate for flows with a high aspect ratio (as is the case for Antarctic and Greenland and all ice streams on long length scales), treats longitudinal and lateral stresses as depth independent, but accounts for vertical velocity gradients in the nonlinear viscosity and in the treatment of basal stress. The equation solved in ``WAVI.jl`` are described in detail in Arthern 2015 and Goldberg 2011.

An accurate representation of floating ice, and the boundary between grounded and floating ice, is particularly important in ice sheet [why?]. WAVI.jl includes a subgrid parametrization of grounding line motion. [Need to explain why melt rate physics is so important to get right!] Owing to the computational expense of coupled ice-ocean models, ice sheet models typically use a parametrization (an analytic approximation to) the melt rate applied to the ice shelf by the adjacent ocean; WAVI.jl includes the most commonly used and physically advanced melt rate parametrizations [cite: Favier, Asay-Davis] including plume model emulator [Lazeroms], PICO, PICOp, and quadratic, linear, mismip etc. WAVI.jl leverages Julia's multiple dispatch paradigm to create a simple, user-friendly interface for embedding further melt rate models into WAVI.jl, as well as models/parametrizations of other important processes such a damage, calving, and solid earth effects. Julia's high level programming interface also facilitates this process.

``WAVI.jl`` was designed to be used by both researchers (including those without prior ice sheet modelling experience) and students. ``WAVI.jl`` employs an adaptive-grid preconditioner (explain what adaptive means here?) that uses a wavelet-based method that originated in solutions of similar elliptic systems of equations(?) which arise in studies of electrical impedence (ref). (Points for friendly:)[!!] The repository includes a number of well-documented examples which can be played with, and a user friendly API facilitates use.

MATLAB version WAVI has been used in research and participated in MISMIP (etc) experiments [figure], version written in Julia is open source, more friendly, similar in speed (with plans to develop parallelisation)

Other ice similar ice sheet models are...all ice sheet models make approximations and include different physics, so comparison between them is where there value arises.



# Acknowledgements

We acknowledge contributions from Rosie Williams, David Bett, Xy Wang, and Daniel Goldberg.

# References
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
Ice sheet models allow us to simulate the behaviour of ice sheets. Such models can be used in a prognostic manner, to make predictions about the future, or in a diagnostic manner, to investigate the behaviour of an ice sheet. The most pertinent prognostic use of ice sheet models relates to contributions to sea level rise: the world's two largest ice sheets -- Antarctica and Greenland -- hold enough ice to raise sea levels by 61 and 8 meters, respectively; prognostic modelling of these ice sheets enables us to make predictions on (for example) how much of this ice will be lost [@Edwards2021Nature], to investiagte the possibility of runaway ice loss  [@DeConto2021], and whether such instabilities have already been initiated [@Favier2014]. Diagnostic modelling allows us to investigate processes controlling the behaviour of an ice sheet, such as how loss of floating ice influences the rest of an ice sheet [@Joughin2021], how different bed conditions affect ice sliding [@DeRydt2021], and probing the conditions under which the aforementioned instabilities might be triggered [@Schoof2007].

On the extemely long lengthscales (thousands of kilometers), flowing ice behaves approximately as a highly viscous fluid with a shear-thinning rheology, meaning that the viscosity reduces with deformation. WAVI.jl (Wavelet-based Adaptive-grid Vertically-integrated Ice-sheet-model) is a Julia package for the numerical solution of a first-order approximation to the Stokes equatiosn, which describe conservation of mass and momentum in such a fluid. This approximation, which is appropriate for fluid flows with a high aspect ratio (as is the case for the vast majority of ice sheets), treats longitudinal and lateral stresses as depth independent, but accounts for vertical velocity gradients in the nonlinear viscosity and in the treatment of basal stress[@Goldberg2011] [@Arthern2015]. 

Physically, ice sheets to do not stand alone, but are forced by their interactions with both other climate systems elements. For the West Antarctic Ice Sheet in particular, the most important of these is with the ocean, which interacts with ice sheets via basal melting of ice shelves. `WAVI.jl` includes a number of the most prominent parametrizations (analytic approximations to) of basal melting [@AsayDavis2017], as well as a developmental coupling with the ocean general circulation model MITgcm [@Marshall1997]. In addition, `WAVI.jl` leverages Julia's multiple dispatch paradigm to create a simple, user-friendly interface for embedding models of other physical processes, such as ice damage, ice shelf calving, and solid earth effects, into `WAVI.jl`.

``WAVI.jl`` was designed to be used by both researchers (including those without prior ice sheet modelling experience) and students. ``WAVI.jl`` employs a number of tools to improve computational speed, making it appropriate for high-resolution simulations of large ice sheets. (Wavelet-based Adaptive-grid refers to one such tools: an adaptive-grid preconditioner that uses a wavelet-based method that originated in solutions of similar elliptic systems of equations which arise in studies of electrical impedence [Vasilyev2005@].)  `WAVI.jl` also includes a simple, user friendly API, which is facilitated by Julia's convenient syntax, and the repository includes a number of well-documented examples, which demonstrate the software's capabilities in a wide variety of situations.

`WAVI.jl` is successor to a similar code that was written in MATLAB, but was retained "in-house" at the British Antarctic Survey. This predecessor code has been used extensively as research software (e.g. [@Arthern2017] [@Arthern2015]), as well as having participated in the most recent ice sheet model intercomparison exercise [@Cornford2020], which acts as a benchmark for ice-flow models. The successor code `WAVI.jl` has been verified independently to agree well with its predecessor; similar speed. `WAVI.jl` also benefits from being more transparent than its predecessor, having been written in an open source language and hosted in a public repository. 

Being in an open repository, and written in an open source language, the successor `WAVI.jl` is more transparent

The new code is written in an open source language and also 

MATLAB version WAVI has been used in research and participated in MISMIP (etc) experiments [figure], version written in Julia is open source, more friendly, similar in speed (with plans to develop parallelisation)

Other ice similar ice sheet models are...all ice sheet models make approximations and include different physics, so comparison between them is where there value arises.



# Acknowledgements

We acknowledge contributions from Rosie Williams, David Bett, Xy Wang, and Daniel Goldberg.

# References
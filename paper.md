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
  - name: Alexander T. Bradley, Robert Arthern, David Bett, C. Rosie Williams
    orcid: 0000-0001-8381-5317
    affiliation: 1
affiliations:
 - name: British Antarctic Survey, Cambridge, UK
   index: 1
date: 8 October 2021
bibliography: paper.bib
---

# Summary
Ice sheet models allow us to simulate the behaviour and evolution of ice sheets, which are large (>50000 km<sub>2</sub>) masses of glacial land ice. There are two main flavours of ice sheet model usage: prognostic, which involves using ice sheet models to make predictions about the future, and diagnostic, to investigate the behaviour of an ice sheet. Today, the most pertinent prognostic use of ice sheet models relates to contributions to sea level rise: the world's two largest ice sheets, located in Antarctica and Greenland, hold enough ice to raise sea levels by approximately 58 and 7 meters, respectively [@Bamber2018]; prognostic modelling of these ice sheets enables us to make predictions on (for example) how much of this ice will be lost [@Edwards2021], to investigate the possibility of runaway ice loss  [@DeConto2021], and analyze whether or not such instabilities have already been initiated [@Favier2014]. Diagnostic modelling allows us to investigate processes controlling the behaviour of an ice sheet, such as how loss of ice shelves -- the floating extensions of ice sheets -- ice influences ice flow speed [@Joughin2021], how different bed conditions affect ice sliding [@DeRydt2021], and probing the conditions under which the aforementioned instabilities might be triggered [@Schoof2007].

On the extemely long lengthscales that are relevant to ice sheets, flowing ice behaves approximately as a highly viscous fluid with a shear-thinning rheology (as the ice deforms, it becomes thinner). WAVI.jl (Wavelet-based Adaptive-grid Vertically-integrated Ice-sheet-model) is a Julia package for the numerical solution of a leading-order approximation to the Stokes equations, which describe conservation of mass and momentum in such a fluid. This approximation, which is appropriate for fluid flows with a high aspect ratio (as is the case for the vast majority of ice sheets), treats longitudinal and lateral stresses as depth independent, but accounts for vertical velocity gradients in the nonlinear viscosity and in the treatment of basal stress[@Goldberg2011] [@Arthern2015]. 

Physically, ice sheets do not stand alone, but are forced by their interactions with both other parts of the climate system. For example, the rapid changes that have occured in the West Antarctic Ice Sheet in the previous decades are thought to be driven by an increase in oceanic heat content reaching [@Pritchard2012], and thus basal melting of, the ice shelves in this region. For basal melting in particular, `WAVI.jl` includes a number of the most prominent parametrizations (analytic approximations to) of basal melting [@AsayDavis2017], as well as a developmental coupling with the ocean general circulation model MITgcm [@Marshall1997]. More generally, `WAVI.jl` leverages Julia's multiple dispatch paradigm to create a simple, user-friendly interface for embedding models of other physical processes, such as ice damage, ice shelf calving, and solid earth effects, into `WAVI.jl`.

WAVI.jl was designed to be used by both researchers (including those without prior ice sheet modelling experience) and students. To that end, `WAVI.jl` employs a number of tools to improve computational speed, making it appropriate for high-resolution simulations of large ice sheets. These include multithreading capabilities, and a wavelet-based adaptive-grid preconditioner (to which the name WAVI.jl refers), which originated in solutions of similar systems of equations that arise in electrical impedence [@Vasilyev2005].)  `WAVI.jl` also includes a simple, user friendly API, which is facilitated by Julia's convenient syntax. In addition, the repository includes a number of well-documented examples, which demonstrate the software's capabilities in a wide variety of situations.

WAVI.jl is the successor to a similar code that was written in the proprietary programming language MATLAB, and was not publicly released. This predecessor code has been used extensively as research software [e.g. [@Arthern2017]], as well as having participated in the most recent ice sheet model intercomparison exercise [@Cornford2020], which acts as a benchmark for ice-flow models. WAVI.jl, has been verified independently against these benchmark experiments.

There currently exists a wide variety of ice sheet models of varying complexity. Examples include (but are certainly not limited to) the Ice Sheet System Model [@Larour2012], the Parallel Ice Sheet Model [@Bueler2009], BISICLES [@Cornford2013],  Elmer/Ice [@Gagliardini2013], and Úa [@Gudmundsson2019]. Since every ice sheet model makes approximations in order to facilitate the numerical solution of the appropriate governing equations, and the fact that these equations lack analytic solutions, means that the _intercomparison_ between ice sheet models is of paramount importance when assessing the trustworthiness of models; `WAVI.jl` contributes to this community of ice sheet models.

# Acknowledgements

We acknowledge contributions from Xy Wang, and Daniel Goldberg.

# References
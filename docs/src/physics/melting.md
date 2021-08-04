# Melt Rates
WAVI.jl includes a number of community melt rate parametrizations of melt rate:
- ['Analytic' parametrizations](#Analytic-Melt-Rate-Parametrizations) 
- [Input file melt rates](#Input-File-Melt-Rates) 
- [Quadratic temperature parametrization](#Quadratic-Temperature-Melt-Rate-Parametrization)
- [Plume emulator parametrization](#Plume-Emulator-Melt-Rate-Parametrization)
- [PICO parametrization](#PICO-Melt-Rate-Parametrization)
- [PICOp parametrization](#PICOp-Melt-Rate-Parametrization)

Further details of these parametrizations can be found on this page.  For details of their use in WAVI.jl, see the [melt rate models](../data_structure/melt_rate_models.md) section. We strongly encourage those who have developed melt rate parametrizations to consider implementing them in WAVI.jl; if you are interested in doing so, see the [melt rate models](../data_structure/melt_rate_models.md) and contributors guide(../contributing.md), or [get in touch](mailto:aleey@bas.ac.uk).

WAVI.jl also supports coupling to the ocean model MITgcm. Please [get in touch](mailto:aleey@bas.ac.uk), or see the [MITgcm coupling](../mitgcm_coupling.md) tab if you are interested in running coupled WAVI.jl-MITgcm simulations.

## Analytic Melt Rate Parametrizations

## Input File Melt Rates

## Quadratic Temperature Melt Rate Parametrization
The quadratic temperature melt rate parametrization [[Holland2008](@cite)] parametrizes the melt rate as a quadratic function of the thermal driving:
```math
\begin{equation}\label{E:quadratic_parametrization}
M = \gamma_T \left( \frac{\rho_w c_p}{\rho_i L} \right) \left(T_0 - T_f\right)^2.
\end{equation}
```
Here $\gamma_T$ is a heat exchange velocity, $\rho_w$ is the density of water, $\rho_i$ is the density of ice, $c_p$ is the heat capacity of the ocean, $L$ is the latent heat of fusion of ice. In addition, $T_f = \lambda_1 S_0 + \lambda_2 + \lambda_3 z_b$ is the local freezing point, with $\lambda_1$, $\lambda_2$, and $\lambda_3$ the liquidus slope, intercept, and pressure coefficient, respectively, and $z_b$ the height of the ice shelf draft above sea level (i.e. $z_b$ is negative).  $T_0$ and $S_0$ are the depth-dependent potential temperatre and practical salinity taken from the far field.

The quadratic formulation \eqref{E:quadratic_parametrization} attempts to account for heat providing both more heat for melting and the feedback between sub-shelf melting and circulation in the cavity (higher temperatures result in a more vigorous circulation, promoting enhanced melt rates).

## Plume Emulator Melt Rate Parametrization

## PICO Melt Rate Parametrization

## PICOp Melt Rate Parametrization
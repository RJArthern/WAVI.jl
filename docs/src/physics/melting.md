# Melt Rates
WAVI.jl includes a number of community melt rate parametrizations of melt rate:
- [Input file melt rates](#Input-File-Melt-Rates) 
- [Quadratic temperature parametrization](#Quadratic-Temperature-Melt-Rate-Parametrization)
- [Plume emulator parametrization](#Plume-Emulator-Melt-Rate-Parametrization)
- [PICO parametrization](#PICO-Melt-Rate-Parametrization)
- [PICOP parametrization](#PICOp-Melt-Rate-Parametrization)

Further details of these parametrizations can be found on this page. Much of these descriptions is taken from [Favier2019(@cite), which describes a thorough assessment of different melt rate parametrizations. For details of the use of melt rate parametrizations in WAVI.jl, see the [melt rate models](../data_structure/melt_rate_models.md) section. We strongly encourage those who have developed melt rate parametrizations to consider implementing them in WAVI.jl; if you are interested in doing so, see the [melt rate models](../data_structure/melt_rate_models.md) and contributors guide(../contributing.md), or [get in touch](mailto:aleey@bas.ac.uk).

WAVI.jl also supports coupling to the ocean model MITgcm. Please [get in touch](mailto:aleey@bas.ac.uk), or see the [MITgcm coupling](../mitgcm_coupling.md) tab if you are interested in running coupled WAVI.jl-MITgcm simulations.

## Input File Melt Rate Parametrizations 

## Quadratic Temperature Melt Rate Parametrization
The quadratic temperature melt rate parametrization [[Holland2008](@cite)] parametrizes the melt rate as a quadratic function of the thermal driving:
```math
\begin{equation}\label{E:quadratic_parametrization}
M = \gamma_T \left( \frac{\rho_w c_p}{\rho_i L} \right) \left(T_0 - T_F \right)^2.
\end{equation}
```
Here $\gamma_T$ is a heat exchange velocity, $\rho_w$ is the density of water, $\rho_i$ is the density of ice, $c_p$ is the heat capacity of the ocean, $L$ is the latent heat of fusion of ice. In addition, $T_F = \lambda_1 S_0 + \lambda_2 + \lambda_3 z_b$ is the local freezing point, with $\lambda_1$, $\lambda_2$, and $\lambda_3$ the liquidus slope, intercept, and pressure coefficient, respectively, and $z_b$ the height of the ice shelf draft above sea level (i.e. $z_b$ is negative).  $T_0$ and $S_0$ are the depth-dependent potential temperatre and practical salinity taken from the far field.

The quadratic formulation \eqref{E:quadratic_parametrization} attempts to account for heat providing both more heat for melting and the feedback between sub-shelf melting and circulation in the cavity (higher temperatures result in a more vigorous circulation, promoting enhanced melt rates).

## Plume Emulator Melt Rate Parametrization
The plume emulator melt rate parametrization from [Lazeroms2018](@cite) emulates the 2-D behaviour of the 1-D plume model of [Jenkins1991](@cite). This model describes the evolution of a buoyant plume originating from the grounding line
with zero thickness and velocity, and temperature and salinity taken from the ambient ocean. Away from the grounding line, the thickness, velocity, temperature, and salinity of the plume evolve through advection, turbulent exchange across
the ocean boundary layer underneath the ice shelf, and entrainment of deep water. 

The melt rate in the plume model emulator can be expressed as 
```math
\begin{equation}\label{E:plume_parametrization}
M = \alpha M_0 g(\theta)(T_0 - T_{F,gl})^2 \hat{M}(\hat{X})
\end{equation}
```
where $M_0$ is a the melt rate prefactor, $\alpha$ is a calibration coefficient, $g(\theta)$ is an expression of the mean basal slope at shelf points and $\hat{M}$ is a universal, dimensionless function of $\hat{X}$, an expression for the dimensionless distance from the grounding line. The mean basal slope $\theta$ and dimensionless distance $\hat{X}$ are determined using a path-finding algorithm that is described in detail in [Lazeroms2018](@cite).

## PICO Melt Rate Parametrization
The Potsdam Ice-shelf Cavity mOdel (PICO) melt rate from [Reese2018](@cite) is based on a one-dimensional ocean box-model which coarsely resolves ice shelf cavities. The box model represents the buoyancy-driven advection of ambient ocean water into the ice-shelf cavity at depth up to the grounding line, then upward along the ice
draft in consecutive boxes. The melt rates in the box model are given by 
```math
\begin{equation}\label{E:pico_parametrization}
M = \gamma_T \left( \frac{\rho_w c_p}{\rho_i L} \right) \left(T_0 - T_{F,k} \right)
\end{equation}
```
where the subscript $k$ indicates properties evaluated in box $k$.  Those properties account for the transformation of ocean temperature and salinity in consecutive boxes through heat and salt turbulent exchange across the ocean boundary layer underneath ice shelves.

## PICOP Melt Rate Parametrization
The PICOp melt rate parametrization from [Pelle2019] is a coupling between the PICO and Plume Emulator melt rate parametrizations. This parametrization uses the box model formulation of PICO, but the melt rate within each box is determined using the plume parametrization melt rate \eqref{E:plume_parametrization}, rather than \eqref{E:pico_parametrization}.
# Governing Equations

## Preliminaries
WAVI.jl uses a Cartesian co-ordinate system $\mathbf{x} = (x,y,z)$, with $z$ positive upwards; the corresponding velocity components are $\mathbf{u} = (u,v,w)$. We use bar notation to denote depth averaged quantities, for example:
```math
    \begin{equation}
        \bar{f} =\frac{1}{h} \int_{z = b(x,y)}^{z = s(x,y,t)}f~\mathrm{d}z
    \end{equation}
```
is the depth average of the quantity $f$. Here, $t$ is the denotes time, $b(x,y)$ is the (known) bed elevation (measured positive upwards), and $s(x,y,t)$ is the surface elevation. 

 We assume that the ice is in hydrostratic equilibrium, so that regions are with $h < -(\rho_i/\rho_w) b$ are floating, and regions with $h \geq -(\rho_i/\rho_w) b$ are grounded, where $\rho_i$ and $\rho_w$ the ice and ocean density, respectively.  Where the ice is grounded, we have $s = h + b$, while where the ice is floating, the hydrostratic assumption enforces $s = (1 - \rho_i/\rho_w)h$.

WAVI.jl solves equations describing conservation of momentum and conservation of mass for $\mathbf{\bar{u}}(x,y,t) = (\bar{u}(x,y,t), \bar{v}(x,y,t))$, the depth averaged velocity components in the $(x,y)$ directions, respectively, and the ice thickness $h(x,y,t)$

## Conservation of Momentum
Conservation of momentum requires that the $\bar{u}$ and $\bar{v}$ satisfy ([Goldberg2011](@cite)):
```math
\begin{align}
    \frac{\partial}{\partial x}\left(4 \bar{\eta}h \frac{\partial \bar{u}}{\partial x} + 2 \bar{\eta}h \frac{\partial \bar{v}}{\partial y})\right) +    \frac{\partial}{\partial y}\left(\bar{\eta}h \frac{\partial \bar{v}}{\partial x} +  \bar{\eta}h \frac{\partial \bar{u}}{\partial y}\right) - \tau_{b,x} &= \rho_i g h \frac{\partial s}{\partial x}, \label{E:x-momentum}\\
    \frac{\partial}{\partial y}\left(4 \bar{\eta}h \frac{\partial \bar{v}}{\partial y} + 2 \bar{\eta}h \frac{\partial \bar{u}}{\partial x})\right) +    \frac{\partial}{\partial x}\left(\bar{\eta}h \frac{\partial \bar{u}}{\partial y} +  \bar{\eta}h \frac{\partial \bar{v}}{\partial x}\right) - \tau_{b,y} &= \rho_i g h \frac{\partial s}{\partial y},\label{E:y-momentum}
\end{align}
```
where $\rho_i$ is the ice density, $g$ is the gravitational acceleration, $\mathbf{\tau}_b = (\tau_{b,x}, \tau_{b,y})$ is the basal drag in the $(x,y)$ directions, and $\eta$ is the ice viscosity, defined implicity in terms of the velocity components (the strain components are themselves functions of $\eta$, see below):
```math
\begin{equation}\label{E:viscosity}
    \eta = \frac{B}{2} \left[\left(\frac{\partial \bar{u}}{\partial x}\right)^2  + \left(\frac{\partial \bar{v}}{\partial y}\right)^2 + \frac{\partial \bar{u}}{\partial x}\frac{\partial \bar{v}}{\partial y} + \frac{1}{4}\left( \frac{\partial \bar{u}}{\partial y} + \frac{\partial \bar{u}}{\partial x}\right)^2 + \frac{1}{4}\left(\frac{\partial \bar{u}}{\partial z}\right)^2 + \left(\frac{\partial \bar{v}}{\partial z}\right)^2 + \epsilon^2\right]^\frac{1-n}{2n}.
\end{equation}
```
Here $n$ is the exponent in a nonlinear Glen flow law, $\epsilon$ is a regularization parameter that prevents the viscosity becoming unbounded at small strain rates (for small strain rates, $\eta$ is constant, corresponding to a linear rheology), and $B(x,y,z)$ is a temperature-dependent coefficient that determines the stiffness of the ice. 

## Boundary Conditions
The momentum equations are solved alongside boundary conditions at the lateral boundary of the ice sheet,
```math
    \begin{align}
        -\frac{1}{2}\rho_w  h_w^2 \hat{n}_x = 2\bar{\eta}h\left(2 \frac{\partial \bar{u}}{\partial x} + \frac{\partial \bar{v}}{\partial y}\right)\hat{n}_y - \frac{1}{2}\rho_i g h^2 \hat{n}_x + \bar{\eta}h \left(\frac{\partial \bar{u}}{\partial y} + \frac{\partial \bar{v}}{\partial x}\right)\hat{n}_y,\label{E:bc1}\\
        -\frac{1}{2}\rho_w  h_w^2 \hat{n}_y = 2\bar{\eta}h\left(2 \frac{\partial \bar{v}}{\partial y} + \frac{\partial \bar{u}}{\partial x}\right)\hat{n}_y - \frac{1}{2}\rho_i g h^2 \hat{n}_y + \bar{\eta}h \left(\frac{\partial \bar{u}}{\partial y} + \frac{\partial \bar{v}}{\partial x}\right)\hat{n}_x, \label{E:bc2}
     \end{align}
```
which impose continuity of depth-integrated momentum there. In \eqref{E:bc1}--\eqref{E:bc2}, $h_w = \max(h - s + \zeta, 0)$ is the thickness of ice below the water level, where $\zeta$ is the sea level with respect to $z = 0$, and $\hat{\mathbf{n}} = (\hat{n}_x, \hat{n}_y)$ is the normal to the lateral boundary. 

In addition, a Robin boundary condition at the bed linearly relates the basal stress $\tau_b = (\tau_{b,x},~\tau_{b,y})$ to the basal velocity $\mathbf{u}_b$ via a multiplicative drag coefficient $\beta$:
```math
    \begin{align}
        \tau_b = \beta  \mathbf{u}_b.
    \end{align}
```

(A no-stress condition at the surface is also implicit in the derivation of \eqref{E:x-momentum}--\eqref{E:y-momentum}.)

## Conservation of Mass
For a given depth-averaged velocity $\mathbf{\bar{u}}$, accumulation rate $a(x,y,t)$ (positive for ice gain), and basal melt rate $m(x,y,t)$ (positive for ice loss), conservation of ice mass requires that the ice thickness $h$ satisfies
```math
    \begin{equation}\label{E:mass_cons}
        \frac{\partial h}{\partial t} = a - m - \nabla. \left(h \mathbf{u}\right).
    \end{equation}
```




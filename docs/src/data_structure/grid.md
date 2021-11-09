# Grid
A `Grid` is a WAVI.jl object that stores information relating to the [numerical grid](../numerical_procedure/numerical_procedure.md#Numerical-Grid). The key parameters, which are passed via keyword arguments in the `OutputParams` constructor, are as follows:
- `nx`: number of x grid points
- `ny`: number of y grid points
- `nσ`: number of levels in the vertical
- `dx`: grid spacing in x 
- `dy`: grid spacing in y 
- `x0`: grid origin x co-ordinate 
- `y0`: grid origin y co-ordinate
- `h_mask`: Array of size `nx` x `ny` that defines domain points: ones in the `h_mask` indicate loations within ice domain, while zeros in the `h_mask` indicate locations outside of the ice domain. Defaults to `trues(nx,ny)`, corresponding to every grid point being in the ice domain.
- `u_iszero`: Array of size `nx` x `ny` that defines locations of zero velocity in the x-direction: ones in the `u_iszero` correspond to points where the velocity in the x-direction is forced to be zero, while zeros in the `u_iszero` array corresponds to points where the velocity is free. Defaults to `falses(nx,ny)`, corresponding to no restrictions on u velocity anywhere. `u_iszero` can also be specified by an array with entries from "North", "South", "East", and "West", which set the corresponding edges of the domain to have zero velocity boundary conditions there.
- `v_iszero`: As in `u_iszero` for the velocity in the y-direction.
- `σ`: Dimensionless locations of sigma (vertical) levels. Values must be increasing between 0 and 1.
- `quadrature_weights`: weights associated with sigma (vertical) levels used in quadrature scheme to determine depth averaged quantities (e.g. viscosity).

WAVI.jl has several grids, on which different quantities are defined (see [numerical grid](../numerical_procedure/numerical_procedure.md#Numerical-Grid)). A `Grid` object contains explicit definitions of these grids; for example, the x-coordinates of the "h grid" are stored in `grid.xxh` (where `grid` is the name of the `Grid` instance).
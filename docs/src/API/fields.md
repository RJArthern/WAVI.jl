# Fields
WAVI.jl stores information relating to the solutions on fields, which are stored in the `model`. Fields are organised based according to the different grids (HGrid, CGrid etc -- see [Grids](grid.md)) on which quantities are stored. Note that all quantities here are set internally (i.e. they cannot be modified by the user). This page contains a full directory of the quantities defined on each of the Grids: [HGrid](#HGrid), [UGrid](#UGrid), [VGrid](#VGrid), [CGrid](#CGrid), [SigmaGrid](#SigmaGrid).

## HGrid
The following quantities are stored in the `HGrid` structure within `Fields` (accessible via `model.fields.gh.<field_name>`):
- `h`: Ice thickness
- `b`: Bed elevation
- `s`: Ice surface elevation ($s = b + h$)
- `dhdt`: Time rate of change of surface elevation
- `accumulation`: Ice accumulation rate (postive for mass added to the surface, negative for mass removed)
- `basal_melt`: Melt rate applied to the base of the ice sheet (positive for mass removed from the base, negative for mass added).
- `grounded_fraction`: Grid cell grounded fraction (zero in the shelf, one in the shelf and interpolated across the grounding line, see [Seroussi et al. (2014)](https://tc.copernicus.org/articles/8/1699/2014/))
- `u`: Depth average ice velocity in the $x$ direction.
- `v`: Depth average ice velocity in the $y$ direction.
- `av_speed`: Depth averaged ice speed 
- `us`: Surface ice velocity in the $x$ direction.
- `vs`: Surface ice velocity in the $y$ direction.
- `ub`: Base ice velocity in the $x$ direction.
- `vb`: Base ice velocity in the $y$ direction.
- `bed_speed`: Ice speed at the base
- `weertman_c`: Weertman C drag cofficient
- `haf`: The height above floatation
- `dsdh`: Rate of change of surface elevation with respect to thickness change
- `shelf_strain_rate`: Strain rate computed using only longitudinal strain components
- `β`: Raw multiplicative drag coefficient
- `βeff`: Effective multiplicative drag coefficient
- `τbed`: Ice stress at the bed.
- `ηav`: Depth averaged viscosity
- `quad_f1`: $F_1$ computed from numerical quadrature
- `quad_f2`: $F_2$ computed from numerical quadrature
- `mask`: boolean matrix that defines the solution space

For convenience, various utility matrices and rheological operators are also stored on the `HGrid`:
- `crop`: diagonal matrix with mask entries on the diagonal 
- `samp`: boolean matrix that takes full domain to the model domain
- `spread`: sparse form of the sampling matrix
- `dneghηav`: $-h \times \bar{\eta}$
- `dimplicit`: $-\rho_i \times g \times \mathrm{d}t \times \mathrm{d}h/\mathrm{d}s$

## UGrid
The following quantities are stored on the `UGrid` structure within `Fields` (accessible via `model.fields.gu.<field_name>`):
- `nxu`: Number of grid cells in x-direction in the `UGrid` (equals `nx + 1`, where `nx` is the number of grid cells in the x direction in the `HGrid`).
- `nyu`: Number of grid cells in y-direction in the `UGrid` (equals `ny`, where `ny` is the number of grid cells in the x direction in the `HGrid`).
- `mask`: Mask specifying model domain with respect to the `UGrid`.
- `n`: Total number of cells in `UGrid`.
- `levels`: Number of levels in the preconditioner
-  `dwt`: Wavelet matrix product at `UGrid` points 
- `s`: Ice surface elevation at `UGrid` points
- `h`: Ice thickness at `UGrid` points
- `grounded_fraction`: Grounded fraction at `UGrid` points
- `βeff`: Effective sliding coefficient at `UGrid` points.
- `u`: Ice velocities in the x-direction

For convenience, various utility matrices and rheological operators are also stored on the `UGrid`:
- `crop`: diagonal matrix with mask entries on the diagonal 
- `samp`: boolean matrix that takes full domain to the model domain
- `spread`: sparse form of the sampling matrix
- `cent`: maps quantities from the `UGrid` to the `HGrid`
- `∂x`: matrix representation of differentiation with respect to x 
- `∂y`: matrix representation of differentiation with respect to y

## VGrid
The following quantities are stored on the `VGrid` structure within `Fields` (accessible via `model.fields.gv.<field_name>`):
- `nxv`: Number of grid cells in x-direction in the `VGrid` (equals `nx`, where `nx` is the number of grid cells in the x direction in the `HGrid`).
- `nyu`: Number of grid cells in y-direction in the `VGrid` (equals `ny+1`, where `ny` is the number of grid cells in the x direction in the `HGrid`).
- `mask`: Mask specifying model domain with respect to the `VGrid`.
- `n`: Total number of cells in `VGrid`.
- `levels`: Number of levels in the preconditioner
-  `dwt`: Wavelet matrix product at `VGrid` points 
- `s`: Ice surface elevation at `VGrid` points
- `h`: Ice thickness at `VGrid` points
- `grounded_fraction`: Grounded fraction at `VGrid` points
- `βeff`: Effective sliding coefficient at `VGrid` points.
- `u`: Ice velocities in the x-direction

For convenience, various utility matrices and rheological operators are also stored on the `VGrid`:
- `crop`: diagonal matrix with mask entries on the diagonal 
- `samp`: boolean matrix that takes full domain to the model domain
- `spread`: sparse form of the sampling matrix
- `cent`: maps quantities from the `VGrid` to the `HGrid`
- `∂x`: matrix representation of differentiation with respect to x 
- `∂y`: matrix representation of differentiation with respect to y

## CGrid
The following quantities are stored in the `CGrid` structure within `Fields` (accessible via `model.fields.gc.<field_name>`):
- `nxc`: Number of grid cells in x-direction in the `CGrid` (equals `nx-1`, where `nx` is the number of grid cells in the x direction in the `HGrid`).
- `nyc`: Number of grid cells in y-direction in the `CGrid` (equals `ny-1`, where `ny` is the number of grid cells in the x direction in the `HGrid`).
- `mask`: Mask specifying model domain with respect to the `CGrid`.
- `n`: Total number of cells in `CGrid`.

The following utility matrices are also stored on the `CGrid`
- `crop`: diagonal matrix with mask entries on the diagonal 
- `samp`: boolean matrix that takes full domain to the model domain
- `spread`: sparse form of the sampling matrix
- `cent`: maps quantities from the `CGrid` to the `HGrid`

## SigmaGrid
The following quantities are stored in the `SigmaGrid` structure within `Fields` (accessible via `model.fields.g3d.<field_name>`):
- `nxs`: Number of grid cells in the x-direction in the SigmaGrid (equals `nx`, where `nx` is the number of grid cells in the x direction in the `HGrid`)
- `nvs`: Number of grid cells in the y-direction in the SigmaGrid (equals `ny`, where `ny` is the number of grid cells in the x direction in the `HGrid`)
- `nxs`: Number of levels in the vertical in the SigmaGrid 
- `σ`: Dimensionless levels in the vertical
- `ζ`: Reverse dimensionless vertical levels
- `quadrature_weights`: weights associated with sigma levels, used in computation of integrals over thickness
- `η`: three dimensional viscosity field
- `θ`: three dimensional temperature field
- `Φ`:  three dimensional damage field
- `glen_b`: three dimensional field of glen_b values in viscosity calcluations


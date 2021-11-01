# Fields
WAVI.jl stores information relating to the solutions on fields, which are stored in the `model`. Fields are organised based according to the different grids (HGrid, CGrid etc -- see [Grids](grid.md)) that quantities are stored. Note that all quantities here are set internally (i.e. they cannot be modified by the user). This page contains a full directory of the quantities defined on each of the Grids: [HGrid](#HGrid), [UGrid](#UGrid), [VGrid](#VGrid), [CGrid](#CGrid), [SigmaGrid](#SigmaGrid).

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

## VGrid

## CGrid

## SigmaGrid
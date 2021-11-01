# Solver Parameters
A `SolverParameters` object is a WAVI.jl object that stores parameters related to the [numerical solution](../numerical_procedure/numerical_procedure.md) of the [governing equations](../physics/governing_equations.md). The following table summarises these parameters


| Keyword Argument   | Description                   | Default Value         |
| ------------------ | ----------------------------- | ------------------- |
| `n_iter_viscosity` | Maximum number of iterations used in determining implicit viscosity        | 2       |
| `maxiter_picard`               | Number of Picard iterations in velocity solve            | 30                  |
| `tol_picard`               | Tolerance at which Picard iteration considered to have converged               | 1e-5                 |
| `tol_coarse`               | Relative tolerance at which the adaptive mesh is coarsened | ∘C/m                |
| `maxiter_coarse`               | Maximum number of coarsening iterations       | 1000       |
| `levels`               | Number of wavelet levels            | 3       |
| `wavelet_threshold`             |    | 10.0       |
| `nsmooth`                       |  | 5          |
| `smoother_omega`               | | 1.0             |
| `stencil_margin`               |  | 3           |



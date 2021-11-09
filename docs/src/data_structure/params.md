# Parameters
A `Params` object is a WAVI.jl object that stores parameters related to physical parameters involved in the problem. The following table summarises these parameters:


| Keyword Argument   | Description                   | Default Value         |
| ------------------ | ----------------------------- | ------------------- |
| `g`                | Gravitional acceleration      | 9.81 m/s^2          |
| `density_ice`      | Ice density                   | 918.0 kg/m^3        |
| `density_ocean`    | Ocean density                 | 1028.0 kg/m^3       |
| `gas_const`        | Gas constant used in Arrhenius relation in determining Glen flow parameter                | 8.314      |
| `sec_per_year`     | Seconds per year                 | 3.16e7 s      |
| `default_thickness`| Default ice thickness choice, value set everywhere if ice thickness is not set           | 100 m     |
| `default_viscosity`| Default ice viscosity choice, value set everywhere if ice thickness is not set           | 1e7 Pa a^{1/3}     |
| `default_temperature`| Default ice temperature choice, value set everywhere if ice thickness is not set           | 265.7 K     |
| `default_damage`| Default dimensionless ice damage choice, value set everywhere if ice thickness is not set           | 0     |
| `accumulation_rate`| Accumulation rate. Can be a scalar (same value everywhere), or an array of the same size as the numerical grid.                 | 0.3 m/a      |
| `glen_a_activation_energy`     | Activation energy used in Glen flow law    |    5.8631e4 J   |
| `glen_a_ref`     |  Reference value used in Glen flow law              | 1.5453e-17      |
| `glen_temperature_ref`     |   Reference temperature used in Glen flow law     | 263.15 K      |
| `glen_n`     | Exponent in Glen flow law                | 3     |
| `glen_reg_strain_rate`     | Strain rate regularization value, sets a lower bound on the strain rate               | 1e-5    |
| `weertman_c`     |  Sliding coefficient in Weertman sliding law. Can be a scalar (same value everywhere), or an array of the same size as the numerical grid.              | 1e4  Pa a^{1/3}/m^{1/3}     |
| `weertman_m`     | Exponent in Weertman sliding law               | 3      |
| `weertman_reg_speed`     | Minimum basal sliding speed                 | 1e-5 m/a     |
| `sea_level_wrt_geoid`     | Reference sea level        | 0 m     |
| `minimum_thickness`     | Minimum ice thickness (on grid points in ice domain)             | 50 m     |




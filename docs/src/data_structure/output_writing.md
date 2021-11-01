# Output Writing
Parameters relating to writing solution data are stored in an `OutputParams` structure. Data is output as individual files at a user-specified interval, and may be zipped to a single file once the simulation has completed.

The key parameters, which are passed via keyword arguments in the `OutputParams` constructor, are as follows:
- `outputs`: a named tuple that points to data fields to be output. Variable names used in the output are those specified in the named tuple. For example, to set the outputs to be the ice sheet thickness, named `h`, and the ice sheet surface elevation, named `s`, we would first instantiate a `Model` (here named `model`) and set `outputs = (h = model.fields.gh.h, s = model.fields.gh.s)`.
!!! note
    If only outputting a single field, you must include a comma so that Julia interprets it as a named tuple, e.g. `outputs = (h = model.fields.gh.h,)`.

!!! note
    Variables named `x`, `y`, and `TIME` are reserved for the spatial and temporal dimensions and are therefore not permitted as variable names.

- `output_freq`: the model time frequency with which the simulation is output. Note that WAVI works with the number of timesteps per output, calculated as `floor(output_freq/dt)` where `dt` is the simulation timestep. If `output_freq` is not an integer multiple of `dt`, this will introduce differences in the actual output frequency, and that specified by the user.

- `output_format` specifies the file format that the solution are to be saved in. Currently only `.jld2` and `.mat` file formats are supported. 

- `prefix` specifies the prefix of solution files. Output files are saved as `<prefix>_xxxxxxxxxx.jld` where `xxxxxxxxxx` is the timestep number at which this solution was outputted, prepended with zeros to be ten digits in length. 

- `output_path` specifies the path (relative or absolute) into which the solution is outputted. Note that if the specified path does not exist, outputs will be written to the local directory. 

- `zip_format` specifies the format that solution files are zipped into. Currently WAVI.jl only supports zipping to an `nc` format. 

- `dump_vel` is a shortcut to dump the velocity at the end of a simulation. This is an experimental feature, necessary for coupling WAVI to the MITgcm (see the [MITgcm coupling](../mitgcm_coupling.md) tab).

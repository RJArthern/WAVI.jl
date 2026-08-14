using CSV
using NCDatasets
using Plots
using Printf

"""Return a NetCDF variable, matching `name` exactly or case-insensitively."""
function nc_var(ds, name::AbstractString)
    haskey(ds, name) && return ds[name]
    lname = lowercase(name)
    for k in keys(ds)
        lowercase(k) == lname && return ds[k]
    end
    throw(KeyError(name))
end

function create_heatmap_animation(netcdf_file::String, variable_name::String, output_dir::String)
    @info "Creating heatmap animation for variable: $variable_name"
    
    # Create output directory if it doesn't exist
    mkpath(output_dir)
    
    # Open NetCDF file
    nc = NCDataset(netcdf_file)
    
    try
        # Read coordinates (WAVI output may use mixed-case names, e.g. TIME)
        x = nc_var(nc, "x")
        y = nc_var(nc, "y")
        time = nc_var(nc, "time")
        data = nc_var(nc, variable_name)

        @info "Data dimensions: $(size(data))"
        @info "Time steps: $(length(time))"
        
        # Create heatmaps for each time step
        for (i, t) in enumerate(time)
            # Extract 2D slice for this time step
            slice_2d = data[:, :, i]

            # Handle case where min == max (prevents GKS invalid range errors)
            d_min, d_max = minimum(data), maximum(data)
            if d_min ≈ d_max
                d_min -= 0.1
                d_max += 0.1
            end

            # Create heatmap
            p = heatmap(
                x, y, slice_2d',
                title="$variable_name at time = $(@sprintf("%.2f", t))",
                xlabel="X",
                ylabel="Y",
                color=:viridis,
                clim=(d_min, d_max),
                aspect_ratio=:equal,
                size=(600, 500),
            )
            
            # Save frame
            frame_filename = joinpath(output_dir, "$(variable_name)_frame_$(lpad(i, 4, '0')).png")
            savefig(p, frame_filename)
            
            if i % 10 == 0 || i == length(time)
                @info "Processed frame $i/$(length(time))"
            end
        end
        
        @info "Animation frames saved to: $output_dir"
        
    finally
        close(nc)
    end
end

"""
    plot_multiple_variables(netcdf_file::String, variables::Vector{String}, output_dir::String)

Create heatmap animations for multiple variables from NetCDF file.
"""
function plot_multiple_variables(netcdf_file::String, 
                                 variables::Vector{String}, 
                                 output_dir::String)
    for var in variables
        var_output_dir = joinpath(output_dir, var)
        try
            create_heatmap_animation(netcdf_file, var, var_output_dir)
        catch e
            @info "Error processing variable $var: $e"
        end
    end
end

"""
    plot_resource_timeseries(csv_paths; labels=String[], reference_cores=nothing, output=...)

Overlay RSS and CPU time series from one or more `resource_timeseries.csv` files.

CPU is plotted as a fraction: uses a `cpu_fraction` column if present, otherwise
`cpu_cores_used / reference_cores`. When `reference_cores` is omitted, each CSV
picks it up from the run’s `benchmark_results.json` (`metadata.reference_cores`),
falling back to `1`. Saves a two-subplot PNG.
"""
function plot_resource_timeseries(
    csv_paths::Vector{<:AbstractString};
    labels::Vector{<:AbstractString} = String[],
    reference_cores::Union{Nothing, Real} = nothing,
    output::AbstractString = joinpath(BENCHMARK_OUTPUT_DIR, "resource_comparison.png"),
)
    isempty(csv_paths) && error("Pass at least one resource_timeseries.csv path.")
    reference_cores === nothing || reference_cores > 0 || error("reference_cores must be positive.")

    labs = if isempty(labels)
        [basename(dirname(abspath(p))) for p in csv_paths]
    else
        collect(String, labels)
    end
    length(labs) == length(csv_paths) ||
        error("Number of labels ($(length(labs))) must match number of CSV paths ($(length(csv_paths))).")

    colours = Plots.palette(:tab10, 10)
    colour_i = 0
    next_colour() = (colour_i += 1; colours[mod1(colour_i, length(colours))])

    plt_rss = Plots.plot(;
        title = "Resident set size",
        ylabel = "RSS (MB)",
        legend = false,
        grid = true,
        gridalpha = 0.25,
        gridcolor = :grey50,
        framestyle = :box,
        bottom_margin = 0Plots.mm,
    )
    plt_cpu = Plots.plot(;
        title = "CPU utilisation",
        xlabel = "Elapsed time (s)",
        ylabel = "CPU fraction",
        legend = false,
        grid = true,
        gridalpha = 0.25,
        gridcolor = :grey50,
        framestyle = :box,
        top_margin = 2Plots.mm,
        bottom_margin = 2Plots.mm,
    )
    # Legend band inset with side margins so labels stay within the figure width,
    # otherwise, it was going outside the figure.
    plt_leg = Plots.plot(;
        legend = :bottom,
        legendcolumns = 3,
        framestyle = :none,
        axis = false,
        grid = false,
        ticks = false,
        lims = (0, 1),
        background_color_inside = :white,
        foreground_color_legend = :grey40,
        background_color_legend = :white,
        top_margin = 0Plots.mm,
        bottom_margin = 2Plots.mm,
        left_margin = 12Plots.mm,
        right_margin = 12Plots.mm,
    )

    t_max = 0.0
    rss_max = 0.0
    cpu_max = 1.0
    legend_entries = Tuple{String, Any}[]

    for (path, lab) in zip(csv_paths, labs)
        isfile(path) || error("CSV not found: $path")
        rows = CSV.File(path)
        cols = propertynames(rows)
        t = collect(rows.elapsed_s)
        t_max = max(t_max, maximum(t; init = 0.0))

        run_colour = next_colour()
        for (rss_bytes, rss_lab, role) in _rss_series_for_plot(rows, cols, lab)
            colour = role === :primary ? run_colour : next_colour()
            y = rss_bytes ./ 1024^2
            rss_max = max(rss_max, maximum(y; init = 0.0))
            Plots.plot!(plt_rss, t, y; label = "", linewidth = 1.5, color = colour)
            push!(legend_entries, (rss_lab, colour))
        end

        frac = if :cpu_fraction in cols
            collect(rows.cpu_fraction)
        else
            ref = something(reference_cores, _reference_cores_for_csv(path))
            collect(rows.cpu_cores_used) ./ Float64(ref)
        end
        cpu_max = max(cpu_max, maximum(frac; init = 0.0))
        # Same colour as the primary run legend entry (no duplicate CPU legend item).
        Plots.plot!(plt_cpu, t, frac; label = "", linewidth = 1.5, color = run_colour)
    end

    x_end = t_max > 0 ? t_max * 1.02 : 1.0
    Plots.xlims!(plt_rss, 0, x_end)
    Plots.xlims!(plt_cpu, 0, x_end)
    Plots.ylims!(plt_rss, 0, max(rss_max * 1.05, eps()))
    Plots.ylims!(plt_cpu, 0, max(1.05, cpu_max * 1.05))

    for (lab, colour) in legend_entries
        Plots.plot!(plt_leg, [NaN], [NaN]; label = lab, linewidth = 1.5, color = colour)
    end
    n_leg = max(length(legend_entries), 1)
    # Prefer wrapping over spilling past the figure sides.
    ncols = n_leg <= 3 ? n_leg : (n_leg <= 6 ? 3 : 4)
    Plots.plot!(plt_leg; legendcolumns = ncols)

    # Compact vertical stack: reduce gap above the legend band.
    lay = Plots.@layout [a{0.46h}; b{0.46h}; c{0.08h}]
    fig = Plots.plot(
        plt_rss,
        plt_cpu,
        plt_leg;
        layout = lay,
        size = (1400, 860),
        plot_title = "Resource usage over time",
        left_margin = 12Plots.mm,
        right_margin = 12Plots.mm,
        top_margin = 4Plots.mm,
        bottom_margin = 4Plots.mm,
        foreground_color_legend = :grey40,
        background_color_legend = Plots.RGBA(1, 1, 1, 0.97),
        legendfontsize = 8,
        titlefontsize = 13,
        guidefontsize = 10,
        tickfontsize = 9,
    )

    out_dir = dirname(output)
    !isempty(out_dir) && mkpath(out_dir)
    Plots.savefig(fig, output)
    @info "Saved resource comparison to: $output"
    return output
end

"""Return sorted `rss_rank*` column symbols present in `cols`."""
function _rss_rank_columns(cols)
    ranks = Tuple{Int, Symbol}[]
    for c in cols
        m = match(r"^rss_rank(\d+)$", String(c))
        m === nothing && continue
        push!(ranks, (parse(Int, m.captures[1]), c))
    end
    sort!(ranks; by = first)
    return [sym for (_, sym) in ranks]
end

"""Choose which memory series to draw for one CSV, with legend labels."""
function _rss_series_for_plot(rows, cols, lab::AbstractString)
    rank_cols = _rss_rank_columns(cols)
    if length(rank_cols) >= 2
        r0 = collect(getproperty(rows, rank_cols[1]))
        n_other = length(rank_cols) - 1
        others = zeros(Float64, length(r0))
        for c in rank_cols[2:end]
            others .+= collect(getproperty(rows, c))
        end
        others ./= n_other
        return (
            (r0, "$lab (rank 0)", :rank0),
            (others, String(lab), :primary),
        )
    elseif length(rank_cols) == 1
        return ((collect(getproperty(rows, rank_cols[1])), String(lab), :primary),)
    else
        return ((_rss_bytes_column(rows, cols), String(lab), :primary),)
    end
end

"""Pick a memory column from the CSV when per-rank columns are not present."""
function _rss_bytes_column(rows, cols)
    if :rss_max_bytes in cols
        return collect(rows.rss_max_bytes)
    elseif :rss_sum_bytes in cols
        return collect(rows.rss_sum_bytes)
    elseif :rss_bytes in cols
        return collect(rows.rss_bytes)
    end
    error("CSV has no RSS column (expected rss_rank*, rss_max_bytes, rss_sum_bytes, or rss_bytes).")
end

"""How many cores this run was meant to use, read from the run’s JSON (or `1`)."""
function _reference_cores_for_csv(csv_path::AbstractString)
    json_path = joinpath(dirname(abspath(csv_path)), "benchmark_results.json")
    isfile(json_path) || return 1.0
    try
        data = JSON3.read(read(json_path, String))
        meta = get(data, :metadata, nothing)
        if meta !== nothing && haskey(meta, :reference_cores)
            return Float64(meta.reference_cores)
        end
        haskey(data, :reference_cores) && return Float64(data.reference_cores)
    catch e
        @debug "Could not read reference_cores from $json_path" exception = e
    end
    return 1.0
end

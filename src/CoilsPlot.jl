module CoilsPlot

using Statistics
using LightGraphs
using ProgressMeter
using PyPlot
using PyCall

using Coils

export plot_vertex, plot_vertices, plot_edges, plot_system,
    plot_cell, plot_edge_current, plot_loop, plot_simpleloop,
    save_plots_simpleloops, save_report, plot_edge_currents,
    plot_cells, plot_loops_field, plot_deviation_histogram


# Would like to do:
# @pyimport matplotlib.patheffects as patheffects
# @pyimport arrow3d
# but annot use precompilation with @pyimport
# ref. https://github.com/JuliaPy/PyCall.jl#using-pycall-from-julia-modules
function __init__()
    global patheffects, arrow3d
    patheffects = pywrap(pyimport("matplotlib.patheffects"))
    pushfirst!(PyVector(pyimport("sys")["path"]), joinpath(dirname(@__FILE__)))
    arrow3d = pywrap(pyimport("arrow3d"))

    using3D()
end

# To be able to write linspace.(starts, stops, lengths)
linspace(start, stop, length) = range(start, stop = stop, length = length)
signif(x, n) = round(x, sigdigits = n)


function plot_poi(poi; standalone = true)
    standalone && subplot(111, projection = "3d")

    # plot needs array of xs, then ys and then zs
    gca().plot(getindex.(poi, 1), getindex.(poi, 2), getindex.(poi, 3), ".")

    if standalone
        xlabel("x")
        ylabel("y")
        # tight_layout()
        equalize_3d_aspect()
    end
end


function plot_vertex(vertex_position; standalone = true, label = "")
    standalone && subplot(111, projection = "3d")

    gca().plot([vertex_position[1]], [vertex_position[2]], [vertex_position[3]],
        "o", color = "black", markersize = 1)
    if label != ""
        txt = gca().text(vertex_position...,  label,
            horizontalalignment = "center", verticalalignment = "center")
        txt.set_path_effects([patheffects.Stroke(linewidth = 3, foreground = "white"),
                           patheffects.Normal()])
    end

    if standalone
        xlabel("x")
        ylabel("y")
        # tight_layout()
        equalize_3d_aspect()
    end
end


function plot_vertices(vertex_positions; standalone = true, labels = true)
    standalone && subplot(111, projection = "3d")

    for (i, vertex_position) in enumerate(vertex_positions)
        plot_vertex(vertex_position, standalone = false, label = labels ? "$i" : "")
    end

    if standalone
        xlabel("x")
        ylabel("y")
        # tight_layout()
        equalize_3d_aspect()
    end
end


function plot_edges(g, vertex_positions; standalone = true, alpha = 1)
    standalone && subplot(111, projection = "3d")

    for edge in edges(g)
        pos1 = vertex_positions[src(edge)]
        pos2 = vertex_positions[dst(edge)]
        gca().plot([pos1[1], pos2[1]], [pos1[2], pos2[2]], [pos1[3], pos2[3]],
            color = "black", alpha = alpha)
    end

    if standalone
        xlabel("x")
        ylabel("y")
        # tight_layout()
        equalize_3d_aspect()
    end
end


function plot_system(g, vertex_positions, poi; standalone = true)
    standalone && subplot(111, projection = "3d")

    plot_poi(poi, standalone = false)
    plot_vertices(vertex_positions, standalone = false)
    plot_edges(g, vertex_positions, standalone = false)

    if standalone
        xlabel("x")
        ylabel("y")
        # tight_layout()
        equalize_3d_aspect()
    end
end


function plot_cell(cell, vertex_positions; shrink = 0.8, label = "")
    vrts = deepcopy(vertex_positions[cell])

    xmean = mean([ vrts[ivrts][1] for ivrts in eachindex(vrts) ])
    ymean = mean([ vrts[ivrts][2] for ivrts in eachindex(vrts) ])
    zmean = mean([ vrts[ivrts][3] for ivrts in eachindex(vrts) ])

    shrunkvrts = deepcopy(vrts)

    for ivrts in eachindex(vrts)
        shrunkvrts[ivrts][1] = (vrts[ivrts][1] - xmean) * shrink + xmean
        shrunkvrts[ivrts][2] = (vrts[ivrts][2] - ymean) * shrink + ymean
        shrunkvrts[ivrts][3] = (vrts[ivrts][3] - zmean) * shrink + zmean
    end


    for k in 1:length(shrunkvrts)-1
        for k in 1:length(shrunkvrts)-1
            x = [ shrunkvrts[k][1], shrunkvrts[k + 1][1] ]
            y = [ shrunkvrts[k][2], shrunkvrts[k + 1][2] ]
            z = [ shrunkvrts[k][3], shrunkvrts[k + 1][3] ]

            arrow3d.add_arrow3D(x, y, z, color = "C0")
        end
    end

    if label != ""
        txt = gca().text(xmean, ymean, zmean, label,
            horizontalalignment = "center", verticalalignment = "center", color = "C0")
        txt.set_path_effects([patheffects.Stroke(linewidth = 3, foreground = "white"),
                           patheffects.Normal()])
    end
end


function plot_cells(cells, vertex_positions, currents = [])
    for icell in eachindex(cells)
        label = "$icell"
        if length(currents) > 0
            label *= ": $(signif(currents[icell], 6))"
        end
        plot_cell(cells[icell], vertex_positions, label = label, shrink = 0.9)
    end
end


function plot_edge_current(edge, current, vertex_positions)
    pos1 = vertex_positions[src(edge)]
    pos2 = vertex_positions[dst(edge)]
    arrow3d.add_arrow3D([pos1[1], pos2[1]], [pos1[2], pos2[2]], [pos1[3], pos2[3]],
        color = "black")

    txt = gca().text(mean([pos1, pos2])...,  "$(signif(current, 4))",
        horizontalalignment = "center", verticalalignment = "center")
    txt.set_path_effects([patheffects.Stroke(linewidth = 3, foreground = "white"),
                       patheffects.Normal()])
end


function plot_edge_currents(g, currents, vertex_positions)
    for edge in edges(g)
        plot_edge_current(edge, currents[edge], vertex_positions)
    end
end


function plot_loop(loop, vertex_positions; color = nothing, standalone = true)
    for k in 1:length(loop)-1
        v1 = loop[k]
        v2 = loop[k + 1]

        x = [ vertex_positions[v1][1], vertex_positions[v2][1] ]
        y = [ vertex_positions[v1][2], vertex_positions[v2][2] ]
        z = [ vertex_positions[v1][3], vertex_positions[v2][3] ]

        arrow3d.add_arrow3D(x, y, z, color = color)
    end
end


function equalize_3d_aspect()
    limmin = minimum([ l[1] for l in [xlim(), ylim(), zlim()] ])
    limmax = maximum([ l[2] for l in [xlim(), ylim(), zlim()] ])
    xlim(limmin, limmax)
    ylim(limmin, limmax)
    zlim(limmin, limmax)
end


function decompose_string(x, coeffs, elements)
    "$(signif(Float64(x), 6)) ≈ $(join([ "$(c)×$e" for (c, e) in zip(coeffs, elements)], " + "))"
end


function plot_simpleloop(loop, current, vertex_positions,
        elemcurrents, simpleloopscurrent_decomp)
    figure(figsize = (15,15))

    plot_vertices(vertex_positions)
    plot_loop(loop, vertex_positions, color = "C0")

    dstr = decompose_string(current, simpleloopscurrent_decomp, elemcurrents)
    title("""current: $(signif(Float64(current), 6))
             vertices: $loop
             decomposition: $dstr
          """, fontsize = "x-large")
end


function save_plots_simpleloops(folder, simpleloops, simpleloopscurrents,
        vertex_positions, elemcurrents, simpleloopscurrents_decomp)
    @showprogress 1 for i in eachindex(simpleloops)
        plot_simpleloop(simpleloops[i], simpleloopscurrents[i],
            vertex_positions, elemcurrents, simpleloopscurrents_decomp[i])
        savefig(joinpath(folder, "simpleloop$(lpad(i, 3, 0)).png"), dpi = 150)
        close("all")
    end
end


"dirs - what to plot, eg. [3, 1, 2] will put z (the 3rd dimension) on the
    x-axis of the plot, x (the 1st dimension) on the y axis and will keep
    y (the 2nd dimension) fixed.
    spanA is the span of the first of the plotted dimension (along the x axis
    of the resulting plot).
    spanB is the span of the second of the plotted dimension (along the y axis
    of the resulting plot)."
function plot_loops_field(loops, currents, vertex_positions,
        dirs, valC;
        n = 50, Bref = x -> [0, 0, 0],
        spanA = extrema([ p[dirs[1]] for p in vertex_positions]),
        spanB = extrema([ p[dirs[2]] for p in vertex_positions]),
        levels = -5:5)

    As = linspace(spanA..., n)
    Bs = linspace(spanB..., n)

    B = zeros(length(As), length(Bs), 3)
    for iA in eachindex(As), iB in eachindex(Bs)
        p = [As[iA], Bs[iB], valC]
        Bp = field_loops(p, loops, currents, vertex_positions) .- Bref(p)
        B[iA, iB, :] .= Bp
    end

    f, axs = subplots(1, 4, figsize = (10,4), gridspec_kw = Dict("width_ratios" => [1, 1, 1, 0.1]))
    for i in 1:3
        sca(axs[i])
        contourf(As, Bs, B[:,:,i]' * 1e6, levels = levels, extend = "both", cmap="RdBu_r")
        gca().set_aspect("equal")
        xlabel("xyz"[dirs[1:1]])
        ylabel("xyz"[dirs[2:2]])
        dirletter = "xyz"[i]
        title("\$B_$dirletter\$")
        xlim(spanA...)
        ylim(spanB...)
    end

    colorbar(cax = axs[4], label = "field deviation from goal (μT)")
    # tight_layout()
end


function plot_deviation_histogram(poi, loops, loopscurrents, vertex_positions,
        Bgoal)

    deviations = vcat([
        field_loops(p, loops, loopscurrents, vertex_positions) .- Bgoal(p)
        for p in poi ]...)
    plt.hist(deviations * 1e6, bins = 100)
    yscale("log", nonposy="clip")
    xlabel("deviation from goal (μT)")
    ylabel("# POI")
end


function printlnflush(s)
    println(s)
    flush(STDOUT)
end


function save_report(folder, vertex_positions, g, poi, Bgoal, simpleloops,
        simpleloopscurrents, elemcurrents, simpleloopscurrents_decomp;
        extension = "png", dpi = 150,
        levels = -5:5)
    isdir(folder) || mkdir(folder)

    printlnflush("Plotting the geometry of the system...")

    figure(figsize = (15,15))
    plot_system(g, vertex_positions, poi)
    savefig(joinpath(folder, "system.$extension"), dpi = dpi)
    close("all")

    figure(figsize = (15,15))
    plot_system(g, vertex_positions, poi)
    gca().azim = -90
    gca().elev = 90
    savefig(joinpath(folder, "system_top.$extension"), dpi = dpi)
    close("all")

    figure(figsize = (15,15))
    plot_system(g, vertex_positions, poi)
    gca().azim = -90
    gca().elev = 0
    savefig(joinpath(folder, "system_front.$extension"), dpi = dpi)
    close("all")

    figure(figsize = (15,15))
    plot_system(g, vertex_positions, poi)
    gca().azim = 0
    gca().elev = 0
    savefig(joinpath(folder, "system_right.$extension"), dpi = dpi)
    close("all")

    simpleloopscurrents_real = real_decomposed_currents(
        simpleloopscurrents_decomp, elemcurrents)

    printlnflush("Plotting the deviation histogram...")
    figure()
    plot_deviation_histogram(poi, simpleloops, simpleloopscurrents_real,
        vertex_positions, Bgoal)
    savefig(joinpath(folder, "deviation_histogram.$extension"), dpi = dpi)
    close("all")

    printlnflush("Plotting the field of the solution...")

    zspan = extrema([ p[3] for p in vertex_positions])
    for z in linspace(zspan..., 5)[2 : end-1]
        plot_loops_field(simpleloops, simpleloopscurrents_real, vertex_positions,
            [1,2,3], z, n = 50, Bref = Bgoal, levels = levels)
        savefig(joinpath(folder, "field_XY_z$(signif(z,2)).$extension"),
            dpi = dpi)
        close("all")
    end

    xspan = extrema([ p[1] for p in vertex_positions])
    for x in linspace(xspan..., 5)[2 : end-1]
        plot_loops_field(simpleloops, simpleloopscurrents_real, vertex_positions,
            [3,2,1], x, n = 50, Bref = Bgoal, levels = levels)
        savefig(joinpath(folder, "field_ZY_x$(signif(x,2)).$extension"),
            dpi = dpi)
        close("all")
    end

    printlnflush("Plotting the simple loops...")

    save_plots_simpleloops(folder, simpleloops, simpleloopscurrents,
        vertex_positions, elemcurrents, simpleloopscurrents_decomp)
end


end

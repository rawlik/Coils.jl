__precompile__()

module Coils

using LightGraphs
using ProgressMeter


export cuboid_system, getedgei, getedge, find_cells, biotsavart,
    biotsavart_cell, system_matrix, find_edgecurrents, simplify_current_graph,
    find_simpleloop, find_all_simpleloops, solve_system, cuboid_poi,
    decompose, decompose_currents, field_loops, μ0, order_cell,
    real_decomposed_currents, save_result



"skipfaces is in the order [x-, x+, y-, y+, z-, z+]"
function cuboid_system(totalsize, ntiles; skipfaces = falses(6))
    gridplanes = linspace.(-totalsize / 2, totalsize / 2, ntiles + 1)
    vertex_positions = []

    # construct the vertices
    for z in gridplanes[3], y in gridplanes[2], x in gridplanes[1]
        any(in.([x,y,z], extrema.(gridplanes))) || continue

        # skip the inside of the X- face
        skipfaces[1] &&
            x == minimum(gridplanes[1]) &&
            !(y in extrema(gridplanes[2])) &&
            !(z in extrema(gridplanes[3])) &&
            continue

        # skip the inside of the X+ face
        skipfaces[2] &&
            x == maximum(gridplanes[1]) &&
            !(y in extrema(gridplanes[2])) &&
            !(z in extrema(gridplanes[3])) &&
            continue

        # skip the inside of the Y- face
        skipfaces[3] &&
            y == minimum(gridplanes[2]) &&
            !(x in extrema(gridplanes[1])) &&
            !(z in extrema(gridplanes[3])) &&
            continue

        # skip the inside of the Y+ face
        skipfaces[4] &&
            y == maximum(gridplanes[2]) &&
            !(x in extrema(gridplanes[1])) &&
            !(z in extrema(gridplanes[3])) &&
            continue

        # skip the inside of the Z- face
        skipfaces[5] &&
            z == minimum(gridplanes[3]) &&
            !(x in extrema(gridplanes[1])) &&
            !(y in extrema(gridplanes[2])) &&
            continue

        # skip the inside of the Z+ face
        skipfaces[6] &&
            z == maximum(gridplanes[3]) &&
            !(x in extrema(gridplanes[1])) &&
            !(y in extrema(gridplanes[2])) &&
            continue

        push!(vertex_positions, [x, y, z])
    end

    # add the edges
    g = DiGraph()
    add_vertices!(g, length(vertex_positions))

    tile_size = mean.(diff.(gridplanes))
    for i in eachindex(vertex_positions)
        for j in (i+1):length(vertex_positions)
            # the "mean(tile_size) is not the best, at the moment. Only works for equal tiles
            if abs(norm(vertex_positions[i] .- vertex_positions[j]) - mean(tile_size)) < 0.01
                add_edge!(g, i, j)
            end
        end
    end
    g, vertex_positions
end


"Identify the edge between vertices a and b (in any direction)."
function getedgei(g, a, b)
    for (i, e) in enumerate(edges(g))
        if (dst(e) == a) & (src(e) == b)
            return i
        elseif (dst(e) == b) & (src(e) == a)
            return i
        end
    end
    return 0
end


"Identify the edge between vertices a and b (in any direction)."
function getedge(g, a, b)
    for e in edges(g)
        if (dst(e) == a) & (src(e) == b)
            return e
        elseif (dst(e) == b) & (src(e) == a)
            return e
        end
    end
    return 0
end


"""Find elementary cells in a graph.
Takes a very long time to find large cells (longer then 13-15 edges).
In those cases large cells can be supplied manually, which drastically speeds up the process.
"""
function find_cells(g; maxpathlength = Inf, cells = [])
    # get all cells in the raster
    edgetaken = zeros(Int, length(edges(g)))

    # mark the manually supplied cells
    for cell in cells
        for k in 1:(length(cell) - 1)
            # mark that the edge is part of one additional cell
            edgetaken[getedgei(g, cell[k], cell[k + 1])] += 1
        end
    end

    # an algorithm to decompose a graph into "cells"
    # we traverse the graph breadth-first, so that the first cell we find
    # has to be the shortest one (smallest cell, so elementary)

    @showprogress 1 for i in eachindex(vertices(g))
        #println("Starting vertex: $i")
        paths = [[i]]
        foundnow = 0

        while !isempty(paths)
            # get all cells including that vertex
            # take the last added path
            path = shift!(paths)
            length(path) > maxpathlength && break

            for n in [ collect(out_neighbors(g, path[end])); collect(in_neighbors(g, path[end])) ]
                # try to extend the path
                if edgetaken[getedgei(g, path[end], n)] >= 2
                    # edge already a member of two cells
                    continue
                elseif (length(path) > 2) && (n == path[1]) && !(sort([path; n]) in sort.(cells))
                    # we are back to the starting point
                    # and the path is not a cell yet
                    for k in 1:(length(path) - 1)
                        # mark that the edge is part of one additional cell
                        edgetaken[getedgei(g, path[k], path[k + 1])] += 1
                    end
                    edgetaken[getedgei(g, path[end], n)] += 1

                    # save the cell
                    push!(cells, [path; n])
                    foundnow += 1
                    break
                elseif !(n in path)
                    # we have just extended the path, maybe it will close to be
                    # a cell further on
                    push!(paths, [path; n])
                end
            end

            # a very generous stop condition
            # if the number of cells found for this vertex is at least number of its neighbors
            # typically can stop earlier
            total_neighbors = (length(out_neighbors(g, i)) + length(in_neighbors(g, i)))
            already_in = count(cell -> i in cell, cells)
            already_in >= total_neighbors && break
        end
    end
    cells
end


"Sine of the angle between vectors"
sinvectors(x1, x2) = norm(x1 × x2) / (norm(x1) * norm(x2))

const μ0 = 4π * 1e-7

function biotsavart(x1, x2, p)
    # n is a normal vector along the wire
    n = normalize(x2 .- x1)

    # d is a vector perpendicular to the wire through the point p
    # ref. https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line#Vector_formulation
    # ref. https://de.wikipedia.org/wiki/Biot-Savart-Gesetz#Gerader_Linienleiter
    ρ = (x1 .- p) .- ((x1 .- p) ⋅ n) * n

    # In the formula on wikipedia the angles are directed.
    # We need to take it into account.
    # s1 = norm(p .+ ρ .- x2) > norm(x1 .- x2) ? 1 : -1
    # s2 = norm(p .+ ρ .- x1) > norm(x1 .- x2) ? 1 : -1
    # normB = μ0 * abs(s1 * sinvectors(ρ, x2 - p) - s2 * sinvectors(ρ, x1 - p)) / (4π * norm(ρ))

    # check with s whether poit p+ρ lies on the wire
    s = norm(p .+ ρ .- (x1 .+ x2) / 2) > norm(x1 .- x2) / 2 ? -1 : 1
    normB = μ0 * abs(sinvectors(ρ, x2 - p) + s * sinvectors(ρ, x1 - p)) / (4π * norm(ρ))

    B = normalize((ρ × n)) * normB
end


function numbiotsavart(x1, x2, p; n = 1000)
    B = zeros(3)
    dl = (x2 .- x1) / n

    for x in [ x1 .+ m * dl for m in 0.5:n ]
        r = p .- x
        B += μ0 / 4π * (dl × r) / norm(r)^3
    end
    B
end


"Magnetic field in point p produced by ith cell"
function biotsavart_cell(p, vertex_positions, cell)
    vrts = [ vertex_positions[k] for k in cell ]
    sum([ biotsavart(vrts[i], vrts[i + 1], p) for i in 1:length(vrts)-1 ])
end


function field_loops(p, loops, currents, vertex_positions)
    field = zeros(3)
    for (loop, current) in zip(loops, currents)
        for iv in 1 : length(loop)-1
            x1 = vertex_positions[loop[iv]]
            x2 = vertex_positions[loop[iv+1]]
            field .+= biotsavart(x1, x2, p) * current
        end
    end
    field
end


function system_matrix(poi, vertex_positions, cells)
    M = hcat([
        # for the ith cell the proportionality constant for the field it produces
        # in every point of interest. x, y and z are simply concatenated in an xyzxyzxyz way
        vcat([ biotsavart_cell(p, vertex_positions, cell) for p in poi ]...)
        # this is done for each cell
        for cell in cells ]...)
end


function find_edgecurrents(g, cells, optI; bigfloat = false)
    edgecurrents = Dict([
        (edge, bigfloat ? BigFloat(0.0) : 0.0)
        for edge in edges(g) ])

    for edge in edges(g)
        for icell in eachindex(cells)
            # for each cell
            vrts = cells[icell]
            for i in 1:length(vrts)-1
                # for each edge of the cell
                if vrts[i] == src(edge) && vrts[i + 1] == dst(edge)
                    # the cell is aligned with the edge
                    edgecurrents[edge] += optI[icell]
                end
                if vrts[i] == dst(edge) && vrts[i + 1] == src(edge)
                    # the cell is anti-aligned with the edge
                    edgecurrents[edge] -= optI[icell]
                end
            end
        end
    end

    edgecurrents
end


function simplify_current_graph(g, edgecurrents; tolerance = 1e-13)
    gs = DiGraph()
    add_vertices!(gs, nv(g))

    edgecurrentss = empty!(Dict(edgecurrents))

    for edge in edges(g)
        if abs(edgecurrents[edge]) < tolerance
            # skip edges below the machine precision
            continue
        elseif edgecurrents[edge] > 0
            # put this edge on the new graph
            add_edge!(gs, src(edge), dst(edge))
            edgecurrentss[getedge(gs, src(edge), dst(edge))] = abs(edgecurrents[edge])
        elseif edgecurrents[edge] < 0
            # put this edge on the new graph, but invert it first
            add_edge!(gs, dst(edge), src(edge))
            edgecurrentss[getedge(gs, dst(edge), src(edge))] = abs(edgecurrents[edge])
        end
    end

    gs, edgecurrentss
end


function findsimpleloop(g, edgecurrents; maxcurrent = Inf, maxpathlength = Inf,
        verbose = false, tolerance = 1e-13)

    for current in sort(sort(collect(values(edgecurrents))), rev = true)
        current <= nextfloat(maxcurrent) || continue
        verbose && println("    Trying current $current, number of edges: $(length(edges(g)))")
        # start from the highest current and gradually reduce it
        for i in vertices(g)
            # try to find a loop with at least this high current
            paths = [[i]]

            safetycounter = 0
            while !isempty(paths)
                # get all cells including that vertex
                # take the last added path
                path = shift!(paths)
                length(path) > maxpathlength && break
                verbose && (safetycounter += 1) % 10_000 == 0 && println("  Pathlength: $(length(path))")

                # only go in the direction of positive current, so outcoming edges
                for n in out_neighbors(g, path[end])
                    # add a safety margin for threshold
                    if edgecurrents[getedge(g, path[end], n)] <= (current - tolerance)
                        # cannot go here, too high current
                        continue
                    end
                    # try to extend the path
                    if (length(path) > 2) && (n == path[1])
                        # we are back to the starting point
                        # and the path is not a cell yet
                        # Return the current and the simple loop
                        verbose && println("  Found a simple loop with current $current : $([path; n])")
                        return (current, [path; n])
                    elseif !(n in path)
                        # we have just extended the path, maybe it will close to be
                        # a cell further on
                        push!(paths, [path; n])
                    end
                end
            end
        end
    end

    # could not find a loop
    NaN, []
end


function find_all_simpleloops(g, edgecurrents; tolerance = 1e-10,
        mincurrent = 0, verbose = false)

    simpleloops = []
    simpleloopscurrents = []

    gs, edgecurrentss = simplify_current_graph(g, edgecurrents,
        tolerance = tolerance)

    # Only look for loops with currents not smaller then the current
    # of the most recently found loop.
    # current = Inf
    while length(edges(gs)) > 0
        current, simpleloop = findsimpleloop(gs, edgecurrentss,
            # maxcurrent = current,
            verbose = verbose,
            tolerance = tolerance)
        if isnan(current)
            println("Could not find a loop, stopping.")
            println("   The graph has $(length(edges(g))) edges left,")
            println("   the maximum current left is $(maximum(abs, values(edgecurrentss)))")
            break
        end

        abs(current) < mincurrent && break

        push!(simpleloopscurrents, current)
        push!(simpleloops, simpleloop)

        # subtract the simple loop from the graph
        for i in 1:length(simpleloop)-1
            edge = getedge(gs, simpleloop[i], simpleloop[i+1])

            if src(edge) == simpleloop[i]
                # the edge is aligned with the simple loop
                edgecurrentss[ edge ] -= current
            elseif dst(edge) == simpleloop[i]
                # the edge is anti-aligned with the simple loop
                edgecurrentss[ edge ] += current
            end
        end

        gs, edgecurrentss = simplify_current_graph(gs, edgecurrentss,
            tolerance = tolerance)
    end

    simpleloops, simpleloopscurrents
end


function printlnflush(s)
    println(s)
    flush(STDOUT)
end

function solve_system(g, vertex_positions, poi, B; initialcells = [],
        verbose = false, tolerance = 1e-10, bigfloat = false, λ = 0,
        mincurrent = 0)

    verbose && printlnflush("Looking for cells in the graph...")
    cells = find_cells(g, cells = initialcells)

    # here the system's matrix is constructed
    verbose && printlnflush("Calculating the system's matrix...")
    M = system_matrix(poi, vertex_positions, cells)

    # the goal field values at each of poi
    # x, y and z are simply concatenated in an xyzxyzxyz way
    Bpoi = vcat(B.(poi)...)

    # calculate the optimal currents for the given goal field
    # \ returns the least-norm solution already (so least current) in case of ambiguity
    # We are using Tikchonov regularization
    verbose && printlnflush("Solving for currents in the cells...")
    optI = [M; λ * eye(size(M,2))] \ [Bpoi; zeros(size(M,2))]

    verbose && printlnflush("Finding the total currents along the edges...")
    edgecurrents = find_edgecurrents(g, cells, optI, bigfloat = bigfloat)

    verbose && printlnflush("Decomposing the solution into simple loops...")
    simpleloops, simpleloopscurrents = find_all_simpleloops(g, edgecurrents,
        tolerance = tolerance, mincurrent = mincurrent)
end


function cuboid_poi(voi_size = [1.0, 1.0, 1.0], voi_centre = [0.0, 0.0, 0.0],
        n = [9, 9, 9]; filled = true)
    voi_walls_positions = voi_size / 2
    voi_pos = linspace.(-voi_walls_positions, voi_walls_positions, n)
    thr = voi_walls_positions .- mean.(diff.(voi_pos)) / 2

    poi = [ [x, y, z] .+ voi_centre for x in voi_pos[1], y in voi_pos[2], z in voi_pos[3] if
            filled ||
            # faces
            ( ~(-thr[1] ≤ x ≤ thr[1]) && (-thr[2] ≤ y ≤ thr[2]) && (-thr[3] ≤ z ≤ thr[3]) ) ||
            ( ~(-thr[2] ≤ y ≤ thr[2]) && (-thr[1] ≤ x ≤ thr[1]) && (-thr[3] ≤ z ≤ thr[3]) ) ||
            ( ~(-thr[3] ≤ z ≤ thr[3]) && (-thr[1] ≤ x ≤ thr[1]) && (-thr[2] ≤ y ≤ thr[2]) ) ||
            # edges
            ( (-thr[1] ≤ x ≤ thr[1]) && ~(-thr[2] ≤ y ≤ thr[2]) && ~(-thr[3] ≤ z ≤ thr[3]) ) ||
            ( (-thr[2] ≤ y ≤ thr[2]) && ~(-thr[1] ≤ x ≤ thr[1]) && ~(-thr[3] ≤ z ≤ thr[3]) ) ||
            ( (-thr[3] ≤ z ≤ thr[3]) && ~(-thr[1] ≤ x ≤ thr[1]) && ~(-thr[2] ≤ y ≤ thr[2]) ) ||
            # corners
            ( ~(-thr[3] ≤ z ≤ thr[3]) && ~(-thr[1] ≤ x ≤ thr[1]) && ~(-thr[2] ≤ y ≤ thr[2]) ) ]
end

"""Decompose a number into given factors.
decompose(823.2, [100., 10., 1.]) = [8, 2, 3]
if roundlast is true, the last element is rounded:
decompose(84.86, [10., 1., 0.1], roundlast = true) = [8, 4, 9] """
function decompose{T<:AbstractFloat}(x::T, elements::Vector{T}; roundlast = false)
    coeffs = similar(elements)

    for i in eachindex(elements)
        coeffs[i] = x
        for j in 1 : i-1
            # use prevfloat for numeric reasons (eg. 4.0 ÷ 0.2 = 19.0)
            # but 4.0 ÷ prevfloat(0.2) = 20.0
            coeffs[i] %= prevfloat(elements[j])
        end
        if roundlast && i == length(elements)
            coeffs[i] /= elements[i]
        else
            coeffs[i] ÷= prevfloat(elements[i])
        end
    end

    round.(Int, coeffs)
end

function decompose(x::Real, elements::AbstractArray; kwargs...)
    xt, elementst = promote([x], elements)
    decompose(xt[], elementst; kwargs...)
end


function decompose_currents(currents, elemcurrents)
    [ decompose(c, elemcurrents, roundlast = true) for c in currents ]
end


function real_decomposed_currents(simpleloopscurrents_decomp, elemcurrents)
    [ sum(decomp .* elemcurrents) for decomp in simpleloopscurrents_decomp ]
end


"Given the set of points that would make up a cell, order them into a proper cell."
function order_cell(g, face_unordered)
    face = [ face_unordered[1] ]

    while true
        for n in [in_neighbors(g, face[end]); out_neighbors(g, face[end])]
            n in face_unordered || continue
            n in face && continue
            push!(face, n)
            break
        end
        length(face) >= length(face_unordered) && break
    end

    push!(face, face[1])
    face
end


function decompose_string(x, coeffs, elements)
    "$(signif(Float64(x), 6)) = $(join([ "$(c)*$e" for (c, e) in zip(coeffs, elements)], " + "))"
end


function save_result(filename, simpleloops, simpleloopscurrents,
                     simpleloopscurrents_decomp, elemcurrents)
    # truncate the contents of the output file
    open(filename, "w") do file
        for i in eachindex(simpleloopscurrents)
            write(file, i, "  :  ")
            write(file, decompose_string(simpleloopscurrents[i],
                                         simpleloopscurrents_decomp[i],
                                         elemcurrents))
            write(file, "  :  ")
            write(file, repr(simpleloops[i]))
            write(file, "\n")
        end
    end
end


end

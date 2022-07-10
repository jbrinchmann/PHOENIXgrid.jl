#
# Routines for visualisation of the results of fitting.
#

# These need to be imported for Plots.create_grid to work properly
import Plots: GridLayout, EmptyLayout

function _make_triangle_string(N)

    txt = "["
    for i=1:N
        line = repeat("x ", i)*repeat("_ ", (N-i))
        if i < N
            txt = txt*line*"\n"
        else
            txt = txt*line*"]\n"
        end
    end
    return txt
end


"""Show a corner plot for the given variables.
"""
function corner(dims, r; to_show=["v", "Teff", "FeH"], vaxis=nothing)

    Nshow = length(to_show)

    if !(:v in keys(dims))
        if vaxis == nothing
            println("If you provide a dimension dict without a velocity axis, you need to provide vaxis!")
            return nothing
        end
        dims[:v] = vaxis
    end
    
    order = Dict("v" => 1, "logg" => 2, "Teff" => 3, "FeH" => 4, "Z" => 4, "alpha" => 5)
    problem = false
    for s ∈ to_show
        if !(s ∈ keys(order))
            println("Unknown key: $s")
            problem = true
        end
    end
    if problem
        return nothing
    end

    m = marginalisation_keymap()
    
    # Order the corner plot in the way we want.
    sort!(to_show, by=x->order[x])
    
    
    # We use an underlying grid to create the triangle plot
    txt = _make_triangle_string(Nshow)

    ex = Meta.parse(txt)
    l = eval(Plots.create_grid(ex))

    plots = []
    for i=1:Nshow
        for j=1:Nshow
            if j > i
                continue
            end
            if i==j
                # 1D histograms
                what = m[to_show[i]]
                xaxis = dims[Symbol(to_show[i])]
                p = plot(xaxis, r[what], legend=false)
            else
                key = to_show[j]*to_show[i]
                println("I am looking up $key")
                what = m[key]
                xaxis = dims[Symbol(to_show[i])]
                yaxis = dims[Symbol(to_show[j])]
                p = heatmap(yaxis, xaxis, transpose((r[what]).^0.1),
                            legend=false)
            end

            push!(plots, p)
        end
    end

    plot(plots..., layout=l)
    
end

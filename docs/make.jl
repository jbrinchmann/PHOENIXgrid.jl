using Documenter, PHOENIXgrid

makedocs(
    modules = [PHOENIXgrid],
    format = Documenter.HTML(; prettyurls = get(ENV, "CI", nothing) == "true"),
    authors = "Jarle Brinchmann",
    sitename = "PHOENIXgrid.jl",
    pages = Any["index.md"]
    # strict = true,
    # clean = true,
    # checkdocs = :exports,
)

deploydocs(
    repo = "github.com/jbrinchmann/PHOENIXgrid.jl.git",
    push_preview = true
)

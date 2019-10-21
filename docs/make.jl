push!(LOAD_PATH,"../")
using Documenter, InterpolatedPDFs

makedocs(sitename="Interpolated PDFs",
         format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
         modules = [InterpolatedPDFs]
        )

deploydocs(
    repo = "github.com/m-wells/InterpolatedPDFs.jl.git",
)

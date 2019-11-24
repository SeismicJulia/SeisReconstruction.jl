using Pkg; Pkg.add("Documenter")
using Documenter, SeisReconstruction

makedocs(
    sitename = "SeisReconstruction.jl",
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    modules = [SeisReconstruction],
    pages = [
        "Home" => "index.md",
	"Library" => Any[
			"Public" => "lib/public.md",
],
],
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.

deploydocs(
	repo = "github.com/SeismicJulia/SeisReconstruction.jl.git",
        target = "build",
	deps = nothing,
	make = nothing
)

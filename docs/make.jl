using WignerFamilies
using Documenter

makedocs(;
    modules=[WignerFamilies],
    authors="xzackli <xzackli@gmail.com> and contributors",
    repo="https://github.com/xzackli/WignerFamilies.jl/blob/{commit}{path}#L{line}",
    sitename="WignerFamilies.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://xzackli.github.io/WignerFamilies.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/xzackli/WignerFamilies.jl",
)

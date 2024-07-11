using EchelleInstruments
using Documenter

makedocs(;
    modules=[EchelleInstruments],
    authors="Eric Ford",
    repo="https://github.com/RvSpectML/EchelleInstruments.jl/blob/{commit}{path}#L{line}",
    sitename="EchelleInstruments.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://RvSpectML.github.io/EchelleInstruments.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Instrument agnostic" => "agnostic.md",
        "Instruments" => [
            "NEID" => "neid.md",
            "EXPRES" => "expres.md"],
        "Index" => "longlist.md",
    ],
    checkdocs=:none,
    #checkdocs=:exports,
    )

deploydocs(;
    repo="github.com/RvSpectML/EchelleInstruments.jl",
)

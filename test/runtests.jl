using CLSGCTGMT

using Aqua: test_all

using Test: @test

# ----------------------------------------------------------------------------------------------- #

test_all(CLSGCTGMT; deps_compat = false, ambiguities = false)

# ---- #

const DA = pkgdir(CLSGCTGMT, "data")

# ---- #

for (cl, ta, re) in (
    (
        "CCLE_mRNA_20Q2_no_haem_phen.cls",
        "HER2",
        [1.087973, -1.358492, -1.178614, -0.77898, 0.157222, 1.168224, -0.360195, 0.608629],
    ),
    ("GSE76137.cls", "Proliferating_Arrested", [1, 2, 1, 2, 1, 2]),
    ("LPS_phen.cls", "CNTRL_LPS", [1, 1, 1, 2, 2, 2]),
)

    cl = joinpath(DA, cl)

    da = CLSGCTGMT.read_cls(cl)

    @test names(da)[1] === "Target"

    @test da[:, 1] == [ta]

    @test all(startswith("Sample "), names(da)[2:end])

    @test eltype(Matrix(da[:, 2:end])) === eltype(re)

    @test collect(da[1, eachindex(re) .+ 1]) == re

end

# ---- #

@test size(CLSGCTGMT.read_gct(joinpath(DA, "a.gct"))) === (13321, 190)

# ---- #

for (gm, re) in (("h.all.v7.1.symbols.gmt", 50), ("c2.all.v7.1.symbols.gmt", 5529))

    @test length(CLSGCTGMT.read_gmt(joinpath(DA, gm))) === re

end

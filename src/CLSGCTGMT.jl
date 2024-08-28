module CLSGCTGMT

using LeMoIO

function read_cls(cl)

    l1, l2, l3 = readlines(cl)

    np = l2[2:end]

    va_ = split(l3)

    us = lastindex(va_)

    sa_ = map(id -> "Sample $id", 1:us)

    if l1 == "#numeric"

        LeMoIO.tabulate("Target", np, sa_, [parse(Float64, st) for _ = 1:1, st in va_])

    else

        l1_ = split(l1)

        u1 = parse(Int, l1_[1])

        if u1 != us

            error("numbers of samples differ: $u1 and $us.")

        end

        u1 = parse(Int, l1_[2])

        gr_ = split(np)

        ug = lastindex(gr_)

        uu = lastindex(unique(va_))

        if !(u1 == ug == uu)

            error("numbers of groups differ: $u1, $ug, and $uu.")

        end

        gr_id = Dict(gr => id for (id, gr) in enumerate(gr_))

        LeMoIO.tabulate("Target", join(gr_, '_'), sa_, [gr_id[st] for _ = 1:1, st in va_])

    end

end

function read_gct(gc)

    LeMoIO.read_table(gc; header = 3, drop = ["Description"])

end

function read_gmt(gm)

    se_ge_ = Dict{String, Vector{String}}()

    for li in eachline(gm)

        sp_ = split(li, '\t')

        se = sp_[1]

        if haskey(se_ge_, se)

            error()

        end

        se_ge_[se] = filter!(!isempty, sp_[3:lastindex(sp_)])

    end

    se_ge_

end

end

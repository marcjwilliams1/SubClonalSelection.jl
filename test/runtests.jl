using SubClonalSelection
using Base.Test
using Plots

#useful to check what the backend is, seems to cause some problems
println(backend())
tests = ["neutral", "1clone"]

println("Running tests ...")

for t in tests
    fn = "test_$t.jl"
    println("* $fn ...")
    include(fn)
end

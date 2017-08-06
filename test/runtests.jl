using SubClonalSelection
using Base.Test

tests = ["neutral", "1clone"]

println("Running tests ...")

for t in tests
    fn = "test_$t.jl"
    println("* $fn ...")
    include(fn)
end

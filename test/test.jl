push!(LOAD_PATH, "/Users/gangchen");
using MultiResidue
using Test
using SymEngine
@vars x y z

dis=solve([x+y-1,2x-y],[x,y])

ideals=[y - 2x + z, x^2 *y - x*z^2 + z^3, y^4 + x*y*z^2 - y^2*z^2 + z^4];
vars=[x,y,z]


@testset "MultiResidue.jl" begin
    @test MultiResidue.degree(x^2+y^3) == 3
    @test multiResidue((2x + 3y + 4z)/(z - 2),ideals,vars)==Basic(93)/Basic(1664)
    @test solve([x+y-1,2x-y],[x,y])==dis
end
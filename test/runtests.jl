using Test

include("../src/Trojan.jl")
using .Trojan

@testset "generate_pyt_triple" begin
	triple1 = generate_pyt_triple(2, 1)
	triple2 = generate_pyt_triple(8, 1)

	@test triple1.a == 3
	@test triple1.b == 4
	@test triple1.c == 5
	@test triple2.a == 16
	@test triple2.b == 63
	@test triple2.c == 65
end

@testset "generate_trojan_triple_120" begin
	triple1 = generate_trojan_triple_120(3, 1)
	triple2 = generate_trojan_triple_120(5, 2)

	@test triple1.a == 3
	@test triple1.b == 5
	@test triple1.c == 7
	@test triple2.a == 5
	@test triple2.b == 16
	@test triple2.c == 19

end

@testset "calc_angle_by_cos_law" begin
	angle1 = calc_angle_by_cos_law(3, 4, 5)
	angle2 = calc_angle_by_cos_law(5, 12, 13)
	angle3 = calc_angle_by_cos_law(7, 24, 25)

	@test isapprox(angle1, 36.87, atol = 0.01)
	@test isapprox(angle2, 22.62, atol = 0.01)
	@test isapprox(angle3, 16.26, atol = 0.01)
end

@testset "generate_trojan_triple_vector" begin
    triple1 = generate_trojan_triple_vector(3)
    triple2 = generate_trojan_triple_vector(15)

    @test triple1[1].a == 3
    @test triple1[1].b == 5
    @test triple1[1].c == 7
    @test triple2[1].a == 3
    @test triple2[1].b == 5
    @test triple2[1].c == 7
    @test triple2[2].a == 7
    @test triple2[2].b == 8
    @test triple2[2].c == 13
    @test triple2[end].a == 29
    @test triple2[end].b == 195
    @test triple2[end].c == 211
end

@testset "get_coordinates_point_C_of_laying_triangle" begin
    koord1 = get_coordinates_point_C_of_laying_triangle(3, 4, 5)
    @test isapprox(koord1.x, 0, atol = 0.01)
    @test isapprox(koord1.y, 3, atol = 0.01)
end
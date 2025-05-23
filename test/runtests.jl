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
	koord1 = get_coordinates_point_C_of_laying_triangle(3, 5, 4)
	@test isapprox(koord1.x, 4, atol = 0.01)
	@test isapprox(koord1.y, 3, atol = 0.01)
end

@testset "rotate_point_by_angle" begin
	point = (x = 0.0, y = 3.0)

	koord1 = rotate_point_by_angle(point, 90.0)
	koord2 = rotate_point_by_angle(point, -90.0)
	@test isapprox(koord1.x, -3, atol = 0.01)
	@test isapprox(koord1.y, 0, atol = 0.01)
	@test isapprox(koord2.x, 3, atol = 0.01)
	@test isapprox(koord2.y, 0, atol = 0.01)
end

@testset "move_point" begin
	point = (x = 0.0, y = 3.0)
	transvektor = (x = 3.0, y = -4.0)
	koord1 = move_point(point, transvektor)
	@test isapprox(koord1.x, 3, atol = 0.01)
	@test isapprox(koord1.y, -1, atol = 0.01)

end

@testset "stretch_point" begin
	point = (x = 1.0, y = 3.0)
	koord1 = stretch_point(point, 2, 1 / 3)
	@test isapprox(koord1.x, 2, atol = 0.01)
	@test isapprox(koord1.y, 1, atol = 0.01)
end

@testset "get_points_normalized" begin
	triple = (a = 3, b = 5, c = 7)
	points_normalized = get_points_normalized(triple)
	@test isapprox(points_normalized.point1.x, 0, atol = 0.01)
	@test isapprox(points_normalized.point1.y, 0, atol = 0.01)
	@test isapprox(points_normalized.point2.x, 1, atol = 0.01)
	@test isapprox(points_normalized.point2.y, 0, atol = 0.01)
	@test isapprox(points_normalized.point3.x, 1.163265, atol = 0.01)
	@test isapprox(points_normalized.point3.y, 1.131071, atol = 0.01)
	@test isapprox(points_normalized.point4.x, 0.163265, atol = 0.01)
	@test isapprox(points_normalized.point4.y, 1.131071, atol = 0.01)
	@test isapprox(points_normalized.point5.x, 0.5, atol = 0.01)
	@test isapprox(points_normalized.point5.y, 0.866025, atol = 0.01)
end

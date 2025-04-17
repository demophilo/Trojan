module Trojan
using StatsBase # Für countmap
using Printf # Für @printf

export generate_pyt_triple,
	generate_trojan_triple_120,
	calc_angle_by_cos_law,
	generate_trojan_triple_vector,
	get_ext_trojan_triple_vector,
	add_angles_to_triple_vector,
	get_trojan_triples_for_a_number_vector,
	expand_trojan_triple_vector,
	analyze_c_frequencies,
	get_coordinates_point_C_of_laying_triangle,
	rotate_point_by_angle,
	move_point,
	stretch_point,
	get_points_normalized,
	draw_custom_lines

"""
	generate_pyt_triple(big_num::Int, small_num::Int)::NamedTuple{(:a, :b, :c), Tuple{Int, Int, Int}}

generates a named pythagorean triple
Input: two intergers
Output: sorted named pythagorean triple
"""
function generate_pyt_triple(big_num::Int, small_num::Int)::NamedTuple{(:a, :b, :c), Tuple{Int, Int, Int}}
	a::Int = big_num^2 - small_num^2
	b::Int = big_num^2 + small_num^2
	c::Int = 2 * big_num * small_num
	gcd_abc = foldl(gcd, [a, b, c])
	if gcd_abc > 1
		a ÷= gcd_abc
		b ÷= gcd_abc
		c ÷= gcd_abc
	end
	triple = sort([a, b, c])
	return (a = triple[1], b = triple[2], c = triple[3])
end


"""
	generate_trojan_triple_120(big_num::Int, small_num::Int)::NamedTuple{(:a, :b, :c), Tuple{Int, Int, Int}}

generates a trojan triple with one angle of 120°
Input: two intergers
Output: sorted named trojan triple
"""
function generate_trojan_triple_120(big_num::Int, small_num::Int)::NamedTuple{(:a, :b, :c), Tuple{Int, Int, Int}}
	a::Int = big_num^2 + small_num^2 - big_num * small_num
	b::Int = abs(big_num^2 - 2 * big_num * small_num)
	c::Int = abs(small_num^2 - 2 * big_num * small_num)
	gcd_abc = foldl(gcd, [a, b, c])
	if gcd_abc > 1
		a ÷= gcd_abc
		b ÷= gcd_abc
		c ÷= gcd_abc
	end
	triple_new = sort([a, b, c])
	return (a = triple_new[1], b = triple_new[2], c = triple_new[3])
end

"""
	calc_angle_by_cos_law(a, b, c)

input: 3 edges of a triangle
Output: angle opposite to the first edge in degrees
"""
function calc_angle_by_cos_law(a, b, c)
	angle::Real = acosd((b^2 + c^2 - a^2) / (2 * b * c))
	return angle
end

"""
	generate_trojan_triple_vector(num::Int)
Input: number
Output: Vector of named tuples with the trojan triples up to the given number
"""
function generate_trojan_triple_vector(num::Int)
	triple_set = Set{NamedTuple{(:a, :b, :c), Tuple{Int, Int, Int}}}()
	for big_num in 3:num
		for small_num::Int in 1:(big_num-1)÷2
			trip = generate_trojan_triple_120(big_num, small_num)
			push!(triple_set, (a = trip.a, b = trip.b, c = trip.c))
		end
	end

	triple_vector = collect(triple_set)

	sort!(triple_vector, by = x -> (x.c, x.b))
	return triple_vector
end

"""
	get_ext_trojan_triples(num::Int)

Input: number
Output: Vector of named tuples with the trojan triples up to the given number extended with the numbers they where generated from
"""
function get_ext_trojan_triple_vector(num::Int)
	triples::Vector{NamedTuple{(:a, :b, :c), Tuple{Int, Int, Int}}} = []
	ext_triples::Vector{NamedTuple{(:p, :q, :a, :b, :c), Tuple{Int, Int, Int, Int, Int}}} = []
	for big_num::Int in 3:num
		for small_num::Int in 1:floor((big_num - 1) / 2)
			triple = generate_trojan_triple_120(big_num, small_num)
			if triple ∉ triples
				push!(triples, (a = triple.a, b = triple.b, c = triple.c))
				push!(ext_triples, (p = big_num, q = small_num, a = triple.a, b = triple.b, c = triple.c))
			end
		end
	end
	sort!(ext_triples, by = x -> (x.c, x.b))
	return ext_triples
end


"""
	add_angles(triples::Vector{<:NamedTuple})

Input: vector of named tuples containing a,b,c
Output: vector of named tuple with the angles α, β, γ
"""
function add_angles_to_triple_vector(triples::Vector{<:NamedTuple})
	ext_triples = []
	for item in triples

		α = calc_angle_by_cos_law(item.a, item.b, item.c)
		β = calc_angle_by_cos_law(item.b, item.a, item.c)
		γ = calc_angle_by_cos_law(item.c, item.a, item.b)

		new_item = (; item..., α = α, β = β, γ = γ)


		push!(ext_triples, new_item)
	end
	return ext_triples
end

"""
	get_trojan_triples_for_a_number(num::Int)::Vector{NamedTuple{(:a, :b, :c), Tuple{Int, Int, Int}}}

Input: number
Output: vector of all possible named trojan triples, which have an edge of the size of the input number
"""
function get_trojan_triples_for_a_number_vector(num::Int)
	big_num::Int = ceil(sqrt(4 * num / 3 + 1))
	triple_vector = generate_trojan_triple_vector(big_num)
	every_triple_vector = expand_trojan_triple_vector(triple_vector)
	right_triple_set = Set()
	for item in every_triple_vector
		for edge in item
			if num % edge == 0
				a = item.a * num / edge
				b = item.b * num / edge
				c = item.c * num / edge
				triple = (a = a, b = b, c = c)
				push!(right_triple_set, triple)
			end
		end

	end
	right_triple_vector = collect(right_triple_set)
	sort!(right_triple_vector, by = x -> (x.c, x.b))
	return right_triple_vector
end

"""
	get_every_trojan_triple(triples::Vector{NamedTuple{(:a, :b, :c), Tuple{Int, Int, Int}}})

Input: vector of named trojan triples with one angle of 120°
Output: vector of named trojan triples with all possible angles and the same hypothenuse
"""
function expand_trojan_triple_vector(triples)
	every_triple_set = Set()
	for item in triples
		first_triple = (a = item.a, b = item.c, c = item.a + item.b)
		second_triple = (a = item.b, b = item.c, c = item.a + item.b)
		push!(every_triple_set, first_triple)
		push!(every_triple_set, second_triple)
		push!(every_triple_set, item)
	end
	every_triple_vector = collect(every_triple_set)
	sort!(every_triple_vector, by = x -> (x.c, x.b))
	return every_triple_vector
end


"""
	analyze_c_frequencies(ext_triples, max_num)

analyzes the frequency of multiple c values of the triples up to a given number
input: vector of triples, number
Output: dictionary with the key of the frequency and the value of frequency the same number of c values
"""
function analyze_c_frequencies(ext_triples, max_num)
	c_values = [item.c for item in ext_triples if item.c < max_num^2 * 3 ÷ 4 + 1]
	c_frequencies = countmap(c_values)
	frequency_counts = Dict()
	for freq in values(c_frequencies)
		frequency_counts[freq] = get(frequency_counts, freq, 0) + 1
	end

	return frequency_counts
end

#####################################################
# Functions for the visualization of the triangles
#####################################################

function get_coordinates_point_C_of_laying_triangle(opposite_edge, adjacent_edge, base_edge)::NamedTuple{(:x, :y), Tuple{Real, Real}}
	cos_α = (base_edge^2 + adjacent_edge^2 - opposite_edge^2) / (2 * base_edge * adjacent_edge)
	x = adjacent_edge * cos_α
	y = adjacent_edge * sqrt(1 - cos_α^2)
	return (x = x, y = y)
end

function rotate_point_by_angle(point::NamedTuple{(:x, :y), <:Tuple{Real, Real}}, angle::Real)::NamedTuple{(:x, :y), Tuple{Real, Real}}
	x = point.x * cosd(angle) - point.y * sind(angle)
	y = point.x * sind(angle) + point.y * cosd(angle)
	return (x = x, y = y)
end

function move_point(point::NamedTuple{(:x, :y), <:Tuple{Real, Real}}, transvector::NamedTuple{(:x, :y), <:Tuple{Real, Real}})::NamedTuple{(:x, :y), Tuple{Real, Real}}
	x_new = point.x + transvector.x
	y_new = point.y + transvector.y
	return (x = x_new, y = y_new)
end

function stretch_point(point::NamedTuple{(:x, :y), <:Tuple{Real, Real}}, x_factor::Real, y_factor::Real)::NamedTuple{(:x, :y), Tuple{Real, Real}}
	x_new = point.x * x_factor
	y_new = point.y * y_factor
	return (x = x_new, y = y_new)

end

function get_points_normalized(triple::NamedTuple{(:a, :b, :c), <:Tuple{Int, Int, Int}})
	point_1 = (x = 0.0, y = 0.0)
	point_2 = (x = 1.0, y = 0.0)
	point_C_of_small_triangle = get_coordinates_point_C_of_laying_triangle(triple.a, triple.a + triple.b, triple.c)
	point_4 = point_C_of_small_triangle |> p -> stretch_point(p, 1000 / triple.c, 1000 / triple.c) |> p -> rotate_point_by_angle(p, 60.0)
	point_3 = (x = point_4.x + 1000, y = point_4.y)
	point_5 = (x = 0.5, y = sqrt(3) / 2)
	return (point1 = point_1, point2 = point_2, point3 = point_3, point4 = point_4, point5 = point_5)
end


function draw_custom_lines(points::NamedTuple{(:point1, :point2, :point3, :point4, :point5), <:Tuple{NamedTuple{(:x, :y), <:Tuple{Real, Real}}, NamedTuple{(:x, :y), <:Tuple{Real, Real}}, NamedTuple{(:x, :y), <:Tuple{Real, Real}}, NamedTuple{(:x, :y), <:Tuple{Real, Real}}, NamedTuple{(:x, :y), <:Tuple{Real, Real}}}}, filename::String, line_width::Real, line_color::String)
    height = 1450.0
    open(filename, "w") do file
        @printf(file, "<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 1500 1500'>\n")
        @printf(file, "  <line x1='%f' y1='%f' x2='%f' y2='%f' stroke='%s' stroke-width='%f'/>\n", points.point1.x, height - points.point1.y, points.point2.x, height - points.point2.y, line_color, line_width)
        @printf(file, "  <line x1='%f' y1='%f' x2='%f' y2='%f' stroke='%s' stroke-width='%f'/>\n", points.point2.x, height - points.point2.y, points.point3.x, height - points.point3.y, line_color, line_width)
        @printf(file, "  <line x1='%f' y1='%f' x2='%f' y2='%f' stroke='%s' stroke-width='%f'/>\n", points.point3.x, height - points.point3.y, points.point4.x, height - points.point4.y, line_color, line_width)
        @printf(file, "  <line x1='%f' y1='%f' x2='%f' y2='%f' stroke='%s' stroke-width='%f'/>\n", points.point4.x, height - points.point4.y, points.point1.x, height - points.point1.y, line_color, line_width)
        @printf(file, "  <line x1='%f' y1='%f' x2='%f' y2='%f' stroke='%s' stroke-width='%f'/>\n", points.point1.x, height - points.point1.y, points.point5.x, height - points.point5.y, line_color, line_width)
        @printf(file, "  <line x1='%f' y1='%f' x2='%f' y2='%f' stroke='%s' stroke-width='%f'/>\n", points.point2.x, height - points.point2.y, points.point5.x, height - points.point5.y, line_color, line_width)
        @printf(file, "  <line x1='%f' y1='%f' x2='%f' y2='%f' stroke='%s' stroke-width='%f'/>\n", points.point3.x, height - points.point3.y, points.point5.x, height - points.point5.y, line_color, line_width)
        @printf(file, "  <line x1='%f' y1='%f' x2='%f' y2='%f' stroke='%s' stroke-width='%f'/>\n", points.point4.x, height - points.point4.y, points.point5.x, height - points.point5.y, line_color, line_width)
        @printf(file, "</svg>\n")
    end
end
end # module Trojan

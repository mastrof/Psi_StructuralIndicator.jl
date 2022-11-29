using DelimitedFiles
using LightGraphs

"""
    linestovec(fname; type=Int sep='\t')
Read each line of a file into a distinct vector.

Input:
- `fname`: name of the file to read 
- `type`: type of the data to parse (defaults to `Int` to read links)
- `sep`: separator used in the file (defaults to tab)
Output:
- a vector containing a vector of data for each line in the file
"""
linestovec(fname; type=Int, sep='\t') = [parse.(type, split(line, sep)) for line in readlines(fname)]

"""
    makegraph(links::Vector{Vector{Int}})
Turn a list of links into a network representation.

Input:
- `links`: a vector whose each element `i` is a vector of the molecules `j` linked to `i`
Output:
- a `SimpleGraph` representation of `links`
"""
function makegraph(links::Vector{Vector{Int}})
    g = SimpleGraph(length(links))
    for i in eachindex(links)
        for j in links[i]
            add_edge!(g, (i,j))
        end
    end
    return g
end

"""
    shell(g, i, D)
Find all the neighbors at chemical distance `D` from `i` in the network `g`
    
Input:
- `g`: a `SimpleGraph` representing connected molecules in the system
- `i`: the id of the molecule for which neighbors will be calculated
- `D`: chemical distance
Output:
- a vector of molecule ids at chemical distance `D` from `i`
"""
shell(g,i,D) = setdiff(neighborhood(g,i,D), neighborhood(g,i,D-1))

"""
    calc_Ψ(linksD4::Vector{Vector{Int}}, distances::Matrix{Float64})
Evaluate the structural indicator Ψ.

Input:
- `linksD4`: a vector where each element `i` is a list of molecules at chemical distance 4 to `i`
- `distances`: matrix with distances between each pair of molecules in the system
Output:
- a vector where each element `i` is the value of Ψ for molecule `i`

---

    calc_Ψ(fname_links::String, fname_distances::String)
Evaluate the structural indicator Ψ.

Input:
- `fname_links`: path to file containing all links in the system
- `fname_distances`: path to file containing all pair distances in the system
Output:
- a vector where each element `i` is the value of Ψ for molecule `i`
"""
function calc_Ψ(linksD4::Vector{Vector{Int}}, distances::Matrix{Float64})
    Ψ = zeros(length(linksD4))
    for i in eachindex(Ψ)
        js = linksD4[i]
        Ψ[i] = minimum(distances[js,i])
    end
    return Ψ
end

function calc_Ψ(fname_links::String, fname_distances::String)
    links = linestovec(fname_links)
    distances = readdlm(fname_distances)
    g = makegraph(links)
    linksD4 = map(i -> shell(g,i,4), eachindex(links))
    Ψ = calc_Ψ(linksD4, distances)
    return Ψ
end

Ψ = calc_Ψ("links.dat", "distances.dat")


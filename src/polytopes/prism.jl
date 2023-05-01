# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Prism(base, side)

A [prism](https://en.wikipedia.org/wiki/Prism_(geometry)) consisting of a
polygonal `base` and a copy of that base translated by a `side` vector.

"""
struct Prism{N,Dim,T,V<:AbstractVector{Point{Dim,T}}} <: Polyhedron{Dim,T}
  base::Ngon{N,Dim,T,V}
  side::Vec{Dim,T}

  function Prism(base::Ngon{N,Dim,T,V}, side::Vec{Dim,T}; fix=(Dim==3)) where {N,Dim,T,V}
    # orient base such that the normal and the side point in the same direction
    if fix && dot(normal(base), side) < zero(T)
      base = Ngon(reverse(vertices(base)))
    end
    new{N,Dim,T,V}(base, side)
  end
end

Prism(base::AbstractVector, args...; kwargs...) =
  Prism(Ngon(base), args...; kwargs...)

Prism(base::Ngon{N,Dim,T}, length::T) where {N,Dim,T} =
  Prism(base, length * normal(base); fix=false)

isconvex(p::Prism) = isconvex(p.base)

measure(p::Prism) = measure(p.base) * abs(dot(normal(p.base), p.side))

nvertices(::Type{<:Prism{N}}) where N = 2 * N
nvertices(p::Prism) = nvertices(typeof(p))

function vertices(p::Prism{N}) where N
  v = vertices(p.base)
  Vec(ntuple(i -> i <= N ? v[i] : v[i-N] + p.side, 2*N))
end

nedges(p::Type{<:Prism{N}}) where N = 3 * N
nedges(p::Prism) = nedges(typeof(p))

function edges(p::Prism{N}) where N
  v = vertices(p)
  (materialize(edge_connectivity(p, ind), v) for ind in 1:nedges(p))
end

nfacets(p::Type{<:Prism{N}}) where N = 2 + N
nfacets(p::Prism) = nfacets(typeof(p))

function facets(p::Prism{N}) where N
  v = vertices(p)
  (materialize(facet_connectivity(p, ind), v) for ind in 1:nfacets(p))
end

# NOTE: creating a boundary allocates a vector for the connections
boundary(p::Prism) = SimpleMesh(vertices(p), collect(facet_connectivity(p, ind) for ind in 1:nfacets(p)))

function edge_connectivity(p::Prism{N}, ind::Int) where N
  segment = if ind < N # lower ring without last segment
    (ind, ind+1)
  elseif ind == N # last segment of lower ring
    (N, 1)
  elseif N < ind <= 2*N # side segments
    (ind - N, ind)
  elseif 2*N < ind < 3*N # upper ring without last segment
    (ind - N, ind - N + 1)
  elseif ind == 3*N # last segment of upper ring
    (2*N, N+1)
  else
    error("Invalid edge index `$ind` (should be between 1 and $(3*N))")
  end
  connect(segment)
end

# NOTE: this function is not type stable unless the prism base is a quadrangle,
# since the prism sides are always quadrangles
function facet_connectivity(p::Prism{N}, ind::Int) where N
  # oriented such that vertices are counter-clockwise when viewed from the
  # outside, assuming that the side vector and the normal vector of the base
  # point in the same direction
  facet = if ind == 1 # lower facet
    Tuple(N:-1:1)
  elseif 1 < ind <= N # side facets without last side
    (ind-1, ind, ind+N, ind+N-1)
  elseif ind == N+1 # last side facet
    (N, 1, 1+N, 2*N)
  elseif ind == N+2 # upper facet
    Tuple(N+1:2*N)
  else
    error("Invalid facet index `$ind` (should be between 1 and $(N+2))")
  end
  connect(facet)
end

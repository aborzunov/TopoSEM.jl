"""
    BivariatePoly{Order,T}

Higher-order (total degree 2..Order) part of the bivariate polynomial φ(u,v)
from Eq. (6) of Neggers et al. 2022. The constant (`p₀`) and linear (`p₂`)
parts are stored separately on the `DetectorModel`. The argument is rescaled
by `scale` before evaluation so that the design matrix in Step C stays well
conditioned for `Order ≥ 3`.

Monomial ordering inside `coeffs` is by total degree, then by power of `v`:
    d=2 :  u²,  uv,  v²
    d=3 :  u³,  u²v, uv², v³
    d=4 :  u⁴,  u³v, u²v², uv³, v⁴
"""
struct BivariatePoly{Order,T}
    coeffs::Vector{T}
    scale::T
end

n_higher_coeffs(Order::Integer) = Order < 2 ? 0 : sum(d + 1 for d in 2:Order)

function BivariatePoly{Order,T}(; scale::T = one(T)) where {Order,T}
    return BivariatePoly{Order,T}(zeros(T, n_higher_coeffs(Order)), scale)
end

function BivariatePoly(::Type{T}, Order::Int; scale::Real = one(T)) where {T}
    return BivariatePoly{Order,T}(zeros(T, n_higher_coeffs(Order)), T(scale))
end

@inline poly_order(::BivariatePoly{Order}) where {Order} = Order
@inline n_coeffs(p::BivariatePoly) = length(p.coeffs)

"""
    evaluate_higher(p, u, v)

Return Σᵢ cᵢ · (u/scale)^pᵢ · (v/scale)^qᵢ. Linear and constant parts are not
included — they live on the `DetectorModel`.
"""
function evaluate_higher(p::BivariatePoly{Order,T}, u::Real, v::Real) where {Order,T}
    Order < 2 && return zero(T)
    us = T(u) / p.scale
    vs = T(v) / p.scale
    s = zero(T)
    idx = 1
    @inbounds for d in 2:Order
        for j in 0:d
            i = d - j
            s += p.coeffs[idx] * us^i * vs^j
            idx += 1
        end
    end
    return s
end

"""
    gradient_higher(p, u, v) -> (∂φₕ/∂u, ∂φₕ/∂v)

Partial derivatives of the higher-order contribution. Used to build the
tangent model in Step A.
"""
function gradient_higher(p::BivariatePoly{Order,T}, u::Real, v::Real) where {Order,T}
    Order < 2 && return (zero(T), zero(T))
    us = T(u) / p.scale
    vs = T(v) / p.scale
    inv_s = one(T) / p.scale
    du = zero(T); dv = zero(T)
    idx = 1
    @inbounds for d in 2:Order
        for j in 0:d
            i = d - j
            c = p.coeffs[idx]
            if i ≥ 1
                du += c * i * us^(i-1) * vs^j
            end
            if j ≥ 1
                dv += c * j * us^i * vs^(j-1)
            end
            idx += 1
        end
    end
    return (du * inv_s, dv * inv_s)
end

"""
    monomial_basis_higher!(out, p, u, v)

Fill `out[k] = (u/scale)^pₖ · (v/scale)^qₖ` for the higher-order monomials in
the same order as `coeffs`. Used to assemble the design matrix in Step C.
"""
function monomial_basis_higher!(out::AbstractVector, p::BivariatePoly{Order,T},
                                 u::Real, v::Real) where {Order,T}
    Order < 2 && return out
    us = T(u) / p.scale
    vs = T(v) / p.scale
    idx = 1
    @inbounds for d in 2:Order
        for j in 0:d
            i = d - j
            out[idx] = us^i * vs^j
            idx += 1
        end
    end
    return out
end

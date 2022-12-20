# constants
const L0::Float64 = 0.052917721 #nm
const R::Float64 = 13605.8 * 2 #meV


# Parameters
const N::Int64 = 9
const M::Float64 = 0.067 # m_e
const L::Float64 = 10 / L0 # atomic
const V0::Float64 = 10 / R
const OMEGA::Float64 = 1.0 / R
const A::Float64 = 2 * M * OMEGA / 2
const Δx::Float64 = 3 / sqrt(A)
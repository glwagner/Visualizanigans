using Printf
using Oceananigans
using Oceananigans.BoundaryConditions
using Oceananigans.Fields
using Oceananigans.OutputWriters
using Oceananigans.Diagnostics
using Oceananigans.Utils
using Oceananigans.AbstractOperations
using Vizinanigans

## Create Oceananigans Model
arch = CPU()
FT   = Float64

end_time = 100*365day
const Lx = 1000kilometer # 
const Ly = 1000kilometer #
const Lz = 3kilometer    #

# Rough resolution
Δx = Δy = scale * 250meter # 5km
Δz = 100meter

const Nx = round(Int, Lx/ Δx / 16) * 16
const Ny = round(Int, Ly/ Δy / 16) * 16
const Nz = round(Int, Lz/ Δz / 16) * 16

topology = (Periodic, Bounded, Bounded)
grid = RegularCartesianGrid(topology=topology, size=(Nx, Ny, Nz), x=(0, Lx), y=(0, Ly), z=(-Lz, 0))

const f = -1e-4
const β = 1 * 10^(-11)
coriolis = FPlane(FT, f=f)
coriolis = BetaPlane(FT, f₀ = f, β = β)

α = 2e-4  # Thermal expansion coefficient [K⁻¹]
eos = LinearEquationOfState(FT, α=α, β=0)
buoyancy = BuoyancyTracer()

κh = 0.5e-5
νh = 12.0 
κv = 0.5e-5 
νv = 3e-4 
closure = AnisotropicDiffusivity(νx = νh, νy = νh, νz =νv, κx = κh, κy = κh, κz=κv)
ΔB = 10 * 2e-3
h = 1000.0 # [m] relexaction profile scale
# Initial Conditions
U₀(x,y,z) = (z + Lz)/(f + β * y) * (ΔB / Ly)
B₀(x, y, z) = ΔB * (y / Ly +  ( exp(z/h) - exp(-Lz/h) ) / (1 - exp(-Lz/h)) - 1)

## checkpointing
model = IncompressibleModel(
            architecture = arch,
                float_type = FT,
                    grid = grid,
                coriolis = coriolis,
                buoyancy = buoyancy,
                closure = closure,
                tracers = (:b,),
)
set!(model, u = U₀, b=B₀)

## check Makie
visualize(model)
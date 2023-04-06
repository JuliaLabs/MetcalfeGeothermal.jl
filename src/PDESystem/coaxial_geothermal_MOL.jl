using ModelingToolkit, MethodOfLines, OrdinaryDiffEq, DomainSets, Plots

println("1D PDE Depth Model of Coaxial Geothermal - Bob.Metcalfe@UTexas.edu TexasGEO.org")
# ---------- Earth
const earth = (EarthSurfaceTemperature=15.0, EarthTemperatureGradient=0.025) # [C], [C/m]
# Ambient temperature in Celsius at depth in meters
function AmbientTemperature(earth, depth)
    earth.EarthSurfaceTemperature + depth * earth.EarthTemperatureGradient
end
# Rock is Granite
const RockDensity = GraniteDensity = 2750000 # [g/m3]
const RockSpecificHeat = GraniteSpecificHeat = 0.790 # [J/gC]
const RockThermalConductivity = GraniteThermalConductivity = 2.62 # [W/mC] range average
# Fluid is Water
const FluidDensity = WaterDensity = 997000.0 # [g/m3]
const FluidSpecificHeat = WaterSpecificHeat = 4.186 # [J/gC]
const FluidThermalConductivity = WaterThermalConductivity = 0.6 # [W/mC]
# ---------- Drilling
# Volume of inner and outer pipes is equal for now
const NumberOfPipes = 15 # [n]
const LengthOfPipe = 10.0 # [m]
const DepthOfWell = NumberOfPipes * LengthOfPipe # [m]
const InnerRadius = 0.1 # [m]
const OuterRadius = InnerRadius * sqrt(2) # [m]
println("Well depth ", DepthOfWell, " meters, pipe radius ", InnerRadius, " meters.")
const InnerArea = 2 * InnerRadius * pi * LengthOfPipe # [m2]
const OuterArea = 2 * OuterRadius * pi * LengthOfPipe # [m2]
const InnerVolume = pi * InnerRadius^2 * LengthOfPipe # [m3]
const OuterVolume = (pi * OuterRadius^2 * LengthOfPipe) - InnerVolume # [m3]
# Check to be sure inner and outer pipe volumes are equal (for now)
if round(InnerVolume, digits=4) != round(OuterVolume, digits=4)
    println("Inner and outer volumnes are not equal.")
end

const PumpSpeed = InnerVolume /5 # [m3/s]
const PumpTime = InnerVolume / PumpSpeed # [s]
const WaterSpeed = LengthOfPipe/PumpTime # [m/s]

# Let the drilling begin! - With MethodOfLines.jl
############################################################################################
# 1D domain: z
############################################################################################

println("Water speed ", WaterSpeed, " meters per second.")


@parameters t z
@variables Tinner(..), Touter(..)

Dz = Differential(z)
Dt = Differential(t)

# Let the drilling begin

@parameters t z
@variables Tinner(..), Touter(..)

Dz = Differential(z)
Dt = Differential(t)

# Define constants
V = WaterSpeed #[m / s]
c_p = FluidSpecificHeat
ρ = FluidDensity

u_inner = 2π * InnerRadius
u_outer = 2π * OuterRadius

# ! These are far too small, we need to find another way to get the correct values
α_inner = FluidThermalConductivity
α_outer = RockThermalConductivity

A_inner = π * InnerRadius^2 # cross sectional area of inner pipe
A_outer = π * (OuterRadius^2 - InnerRadius^2) # cross sectional area of outer pipe

tmax = 3000.0

# Structure defined in the PDF
eqs = [A_inner * ρ * c_p * Dt(Tinner(t, z)) - V * Dz(Tinner(t, z)) ~ -u_inner * α_inner * (Tinner(t, z) - Touter(t, z)),
    A_outer * ρ * c_p * Dt(Touter(t, z)) + V * Dz(Touter(t, z)) ~ u_inner * α_inner * (Tinner(t, z) - Touter(t, z)) +
                                                              u_outer * α_outer * (AmbientTemperature(earth, z) - Touter(t, z)),
]

bcs = [Tinner(0, z) ~ AmbientTemperature(earth, z),
    Touter(0, z) ~ AmbientTemperature(earth, z),
    Touter(t, 0) ~ AmbientTemperature(earth, 0),
    Tinner(t, DepthOfWell) ~ Touter(t, DepthOfWell),
]

domains = [z in Interval(0.0, DepthOfWell),
    t in Interval(0.0, tmax),
]

@named pde = PDESystem(eqs, bcs, domains, [t, z], [Tinner(t, z), Touter(t, z)])

# Create the discretization

disc = MOLFiniteDifference([z => 100], t)

# Generate the ODE problem
prob = discretize(pde, disc)

# Solve the ODE problem
sol = solve(prob, QBDF(), saveat=10.0)

sol_Tinner = sol[Tinner(t, z)]
sol_Touter = sol[Touter(t, z)]

solt = sol[t]
solz = sol[z]

println("At t = $tmax, the inner pipe temperature at ground level is ", sol_Tinner[end, 1], " and the outer pipe temperature is ", sol_Touter[end, 1], ".")

# Plot the solution
p = plot(solz, [sol_Tinner[end, 1:end], sol_Touter[end, 1:end]], label=["Tinner" "Touter"], xlabel="z", ylabel="T", title="Temperature vs. Depth")
savefig(p, "src/PDESystem/metcalfe_geothermal.png")

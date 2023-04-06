using ModelingToolkit, MethodOfLines, OrdinaryDiffEq, DomainSets, Plots

println("1D PDE Model of Coaxial Geothermal - Bob.Metcalfe@UTexas.edu TexasGEO.org")
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

const PumpSpeed = InnerVolume / 5.0 # [m3/s]
const PumpTime = InnerVolume / PumpSpeed # [s]

############################################################################################
# cylindrical domain **WIP**
############################################################################################

## Now try to do the same thing with radial gradients and a cylindrical domain
@parameters z r t
@variables T(..)

Dt = Differential(t)
Dz = Differential(z)
Dr = Differential(r)
Drr = Differential(r)^2

# r dependent velocity
v(x) = ifelse(x < InnerRadius, PumpSpeed, -PumpSpeed)

# ? TODO: Define constants

# Single equation this time
eqs = [K * Dt(T(t, z, r)) + v(r) * Dz(T(t, z, r)) ~ u_inner * α_inner * 1 / r^2 * Dr(r^2 * Dr(T(t, z, r)))]

bcs = [T(0, z, r) ~ AmbientTemperature(earth, z),
    T(0, z, r) ~ AmbientTemperature(earth, z),
    T(t, z, 0) ~ 6 * Drr(T(u, z, r)), # from the paper https://web.mit.edu/braatzgroup/analysis_of_finite_difference_discretization_schemes_for_diffusion_in_spheres_with_variable_diffusivity.pdf
    T(t, z, OuterRadius) ~ AmbientTemperature(earth, z),
    T(t, 0, r) ~ AmbientTemperature(earth, 0),
    Touter(t, DepthOfWell, r) ~ Tinner(t, DepthOfWell, r),
]

domains = [z in Interval(0.0, DepthOfWell),
    r in Interval(0.0, OuterRadius),
    t in Interval(0.0, 100.0),
]

@named pde = PDESystem(eqs, bcs, domains, [t, z, r], [T(t, z, r)])

# Create the discretization

disc = MOLFiniteDifference([z => 32, r => 16], t)

# Generate the ODE problem
prob = discretize(pde, disc)

# Solve the ODE problem
sol = solve(prob, QBDF(), saveat=10.0)

sol_Tinner = sol[Tinner(t, z, r)]
sol_Touter = sol[Touter(t, z, r)]

solt = sol[t]
solz = sol[z]

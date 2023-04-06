module MetcalfeGeothermal
using Plots
println("Toy Model of Coaxial Geothermal - Bob.Metcalfe@UTexas.edu TexasGEO.org")
# ---------- Earth
const earth = (EarthSurfaceTemperature = 15.0, EarthTemperatureGradient = 0.025) # [C], [C/m]
# Ambient temperature in Celsius at depth in meters
function  AmbientTemperature(earth, depth)
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
println("Well depth ",DepthOfWell," meters, pipe radius ",InnerRadius," meters.")
const InnerArea = 2 * InnerRadius * pi * LengthOfPipe # [m2]
const OuterArea = 2 * OuterRadius * pi * LengthOfPipe # [m2]
const InnerVolume = pi * InnerRadius^2 * LengthOfPipe # [m3]
const OuterVolume = (pi * OuterRadius^2 * LengthOfPipe) - InnerVolume # [m3]
# Check to be sure inner and outer pipe volumes are equal (for now)
if round(InnerVolume, digits=4) != round(OuterVolume, digits=4)
    println("Inner and outer volumnes are not equal.")
end
# Let the drilling begin
mutable struct Pipe{T,V}
    innertemp::T
    outertemp::V
end
Base.show(io::IO, p::Pipe) = print(io, "Pipe: innertemp = $(p.innertemp), outertemp = $(p.outertemp)")
PipeString = [(ambienttemp = AmbientTemperature(earth, (ip-1)*LengthOfPipe); Pipe(ambienttemp, ambienttemp)) for ip = 1:NumberOfPipes]
# ---------- Completion of well
const PumpSpeed = InnerVolume/5.0 # [m3/s]
const PumpTime = InnerVolume/PumpSpeed # [s]
# ---------- Operation
for run in 1:100 # iteration on thermo updates

    for ip in 1:NumberOfPipes
        pipe = PipeString[ip]
        # outer pipe to inner pipe Joule flux
        O2Ic = FluidThermalConductivity
        O2Idt = pipe.outertemp - pipe.innertemp
        O2Ia = InnerArea
        O2Im = InnerVolume * FluidDensity
        O2Ipt = PumpTime

        O2Ijf = O2Ic * O2Idt * O2Ia * O2Im * O2Ipt
        # print("O2I ",run, O2Ijf,O2Ic,O2Idt,O2Ia,O2Im,O2Ipt)
        # ambient to outer Joule flux
        A2Oc = RockThermalConductivity
        A2Odt = AmbientTemperature(earth, (ip-1)*LengthOfPipe) - pipe.outertemp
        A2Oa = OuterArea
        A2Om = OuterVolume * FluidDensity
        A2Opt = PumpTime

        A2Ojf = A2Oc * A2Odt * A2Oa * A2Om * A2Opt
        # print("A2O ",run,A2Ojf,A2Oc,A2Odt,A2Oa,A2Om,A2Opt)
        # Update temperatures
        pipe.innertemp += O2Ijf / (FluidSpecificHeat * O2Im)
        pipe.outertemp += -O2Ijf / (FluidSpecificHeat * A2Om)
        pipe.outertemp += A2Ojf / (RockSpecificHeat * A2Om)
    end
end
# pump - move fluid down and up coaxial pipes through rock
function printPipeString(message)
    println(message)
    for ip in 1:NumberOfPipes
        println(ip," ", PipeString[ip])
    end
end
# printPipeString("---------- start ----------")
savedinnerpipe = PipeString[1].innertemp
for ip in 1:NumberOfPipes-1
    PipeString[ip].innertemp = PipeString[ip+1].innertemp # fluid moving up
end
savedouterpipe = PipeString[NumberOfPipes].outertemp
for ip in NumberOfPipes:-1:2
    PipeString[ip].outertemp = PipeString[ip-1].outertemp
end
PipeString[1].outertemp = PipeString[end].outertemp
PipeString[1].outertemp = savedinnerpipe
PipeString[NumberOfPipes].innertemp = savedouterpipe
printPipeString("----------- stop -----------")

plot([pipe.innertemp for pipe in PipeString], label="inner")
plot!([pipe.outertemp for pipe in PipeString], label="outer")
end # module MetcalfeGeothermal

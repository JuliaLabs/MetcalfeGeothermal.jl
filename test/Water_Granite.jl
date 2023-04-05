cd(@__DIR__)
using Pkg
Pkg.activate("..")

using MetcalfeGeothermal
@show MetcalfeGeothermal.PipeString

MetcalfeGeothermal.printPipeString("----------- stop -----------")

# 01_solvemodel.jl

using DifferentialEquations
using Plots, CSV

PATH_TO_MODELS = pwd()*"/models/"
PATH_TO_OUTPUT = pwd()*"/output/" # Path for output data, plots
FILENAME = "outputMM.csv" # Filename for output data
PLOTNAME = "plottraces.pdf"
MODELNAME = "TXTL_MM.jl"
include(PATH_TO_MODELS*MODELNAME)

###########################
# 1. Set up and solve ODEs
###########################

u0 = valuesVar
params = valuesPar
tspan = 0,Float64(12*60) # time in min
prob = ODEProblem(model!,u0,tspan,params)

# Solve ODEs
sol_ODE = solve(prob,saveat=10)

###########################
# 2. Output results
###########################

# Plot options

fntsm = Plots.font("sans-serif", pointsize=round(8.0))
fntlg = Plots.font("sans-serif", pointsize=round(12.0))
default(titlefont=fntlg, guidefont=fntlg, tickfont=fntlg, legendfont=fntsm)
p = plot(grid=:false,legend=:true,framestyle=:frame,size=(500,300))

# Plot data
plot!(sol_ODE, lw=2, labels=permutedims(keysVar),legend = :outertopright) 
plot!(xaxis = "t (min)",yaxis="concs (uM)")


# Save plot
savefig(PATH_TO_OUTPUT*PLOTNAME)

# Save data
CSV.write(PATH_TO_OUTPUT*FILENAME,sol_ODE)
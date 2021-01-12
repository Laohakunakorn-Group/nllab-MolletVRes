# 04_calibrateGlobal.jl

using DifferentialEquations, DataFrames
using Plots, CSV, XLSX, StatsPlots
using Optim, DiffEqParamEstim, Turing, Distributions

PATH_TO_MODELS = pwd()*"/models/"
PATH_TO_OUTPUT = pwd()*"/output/" # Path for output data, plots
FILENAME = "outputMM.csv" # Filename for output data
PLOTNAME = "plotGlobalFit.pdf"
MODELNAME = "TXTL_MM.jl"
include(PATH_TO_MODELS*MODELNAME)


# Get data (more info about experiment in data folder)

PATH_TO_DATA = "./data/"
xf = XLSX.readxlsx(PATH_TO_DATA*"data_capacity_in_nM.xlsx")
sheetname = XLSX.sheetnames(xf)[1]

sh = xf[sheetname]
data = sh[:]

t = data[2:end,1]*60 # convert time to minutes
y = data[2:end,2:10]
dy = data[2:end,12:20]

conds = [
60.3,
44.5,
22.25,
11.125,
5.5625,
2.78125,
1.390625,
0.6953125,
0.34765625]

# Experimental data for 60.3 nM
ydata=y[:,1]/1000; # in uM
ystdev=dy[:,1]/1000; # in uM

# Model output for 60.3 nM
dsel = 60.3/1000 # DNA concentration in uM

m0=0.0;
P0=2.0;
Pmat0=0.0;
A0=1.0;

Vmaxtx=1.0;
KTX=0.01;
kdeg=0.01;
Vmaxtl=0.1;
KTL=0.1;
kmat=0.03;
lamb=0.02;

idx_d=1;
idx_m=2;
idx_P=3;
idx_Pmat=4;
idx_A=5;

# Define simple L2 cost function
params = [Vmaxtx,KTX,kdeg,Vmaxtl,KTL,kmat,lamb]
u0=[dsel,m0,P0,Pmat0,A0]
tspan = 0,Float64(12*60) # time in min
prob = ODEProblem(model!,u0,tspan,params)
sol = solve(prob,saveat=60);
tsim = sol.t

cost_function = build_loss_objective(prob, Tsit5(),
                L2Loss(tsim,ydata),
    maxiters=10000000,verbose=true,save_idxs=4) # save only observable protein output

# Check to see if cost function evaluates
cost_function(params)
print("cost function OK \n")

#################################################
# Global optimization example (Simulated Annealing with bounds)
#################################################

lowerbounds = params*0.0
upperbounds = params*10.0
result = optimize(cost_function,lowerbounds,upperbounds,params, SAMIN(), 
    Optim.Options(iterations=100000,time_limit=100)) 
paramsfit=Optim.minimizer(result)

mod_prob = remake(prob, p=paramsfit)
mod_sol = solve(mod_prob, Tsit5(), saveat=60, verbose=true,maxiters=1000000,save_idxs=4)
ymodelfit = Array(mod_sol)

# Plot comparison
fntsm = Plots.font("sans-serif", pointsize=round(8.0))
fntlg = Plots.font("sans-serif", pointsize=round(12.0))
default(titlefont=fntlg, guidefont=fntlg, tickfont=fntlg, legendfont=fntsm)
pl = plot(grid=:false,legend=:true,framestyle=:frame,size=(500,300))

scatter!(t,ydata, lw=2, 
    labels="60.3 nM",legend = :outertopright)
plot!(tsim,ymodelfit, lw=2, 
        labels="sim 60.3 nM",legend = :outertopright) 
plot!(xaxis = "t (min)",yaxis="concs (uM)", title="Simulated Annealing fit")

print("Best fit params: ", paramsfit)

# Save plot
savefig(PATH_TO_OUTPUT*PLOTNAME)
# 05_calibrateMCMC.jl

using DifferentialEquations, DataFrames
using Plots, CSV, XLSX, StatsPlots
using Optim, DiffEqParamEstim, Turing, Distributions

PATH_TO_MODELS = pwd()*"/models/"
PATH_TO_OUTPUT = pwd()*"/output/" # Path for output data, plots
FILENAME = "outputMM.csv" # Filename for output data
PLOTNAME = "plotMCMCFit.pdf"
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

u0=[dsel,m0,P0,Pmat0,A0]
tspan = 0,Float64(12*60) # time in min

########################################
# Bayesian example
########################################

Turing.setadbackend(:forwarddiff) # Details here: forwarddiff ok for small models, for 
                                  # larger models require adjoint 

# Define our probabilistic model
@model function fitlv(data)
    σ ~ InverseGamma(2, 3) 
    Vmaxtx ~ truncated(Normal(1.0,0.5),0,10) # mean, sigma, lower bound, upper bound
    KTX ~ truncated(Normal(0.01,0.005),0,0.1)
    kdeg ~ truncated(Normal(0.01,0.005),0,0.1)
    Vmaxtl ~ truncated(Normal(0.1,0.05),0,1)
    KTL ~ truncated(Normal(0.1,0.05),0,1)
    kmat ~ truncated(Normal(0.03,0.015),0,0.3)
    lamb ~ truncated(Normal(0.02,0.01),0,0.2)
    
    params = [Vmaxtx,KTX,kdeg,Vmaxtl,KTL,kmat,lamb] # Draw parameters from prior distribution
    prob = ODEProblem(model!,u0,tspan,params)
    predicted = solve(prob,Tsit5(),saveat=60,save_idxs=4) 

    for i = 1:length(predicted)
        data[i] ~ Normal(predicted[i], σ) # generate synthetic data with σ-normal distribution noise
    end
end

# Start sampling
modelbayes = fitlv(ydata) # Fit model to synthetic data generated earlier
chain = sample(modelbayes, NUTS(.65),1000) # No U turn sampler
# Finish sampling

# Let's convert the chain to a dataframe and use this to 
# obtain samples, to plot an ensemble model result
dfc=DataFrame(chain);
print(names(dfc))
samples = [dfc[!,:Vmaxtx] dfc[!,:KTX]  dfc[!,:kdeg] dfc[!,:Vmaxtl] dfc[!,:KTL] dfc[!,:kmat] dfc[!,:lamb]];

# Experimental data for 60.3 nM
ydata=y[:,1]/1000; # in uM
ystdev=dy[:,1]/1000; # in uM

# Model output for 60.3 nM
dsel = 60.3/1000 # DNA concentration in uM
u0=[dsel,m0,P0,Pmat0,A0]
params = [Vmaxtx,KTX,kdeg,Vmaxtl,KTL,kmat,lamb];
prob = ODEProblem(model!,u0,tspan,params)

# Plot comparison
fntsm = Plots.font("sans-serif", pointsize=round(8.0))
fntlg = Plots.font("sans-serif", pointsize=round(12.0))
default(titlefont=fntlg, guidefont=fntlg, tickfont=fntlg, legendfont=fntsm)
pl = plot(grid=:false,legend=:false,framestyle=:frame,size=(500,300))

# Run simulation 50 times with random parameter values drawn from posterior distribution
for j in 1:50
    modprob = remake(prob, p=samples[rand(1:1000),:], u0=[dsel,m0,P0,Pmat0,A0])
    sol = solve(modprob,saveat=60);
    tsim = sol.t 
    Pmat = [datum for subarr in sol.u for datum in subarr[idx_Pmat]];
    plot!(tsim,Pmat, lw=2, alpha=0.05, color="blue") 
end

scatter!(t,ydata, yerror=ystdev, lw=2)

plot!(xaxis = "t (min)",yaxis="concs (uM)")

# Save plot
savefig(PATH_TO_OUTPUT*PLOTNAME)

# Save MCMC results
plot(chain)
savefig(PATH_TO_OUTPUT*"MCMCchain.pdf")
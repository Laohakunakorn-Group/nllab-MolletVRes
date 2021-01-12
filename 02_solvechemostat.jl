# 02_solvechemostat.jl

using DifferentialEquations, DataFrames
using Plots, CSV, XLSX, StatsPlots
using Optim, DiffEqParamEstim, Turing, Distributions

PATH_TO_MODELS = pwd()*"/models/"
PATH_TO_OUTPUT = pwd()*"/output/" # Path for output data, plots
FILENAME = "outputMM.csv" # Filename for output data
PLOTNAME = "plottraceschemostat.pdf"
MODELNAME = "TXTL_MM.jl"
include(PATH_TO_MODELS*MODELNAME)

######################################################
# 1. Encode experimental protocol in callbacks
######################################################

function marktime!(integrator)   
end

function condition(u,t,integrator) 
    t%INTERVAL_DIL # condition true when t is multiple of dilution interval
end

function dilute!(integrator)
    if integrator.t<SWITCHTIMES[1]
        # Stage 1
        INDEX_REFRESH = [idx_d,idx_A]
        CONC_REFRESH = [d0,A0]
        
        # first dilute everything
        for j in 1:NSPECIES
            integrator.u[j] = integrator.u[j]*(1-DIL_FRAC)
        end
        
        # then refresh appropriate species: here DNA and system activity A refreshed
        for j in 1:size(INDEX_REFRESH)[1]
            integrator.u[INDEX_REFRESH[j]] = integrator.u[INDEX_REFRESH[j]] + DIL_FRAC*CONC_REFRESH[j]
        end
        
    elseif SWITCHTIMES[1]<=integrator.t
        # Stage 2 - right now same as Stage 1
        INDEX_REFRESH = [idx_d,idx_A]
        CONC_REFRESH = [d0,A0]
        
        # first dilute everything
        for j in 1:NSPECIES
            integrator.u[j] = integrator.u[j]*(1-DIL_FRAC)
        end
        
        # then refresh appropriate species: here DNA and system activity A refreshed
        for j in 1:size(INDEX_REFRESH)[1]
            integrator.u[INDEX_REFRESH[j]] = integrator.u[INDEX_REFRESH[j]] + DIL_FRAC*CONC_REFRESH[j]
        end
        
    end
end

function solvemodel(grads,u0,params,TMAX,INTERVAL_DIL,tsave,DIL_FRAC,NSPECIES,SWITCHTIMES)
    tspan = (0.0,TMAX); 
    prob = ODEProblem(grads,u0,tspan,params);

    # Callbacks

    contcb = ContinuousCallback(condition,dilute!;save_positions=(false,false))
    periodcb = PeriodicCallback(marktime!,INTERVAL_DIL;save_positions=(false,false)) # hack to mark time points for dilution
    cb = CallbackSet(periodcb,contcb)

    # Solve
    sol = solve(prob, callback=cb,saveat=tsave);
    return(sol) # Return trajectories
end

######################################################
# 2. Simulate
######################################################

# Global simulation settings
TMAX = 12.0*60 # in minutes
INTERVAL_DIL = 15.0 # in minutes
DIL_FRAC = 0.2;
NSPECIES = 5
SWITCHTIMES = [120] # Time to swap between experimental stages, in minutes

# 1. Define names for species indices
idx_d=1;
idx_m=2;
idx_P=3;
idx_Pmat=4;
idx_A=5;

# 2. Set initial conditions and parameters
d0=0.001;
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

u0 = [d0,m0,P0,Pmat0,A0];
params = [Vmaxtx,KTX,kdeg,Vmaxtl,KTL,kmat,lamb];

# 3. Solve

# Save data frequently e.g. every 2 minutes to observe sawtooth behaviour:
TSAVE = collect(0:2:TMAX) 

# Or save at each dilution step to see smooth trace (similar to experiment):
#TSAVE = collect(0:INTERVAL_DIL:TMAX) 

solU,solDU=solvemodel(model!,u0,params,TMAX,INTERVAL_DIL,TSAVE,DIL_FRAC,NSPECIES,SWITCHTIMES);
t = solU.t;
d = [datum for subarr in solU.u for datum in subarr[idx_d]];
m = [datum for subarr in solU.u for datum in subarr[idx_m]];
P = [datum for subarr in solU.u for datum in subarr[idx_P]];
Pmat = [datum for subarr in solU.u for datum in subarr[idx_Pmat]];
A = [datum for subarr in solU.u for datum in subarr[idx_A]];

# Plot options

fntsm = Plots.font("sans-serif", pointsize=round(8.0))
fntlg = Plots.font("sans-serif", pointsize=round(12.0))
default(titlefont=fntlg, guidefont=fntlg, tickfont=fntlg, legendfont=fntsm)
p = plot(grid=:false,legend=:true,framestyle=:frame,size=(500,300))

# Plot data
plot!(t,d, lw=2, labels="d",legend = :outertopright) 
plot!(t,m, lw=2, labels="m") 
plot!(t,P, lw=2, labels="P") 
plot!(t,Pmat, lw=2, labels="Pmat") 
plot!(t,A, lw=2, labels="A") 
plot!(xaxis = "t (min)",yaxis="concs (uM)")

# Save plot
savefig(PATH_TO_OUTPUT*PLOTNAME)
# TXTL_MM.jl
# Simple saturation kinetics with DNA 
# Julia model ODEs from antimony file

function model!(du, u, params, t)

	d,m,P,Pmat,alpha = u

	Vmaxtx,KTX,kdeg,Vmaxtl,KTL,kmat,lamb = params

	du[1] = 0.0
	du[2] = + 1.0 * Vmaxtx*alpha*d/(d+KTX)-kdeg*m
	du[3] = + 1.0 * Vmaxtl*alpha*m/(m+KTL) - 1.0 * kmat*P
	du[4] = + 1.0 * kmat*P
	du[5] = - 1.0 * lamb*alpha
end

keysVar = ["d","m","P","Pmat","alpha"]
valuesVar = [0.001,0.0,0.0,0.0,1.0]
dictVar = Dict(keysVar .=> valuesVar)

keysPar = ["Vmaxtx","KTX","kdeg","Vmaxtl","KTL","kmat","lamb"]
valuesPar = [1.0,0.01,0.01,0.1,0.1,0.04,0.02]
dictPar = Dict(keysPar .=> valuesPar)


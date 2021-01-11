# TXTL_MA.jl
# Mass action kinetics
# Julia model ODEs from antimony file

function model!(du, u, params, t)

	D,R,RD,m,r,rm,P,Pmat,ATX,ATL = u

	k1,k2,k3,ka,kb,kc,kmat,k4,kx,kl = params

	du[1] = - 1.0 * k1*D*R + 1.0 * k2*RD + 1.0 * k3*RD*ATX
	du[2] = - 1.0 * k1*D*R + 1.0 * k2*RD + 1.0 * k3*RD*ATX
	du[3] = + 1.0 * k1*D*R - 1.0 * k2*RD - 1.0 * k3*RD*ATX
	du[4] = + 1.0 * k3*RD*ATX - 1.0 * ka*m*r + 1.0 * kb*rm + 1.0 * kc*rm*ATL - 1.0 * k4*m
	du[5] = - 1.0 * ka*m*r + 1.0 * kb*rm + 1.0 * kc*rm*ATL
	du[6] = + 1.0 * ka*m*r - 1.0 * kb*rm - 1.0 * kc*rm*ATL
	du[7] = + 1.0 * kc*rm*ATL - 1.0 * kmat*P
	du[8] = + 1.0 * kmat*P
	du[9] = - 1.0 * kx*ATX
	du[10] = - 1.0 * kl*ATL
end

keysVar = ["D","R","RD","m","r","rm","P","Pmat","ATX","ATL"]
valuesVar = [1.0,870.0,0.0,0.0,108.0,0.0,0.0,0.0,1.0,1.0]
dictVar = Dict(keysVar .=> valuesVar)

keysPar = ["k1","k2","k3","ka","kb","kc","kmat","k4","kx","kl"]
valuesPar = [0.001,0.01,0.025,0.001,0.05,0.0112,0.000396,0.00121,3.8e-10,0.000142]
dictPar = Dict(keysPar .=> valuesPar)


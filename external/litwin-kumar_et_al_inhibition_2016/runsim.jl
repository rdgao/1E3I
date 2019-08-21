#this file is part of litwin-kumar_et_al_inhibition_2016
#Copyright (C) 2016 Ashok Litwin-Kumar
#see README for more information

include("sim.jl")
include("genweights.jl")

Ntrials = 1

T = 20000
dt = 0.1
Ncells = [400, 50, 25, 25] #number of cells in each population: [e, pv, som, vip]
Ntot = sum(Ncells)
Npop = length(Ncells)

avcond = zeros(Npop,Ntot)
rates = zeros(Ntot)

for nn = 1:Ntrials
	println("trial ",nn)
	wind,wipost,wstr,wext,rext = genweights(Ncells)
	#control
	times,tinds,rates1,avcond1,vall = sim(T,dt,Ncells,wind,wipost,wstr,wext,rext,istim=2,curstim=0)

	#pv stimulation
	#times,tinds,rates1,avcond1,vall = sim(T,dt,Ncells,wind,wipost,wstr,wext,rext,istim=2,curstim=25)

	#som stimulation
	#times,tinds,rates1,avcond1,vall = sim(T,dt,Ncells,wind,wipost,wstr,wext,rext,istim=3,curstim=45)

	global rates += rates1 
end

rates /= Ntrials;

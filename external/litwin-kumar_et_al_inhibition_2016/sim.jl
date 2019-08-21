#this file is part of litwin-kumar_et_al_inhibition_2016
#Copyright (C) 2016 Ashok Litwin-Kumar
#see README for more information
using Random,Statistics
include("params.jl")

function sim(T,dt,Ncells,wind,wipost,wstr,wext,rext; istim=1,curstim=0)
	#initialization-------------------------------------------------------------
	@params()
	Ntot = sum(Ncells)
	Npop = length(Ncells)
	pinds = [1; 1 .+ cumsum(Ncells)] #start index of each population (ends with Ntot+1)

	#simulation parameters
	NT = round(Int,T/dt)
	recstart = 2000 #time at which to start analysis, to remove transient
	irecstart = round(Int,recstart/dt)

	#external stimulation
	curext = zeros(Ntot) #constant bias current to one population
	curext[pinds[istim]:(pinds[istim+1]-1)] .= curstim

	#state vectors
	v = -60. * ones(Ntot)
	vT = zeros(Ntot)
	lastSpike = -100. * ones(Ntot)
	whichpop = zeros(Int,Ntot)
	for pp = 1:Npop
		whichpop[pinds[pp]:(pinds[pp+1]-1)] .= pp
		vT[pinds[pp]:(pinds[pp+1]-1)] .= vTpop[pp] .+ sigvT[pp]*randn(Ncells[pp])
	end
	nextext = zeros(Ntot)
	for cc = 1:Ntot
		nextext[cc] = randexp()/rext[cc]
	end
	w_adapt = zeros(Ntot)
	xrise = zeros(4,Ntot) #for each cell, one set of ODEs for each presynaptic population 
	xdecay = zeros(4,Ntot)
	D = ones(4,Ntot) #for each cell, one set of ODes for each postsynaptic population
	F = ones(4,Ntot)

	#return values-------------------------------------------------------------
	#spike counts and times
	ns = 0
	maxns = round(Int,T*Ntot*.05)
	times = zeros(maxns)
	tinds = zeros(Int,maxns)
	vall = zeros(1)

	#average conductance each neuron received
	avcond = zeros(Npop,Ntot) 

	#simulation-------------------------------------------------------------
	println("starting sim")
	for tt = 1:NT
		t = dt*tt

		if mod(tt,NT/100) == 1 #print percent complete
			print("\r",round(Int,100*tt/NT),"%")
		end		

		for cc = 1:Ntot #update synaptic/adaptation parameters
			pc = whichpop[cc]

			w_adapt[cc] += (dt/tauw_adapt[pc])*(a_adapt[pc]*(v[cc] - EL[pc]) - w_adapt[cc])
			for qq = 1:Npop
				xrise[qq,cc] -= dt*xrise[qq,cc]/taurise[qq,pc]
				xdecay[qq,cc] -= dt*xdecay[qq,cc]/taudecay[qq,pc]
				D[qq,cc] += dt*(1-D[qq,cc])/tauD[pc,qq]
				F[qq,cc] += dt*(1-F[qq,cc])/tauF[pc,qq]
			end
		end

		for cc = 1:Ntot
			pc = whichpop[cc]

			if t > (lastSpike[cc]+tauref[pc]) #not in refractory period
				dv = curext[cc] - w_adapt[cc] - gL[pc]*(v[cc]-EL[pc]) + gL[pc]*deltaT[pc]*exp((v[cc]-vT[cc])/deltaT[pc])

				for qq = 1:Npop
					dv -= gSyn[qq,pc]*(v[cc] - Esyn[qq,pc])*(xdecay[qq,cc]-xrise[qq,cc])/(taudecay[qq,pc]-taurise[qq,pc])
				end
				v[cc] += dt*dv/Cm[pc]
			end #end if(refractory)
		end

		for cc = 1:Ntot #do spikes
			pc = whichpop[cc]
			if (v[cc] > vth[pc]) && (ns<maxns) #spike occurred
				v[cc] = vre[pc]
				lastSpike[cc] = t
				w_adapt[cc] += b_adapt[pc]
				ns = ns+1
				times[ns] = t
				tinds[ns] = cc

				for kk = wind[cc]:(wind[cc+1]-1) #propagate spike
					ipost = wipost[kk]
					ppost = whichpop[ipost]
					xrise[pc,ipost] += wstr[kk]*F[ppost,cc]*D[ppost,cc]
					xdecay[pc,ipost] += wstr[kk]*F[ppost,cc]*D[ppost,cc]
				end

				for qq = 1:Npop
					F[qq,cc] += UF[pc,qq]*(Fmax[pc,qq]-F[qq,cc])
					D[qq,cc] = D[qq,cc]*UD[pc,qq]
				end
			end
		end

		for cc = 1:Ntot
			while(t > nextext[cc])
				nextext[cc] += randexp()/rext[cc]
				xrise[1,cc] += wext[cc] #increment excitation
				xdecay[1,cc] += wext[cc]
			end
		end

		
		if tt > irecstart
			for cc = 1:Ntot
				pc = whichpop[cc]
				for qq = 1:Npop
					avcond[qq,cc] += gSyn[qq,pc]*(xdecay[qq,cc]-xrise[qq,cc])/(taudecay[qq,pc]-taurise[qq,pc])
				end
			end
		end

		
	end #end loop over time
	println("\r100%")

	println("calculating rates")
	counts = zeros(Ntot)
	for cc = 1:Ntot
		counts[cc] = sum((tinds.==cc) .* (times.>recstart))
	end

	rates = 1000. * counts / (T-recstart)

	for pp = 1:Npop
		println("population ",pp," firing rate: ",mean(rates[pinds[pp]:(pinds[pp+1]-1)]))
	end

	times = times[1:ns]
	tinds = tinds[1:ns]
	avcond = avcond/(NT-irecstart)

	return times,tinds,rates,avcond,vall
end #sim()

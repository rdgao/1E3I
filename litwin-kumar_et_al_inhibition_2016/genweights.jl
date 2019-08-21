#this file is part of litwin-kumar_et_al_inhibition_2016
#Copyright (C) 2016 Ashok Litwin-Kumar
#see README for more information

include("params.jl")

#generates connectivity: recurrent and external weights and external input rates
function genweights(Ncells)
	@params()
	Ntot = sum(Ncells)
	Npop = length(Ncells)

	#set up recurrent weight matrix
	syncount = 1;
	Maxw = round(Int,Ntot*Ntot*0.3) #maximum number of weights in the weight matrix
	wind = zeros(Int,Ntot+1) #column of w corresponding to the start of the ith neuron's projections
	wipost = zeros(Int,Maxw)
	wstr = zeros(Maxw)

	syncount = 1
	for pp = 1:Npop
		theta_p = range(-pi/2,stop=pi/2,length=Ncells[pp])

		for cc = 1:Ncells[pp]
			wind[cc + sum(Ncells[1:(pp-1)])] = syncount

			theta = theta_p[cc]
			
			for qq = 1:Npop
				theta_q = range(-pi/2,stop=pi/2,length=Ncells[qq])

				dtheta = theta .- theta_q
				dtheta = min.(abs.(dtheta),pi.-abs.(dtheta))

				prob = p0[pp,qq]*(1 .+ p2[pp,qq]*cos.(2*dtheta))

				iconns = findall(rand(Ncells[qq]) .< prob) .+ sum(Ncells[1:(qq-1)])

				wipost[syncount:(syncount+length(iconns)-1)] = iconns
				wstr[syncount:(syncount+length(iconns)-1)] .= J[pp,qq]

				syncount = syncount+length(iconns)
			end
		end
	end
	wind[Ntot+1] = syncount
	wipost = wipost[1:syncount]
	wstr = wstr[1:syncount]

	rext = zeros(Ntot)
	pinds = [1; 1 .+ cumsum(Ncells)] #start index of each population (ends with Ntot+1)
	for pp = 1:Npop
		theta_p = range(-pi/2,stop=pi/2,length=Ncells[pp])
		rext[pinds[pp]:(pinds[pp+1]-1)] = r0[pp]*(1 .+ r2[pp]*cos.(2*theta_p))
	end

	#synaptic strength of external inputs.  set to 1 and just edit gSyn
	wext = ones(Ntot)


	return wind,wipost,wstr,wext,rext
end #genweights()

#generates connectivity for testing PSPs (one neuron from specified population projection to one of each other population
function genweights_PSPs(Ncells,popPSPs)
	@params()
	Ntot = sum(Ncells)
	Npop = length(Ncells)

	wext = zeros(Ntot)
	rext = zeros(Ntot)

	ind = 2*popPSPs

	wind = zeros(Int,Ntot+1)
	wipost = zeros(Int,Npop)
	wstr = zeros(Npop)
	wstr[:] .= J[popPSPs,:]
	wipost[:] .= 2*(1:Npop) - 1

	#no one projects to anyone except index ind, who projects one synapse to each pop.
	syncount = 1
	for cc = 1:Ntot
		wind[cc] = syncount
		if cc == ind
			syncount += Npop
		end
	end
	wind[end] = syncount

	return wind,wipost,wstr,wext,rext 
end #genweights_PSPs()

#generates connectivity: recurrent and external weights and external input rates
function genweights_shape!(Ncells,istim)
	@params()
	Ncells[istim] *= 2
	Ntot = sum(Ncells)
	Npop = length(Ncells)

	#set up recurrent weight matrix
	syncount = 1;
	Maxw = round(Int,Ntot*Ntot*0.3) #maximum number of weights in the weight matrix
	wind = zeros(Int,Ntot+1) #column of w corresponding to the start of the ith neuron's projections
	wipost = zeros(Int,Maxw)
	wstr = zeros(Maxw)

	syncount = 1
	for pp = 1:Npop
		theta_p = range(-pi/2,stop=pi/2,length=Ncells[pp])
		if pp == istim
			theta_p = repmat(range(-pi/2,stop=pi/2,length=round(Int,Ncells[pp]/2)),2)
		end

		for cc = 1:Ncells[pp]
			wind[cc + sum(Ncells[1:(pp-1)])] = syncount

			if (pp == istim) && (cc > round(Int,Ncells[pp]/2))
				continue
			end

			theta = theta_p[cc]

			for qq = 1:Npop
				theta_q = range(-pi/2,stop=pi/2,length=Ncells[qq])
				if qq == istim
					theta_q = repmat(range(-pi/2,stop=pi/2,length=round(Int,Ncells[qq]/2)),2)
				end

				dtheta = theta - theta_q
				dtheta = min.(abs.(dtheta),pi-abs.(dtheta))

				prob = p0[pp,qq]*(1 .+ p2[pp,qq]*cos.(2*dtheta))

				iconns = findall(rand(Ncells[qq]) .< prob) + sum(Ncells[1:(qq-1)])

				wipost[syncount:(syncount+length(iconns)-1)] = iconns
				wstr[syncount:(syncount+length(iconns)-1)] = J[pp,qq]

				syncount = syncount+length(iconns)
			end
		end
	end
	wind[Ntot+1] = syncount
	wipost = wipost[1:syncount]
	wstr = wstr[1:syncount]

	rext = zeros(Ntot)
	pinds = [1; 1+cumsum(Ncells)] #start index of each population (ends with Ntot+1)
	for pp = 1:Npop
		theta_p = range(-pi/2,stop=pi/2,length=Ncells[pp])
		if pp == istim
			theta_p = repmat(range(-pi/2,stop=pi/2,length=round(Int,Ncells[pp]/2)),2)
		end
		rext[pinds[pp]:(pinds[pp+1]-1)] = r0[pp]*(1+ r2[pp]*cos.(2*theta_p))
	end

	#synaptic strength of external inputs.  set to 1 and just edit gSyn
	wext = ones(Ntot)


	return wind,wipost,wstr,wext,rext
end #genweights_shape!()


#generates connectivity: recurrent and external weights and external input rates
function genweights_current!(Ncells,istim,thetastim)
	@params()
	Ncells[istim] *= 2
	Ntot = sum(Ncells)
	Npop = length(Ncells)

	#set up recurrent weight matrix
	syncount = 1;
	Maxw = round(Int,Ntot*Ntot*0.3) #maximum number of weights in the weight matrix
	wind = zeros(Int,Ntot+1) #column of w corresponding to the start of the ith neuron's projections
	wipost = zeros(Int,Maxw)
	wstr = zeros(Maxw)

	syncount = 1
	for pp = 1:Npop
		theta_p = range(-pi/2,stop=pi/2,length=Ncells[pp])
		if pp == istim
			theta_p = [range(-pi/2,stop=pi/2,length=round(Int,Ncells[pp]/2)),thetastim*ones(round(Int,Ncells[pp]/2))]
		end

		for cc = 1:Ncells[pp]
			wind[cc + sum(Ncells[1:(pp-1)])] = syncount

			if (pp == istim) && (cc > round(Int,Ncells[pp]/2))
				continue
			end

			theta = theta_p[cc]
			
			for qq = 1:Npop
				theta_q = range(-pi/2,stop=pi/2,length=Ncells[qq])
				if qq == istim
					theta_q = [range(-pi/2,stop=pi/2,length=round(Int,Ncells[qq]/2)),thetastim*ones(round(Int,Ncells[qq]/2))]
				end

				dtheta = theta - theta_q
				dtheta = min.(abs.(dtheta),pi-abs.(dtheta))

				prob = p0[pp,qq]*(1 + p2[pp,qq]*cos.(2*dtheta))

				iconns = findall(rand(Ncells[qq]) .< prob) + sum(Ncells[1:(qq-1)])

				wipost[syncount:(syncount+length(iconns)-1)] = iconns
				wstr[syncount:(syncount+length(iconns)-1)] = J[pp,qq]

				syncount = syncount+length(iconns)
			end
		end
	end
	wind[Ntot+1] = syncount
	wipost = wipost[1:syncount]
	wstr = wstr[1:syncount]

	rext = zeros(Ntot)
	pinds = [1; 1+cumsum(Ncells)] #start index of each population (ends with Ntot+1)
	for pp = 1:Npop
		theta_p = range(-pi/2,stop=pi/2,length=Ncells[pp])
		if pp == istim
			theta_p = [range(-pi/2,stop=pi/2,length=round(Int,Ncells[pp]/2)),thetastim*ones(round(Int,Ncells[pp]/2))]
		end
		rext[pinds[pp]:(pinds[pp+1]-1)] = r0[pp]*(1+ r2[pp]*cos.(2*theta_p))
	end

	#synaptic strength of external inputs.  set to 1 and just edit gSyn
	wext = ones(Ntot)


	return wind,wipost,wstr,wext,rext
end #genweights_current!()



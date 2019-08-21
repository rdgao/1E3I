#this file is part of litwin-kumar_et_al_inhibition_2016
#Copyright (C) 2016 Ashok Litwin-Kumar
#see README for more information

macro params()
quote

println("unpacking parameters")

pvtuned = false

#baseline probability of connection from pre->post
p0 = [0.1  0.6  0.6  0.6;
      0.6  0.6  0.   0.6;
      0.6  0.6  0.   0.6;
      0.   0.   0.4  0.  ]
  
#tuned component of connections from pre->post (normalized to p0)
#p_ij(dtheta) = p0_ij*(1 + p2_ij*cos(2*dtheta) )
p2 = [0.8  0.1  0.  0. ;
      0.8  0.1  0.  0. ;
      0.   0.   0.  0. ;
      0.   0.   0.  0.  ]

#synaptic strength for existing connection from pre->post.  set to 1 and just use gsyn
J = [1. 1. 1. 1.;
     1. 1. 1. 1.;
     1. 1. 1. 1.;
     1. 1. 1. 1. ]

#set up external input rates
#r0 = [3*2.4, 3*0.4, 0.8, 0.3] #high rate
r0 = [2.4, 0.4, 0.8, 0.3]
r2 = [0.2, 0.2, 0., 0.]

if !pvtuned
	p2[2,:] .= 0.
	r2[2] = 0.
end
   
#synaptic parameters, Matrix for pre- and post- pop., nS
gSyn = [  .1/60  .3/60  .05/60   .1/60;
          3.0/22  3.0/22   0.      .6/22;
         1.5/22  1.0/22   0.0/22  2.5/22;
          0.     0.      3.0/22   0.    ]*1000.
#reversal potential, mV
Esyn = [  0.   0.   0.   0.;
        -67. -67. -67. -67.;
        -67. -67. -67. -67.;
        -67. -67. -67. -67. ]

#rise time, ms
taurise = [ 0.5  0.5  0.5  0.5;
            0.5  0.5  0.5  0.5;
            1.0  1.0  1.0  1.0;
            1.0  1.0  1.0  1.0 ]

#decay time, ms
taudecay = [2. 2. 2. 2.;
            3. 3. 3. 3.;
            4. 4. 4. 4.;
            4. 4. 4. 4. ]

#depression time constant, ms
tauD = [800. 800. 800. 800.; 
        800. 800. 800. 800.;
        800. 800. 800. 800.;
        800. 800. 800. 800. ]

#amount of depression
UD = [0.75  0.75  1.   1. ;
      0.9   0.9   0.9  0.9 ;
      1.    1.    1.   1. ;
      1.    1.    1.   1.  ]

#facilitation time constant, ms
tauF = [200. 200. 200. 200.;
        200. 200. 200. 200.;
        200. 200. 200. 200.;
        200. 200. 200. 200. ]

#amount of facilitation
UF=[0.  0.  0.5  0. ;
    0.  0.  0.   0. ;
    0.  0.  0.   0. ;
    0.  0.  0.   0.  ]

#maximum facilitation
Fmax = [2. 2. 2. 2.; #changed here, facilitation decreased for E->SOM
        2. 2. 2. 2.;
        2. 2. 2. 2.;
        2. 2. 2. 2. ]

#cell parameters, one value for each pop
Cm = [180., 80., 80., 80.] #capacitance, pF
gL = [1/0.16, 1/0.1, 1/0.2, 1/0.2] #leak conductance, nS

tau = Cm./gL

EL = [-60., -60., -60., -60.] #leak voltage, mV
deltaT = [1., .25, 1., 1.] #EIF slope parameter, mV
vTpop = [-40., -40., -45., -45.] #EIF threshold voltage, mV
sigvT = [3., 3., 3., 3.] #stdev of the above
vth = [20., 20., 20., 20.] #threshold voltage, mV
vre = [-60., -60., -60., -60.] #reset voltage, mV
tauref = [2., 2., 2., 2.] #refractory period, ms
tauw_adapt = [150., 150., 150., 150.] #adaptation timescale, ms

a_adapt = [4., 0., 4., 4.] #adaptation slope, nS
b_adapt = 10.0*[0.8, 0., 0.8, 0.8]  #adaptation increment, pA

end |>esc
end

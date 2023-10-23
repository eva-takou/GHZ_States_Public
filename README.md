# GHZ_states_Preparation

Scripts used to produce the results of the paper "Generation of genuine all-way entanglement in defect-nuclear spin systems through
dynamical decoupling sequences", by E. Takou, E. Barnes, and S. E. Economou, arXiv:2302.05580.

We use decoupling sequences that drive the central electronic spin so as to manipulate the entanglement within the electron-nuclear spin register.


There are three main classes used throughout the simultations. 

-The class "SuperClass_Sequences" contains various sequences such as the CPMG, UDD and XY2 sequences. 

-The classs "SubClass_U4Operations" inherits the properties from the class "SuperClass_Sequences", and can be used to simulate the evolution of a single nuclear spin coupled with the electron, to calculate the Makhlin invariants or entangling power of two-qubit gate.  

-The class "SubClass_Ent_and_Fid" can be used to:

i) Calculate the Kraus operators related to unwanted nuclear spins.

ii) Find the gate error of a target evolution due to the presence of unwanted spins.

iii) Calculate one-tangles, a metric of correlations defined on the evolution operator level, and quantifies the entanglement between a single spin and the remaining register.

iv) Calculate M-way entangling metrics, such as the unitary M-way entangling power, the non-unitary M-way entangling power, the "dephased" M-way entangling power. All these are defined for CR-type evolution operators. 

v) Calculate the unitary M-way entangling power of an arbitrary evolution operator U for even M, and the one-tangles for an arbitrary multipartite evolution operator U. 


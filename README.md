# hydrogenbond_reaction_flux
 0. This function is used to calculate the rate of relaxation to equilibrium, ie. the reactive flux hydrogen bond (HB) correlation function
k(t). The function -k(t) measures the average rate of change of HB population (initial set of hydrogen bonds) for the bond between a 
tagged pair of molecules is broken at time t later.

 1. To run this function, run the command:
$./rfhb_in_and_n.x < input_hbacf

 2. The content of input file "input_hbacf" is:
0.0005
118_2LiI
traj_pos.xyz
118_2LiI_Chandra_list.dat
50000
358
13806
50

The meaning of each parameters are as follows.
0.0005         # the original time step delta_t0 in the trajectory file (unit: ps) 
118_2LiI       # the name of the system
traj_pos.xyz   # the trajectory file
118_2LiI_Chandra_list.dat    # the "list file", which give all the indices of atoms in each tagged pair of molecules
1000           # the total steps. Usually, the total steps must be at least 20000 step (10-ps trajectory)
358            # the total number of atoms in the system
13806          # the number of tagged pairs of molecules. This value equals to the lines of the valid lines of the "list file".
50             # the number of skipping steps (n). The new time step delta_t can be determined by: delta_t=delta_t0


Ref.: A.Luzar, J. Chem. Phys. 113,23, 10663

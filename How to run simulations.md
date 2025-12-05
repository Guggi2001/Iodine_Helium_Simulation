# Single pulse simulation

to start a simulation run, open the script
*run_simulation.m*

in this script, the parameters of the simulation are set
and input files can be run

the input files set  global variables
(like *single_pulse*, which determines if a single pulse or pump-probe simulation is run)


## comparison to HeDFT (neutrals)

use the input file *inputfiles\pumpprobe\HeDFT_mimic_neutral_and_ion_N2000_centerstart*
and run only the neutral simulation

then use 
*\I2HeN_velocity_simulation\HeDFT_MD_comparison_neutral\compare_neutral_dynamics_to_HeDFT.m*
to compare the neutral velocities to neutral HeDFT results

## Comparison to HeDFT (ions)

![[single_pulse_inputs.png]]

uncomment the line shown in the image.

*single_pulse_N2000.m*
contains the input parameters needed to create a simulation run that 
can be compared to the HeDFT simulation

most critical are the crosssection and ion binding energy (found in *single_pulse_N2000.m*)
![[crosssection_and_binding_energy.png]]


after the simulation has run, execute
*\I2HeN_velocity_simulation\single_pulse_simulation\HeDFT_comparison\simulation_image_only_trajectories.m*
to check how the ion MD velocities compare to the HeDFT results.


for the newer 18 Angström HeDFT data
the script *single_pulse_N2000_18Angst.m*
was used. A scattering crosssection of 1600 $A^2$ and a binding energy of only 50 meV was sufficient to reproduce the HeDFT trajectory. 



To simulate the droplet size used in the experiment, with scaled down Coulomb potential, select the inputfile *single_pulse_droplet_distributions.m*

make sure to adjust binding energy and crosssection to match the result used to match the HeDFT calculation.

once the results of this simulation agrees with the measured velocity distributions,
use 
*\I2HeN_velocity_simulation\single_pulse_simulation\HeDFT_comparison\simulation_image.m*
to generate Figure 6.8 from the thesis. 

Figure 6.3 is generated with the script
*post_process_single_pulse_paper_IplusHe_comparison_cov.m*

# Pumpprobe simulation

to generate 7.18 b and 7.17  of the thesis use the 
inputfile: 
*inputfiles_dft_comparison\single_pulse_droplet_distribution_T_limit_9Angstroem.m*

and the script
*simulation_image_pumpprobe_comp*


full neutral + ion dynamics simulations are run with the input file

*inputfiles\pumpprobe\decreased_sigma0.m* 
for example 

and are visualized immediately after executing the *run_simulation.m* script

comparison to a timescan can be made with the script
*vmi_sim_compare_to_timescan.m*


# current state
the 18 Angström DFT trajectory can be reached with 
sigma0 = 1900 Angström^2
and E_bind = 0.00 meV

however, the same parameters lead to I+He from the 2.666 Angström simulation that are too fast.

on the other hand,
 E_bind = 0 meV
and a sigma0 of 5000 Angström^2
leads to the correct I+He distribution from 2.666 Angström 
![[single_pulse_gs_sigma05000.png]]
and from 9 Angström (compared to pumpprobe distribution at 200 ps )
(inputfile *inputfiles_dft_comparison\single_pulse_droplet_distribution_T_limit_9Angstroem_zero_bind.m*)

this result suggests a sigma(v) = sigma_0 v^-n  
where n>2

(this would lead to less deceleration for simulation starting from 18 Angström  and more simulation for simulation starting from 2.66 Angström



the full neutral and ion dynamic simulation suffers from the central issue that once the I+ are further apart than about 9 Angström, the MD ion velocities are much lower than the expected ion velocities at 200 ps.
furthermore, there is no good way to model the recombination of the molecule at 9 Angström


determining the correct crosssection would be a crucial step towards a full pumpprobe dynamic simulation


NATRC (Non-Adiabatic  Transfer Rates Calculations)

1. Description
NATRC is the program designed for calculation of non-adibatic transition rates between two states in molecules. Calculations are carried out based on the classical and modified Bixon-Jortner-Plotnikov models.

2. Requirements
Linux OS 32/64-bit
gfortran ≥4.8 
Intel MKL ≥10.2.2

3. Installation
Download NATRC from the https://github.com/Arthur140/NATRC.
Unpack the downloaded archive in any folder and enter it.

  cd NATRC

Compile the program using the command below. In this command, you have to replace 'mkllibpath' with the path to the mkl library that is relevant for your system.

  gfortran NATRC.f90 -Wl,--start-group mkllibpath/libmkl_gf_lp64.a mkllibpath/libmkl_sequential.a mkllibpath/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl -o NATRC

If your system is 32-bit, then replace 'libmkl_gf_lp64.a' with 'libmkl_gf.a'

After compilation, a NATRC file will be generated in the folder. This file is a program for calculating non-adiabatic transfer rates. To run the program, write the next command:

  ./NATRC test | tee -a test.log

If the program works correctly, you will obtain 'Kic=3.334E+06 sec**-1'

4. The input data
NATRC use output files of GAMESS as input files for performing calculations. You should have a point energy calculation for equilibrium geometries of the initial and final states, two hessian files for the equilibrium geometry of the ground state. This files can be performed at the RHF, ROHF, UHF, DFT, TD-DFT, or MCSCF levels of theory.
Also you should have either a file with MCSCF/NACME calculation or DFT/ELFLDG calculation for the equilibrium geometry of the initial state.
All these calculations can be performed with solvation models.

5. Preparation the input file
The file name must be written according to the pattern 'name.inp'.
When you run program, you should skip '.inp'

  ./NATRC name

The input file includes a list of parameters and their values.
All parameters in the input file has to be described as:

 parameter: value

One line should contain only one parameter. The comment is indicated with an exclamation mark! Lines containing only comments are allowed.
 
 parameter: value !comment

You can only set a limited set of keywords as a parameter. A colon is used for separation parameter and its value. A value can be either a text or a number depending on a parameter. 

6. Description of input parameters
 
 DOS:
 The approximation of the density of spectra (DOS): Pekar, Hybrid or Gauss.
 (Default is Pekar)

 threshold: !Only for DOS=Hybrid
 The threshold parameter for chosing vibronic modes which will be approximated by the Gaussian function.

 number_states_nacme: 
 The number of states for which the NACME calculation was performed. It can only be ≥2.

 initial_state:
 The number of intial state where the ground state is of 1 and the first excited stated is 2. This parameter can only be ≤number_states_nacme

 final_state:
 The number of the final state where the ground state is of 1 and the first excited stated is 2. This parameter can only be <number_states_nacme and <initial_state.

 Eif:
 The energy difference (in Hartree) between the initial and final states at their equilibrium geometries.

 Temperature: !Optional
 It's temperature of the system (default is 298.15).

 initial_coord_file:
 The path to log-file of one-point calculation which contains the equilibrium geometry of the initial state.

 final_coord_file:
 The path to log-file of one-point calculation which contains the equilibrium geometry of the final state.

 hess_file:
 The path to dat-file which contains the hessian matrix of studied molecules at the ground state

 hess_coord:
 The path to log-file which contains the coordinats for which the hessian matrix were obtained.

 grad_file:
 The path to log-file which contains gradient of the final state which was calculated for the equlibrium geometry of the initial state. 

 naccalc:
 Flag of evaluation of NACME basing on transition integrals of Coulomb field. If you set 'true', the program will evaluate NACME using these integrals. If you set 'false',  the program will read NACME from an input file (nacme_file).
 (default is 'false')

 nacme_file: !Only for naccalc=false
 The path to the file which contains NACME obtained at the MCSCF level of theory

 elf_file: !Only for naccalc=true
 The path to the file which contains transition integrals of Coulomb field obtained at the TD-DFT level of theory

 wdebye: !Optional
 The Debye temperature (in Hartree) for the solvent. If wdebay>0, the program will take into account the solvatation effect. 

 Esolv: !Only for wdebye>0
 The reorganization energy (default is 0)

 symmetry: !Optional
 Flag of molecular symmetry. You should set 'true' for symmetric molecules and 'false' for unsymmetrical molecules (default is 'false')

 total_symmetry: !Only for symmetry=true
 Flag of transition symmetry. You have to set 'true' for totally symmetric transition (for instance, Ag->Ag, B3u->B3u or A"->A") and you have to set 'false' for non-totally symmetric transition (for instance, B2u->Ag, A"->A').
 
 NumRotTrM:
 The number of translational and rotational modes that are excluded from the calculation.
 (Default is 6)

 cutoff: !Only for Hybrid
 This parameter controls the number of modes that will be excluded from the calculation. The smaller it is, the greater the number of modes included in the calculation. Calculation is carried out for all modes at cutoff=0, and calculation isn't performed at cutoff=1.0. (Default 0.05)

 deep:
 All quantum numbers that result in values ​​less than the maximum multiplied by deep are cut off. If deep=1.0, then the quantum numbers will not vary. If deep=0.0, then the calculation will use all quantum numbers that give a speed different from 0 in the calculation.
(Default is 0.001)

7. Authors
 Arthur I. Martynov - Programming, software development; implementation of the computer code and supporting algorithms.
 Alexander S. Belov - Testing of existing code components.

8. License
 This project is licensed under the GNU GPL v3.0 License—see COPYING for details.
 
 
 





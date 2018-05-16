% /***************************************************************************
% Copyright 
% MUSiiC Laboratory
% Haichong Zhang,Emad M Boctor
% Johns Hopkins University
% 
% For commercial use/licensing, please contact Emad Boctor, Ph.D. at eboctor@jhmi.edu.
% ***************************************************************************/

Simulation package for "A Phantom with Multiple Active Points for Ultrasound Calibration"

Main execution code
	main_AXXB_2d_Simulation_JMI.m 
	- Run simulation for AXXB calibration based on multiple active points using a 1D arrary ultrasound probe

Functions
	decompose_rotation_d.m
	- Decompose a rotation matrix to Euler angle representation

	funAXXB.m
	- Function to solve AXXB through minimization

	funPhantomModel.m
	- Build the phantom model

	funPreparingAXXB.m
	- Initializing randomized variables for simulation

	funRecoverA_2invariants.m
	- Recover full superscript A from partially recovered A

	funRecoverPartialA.m
	- Recover superscript A with 1DoF unknown rotation axis

	funRecoverPartialA_OneMidPlane.m
	- Recover full supoerscript A assuming one ultrasound plane cross a point while elevational component is zero

	param_extract2.m
	- Calculating invariants

	solveAXXB.m
	- Least square based AXXB solver 
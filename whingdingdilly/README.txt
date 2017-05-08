The objective of this directory is to locally write, run, and analyze simulations using HOOMD, GSD, and Freud
This task is broken down accordingly:

Writing Files:
	Parameters:
		-x_a  = particle fraction of particle whose activity is varied
		-Pe_a = activity of species being varied
		
	1.)Enter the above parameters and how to change them (stepsize)
	2.)Using these inputs, generate simulation aspect of infile
	
Note: HOOMD scripts allow all of this to be done in one file (so analysis will have to be edited as well)
Note: you can access system values within a script (like particle number) via a snapshot

	3.)Params cascade into analysis methods...
		Note: How should we be saving this output (need hard data)
	4.)Run infiles with SLURM manager (batch submission)
	
Jump for joy as your computer does your job :)


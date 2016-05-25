#ge2-processing
A hexrd, numpy, scipy based collection of Python utilities to process GE2 files.

Currently spot IDing and cleaning is implemented.

To run on W315 workstation:
`python smooth_ge_data.py config.yml` OR `ge-spot-cleanup config.yml`

To run on Stampede:
Create a job file (e.g. `spot_cleanup.job`) and submit it with sbatch

	#!/bin/bash

	#SBATCH -J spot_cleanup          # Job name
	#SBATCH -o spot_cleanup.%j.out   # stdout; %j expands to jobid
	#SBATCH -e spot_cleanup.%j.err   # stderr; skip to combine stdout and stderr
	#SBATCH -p largemem              # queue
	#SBATCH -N 1                     # Number of nodes, not cores (16 cores/node)
	#SBATCH -n 32                    # Total number of MPI tasks (if omitted, n=N)
	#SBATCH -t 00:30:00              # max time

	python /home1/04132/harshad/ge_spot_cleanup/ge2-processing/smooth_ge_data.py config.yml`

This will read the GE2 files specified in the config.yml file and detect centers of the spots.
Then it will write the spots to a new GE2 file named `synth_spots.ge2` or the spot info to a txt file based
on config options.

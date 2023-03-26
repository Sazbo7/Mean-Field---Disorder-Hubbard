for INPUT_FILENAME in DisHB_disint_36_*.fits
do

	cat <<-EOF > script
	#!/bin/bash

	#PBS -l nodes=1:ppn=4
	#PBS -l walltime=20:00:00
	#PBS -l mem=5gb
	#PBS -j oe
	#PBS -N HF_SOLVER_${INPUT_FILENAME%%.fits}

	module load intel/18.0.2
	source activate local

	export MKL_NUM_THREADS=4
	export OMP_NUM_THREADS=4
	
	cd Trivedi_Group/fits_files_toAnalyze
	
	hostname
	pwd
	free
	echo
	
	./HF_U_repulsive_mu.x ${INPUT_FILENAME}

	EOF

	qsub script

done

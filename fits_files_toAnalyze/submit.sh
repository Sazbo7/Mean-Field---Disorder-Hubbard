for INPUT_FILENAME in fs_cm_0_44*.fits
do

	cat <<-EOF > script
	#!/bin/bash

	#PBS -l nodes=1:ppn=6
	#PBS -l walltime=80:00:00
	#PBS -l mem=8gb
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
	
	./HF_U_repulsive_UNITY.x ${INPUT_FILENAME}

	EOF

	qsub script

done

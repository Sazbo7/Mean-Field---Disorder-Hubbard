for INPUT_FILENAME in fs_cm_0_40*.fits
do

	cat <<-EOF > script
	#!/bin/bash

	#PBS -l nodes=1:ppn=4
	#PBS -l walltime=4:00:00
	#PBS -l mem=5gb
	#PBS -j oe
	#PBS -m abe
	#PBS -M szabo.48@osu.edu
	#PBS -N LDOS_${INPUT_FILENAME%%.fits}

	module load intel/18.0.2
	source activate local

	export MKL_NUM_THREADS=4
	export OMP_NUM_THREADS=4
	
	cd Trivedi_Group/fits_files_toAnalyze
	
	hostname
	pwd
	free
	echo
	
	./U_repulsive_append_ldos_to_fits.x ${INPUT_FILENAME} 0.00625 -7.25 7.25 5001 LDOS_0

	EOF

	qsub script

done

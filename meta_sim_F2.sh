#!/bin/bash
# This runs main_workflow.sh to perform several simulations while variating the RD (Read Depth) and MPS (Mapping Mopulation Size) parameters. 

# command structure: 	./meta_sim.sh meta_name in_fasta map_pop
# test command: 		./meta_sim.sh testing check.1.fa M2



######################################################### STRATEGY #########################################################
#
# 1- Simulation of CHR4-Lab -> lab.va with variants of lab strain and lab.fa with the sequence
# 2- lab.fa is used as a base for the second mutagenesis 
# 4- Analysis of output
# 5- Meta-analysis
#
############################################################################################################################


# 1) Initial parameters and accomodations ______________________________________________________________________________________________________________________________________________________________________________________________________
timestamp=$(date "+%F-%T")

meta_name=u_projects/$timestamp"_"$1
in_fasta=$2
map_pop=$3
mkdir $meta_name
mkdir $meta_name/META
meta_folder=$meta_name/META
nbr_mutations=47																			# <------------------------- not used
mut_pos='1,5845220'

my_meta_log=$meta_folder/meta.log; touch $my_meta_log
my_meta_info=$meta_folder/meta_info.txt; touch $my_meta_info
echo "#RD	MPS	CANDIDATES_95	SPAN_95	CANDIDATES_98	SPAN_98" >> $my_meta_info

nbr_background_mutations=109																# <------------------------- Backcross
#nbr_background_mutations=175000															# <------------------------- Outcross

#rd_list=(10)																				# <------------------------- SET
#mps_list=(10)																				# <------------------------- SET
mps_list=(40 80 160 320)
rd_list=(30 60 200)
export location="$PWD" 			#Save path to bowtie2-build and bowtie2 in variable BT2


# 4) Running the individual simulations ____________________________________________________________________________________________________________________________________________________________________________________________________________________

rec_freq_distr='0,24-1,43-2,25-3,6-4,1-5,1'							     		# <------------------------- SET

n_jobs=0
maxjobs=${SLURM_CPUS_ON_NODE} 							 									# <------------------------- SET Number of CPUs

for n in `seq 100`; do 							 								# <------------------------- SET Number of replicates
	for i in ${rd_list[@]}; do
		for j in ${mps_list[@]}; do
			for p in *.m4a ; do
				(
				rd=$i
		        mps=$j

		        project_name=$meta_name/$rd"_"$mps"_"$n

		        ./main_workflow_F2.sh $project_name $in_fasta $nbr_mutations $rec_freq_distr $mut_pos $mps $rd $meta_name $map_pop

		        echo $project_name ' done!'
		        ) &

				jobs_running=$(jobs -p | wc -w)
				while [ $jobs_running -ge "${maxjobs}" ]; do
				     sleep 30
				     jobs_running=$(jobs -p | wc -w)
				done
			done
		done
	done
done


wait


# 5) Analizing the metadata _______________________________________________________________________________________________________________________________________________________________________________________________________________
{
	python2 ./an_scripts/meta-analysis.py -meta_in $my_meta_info -out $meta_folder/averaged_data_95.txt -step 95 2>> $my_meta_log
	python2 ./an_scripts/meta-analysis.py -meta_in $my_meta_info -out $meta_folder/averaged_data_98.txt -step 98 2>> $my_meta_log

} || {
	echo $(date "+%F > %T")': Error during execution of meta-analysis.py' >> $my_meta_log
	exit_code=1; echo $exit_code; exit

}
echo $(date "+%F > %T")': meta-analysis.py finished.' >> $my_meta_log


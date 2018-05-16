# This workflow (1) simulates read files from different bulked segregant mapping-by-sequencing experimental setups, (2) aligns the reads to a reference genome,
# (3) generates .vcf and .va files and processes them into relevant data


# command structure: 	./main_workflow.sh project_name in_fasta nbr_mutations
# test command: 		./main_workflow.sh testing check.1.fa 20



# _________________________________________________________________Some initial arrangements______________________________________________________________________________________
timestamp=$(date "+%F-%T")

# Get command arguments and assign them to variables
project_name=u_projects/$timestamp"_"$1
in_fasta=$2
nbr_mutations=$3

# Set 'exit_code' (flag variable) to it's initial value (0)
exit_code=0

# Store the location of each folder in a variable
f0=u_data
f1=$project_name/1_intermediate_files
f2=$project_name/2_logs; my_log_file=$f2/log.log
f3=$project_name/3_output

# Creating directories
mkdir $project_name
mkdir $f1
mkdir $f2
mkdir $f3

touch $my_log_file

sim_mut_output_folder=$f1/sim_data/sim_mut_output/mutant_strain

sim_recsel_output_folder_recessive=$f1/sim_data/sim_recsel_output_recessive
sim_recsel_output_folder_dominant=$f1/sim_data/sim_recsel_output_dominant



##################################################################################################################################################################################
#																																												 #
#																																												 #
#																		SIMULATOR MODULE																						 #
#																																												 #
#																																												 #
##################################################################################################################################################################################


# 1) Simulation of mutagenesis with sim_mut.py

{
	python2 simulator/sim-mut.py -nbr $nbr_mutations -mod d -con $f0/$in_fasta -out $sim_mut_output_folder 2>> $my_log_file

} || {
	echo $(date "+%F > %T")": Simulation of mutagenesis failed. Quit." >> $my_log_file
	exit_code=1; echo $exit_code; exit
}
echo $(date "+%F > %T")": Simulation of mutagenesis completed." >> $my_log_file



# 2) Simulating recombination with sim_rec.py 
# This workflow (1) simulates read files from different bulked segregant mapping-by-sequencing experimental setups, (2) aligns the reads to a reference genome,
# (3) generates .vcf and .va files and processes them into relevant data


# command structure: 	./main_workflow.sh project_name in_fasta nbr_mutations rec_freq_distr mut_pos nbr_rec_chrs read_depth
# test command: 		./main_workflow.sh testing check.1.fa 40 0,24-1,42-2,25-3,6-4,1-5,2 1,10000 200 20



# _________________________________________________________________Some initial arrangements______________________________________________________________________________________
timestamp=$(date "+%F-%T")

# Get command arguments and assign them to variables
project_name=u_projects/$timestamp"_"$1
in_fasta=$2
nbr_mutations=$3

rec_freq_distr=$4
mut_pos=$5
nbr_rec_chrs=$6

read_depth=$7


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
sim_recsel_output_folder=$f1/sim_data/sim_recsel_output
sim_seq_output_folder=$f1/sim_data/sim_seq_output/sample



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
parmut_sample=$sim_mut_output_folder/mutated_genome/mutated_genome.fa
parpol_sample=$f0/$in_fasta

{
	python2 simulator/sim-recsel.py -outdir $sim_recsel_output_folder -rec_freq_distr $rec_freq_distr -parmut $parmut_sample -parpol $parpol_sample -mutpos $mut_pos -smod r -nrec $nbr_rec_chrs 2>> $my_log_file 

} || {
	echo $(date "+%F > %T")": Simulation of recombination and phenotype selection failed. Quit." >> $my_log_file
	exit_code=1; echo $exit_code; exit
}
echo $(date "+%F > %T")": Simulation of recombination and phenotype selection completed." >> $my_log_file


# 3) Simulating HTS reads
lib_type=pe 									#<------------- Comprobar y establecer parametros por defecto
read_length_mean=100
read_length_sd=0
fragment_length_mean=500
fragment_length_sd=100
basecalling_error_rate=1
gc_bias_strength=50

{
	python2 simulator/sim-seq.py -input_folder $sim_recsel_output_folder -out $sim_seq_output_folder -mod $lib_type -rd $read_depth -rlm $read_length_mean -rls $read_length_sd -flm $fragment_length_mean -fls $fragment_length_sd -ber $basecalling_error_rate -gbs $gc_bias_strength 2>> $my_log_file

} || {
	echo $(date "+%F > %T")": Simulation of high-throughput sequencing on F2 recombinant population failed. Quit." >> $my_log_file
	exit_code=1; echo $exit_code; exit
}
echo $(date "+%F > %T")": Simulation of high-throughput sequencing reads on F2 recombinant population completed." >> $my_log_file



##################################################################################################################################################################################
#																																												 #
#																																												 #
#																		READ PROCESSING MODULE																					 #
#																																												 #
#																																												 #
##################################################################################################################################################################################


# 1) Indexing genome and aligning reads to reference with bowtie2

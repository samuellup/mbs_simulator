#!/bin/bash
# This workflow (1) simulates read files from different bulked segregant mapping-by-sequencing experimental setups, (2) aligns the reads to a reference genome,
# (3) generates .vcf and .va files and processes them into relevant data


# command structure: 	./main_workflow.sh project_name in_fasta nbr_mutations rec_freq_distr mut_pos nbr_rec_chrs read_depth meta_folder map_pop
# test command: 		./main_workflow.sh testing check.1.fa 40 0,24-1,42-2,25-3,6-4,1-5,2 1,10000 200 20 META M2

# _________________________________________________________________Some initial arrangements______________________________________________________________________________________
timestamp=$(date "+%F-%T")
export location="$PWD" 			#Save path to bowtie2-build and bowtie2 in variable BT2

# Get command arguments and assign them to variables
project_name=$1
in_fasta=$2
nbr_mutations=$3
rec_freq_distr=$4
mut_pos=$5
nbr_rec_chrs=$6
read_depth=$7
meta_folder=$8/META
map_pop=$9


causal_mut=$(echo $mut_pos | cut -d'-' -f 1) 


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

my_meta_log=$meta_folder/meta.log
my_meta_info=$meta_folder/meta_info.txt

echo $project_name >> $my_meta_log

sim_mut_output_folder=$f1/sim_data/sim_mut_output/mutant_strain
sim_recsel_output_folder=$f1/sim_data/sim_recsel_output
sim_seq_output_folder=$f1/sim_data/sim_seq_output/sample
sim_mut_output_folder_1=$f1/sim_data/sim_mut_output/mutant_strain_1
sim_mut_output_folder_2=$f1/sim_data/sim_mut_output/mutant_strain_2


##################################################################################################################################################################################
#																																												 #
#																																												 #
#																		SIMULATOR MODULE																						 #
#																																												 #
#																																												 #
##################################################################################################################################################################################


if [ $map_pop == 'F2' ]; then

	# 1) Simulation of mutagenesis with sim_mut.py
	{
		python2 sim_scripts/sim-mut.py -nbr $nbr_mutations -mod e -con $meta_folder/mutated_genome/mutated_genome.fa -causal_mut $causal_mut -out $sim_mut_output_folder 2>> $my_log_file

	} || {
		echo $(date "+%F > %T")": Simulation of mutagenesis failed. Quit." >> $my_log_file
		exit_code=1; echo $exit_code; exit
	}
	echo $(date "+%F > %T")": Simulation of mutagenesis completed." >> $my_log_file


	# 2) Simulating recombination with sim_rec.py 
	parmut_sample=$sim_mut_output_folder/mutated_genome/mutated_genome.fa
	parpol_sample=$meta_folder/mutated_genome/mutated_genome.fa
	mkdir $sim_recsel_output_folder
	{
		python2 sim_scripts/sim-recsel.py -outdir $sim_recsel_output_folder -rec_freq_distr $rec_freq_distr -parmut $parmut_sample -parpol $parpol_sample -mutpos $mut_pos -smod r -nrec $nbr_rec_chrs 2>> $my_log_file 

	} || {
		echo $(date "+%F > %T")": Simulation of recombination and phenotype selection failed. Quit." >> $my_log_file
		exit_code=1; echo $exit_code; exit
	}
	echo $(date "+%F > %T")": Simulation of recombination and phenotype selection completed." >> $my_log_file

fi


if [ $map_pop == 'M2' ]; then


	# We have to define the number of mutations for each chromatide, aproximately 1/2 of the total.


	nbr_mutations_1=`python2 an_scripts/rand.py -nbr_mutations $nbr_mutations 2>> $my_log_file`


	nbr_mutations_2=$(( nbr_mutations-nbr_mutations_1 ))
	
	echo $nbr_mutations_1 >> $my_log_file
	echo $nbr_mutations_2 >> $my_log_file


	# 1) Simulation of mutagenesis of each chromosome with sim_mut.py
	{
		python2 sim_scripts/sim-mut.py -nbr $nbr_mutations_1 -mod e -con $meta_folder/mutated_genome/mutated_genome.fa -causal_mut $causal_mut -out $sim_mut_output_folder_1 2>> $my_log_file
		python2 sim_scripts/sim-mut.py -nbr $nbr_mutations_2 -mod e -con $meta_folder/mutated_genome/mutated_genome.fa -out $sim_mut_output_folder_2 2>> $my_log_file

	} || {
		echo $(date "+%F > %T")": Simulation of mutagenesis failed. Quit." >> $my_log_file
		exit_code=1; echo $exit_code; exit
	}
	echo $(date "+%F > %T")": Simulation of mutagenesis completed." >> $my_log_file



	# 2) Simulating recombination with sim_rec.py 
	parmut_sample=$sim_mut_output_folder_1/mutated_genome/mutated_genome.fa
	parpol_sample=$sim_mut_output_folder_2/mutated_genome/mutated_genome.fa
	mkdir $sim_recsel_output_folder
	{
		python2 sim_scripts/sim-recsel.py -outdir $sim_recsel_output_folder -rec_freq_distr $rec_freq_distr -parmut $parmut_sample -parpol $parpol_sample -mutpos $mut_pos -smod r -nrec $nbr_rec_chrs 2>> $my_log_file 

	} || {
		echo $(date "+%F > %T")": Simulation of recombination and phenotype selection failed. Quit." >> $my_log_file
		exit_code=1; echo $exit_code; exit
	}
	echo $(date "+%F > %T")": Simulation of recombination and phenotype selection completed." >> $my_log_file


fi



# 3) Simulating HTS reads
lib_type=se 									#<------------- Comprobar y establecer parametros por defecto, establecer RD para la muestra control
read_length_mean=200
read_length_sd=40
fragment_length_mean=500
fragment_length_sd=100
basecalling_error_rate=1
gc_bias_strength=100

{
	python2 sim_scripts/sim-seq.py -input_folder $sim_recsel_output_folder -out $sim_seq_output_folder -mod $lib_type -rd $read_depth -rlm $read_length_mean -rls $read_length_sd -flm $fragment_length_mean -fls $fragment_length_sd -ber $basecalling_error_rate -gbs $gc_bias_strength 2>> $my_log_file

} || {
	echo $(date "+%F > %T")": Simulation of high-throughput sequencing failed. Quit." >> $my_log_file
	exit_code=1; echo $exit_code; exit
}
echo $(date "+%F > %T")": Simulation of high-throughput sequencing completed." >> $my_log_file



##################################################################################################################################################################################
#																																												 #
#																																												 #
#																		READ PROCESSING MODULE																					 #
#																																												 #
#																																												 #
##################################################################################################################################################################################


# 1) Indexing genome and aligning reads to reference with bowtie2

#Run bowtie2-build on genome sequence 
{
	$location/toolshed/bowtie2/bowtie2-build $f0/$in_fasta $f1/genome_index 1> $f2/bowtie2-build_std1.txt 2> $f2/bowtie2-build_std2.txt

} || {
	echo $(date "+%F > %T")': Bowtie2-build on genome sequence returned an error. See log files.' >> $my_log_file
	exit_code=1; echo $exit_code; exit
}
echo $(date "+%F > %T")': Bowtie2-build finished.' >> $my_log_file

#Run bowtie2 paired to align raw F2 reads to genome 						<------- Ajustar los parametros del alineador / variant calling / etc para el articulo
my_rd=$sim_seq_output_folder/se_reads.fq
{
	$location/toolshed/bowtie2/bowtie2 --very-sensitive  --mp 3,2 -X 1000  -x $f1/genome_index -U $my_rd -S $f1/alignment1.sam 2> $f2/bowtie2_problem-sample_std2.txt

} || {
	echo $(date "+%F > %T")': Bowtie2 returned an error during the aligment of reads. See log files.' >> $my_log_file
	exit_code=1; echo $exit_code; exit

}
echo $(date "+%F > %T")': Bowtie2 finished the alignment of reads to genome.' >> $my_log_file

# 2) Variant calling

# SAM-TO-BAM
{
	$location/toolshed/samtools1/samtools sort $f1/alignment1.sam > $f1/alignment1.bam 2> $f2/sam-to-bam_problem-sample_std2.txt
	#rm -rf $1/alignment1.sam

} || {
	echo 'Error transforming SAM to BAM.' >> $my_log_file
	exit_code=1; echo $exit_code; exit

}
echo $(date "+%F > %T")': SAM to BAM finished.' >> $my_log_file

# VC pipeline
{
	$location/toolshed/samtools1/samtools mpileup  -B -t DP,ADF,ADR -C50 -uf $f0/$in_fasta $f1/alignment1.bam 2> $f2/mpileup_problem-sample_std.txt | $location/toolshed/bcftools-1.3.1/bcftools call -mv -Ov > $f1/raw_variants.vcf 2> $f2/call_problem-sample_std.txt

} || {
	echo $(date "+%F > %T")': Error during variant-calling' >> $my_log_file
	exit_code=1; echo $exit_code; exit

}
echo $(date "+%F > %T")': Variant calling finished.' >> $my_log_file

#Groom vcf
{
	python2 $location/an_scripts/vcf-groomer.py -a $f1/raw_variants.vcf -b $f1/raw_variants.va  2>> $my_log_file

} || {
	echo $(date "+%F > %T")': Error during execution of vcf-groomer.py' >> $my_log_file
	exit_code=1; echo $exit_code; exit

}
echo $(date "+%F > %T")': VCF grooming finished.' >> $my_log_file

# Filter .va file
{
	python2 $location/an_scripts/variants-filter.py -a $f1/raw_variants.va -b $f1/filtered_variants.va -step 3 -fasta $f0/$in_fasta -dp_min 7 -qual_min 20  2>> $my_log_file

} || {
	echo 'Error during execution of variants-filter.py with F2 data.' >> $my_log_file
	exit_code=1; echo $exit_code; exit

}
echo $(date "+%F > %T")': VCF filtering step of sample data finished.' >> $my_log_file

# Substraction of background variants
my_operation_mode=A
{
	python2 $location/an_scripts/variants-operations.py -a $f1/filtered_variants.va -b $meta_folder/lab.va -c $f1/variants.va -mode $my_operation_mode -primary 1  2>> $my_log_file

} || {
	echo $(date "+%F > %T")': Error during execution of variants-operations.py .' >> $my_log_file
	exit_code=1; echo $exit_code; exit

}
echo $(date "+%F > %T")': VCF operations finished.' >> $my_log_file


##################################################################################################################################################################################
#																																												 #
#																																												 #
#																		OUTPUT GENERATION MODULE																				 #
#																																												 #
#																																												 #
##################################################################################################################################################################################

# Scatter plot
{
	python2 $location/an_scripts/plot.py -in_va $f1/variants.va -out $f3/ye.png  2>> $my_log_file

} || {
	echo $(date "+%F > %T")': Error during execution of plot.py' >> $my_log_file
	exit_code=1; echo $exit_code; exit

}
echo $(date "+%F > %T")': Data plotting finished.' >> $my_log_file

#xdg-open $f3/ye.png

# Candidates analysis
{
	python2 $location/an_scripts/cr-analysis.py -in_va $f1/variants.va -out $my_meta_info -mps $nbr_rec_chrs -rd $read_depth  2>> $my_log_file

} || {
	echo $(date "+%F > %T")': Error during execution of cr-analysis.py' >> $my_log_file
	exit_code=1; echo $exit_code; exit

}
echo $(date "+%F > %T")': Data analysis finished.' >> $my_log_file


# Intermediate files cleanup and re-organization
mv $f2/log.log $project_name
mv $f3/ye.png $project_name
mv $f1/variants.va $project_name
rm -rf $f1
rm -rf $f2
rm -rf $f3

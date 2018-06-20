# This runs main_workflow.sh to perform several simulations while variating the RD (Read Depth) and MPS (Mapping Mopulation Size) parameters. 

# command structure: 	./meta_sim.sh meta_name in_fasta map_pop
# test command: 		./meta_sim.sh testing check.1.fa M2



######################################################### STRATEGY #########################################################
#
# 1- Simulation of CHR4-Lab -> lab.va with variants of lab strain and lab.fa with the sequence
# 2- lab.fa is used as a base for the second mutagenesis within main_workflow.sh; -> variants.va, substraction of variants.va - lab.va 
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

my_meta_log=$meta_folder/meta.log; touch $my_meta_log
my_meta_info=$meta_folder/meta_info.txt; touch $my_meta_info
echo "#RD	MPS	CANDIDATES_95	SPAN_95	CANDIDATES_98	SPAN_98" >> $my_meta_info

#nbr_background_mutations=109															# <------------------------- Backcross
nbr_background_mutations=175000															# <------------------------- Outcross

rd_list=(30)																		# <------------------------- SET
mps_list=(40)																	# <------------------------- SET
#mps_list=(40 80 160 320)
#rd_list=(30 60 200)
export location="$PWD" 			#Save path to bowtie2-build and bowtie2 in variable BT2


# 2) Generation of lab strain files for subsequent analysis ____________________________________________________________________________________________________________________________________________________________________________________

# Simulation of mutagenesis with sim_mut.py
{
	python2 sim_scripts/sim-mut.py -nbr $nbr_background_mutations -mod d -con ./u_data/$in_fasta -out $meta_folder 2>> $my_meta_log

} || {
	echo $(date "+%F > %T")": Simulation of mutagenesis failed. Quit." >> $my_meta_log
	exit_code=1; echo $exit_code; exit
}
echo $(date "+%F > %T")": Simulation of mutagenesis completed." >> $my_meta_log


# Simulating HTS reads
lib_type=se 									#<------------- Comprobar y establecer parametros por defecto, establecer RD para la muestra control
read_length_mean=200
read_length_sd=40
fragment_length_mean=500
fragment_length_sd=100
basecalling_error_rate=1
gc_bias_strength=100
control_rd=120 									#<------------- SET

{
	python2 sim_scripts/sim-seq.py -input_folder $meta_folder/mutated_genome/ -out $meta_folder/seq_out -mod $lib_type -rd $control_rd -rlm $read_length_mean -rls $read_length_sd -flm $fragment_length_mean -fls $fragment_length_sd -ber $basecalling_error_rate -gbs $gc_bias_strength 2>> $my_meta_log

} || {
	echo $(date "+%F > %T")": Simulation of high-throughput sequencing failed. Quit." >> $my_meta_log
	exit_code=1; echo $exit_code; exit
}
echo $(date "+%F > %T")": Simulation of high-throughput sequencing completed." >> $my_meta_log



# bowtie2-build 
{
	$location/toolshed/bowtie2/bowtie2-build ./u_data/$in_fasta $meta_folder/genome_index 1> $meta_folder/bowtie2-build_std1.txt 2> $meta_folder/bowtie2-build_std2.txt

} || {
	echo $(date "+%F > %T")': Bowtie2-build on genome sequence returned an error. See log files.' >> $my_meta_log
	exit_code=1; echo $exit_code; exit
}
echo $(date "+%F > %T")': Bowtie2-build finished.' >> $my_meta_log

# bowtie2 												<------- Ajustar los parametros del alineador / variant calling / etc para el articulo
my_rd=$meta_folder/seq_out/se_reads.fq
{
	$location/toolshed/bowtie2/bowtie2 --very-sensitive  --mp 3,2 -X 1000  -x $meta_folder/genome_index -U $my_rd -S $meta_folder/alignment1.sam 2> $meta_folder/bowtie2_problem-sample_std2.txt

} || {
	echo $(date "+%F > %T")': Bowtie2 returned an error during the aligment of reads. See log files.' >> $my_meta_log
	exit_code=1; echo $exit_code; exit

}
echo $(date "+%F > %T")': Bowtie2 finished the alignment of reads to genome.' >> $my_meta_log


# SAM-TO-BAM
{
	$location/toolshed/samtools1/samtools sort $meta_folder/alignment1.sam > $meta_folder/alignment1.bam 2> $meta_folder/sam-to-bam_problem-sample_std2.txt
	rm -rf $meta_folder/alignment1.sam

} || {
	echo 'Error transforming SAM to BAM.' >> $my_meta_log
	exit_code=1; echo $exit_code; exit

}
echo $(date "+%F > %T")': SAM to BAM finished.' >> $my_meta_log


# VC pipeline
{
	$location/toolshed/samtools1/samtools mpileup  -B -t DP,ADF,ADR -C50 -uf ./u_data/$in_fasta $meta_folder/alignment1.bam 2> $meta_folder/mpileup_problem-sample_std.txt | $location/toolshed/bcftools-1.3.1/bcftools call -mv -Ov > $meta_folder/raw_variants.vcf 2> $meta_folder/call_problem-sample_std.txt

} || {
	echo $(date "+%F > %T")': Error during variant-calling' >> $my_meta_log
	exit_code=1; echo $exit_code; exit

}
echo $(date "+%F > %T")': Variant calling finished.' >> $my_meta_log

# Groom vcf
{
	python2 $location/an_scripts/vcf-groomer.py -a $meta_folder/raw_variants.vcf -b $meta_folder/lab.va  2>> $my_meta_log

} || {
	echo $(date "+%F > %T")': Error during execution of vcf-groomer.py' >> $my_meta_log
	exit_code=1; echo $exit_code; exit

}
echo $(date "+%F > %T")': VCF grooming finished.' >> $my_meta_log

# Cleanup
rm -rf $meta_folder/*.bt2 $meta_folder/*.txt $meta_folder/*.vcf $meta_folder/*.bam $meta_folder/seq_out



# 3) Running the simulations ____________________________________________________________________________________________________________________________________________________________________________________________________________________

rec_freq_distr='0,24-1,43-2,25-3,6-4,1-5,1'							     		# <------------------------- SET
nbr_mutations=232 															    # <------------------------- SET
mut_pos='1,5845220'
#mut_pos='1,500000'

for n in `seq 5`; do 							 								# <------------------------- SET Number of replicates
	for i in ${rd_list[@]}; do
		for j in ${mps_list[@]}; do
				rd=$i
		        mps=$j

		        project_name=$meta_name/$rd"_"$mps"_"$n

		        ./main_workflow.sh $project_name $in_fasta $nbr_mutations $rec_freq_distr $mut_pos $mps $rd $meta_name $map_pop

		        echo $project_name ' done!'
		done
	done
done

# 4) Analizing the obtained data _______________________________________________________________________________________________________________________________________________________________________________________________________________
{
	python2 ./an_scripts/meta-analysis.py -meta_in $my_meta_info -out $meta_folder/averaged_data_95.txt -step 95 2>> $my_meta_log
	python2 ./an_scripts/meta-analysis.py -meta_in $my_meta_info -out $meta_folder/averaged_data_98.txt -step 98 2>> $my_meta_log

} || {
	echo $(date "+%F > %T")': Error during execution of meta-analysis.py' >> $my_meta_log
	exit_code=1; echo $exit_code; exit

}
echo $(date "+%F > %T")': meta-analysis.py finished.' >> $my_meta_log


# This runs main_workflow.sh to perform several simulations while variating the RD (Read Depth) and MPS (Mapping Mopulation Size) parameters. 

# command structure: 	./meta_sim.sh meta_name in_fasta 
# test command: 		./meta_sim.sh testing check.1.fa



# 1) Initial parameters and accomodations 
timestamp=$(date "+%F-%T")

meta_name=u_projects/$timestamp"_"$1
in_fasta=$2
mkdir $meta_name
mkdir $meta_name/META
meta_folder=$meta_name/META

my_meta_log=$meta_folder/meta.log; touch $my_meta_log
my_meta_info=$meta_folder/meta_info.txt; touch $my_meta_info
echo "#RD	MPS	CANDIDATES	SPAN" >> $my_meta_info

rd_list=(2 3)
mps_list=(2 10)


# 2) Running the simulations
rec_freq_distr='0,24-1,42-2,25-3,6-4,1-5,2'		  # <------------------------- SET
nbr_mutations=20 								  # <------------------------- SET
mut_pos='1,10000'

for n in `seq 10`; do 							  # Number of replicates
	for i in ${rd_list[@]}; do
		for j in ${mps_list[@]}; do
				rd=$i
		        mps=$j

		        project_name=$meta_name/$rd"_"$mps"_"$n

		        ./main_workflow.sh $project_name $in_fasta $nbr_mutations $rec_freq_distr $mut_pos $mps $rd $meta_name

		        echo $project_name ' done!'

		done
	done
done

# 3) Analizing the obtained data
{
	python2 ./an_scripts/meta-analysis.py -meta_in $my_meta_info -out $meta_folder/averaged_data.txt  2>> $my_meta_log

} || {
	echo $(date "+%F > %T")': Error during execution of meta-analysis.py' >> $my_meta_log
	exit_code=1; echo $exit_code; exit

}
echo $(date "+%F > %T")': meta-analysis.py finished.' >> $my_meta_log

# 4) Creating a heatmap

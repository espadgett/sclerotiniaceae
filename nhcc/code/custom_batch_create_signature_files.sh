#!/bin/bash
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=5G
#SBATCH -t 00:15:00
#SBATCH --output=%x.e%j
#SBATCH --account=investigator      ### DIFFERENCE: added account line (not present in base script)



#----------------EDIT CELL INFO------------------#
batch=1         # number  ex.  [4]           ### DIFFERENCE: filled value (1)

cell1=Smin            # cell_type1_dataset_name  ex. [Astrocyte_spine]       ### DIFFERENCE: filled with “Sminor”
cell2=Smin          # cell_type2_dataset_name  ex. [Astrocyte_cerebellum]   ### DIFFERENCE: filled with “Sminor”
cell3=Smin          # cell_type3_dataset_name  ex. [Aorta_Leung]            ### DIFFERENCE: filled with “Sminor”
#------------------------------------------------#

#---------------EDIT ANALYSIS INFO---------------#
path=/home/user/scratch/streptotinia_caulophylli/hic_asm/interchr_interact/sminor_rerun            
coolerpath=/home/user/scratch/streptotinia_caulophylli/hic_asm/interchr_interact/sminor_rerun #where cooler outputs are

analysis=trans        # either [cis] or [trans]     ### DIFFERENCE: set to “trans”
trans=1vsAll           # if trans selected above: either [1vsAll] or [pairwise]  
resolution=10000         # resolution: [1000000], [50000], etc.                  
res=10KB             # resolution in human readable format; [1MB], [50KB], etc.  

module load StdEnv/2023       ### DIFFERENCE: added module loading block
module load bedtools/2.31.0   ### DIFFERENCE
module load r/4.4.0           ### DIFFERENCE
module load python/3.11.5     ### DIFFERENCE
#------------------------------------------------#






#-----------------DONT EDIT----------------------#
scr=/home/user/scratch/streptotinia_caulophylli/hic_asm/interchr_interact/LWPR/scripts_Signature #where Signature scripts are
path=$path/smin_results
mkdir $path

# LIST OF CELLS
echo $cell1 > $path/batch$batch\_cells.txt
echo $cell2 >> $path/batch$batch\_cells.txt
echo $cell3 >> $path/batch$batch\_cells.txt

mkdir -p $path/schedulers/$cell1.$cell2.$cell3
mkdir -p $path/output/$cell1.$cell2.$cell3
echo "Pathway to configuration and scheduler = $path/schedulers/$cell1.$cell2.$cell3"



# CONFIGURATION FILE

true > $path/schedulers/$cell1.$cell2.$cell3/$cell1.$cell2.$cell3.signature.$analysis$trans.configuration.txt

echo type.class=$analysis$trans >> $path/schedulers/$cell1.$cell2.$cell3/$cell1.$cell2.$cell3.signature.$analysis$trans.configuration.txt
echo cell=$path/batch$batch\_cells.txt >> $path/schedulers/$cell1.$cell2.$cell3/$cell1.$cell2.$cell3.signature.$analysis$trans.configuration.txt
echo res=10000 >> $path/schedulers/$cell1.$cell2.$cell3/$cell1.$cell2.$cell3.signature.$analysis$trans.configuration.txt     
echo name=$cell1.$cell2.$cell3.$analysis$trans.$res >> $path/schedulers/$cell1.$cell2.$cell3/$cell1.$cell2.$cell3.signature.$analysis$trans.configuration.txt
echo out=$path/output/$cell1.$cell2.$cell3 >> $path/schedulers/$cell1.$cell2.$cell3/$cell1.$cell2.$cell3.signature.$analysis$trans.configuration.txt
echo pval=0.05 >> $path/schedulers/$cell1.$cell2.$cell3/$cell1.$cell2.$cell3.signature.$analysis$trans.configuration.txt
echo dir=$coolerpath >> $path/schedulers/$cell1.$cell2.$cell3/$cell1.$cell2.$cell3.signature.$analysis$trans.configuration.txt
echo in=$scr >> $path/schedulers/$cell1.$cell2.$cell3/$cell1.$cell2.$cell3.signature.$analysis$trans.configuration.txt

echo "CONFIGURATION FILE DONE"

# SCHEDULER FILE

true > $path/schedulers/$cell1.$cell2.$cell3/$cell1.$cell2.$cell3.signature.$analysis$trans.scheduler.sh

echo '#!/bin/bash' >> $path/schedulers/$cell1.$cell2.$cell3/$cell1.$cell2.$cell3.signature.$analysis$trans.scheduler.sh
echo "#SBATCH -N 1" >> $path/schedulers/$cell1.$cell2.$cell3/$cell1.$cell2.$cell3.signature.$analysis$trans.scheduler.sh
echo "#SBATCH -c 20" >> $path/schedulers/$cell1.$cell2.$cell3/$cell1.$cell2.$cell3.signature.$analysis$trans.scheduler.sh
echo "#SBATCH --mem=100G" >> $path/schedulers/$cell1.$cell2.$cell3/$cell1.$cell2.$cell3.signature.$analysis$trans.scheduler.sh
echo "#SBATCH --tmp=300G" >> $path/schedulers/$cell1.$cell2.$cell3/$cell1.$cell2.$cell3.signature.$analysis$trans.scheduler.sh
echo "#SBATCH --account=def-mtodesco" >> $path/schedulers/$cell1.$cell2.$cell3/$cell1.$cell2.$cell3.signature.$analysis$trans.scheduler.sh    ### DIFFERENCE: replaced “-p general” with “--account=def-mtodesco”
echo "#SBATCH -t 7-00:00:00" >> $path/schedulers/$cell1.$cell2.$cell3/$cell1.$cell2.$cell3.signature.$analysis$trans.scheduler.sh
echo "#SBATCH --output=%x.e%j" >> $path/schedulers/$cell1.$cell2.$cell3/$cell1.$cell2.$cell3.signature.$analysis$trans.scheduler.sh

echo -e "\n" >> $path/schedulers/$cell1.$cell2.$cell3/$cell1.$cell2.$cell3.signature.$analysis$trans.scheduler.sh

echo "module load StdEnv/2023" >> $path/schedulers/$cell1.$cell2.$cell3/$cell1.$cell2.$cell3.signature.$analysis$trans.scheduler.sh     ### DIFFERENCE: changed module section entirely
echo "module load bedtools/2.31.0" >> $path/schedulers/$cell1.$cell2.$cell3/$cell1.$cell2.$cell3.signature.$analysis$trans.scheduler.sh   
echo "module load r/4.4.0" >> $path/schedulers/$cell1.$cell2.$cell3/$cell1.$cell2.$cell3.signature.$analysis$trans.scheduler.sh            
echo "module load python/3.11.5" >> $path/schedulers/$cell1.$cell2.$cell3/$cell1.$cell2.$cell3.signature.$analysis$trans.scheduler.sh     

echo "" >> $path/schedulers/$cell1.$cell2.$cell3/$cell1.$cell2.$cell3.signature.$analysis$trans.scheduler.sh
echo "toSign=$scr" >> $path/schedulers/$cell1.$cell2.$cell3/$cell1.$cell2.$cell3.signature.$analysis$trans.scheduler.sh
echo "configuration=$path/schedulers/$cell1.$cell2.$cell3/$cell1.$cell2.$cell3.signature.$analysis$trans.configuration.txt" >> $path/schedulers/$cell1.$cell2.$cell3/$cell1.$cell2.$cell3.signature.$analysis$trans.scheduler.sh
echo "" >> $path/schedulers/$cell1.$cell2.$cell3/$cell1.$cell2.$cell3.signature.$analysis$trans.scheduler.sh
echo '$toSign/run_Signature.sh $(cat $configuration)' >> $path/schedulers/$cell1.$cell2.$cell3/$cell1.$cell2.$cell3.signature.$analysis$trans.scheduler.sh

echo "SCHEDULER FILE DONE"

echo "SUBMITTING SCHEDULER ..."
cd $path/schedulers/$cell1.$cell2.$cell3
sbatch --account=def-mtodesco $path/schedulers/$cell1.$cell2.$cell3/$cell1.$cell2.$cell3.signature.$analysis$trans.scheduler.sh   ### DIFFERENCE: added “--account=def-mtodesco”
#-------------------------------------#
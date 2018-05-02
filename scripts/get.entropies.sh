#!/bin/bash
#Mebs v2 

#Export mebs data
scripts=/home/valdeanda/src/metagenome_Pfam_score/scripts/
datadir=home/valdeanda/src/metagenome_Pfam_score/cycles/
outputdir=$pathway_entropy_dir
#tabdir=/home/valdeanda/MEBSv2/GenF_Pfam
pathway=(carbon nitrogen iron sulfur oxygen)

#Compute entropies only for oxygen cycle in Gen and GenF  
if [ $# -eq 0 ]
  then
    echo "# Need folder containing tab files against Pfam database 
            using hmmserach\n";
    exit;
fi

inputdir=$1

echo "# parameters:"
echo "# datadir=$datadir"
echo "# cycle = $pathway"
echo "# inputdir =$inputdir"
echo "Selected pathways are ${#pathway[@]}"

for path in "${pathway[@]}"; do 

  for i in $inputdir/*.tab; do \

perl $scripts/entropy.pl $i \
$datadir/$pathway/$pathway.txt > $i.$pathway.csv 
  done 


done 



#python3 $scripts/extract_entropies.py $output_dir/

#echo "The final output of $pathway entropies is $output_dir_entropies.tab" 
echo "Done" 



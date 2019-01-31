#!/bin/bash -x

# Usage: bash sparcc.sh

### One at a time ###
#source activate sparcc

#cd /Users/abby/Documents/Projects/dietstudy/analysis/tax_network_cor/sparcc

#SparCC.py ../tax_username/MCTs10_tax.txt -i 20 -c corr_sparcc_MCTs10.out
#SparCC.py ../tax_username/MCTs10_tax.txt -i 20 -c corr_pearson_MCTs10.out -a pearson
#SparCC.py ../tax_username/MCTs10_tax.txt -i 20 -c corr_spearman_MCTs10.out -a spearman

# Pseudo p-value calculation
#MakeBootstraps.py ../tax_username/MCTs10_tax.txt -n 100 -t permutation_#.txt -p pvals/

#for i in {0..99};
#do SparCC.py pvals/permutation_$i.txt -i 100 --cor_file=pvals/perm_cor_$i.txt &
#done


#PseudoPvals.py corr_sparcc_MCTs10.out pvals/perm_cor_#.txt 100 -o pvals/pvals.one_sided.txt -t one_sided
#PseudoPvals.py corr_sparcc_MCTs10.out pvals/perm_cor_#.txt 100 -o pvals/pvals.two_sided.txt -t two_sided



####### scale up! ###########
source activate sparcc

cd /Users/abby/Documents/Projects/dietstudy/analysis/tax_network_cor/tax_username

for i in *_tax.txt; 
do SparCC.py $i -i 20 -c ../sparcc_results/corr_sparcc_${i//_tax.txt}_out.txt & 
done

for i in *_tax.txt; 
do MakeBootstraps.py $i -n 100 -t permutation_#.txt -p ../sparcc_results/pvals_${i//_tax.txt}/ & 
done

cd ../sparcc_results/

# the slow and steady version - at least this doesn't break my computer...
# not great, but can run in the background on my computer
for j in pvals_*; 
do for k in {0..99}; 
do SparCC.py $j/permutation_$k.txt -i 100 --cor_file=$j/perm_cor_$k.txt; 
done; 
done

## here is the faster version that probably needs more computing power than I have on my laptop to run.
#for j in pvals_*;  
#do echo "for k in {0..99}; do SparCC.py $j/permutation_\$k.txt -i 100 --cor_file=$j/perm_cor_\$k.txt & done" > job-$j.sh; 
#done

#chmod +x *.sh

#for j in pvals_*; 
#./job-$j.sh;
#done

# final step!
for j in pvals_*;
do PseudoPvals.py corr_sparcc_${j//pvals_}_out.txt $j/perm_cor_#.txt 100 -o pvals_${j//pvals_}_two_sided.txt -t two_sided
done

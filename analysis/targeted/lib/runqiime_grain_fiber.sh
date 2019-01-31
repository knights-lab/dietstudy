# Usage bash div.preprocess.sh

# To run on my own computer
source activate qiime

# change directory
cd ../../../data/processed_food/

# convert food otu tables to biom
biom convert -i grains_fiber.txt -o grains_fiber.biom --to-json --table-type "OTU table" --process-obs-metadata=taxonomy

# calculate beta diversity tables for food otus
beta_diversity.py -i grains_fiber.biom -o grains_beta -t ../../raw/diet.tree.txt

# calculate pcoa for food
principal_coordinates.py -i grains_beta/weighted_unifrac_grains_fiber.txt -o grains_beta/w_grains_fiber_pcoa.txt
principal_coordinates.py -i grains_beta/unweighted_unifrac_grains_fiber.txt -o grains_beta/unw_grains_fiber_pcoa.txt


# change directory
cd ../processed_UN_tax/


# convert taxa counts table to biom
biom convert -i UN_taxonomy_counts_s.txt -o UN_taxonomy_counts_s.biom --to-json --table-type "OTU table"
biom convert -i UN_taxonomy_clr_s.txt -o UN_taxonomy_clr_s.biom --to-json --table-type "OTU table"

# calculate beta diversity tables for taxa
beta_diversity.py -i UN_taxonomy_counts_s.biom -o UN_tax_beta -m chisq,bray_curtis
beta_diversity.py -i UN_taxonomy_clr_s.biom -o UN_tax_beta -m euclidean

# calculate pcoa for taxa
principal_coordinates.py -i UN_tax_beta/chisq_UN_taxonomy_counts_s.txt -o UN_tax_beta/chi_tax_counts_pcoa.txt
principal_coordinates.py -i UN_tax_beta/bray_curtis_UN_taxonomy_counts_s.txt -o UN_tax_beta/bray_tax_counts_pcoa.txt
principal_coordinates.py -i UN_tax_beta/euclidean_UN_taxonomy_clr_s.txt -o UN_tax_beta/euclidean_tax_clr_pcoa.txt

# change directory
cd ../../analysis/targeted/grains/

# run procrustes on different combos
#transform_coordinate_matrices.py -i ../../../data/processed_UN_tax/UN_tax_beta/bray_tax_counts_pcoa.txt,../../../data/processed_food/grains_beta/w_grains_fiber_pcoa.txt  -o grain_fiber_procrustes_results/ -r 999;
#transform_coordinate_matrices.py -i ../../../data/processed_UN_tax/UN_tax_beta/chi_tax_counts_pcoa.txt,../../../data/processed_food/grains_beta/unw_grains_fiber_pcoa.txt  -o grain_fiber_procrustes_results/ -r 999;
#transform_coordinate_matrices.py -i ../../../data/processed_UN_tax/UN_tax_beta/euclidean_tax_clr_pcoa.txt,../../../data/processed_food/grains_beta/w_grains_fiber_pcoa.txt  -o grain_fiber_procrustes_results/ -r 999;
transform_coordinate_matrices.py -i ../../../data/processed_UN_tax/UN_tax_beta/euclidean_tax_clr_pcoa.txt,../../../data/processed_food/grains_beta/unw_grains_fiber_pcoa.txt  -o grain_fiber_procrustes_results/ -r 999;

cd grain_fiber_procrustes_results/

# plot 
make_emperor.py -c -i . -m ../../data/map_smry.txt		 -o plots/ 
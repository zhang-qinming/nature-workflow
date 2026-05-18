GPU=0
TRAIT=$1

CUDA_VISIBLE_DEVICES=$GPU python genebayes.lof.gene_level.py --response traits/Backman_2021_${TRAIT}.summary_statistics.csv --features feature_selection/features.Backman_2021_${TRAIT}.tsv --train_genes all.txt --val_genes all.txt --out Backman_2021_${TRAIT} --integration_lb -3.99 --integration_ub 3.99 --batch_size 100 --n_integration_pts 50000 --lr 0.3 --total_iterations 1000 --n_trees_per_iteration 2 --reg_alpha 2 --reg_lambda 2 --subsample 0.8 --min_child_weight 3 --max_depth 3

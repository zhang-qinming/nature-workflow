GPU=0
source cuda10.2

TRAIT=$1

CUDA_VISIBLE_DEVICES=$GPU  python prepare_features.py Backman_2021_${TRAIT}
CUDA_VISIBLE_DEVICES=$GPU  python prepare_features_signed.py Backman_2021_${TRAIT}

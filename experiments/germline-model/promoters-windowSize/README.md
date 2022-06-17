## Remove test promoters from neutral regions used to estimate substitution probabilities

- [x] Run `experiments/germline-model/promoters-windowSize/train_test_split.ipynb` to split promoters into train and test subsets
- [ ] Run `download-process-data/compute-neutral-regions.sh` to exclude test promoters (`download-process-data/promoters/promoters.grch38.test.csv`) from set of neutral regions (`dist/neutral-regions-germline-grch38-exclude-test-promoters.bed`) on which substitution probabilities, etc, are to be estimated 

## Retrain constraint model not just on new neutral regions, but also using various window sizes 

- [ ] Use line 82 in `train-germline-model.sh` as a starting point in a new batch script (call it `train-germline-model-production.exclude-test-promoters.sh`) that trains three models (windowSize = 101, 501, 1001) with a new `work` directory, a new `neutral_regions` file (`dist/neutral-regions-germline-grch38-exclude-test-promoters.bed`), producing three new `model` files (e.g., `dist/model-germline-grch38-exclude-test-promoters.windowSize-101.json`). Also see `train-germline-model-production.sh`. Also see `train/germline-model/train-map-reduce`. 
- [ ] copy `experiments/germline-model/promoters/promoters-compute-zscores.ipynb` to `experiments/germline-model/promoters-windowSize/compute-zscores-on-test-promoters.ipynb`. Then do `model = read_model('/scratch/ucgd/lustre-work/quinlan/u6018199/constraint-tools/dist/model-germline-grch38-exclude-test-promoters.windowSize-101.json')`, and do `promoters_filename = f'{CONSTRAINT_TOOLS}/download-process-data/promoters/promoters.grch38.test.csv’` and then use that model to make predictions on those promoters

 
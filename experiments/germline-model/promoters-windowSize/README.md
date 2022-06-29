## Remove test promoters from neutral regions used to estimate substitution probabilities

- [x] Run scripts in `experiments/germline-model/promoters-windowSize/` to split promoters into train and test subsets
- [x] Run `download-process-data/compute-neutral-regions.sh` to exclude test promoters (`download-process-data/promoters/promoters.grch38.test.sorted.bed.gz`) from set of neutral regions (`dist/neutral-regions-germline-grch38-exclude-test-promoters.bed`) on which substitution probabilities, etc, are to be estimated 

## Retrain constraint model not just on new neutral regions, but also using various window sizes 

- [x] Use line 82 in `train-germline-model.sh` as a starting point in a new batch script (call it `train-germline-model.exclude-test-promoters.sh`) that trains multiple models (windowSize = 101, 501, 1001, etc) with a new `work` directory, a new `neutral_regions` file (`dist/neutral-regions-germline-grch38-exclude-test-promoters.bed.gz`), producing multiple new `model` files (e.g., `dist/model-germline-grch38-exclude-test-promoters.windowSize-101.json`). 

## Re-do sanity checks

- [x] re-run `experiments/germline-model/sanity-check-*.ipynb` 

## Recompute promoter z-scores and make predictions from them 

- [x] copy `experiments/germline-model/promoters/promoters-compute-zscores.ipynb` to `experiments/germline-model/promoters-windowSize/compute-zscores-on-test-promoters.ipynb`. Then do `model = read_model('/scratch/ucgd/lustre-work/quinlan/u6018199/constraint-tools/dist/model-germline-grch38-exclude-test-promoters.windowSize-101.json')`, and do `promoters_filename = f'{CONSTRAINT_TOOLS}/download-process-data/promoters/promoters.grch38.test.csv’` and then use that model to make predictions on those promoters. 

- [x] use `experiments/germline-model/promoters/promoters-batch.sh` and `experiments/germline-model/promoters/promoters-execute-notebook.sh` as templates to run each model (one for each windowSize-windowStride combination) on the test promoters. 

- [x] use `experiments/germline-model/promoters/promoter-constraint-does-not-predict-gene-constraint.ipynb` as a template to hopefully show that promoter constraint does predict gene constraint, when window-size is large enough 

- [ ] use `experiments/germline-model/promoters/promoter-constraint-does-not-predict-gene-expression-variance.ipynb` as a template to hopefully show that promoter constraint does predict gene expression variance, when window-size is large enough 




 
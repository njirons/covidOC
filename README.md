# covidOC

Code to reproduce the results in
> Irons NJ, Raftery AE (2024). "US COVID-19 school closure was not cost-effective, but other measures were."
> [[arxiv]](https://arxiv.org/abs/2411.12016)

## Workflow

The scripts in the main directory can be run in the following sequence:

+ `10_data-cleaning.ipynb`: Jupyter notebook to download and clean the COVID case, death, and non-pharmaceutical intervention (NPI) data.
+ `20_set-params.R`: load the data into R as a JSON file, set epidemiological parameters for the model, and save the resulting dataframes to the `data` folder.
+ `30_cost-calc.R`: R script to calculate parameters for the NPI cost function.
+ `40_seir.sbatch`: Slurm batch script to fit the Bayesian SEIRD model across all US states in parallel.
    - `41_seir-local.R`: Example R script to fit the Bayesian SEIRD model for a single state locally.
+ `50_seir-output.R`: R script to process the output of the SEIRD model and feed into the NPI Bayesian hierarchical (regression) model (BHM).
+ `60_npi-bhm.sbatch`: Slurm batch script to fit the NPI BHM across SEIRD posterior trajectories in parallel.
    - `61_npi-bhm-local.R`: Example R script to fit the NPI BHM for a single SEIRD posterior trajectory locally.
+ `70_npi-bhm-output.R`: R script to process the output of the NPI BHM and feed into the policy optimization protocol.
+ `80_oc.sbatch`: Slurm batch script to solve for optimal NPI policies across US states and cost function and NPI model specifications in parallel.
    - `81_oc-local.R`: Example R script to run policy optimization for a single (US state, cost function, NPI model) specification.
+ `90_oc-plots.R`: R script to generate NPI policy plots in the paper.
+ `90_posterior-plots.R`: R script to generate parameter posterior plots in the paper.
+ `90_us-map-plots.R`: R script to generate US map plot in the paper.

Any supplementary scripts or data (or code needed to generate them) can be found in the `scripts` and `data` folders, respectively.

## Cite

If you make use of this code in your own work, please cite our paper:
```
@article{irons2024us,
  title={US COVID-19 school closure was not cost-effective, but other measures were},
  author={Irons, Nicholas J and Raftery, Adrian E},
  journal={arXiv preprint arXiv:2411.12016},
  year={2024}
}
```

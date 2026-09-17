# SBRT and Sublobar Resection for Early-Stage NSCLC

Code to generate the results in the paper comparing SBRT and sublobar resection using the SEER Medicare dataset. 

- `preprocess.2020.final.R`: Preprocessing code
- `load.data.R`: function for loading the data obtained from preprocess above
- `codes.R`: Billing codes for variable definitions
- `proximal.nc.2.mc.R`: Analysis code
- `two.step.variable.selection.R`: Wrapper for the pci2s functions
- `survival.curves.R`: Generate the survival/CIF curves
- `utilities.R`: Helper functions
- `file.paths.R`: Paths to raw data for the preprocessing code to work
- `martingale.residual.print.R`: Generate martingale residual plots for the two stage model
- `sensitivity.analysis.R` and `sensitivity.analysis.fig.R`: Generate sensitivity analysis in which each Z variable is excluded 


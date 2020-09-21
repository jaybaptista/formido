# Formido - Stellar Membership Probability Package

`formido` (another word for Deimos, namesake of the Keck Telescope instrument) is a package used to streamline membership probability determinations of a stars within galaxies and clusters. The current methods available can calculate the membership probability based on isochrone distance and distance relative to the half-light radius. These methods were obtained from a paper published by [https://arxiv.org/abs/1910.12879](Collins et al.). Eventually we hope to add another method that uses MCMC sampling on velocity dispersions to determine probabilities based on target and contaminant velocity dispersions.

Currently, `formido` only requires two dependencies: `astropy` and `numpy`.

The data provided in the examples folder is provided by Marla Geha at the Geha Lab at Yale University.

Currently we are trying to research a good probability cut-off that can best select populations in ultra-faint dwarf galaxies amidst contamination from the Milky Way.


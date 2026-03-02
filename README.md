<h2><b>Inferring patterns of purifying, positive and balancing selection in the coppery titi monkey, (<i>Plecturocebus cupreus</i>) utilizing a well-fit evolutionary baseline model</b></h2>
<h3>Vivak Soni, Cyril J. Versoza, John W. Terbot II, Gabriella Spatola, Karen L. Bales, Susanne P. Pfeifer, and Jeffrey D. Jensen</h3>

Repository contains data and scripts to run analyses.

<b>DFE_inference.tar.gz</b>: Contains vcf files for each empirical coppery titi monkey exon (empirical_data.tar.gz), and for simulated exons under various DFEs (simulations.tar.gz), as well as summary statistics calculated across simulated exons under various DFEs (bed files).

<b>genome_scans/power_analysis</b>: Contains selective sweep and balancing selection simulation vcfs (sweep_sims.tar.gz and bs_sims.tar.gz respectively), input files for genome scans (SF2_input.tar.gz and BM_input.tar.gz) and results of sweep scans (SF2_results.tar.gz and BM_results.tar.gz). Finally, mutation and recombination rate maps for simulations are contained in mu_maps.tar.gz and rr_maps.tar.gz respectively.

<b>genome_scans/empirical</b>: Contains input and results files of selective sweep and balancing selection inference (SF2_input.tar.gz, SF2_results.tar.gz, BM_input.tar.gz, BM_results.tar.gz). Also included are bed files of candidate loci, and candidate genes for each form of inference.

<b>scripts</b>: Contains SLiM script for running simulations for DFE inference (titi_DFE_highMut.slim), SLiM scripts for power analysis simulations (titi_demog_sweep_scaled.slim and titi_demog_bs_scaled.slim), python script for running coalescent simulations for generating null thresholds using MSPrime (titi_null_msprime.py), and finally a plotting script to generate main figures from the manuscript (titi_selection_plots.ipynb).

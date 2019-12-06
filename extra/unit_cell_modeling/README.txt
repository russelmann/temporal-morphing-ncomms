Requirements

* Matlab

* cbrewer - for plotting
	https://www.mathworks.com/matlabcentral/fileexchange/34087-cbrewer-colorbrewer-schemes-for-matlab


Optional

* TOMLAB /MAD - commercial algorithmic differentiation tool (processed data is included if no license)
	https://tomopt.com/tomlab/products/mad/
	
	
Description

Launch all_run.m to compute bracket model (required for plots below).

File <plot_dry.m>:         Plot force-displacement under normal conditions (Supplementary Figure 5)
File <plot_speed.m>:       Plot dependence of plasticity on deformation rate (Supplementary Figure 7)
File <plot_thick_rates.m>: Plot deformation rate depending on bracket thickness and load (Figure 2b)
File <plot_wet.m>:         Plot displacement over time in water (Supplementary Figure 6)
File <plot_wet_slice.m>:   Plot an example of displacement over time for a constant load (Figure 2c)

Folder <experiments>: Raw tensile testing data for unit cell specimens for the range of parameters,
length 5-9 mm, thickness 0.3-0.65 mm (shown in Supplementary Figure 5).

	File <dry_20181102.xls>: Force-displacement measurements under normal conditions.

	File <wet_1N_20181119.xls>: Force-displacement measurements in water with specimens
	subject to constant load of 1 N (shown in Supplementary Figure 6).

	File <wet_2N_20181119.xls>: Force-displacement measurements in water with specimens
	subject to constant load of 2 N (shown in Supplementary Figure 6).

	File <wet_3N_20181119.xls>: Force-displacement measurements in water with specimens
	subject to constant load of 3 N (shown in Supplementary Figure 6).

	File <wet_4Nand5N_20181119.xls>: Force-displacement measurements in water with
	specimens subject to constant loads of 4 N and 5 N (shown in Supplementary
	Figure 6).

	File <plastic_20181119.xls>: Plasticity tests. Specimens are compressed up to a
	fixed displacement and then slowly unloaded to study residual deformation.

	File <wet_speed_20181102.xls>: Tests for dependence of plastic deformation on
	deformation rates. Three measurements with different deformation rates for
	a specimen of lengths 8 mm and thickness 0.4 mm (shown in Supplementary
	Figure 7).
	
	File <dry_fit.m>: Fitted model for force-displacement under normal conditions (if no TOMLAB license)

	File <wet_fit.m>: Fitted model for force-displacement in water (if no TOMLAB license)

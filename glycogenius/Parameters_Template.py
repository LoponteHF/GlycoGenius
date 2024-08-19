[running_modes]
mode = analysis
;	Select in which mode you want to run Glycogenius:
;	- Use 'analysis' if you want to start a new analysis.
;	- Use 'library' to just generate a new library and
;	stop the execution afterwards.
;	- Use 'reanalysis' to reanalyze an existing .gg
;	analysis file.
;
use_multiple_CPU_cores = yes
number_cores = all
; 	Allows to use multiple cores for the processing 
;	of the data. If number_cores = all, uses 
;	total_cores-2 (ie. if you have a CPU with 20 
;	cores, it will use 18 cores).
;
working_directory =
; 	Directory to load and save files from script.
;
samples_directory =
; 	Directory where sample files are located. 
;	-> Only needed in analysis mode.
;
exported_library_name = 
;	Choose the name of the library file in library_name. 
;	If left blank, library will be exported with the file 
;	name 'date_time_glycans_library.ggl'.
;	-> Only needed in library mode.
;
file_for_reanalysis = 
; 	Indicate the path to the .gg analysis file. 
;	-> Only needed in reanalysis mode.
;

[library_building_modes]
mode = generate_library
;	Select in what mode you want to build a library:
;	- Use 'import_library' if you want to use an existing 
;	.ggl library file. Ignore common_library_building_
;	settings if importing a library.
;	- Use 'custom_library' if you want to build a library
;	from a list of glycans.
;	- Use 'generate_library' if you want to generate a 
;	new library using combinatorial analysis.
;
export_library = no
exported_library_name = 
; 	Exports the library you generated to the working
;	directory to use in future analysis without having
;	to build a new library. Also creates an excel file 
;	containing a human-readable version of the library 
;	generated and a file compatible with Skyline's 
;	transition list model. Choose the name of the 
;	library file in library_name. If left blank, 
;	library will be exported with the file name
;	'date_time_glycans_library.ggl'.
;
import_library_path = 
; 	Indicate the path to an existing .ggl library file
;	or to where and with which name you'd like to save
;	your exported library. If left blank, library will
;	be exported as 'date_time_glycans_library.ggl' at
;	the working directory folder.
;	-> Only needed in import_library mode.
;
custom_glycans_list = H3N2, H5N2, H5N4S2F1
;	Input the glycans directly, comma separated, or
;	indicate the path to a text file containing
;	the list (line or comma separated).
;	Monosaccharides accepted: Hexoses (H), 
;	HexNAc (N), Acetyl Sialic Acid (S or Am and E
;	if you have lactonized-ethyl sterified glycans),
;	Glycolyl Sialic Acid (G), Deoxyhexose (F). Case
;	sensitive.
;	-> Only needed in custom_library mode. 
;
total_monosaccharides = 5, 22
hexoses = 3, 10
hexosamines = 0, 0
hexnacs = 2, 8
xylose = 0, 0
sialic_acids = 0, 4
uronic_acids = 0, 0
fucoses = 0, 2
neu5ac = 0, 4
neu5gc = 0, 0
; 	Specify the minimum and maximum monosaccharides
;	amounts for generating a library.
;	-> Only needed in generate_library mode.
;

[common_library_building_settings]
force_nglycan = yes
; 	Used to force some monosaccharides compositions
;	associated with N-Glycans biologically known 
;	features. If not used, gives a much broader 
;	library that can be used for analysis of 
;	O-Glycans, for example.
;
max_adducts = H3
adducts_exclusion = 
; 	Indicates the desired adducts and their maximum
;	amount. H3Na1 means a maximum of 3 Hydrogens and
;	a maximum of 1 Sodium per adduct combination. 
;	Case sensitive. Doesn't work with complex adducts,
;	such as NH4. Set adducts in 'adducts_exclusion' to
;	avoid using specific adducts. Comma separated list.
;
max_charges = 3
; 	Limits the maximum amount of calculated charges
;	for each glycan. Set to a negative value if you
;	want to do negative mode analysis [EXPERIMENTAL].
;
reducing_end_tag = 133.0644
; 	If a reducing end tag is added to the glycans, 
; 	insert its added mass or molecular formula here.
;	If no reducing end tag is added to the glycans, 
;	set this value to 0. Procainamide: 219.1735 or
;	C13H21N3; Girard Reagent P: 133.0644 or C7H7N3 
;	(deprotonated, neutral). Can also be used with 
;	a peptide by inputing 'pep-' + the petide 
;	sequence (ie. 'pep-NK' for the peptide NK).
;
permethylated = no
; 	If the sample was permethylated, set this 
;	parameter to "yes". Doesn't take into account 
;	partial permethylations.
;
reduced = no
; 	If the sample doesn't have a tag and the glycans
;	had their reducing end reduced, set this to "yes".
;
aminated_ethyl_esterified = no
; 	Use if the sialic acids were derivitized with 
;	amination (alpha2,3) and ethyl esterification 
;	(alpha2,6). Lactonized Acetyl Sialic Acid will
;	be identified as 'Am' and Ethyl Esterified Acetyl
;	Sialic Acid will be identified as 'E'.
;
min_max_sulfation_per_glycan = 0, 0
;	Minimum and maximum amount of sulfation
;	substituents per glycan.
;
min_max_phosphorylation_per_glycan = 0, 0
;	Minimum and maximum amount of phosphorylation
;	substituents per glycan.
;
fast_iso = yes
; 	Allows you to calculate the isotopic distribution
;	of glycans based only on carbon isotopes (fast, 
;	innacurate) and corrected artificially or on all 
;	the atoms isotopes (VERY SLOW, very accurate). 
;	If not used, library building may take many hours
;	depending on size.
;
high_resolution_isotopic_dist = no
; 	If not used, doesn't clump isotope peaks 
;	together, meaning that you'll have more than one 
;	isotopic peak in a 1 Da interval. Only use this 
;	if you have fast_iso off. Useful when analyzing 
;	very high resolution data, such as data acquired 
;	on FT mass spectrometers.
;
internal_standard_mass = 0.0
; 	If using an internal standard, insert its mass 
;	here for the script to calculate its area. This 
;	also allows the script to output a normalized 
;	metaboanalyst compatible file.
;

[analysis_parameters]
analyze_ms2 = yes
force_fragments_to_glycans = yes
unrestricted_fragments = no
; 	Allows to analyze ms2 data, as well. Fragments 
;	identified will be associated with each glycan. 
;	You can choose to filter identified fragments by 
;	monosaccharides compositions, in order to avoid 
;	reporting fragments that aren't compatible with 
;	detected precursor. If unrestricted_fragments is 
;	used, it searches for glycans in every ms2 scan, 
;	regardless if the glycan was found in full scan. 
;	This will take a bit longer. 
;
accuracy_unit = ppm
; 	Determines the units of mz tolerance to be used 
;	by the script. Options: 'ppm' or 'pw'. 
;	'ppm' = Particles per Million, where 10 ppm is 
;	around 0.01 mz tolerance at mz 1000, 'mz' = Fixed 
;	mz tolerance from centroid, 0.01 mz means it 
;	tolerates a 0.01 variance in mz
;
accuracy_value = 10
; 	The value for the accuracy_unit parameter. You 
;	can use a broader accuracy value and then filter 
;	raw data using max_ppm, but this may lead to 
;	false positives.
;
ret_time_interval = 0, 999
; 	The minimum and maximum retention time, in MINUTES, 
;	used for various portions of the script. A shorter 
;	interval of ret_time makes the script run faster, 
;	so try to trim your sample as much as possible, 
;	if you know when your analytes are leaving the 
;	column.
;
custom_min_points_per_peak = no
number_points_per_peak = 5
; 	If used, set the minimum number of datapoints to 
;	consider a chromatogram peak part of the raw 
;	dataset. If left on False it calculates 
;	automatically.
;
limit_peaks_picked = yes
max_number_peaks = 5
; 	If used, picks only the most intense peak on the 
;	EIC and up to [max_number_peaks]-1 other peaks 
;	closest to it. Warning: This may reduce the range 
;	of your results.
;

[post-analysis/reanalysis]
filter_ms2_by_reporter_ions = N1T1, 366.14
;	Set reporter ions to hide MS2 spectra that don't
;	contain them. Can be based with glycans formula 
;	(with T being the reducing end of the glycan, 
;	including possibly the tag, if used) or an mz.
;	
align_chromatograms = yes
; 	If enabled, will align the assignments and drawn 
;	processed EICs of the different samples. The 
;	alignment is highly dependent on the features 
;	identified and their inherent quality, so things 
;	will change with different quality thresholds. 
;
auc_percentage_threshold = 1
; 	Allows you to supress from the analysis peaks 
;	that are of the specified percentage (from 0% to
;	100%) area under curve related to the most intense
;	peak area within the same adduct (ie. if biggest 
;	peak has a area under curve of 100 and 
;	auc_percentage_threshold is set to 1%, every peak 
;	with an auc of 1 and below will be supressed)
;
max_ppm_threshold = 10
; 	Maximum PPM error for data curation. If value is 
;	greater than equivalent accuracy_value, data won't 
;	be filtered by this criteria, as it was already 
;	filtered during processing by accuracy_value. 
;
isotopic_fitting_score_threshold = 0.9
; 	Minimum score of the isotopic distribution fitting 
;	in order to consider a mz peak viable.
;
curve_fitting_score_threshold = 0.9
; 	Minimum score for the chromatogram peak curve 
;	fitting to a gaussian to consider a viable peak. 
;
signal_to_noise_threshold = 3
; 	Minimum signal-to-noise ratio to consider a 
;	chromatogram peak viable. Can be reapplied on 
;	raw data reanalysis.
;
output_compositions_analysis = yes
; 	If used, also plots data related to the whole 
;	composition of each identified glycan in the 
;	analysis, in addition to the peak-separated data.
;
output_metaboanalyst_file = no
metaboanalyst_groups = CONTROL, TREATED
; 	Here you set up whether or not you want to output 
;	a .csv file to be used for plotting data using 
;	metaboanalyst. If you want that, you must specify 
;	your sample groups, comma separated. Sample 
;	groups specified must be present in sample 
;	filenames for proper identification. If none is 
;	set, samples are defaulted to "ungrouped". Case sensitive. 
;
output_fittings_data = no
; 	Allows to output files with the fittings data to 
;	check scoring criterias. Defaultted to 'no' as 
;	these files will be big. Only use it if you really 
;	need.
;
output_plot_data = no
; 	Allows to output data plotting files for all the 
;	EICs drawn by the program. If set to 'no', it will 
;	still output the found glycans EIC.
;
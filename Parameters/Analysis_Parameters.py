##---------------------------------------------------------------------------------------
##Parameters for the data analysis

multithreaded_analysis = (True, 20)
'''Set to true and specify the amount of threads you want to split your analysis
   in, if you want to do a multithreaded analysis. After running GlycoGenius,
   there will be several .py files in GlycoGenius' folder.
   Run the ones named Multithreaded_0.py-Multithreaded_n.py in whatever time
   you want, in how many amounts at a time you want, in how many computers you
   want. It's only important to run the last one in the same folder where the 
   results_n.py files are.
'''
analyze_ms2 = (False, True)
'''Allows to analyze ms2 data, as well. Fragments identified will be associated with each
   glycan. First boolean determines whether or not to look for ms2 data. Second boolean
   determines if the raw_data acquired will be filtered by monosaccharides compositions, in
   order to avoid reporting fragments that aren't compatible with detected precursor.
'''
accuracy_unit = "pw"
'''Determines the units of mz tolerance to be used by the script. Options: 'ppm' or 'pw'.
   'ppm' = Particles per Million, where 10 ppm is around 0.01 mz tolerance
   'pw' = Peak width, 0.01 pw means it tolerates a 0.01 mz variance
'''
accuracy_value = 0.01
'''The value for the accuracy_unit parameter. You can use a broader accuracy value and then
   filter raw data using max_ppm, but this may lead to false positives.
'''
ret_time_interval = (40, 80)
'''The minimum and maximum retention time used for various portions of the script. A
   shorter interval of ret_time makes the script run faster, so try to trim your sample as
   much as possible, if you know when your analytes are leaving the column.
'''
min_isotopologue_peaks = 3
'''Minimum amount of isotopologue peaks that an identified glycan mz must have to actually
   be taken into account by the script. Minimum amount is 2 (second one necessary to confirm
   charge). May affect isotopic distribution fitting score and can't be recalculated on data 
   reanalysis.
'''
min_points_per_peak = (False, 5)
'''If first parameter is set to True, set the minimum number of datapoints to consider a 
   chromatogram peak part of the raw dataset. If left on False it calculates automatically.
'''
max_ppm = 10
'''Maximum PPM for data curation. If value is greater than equivalent accuracy_value, data won't
   be filtered by this criteria, as it was already filtered during processing by accuracy_value.
   > Can be reapplied on raw data reanalysis.
'''
isotopic_fitting_score = 0.6
'''Minimum score of the isotopic distribution fitting in order to consider a mz peak viable.
   > Can be reapplied on raw data reanalysis.
'''
curve_fitting_score = 0.9
'''Minimum score for the chromatogram peak curve fitting to a gaussian to consider a viable peak.
   > Can be reapplied on raw data reanalysis.
'''
signal_to_noise = 3
'''Minimum signal-to-noise ratio to consider a chromatogram peak viable.
   > Can be reapplied on raw data reanalysis.
'''
custom_noise_level = (False, [])
'''Sets a custom level of noise for each sample, ignoring the automatic noise calculation. Noise for
   each sample is comma separated in second parameter.
   > Warning: NOT RECOMMENDED unless you're really sure about what you're doing.
'''
samples_list = ["D:/Arquivos/Desktop/Data to analyze/221013_PatuS_0cell_01_S1-D4_01_4594.mzXML", "D:/Arquivos/Desktop/Data to analyze/221027_patuS_1cell_01_S1-D4_01_4605.mzXML", "D:/Arquivos/Desktop/Data to analyze/221027_patus_5cell_01_S1-D4_01_4607.mzXML"]
'''A list where each index contains the path to a file to be analyzed together.
'''
results_save_path = ""
'''Directory to save results. Leave blank to save it in GlycoGenius root folder. Must end with '/' if
   not blank.
   > Warning: Must point to an EXISTING directory.
'''
reanalysis = (False, True)
'''Reanalyzes raw data with new max_ppm, isotopic_fitting_score, curve_fitting_score and 
   signal_to_noise criteria. Overrides any other setting besides these mentioned. First parameter 
   produces a new Results file, second parameter also produces a new Plotting Data file (in case you
   deleted your original one. The data in it will not be any different than the former one).
   > Warning: If setting a stricter max_ppm criteria on reanalysis without remaking the whole 
   execution with a new accuracy_value, data may still contain false positives.
'''

##---------------------------------------------------------------------------------------
##Parameters for the data analysis

multithreaded_analysis = (True, 20)
'''Set to true and specify the amount of threads you want to split your analysis
   in, if you want to do a multithreaded analysis. After running GlycoGenius,
   there will be several .py files in GlycoGenius' folder.
   Run the ones named Multithreaded_0.py-Multithreaded_n.py in whatever time
   you want, in how many amounts at a time you want, in how many computers you
   want.
   After all of them are done, run join_results.py to have a single result file
   called results.txt
'''
accuracy_unit = "pw"
'''Determines the units of mz tolerance to be used by the script. Options: 'ppm' or 'pw'.
   'ppm' = Particles per Million, where 10 ppm is around 0.01 mz tolerance
   'pw' = Peak width, 0.01 pw means it tolerates a 0.01 mz variance
'''
accuracy_value = 0.01
'''The value for the accuracy_unit parameter.
'''
threshold = 5000
'''Intensity cutoff to be used at various portions of the script.
'''
ret_time_interval = (10, 20)
'''The minimum and maximum retention time used for various portions of the script. A
   shorter interval of ret_time makes the script run faster, so try trim your sample as
   much as possible, if you know when your analytes are leaving the column.
'''
peak_width = 0.6
'''Expected peaks width, in minutes. Smaller values makes it easy to detect different
   peaks that are close together, but may create false positive peaks.
'''
peak_area_threshold = (threshold*(peak_width*60))/2
'''Minimum area of deconvoluted EIC peak to consider an actual peak. Default formula is
   a triangle area based on threshold and peak_width ((threshold*(peak_width*60))/2).
'''
max_peaks = 3
'''Maximum amount of peaks to be picked per mz.
'''
peaks_relative_intensity = 0.6
'''The less intense peaks in multi-peak picking must be at least this relative intensity
when compared with the most intense peak for them to be considered.
'''
samples_list = ["D:/Arquivos/Desktop/Data to analyze/20231101_BiaRajsfus_CHIKV_HG_1-34_1_1578.mzXML",
                "D:/Arquivos/Desktop/Data to analyze/20231101_BiaRajsfus_CHIKV_NG_1-33_1_1577.mzXML"]
'''A list where each index contains the path to a file to be analyzed together.
'''
analyze_ms2 = False
'''Allows to analyze ms2 data, as well. Fragments identified will be associated with each
   glycan. (Work in Progress)
'''
output_file = 'results.txt'
'''File name to print results after full analysis in singlethreaded.
'''

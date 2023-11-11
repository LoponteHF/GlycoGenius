from pyteomics import mzxml, mzml, mass, auxiliary
from itertools import combinations_with_replacement
from .General_Functions import form_to_comp, form_to_charge
from .Library_Tools import generate_glycans_library, full_glycans_library
from .File_Accessing import eic_from_glycan, peaks_auc_from_eic, peaks_from_eic, check_monoisotopic_charge, deconvoluted_glycan_eic
from re import split
from math import inf
import sys
import datetime

##---------------------------------------------------------------------------------------
##Functions to be used for execution and organizing results data

def list_of_data(samples_list): ##mzML must be worked on or discarded
    '''Detects if a file is mzXML or mzML and processes it into a generator using
    pyteomics.

    Parameters
    ----------
    No parameters needed, but must be executed after parameters section of script.

    Uses
    ----
    pyteomics.mzxml.MzXML() : generator
        Indexes the mzXML file into a generator, allowing you to parse the file for
        analysis.

    pyteomics.mzml.MzML() : generator
        Indexes the mzML file into a generator, allowing you to parse the file for
        analysis.

    Returns
    -------
    data : list
        A list containing the generator of each file name at each index.
    '''
    data = []
    for i in samples_list:
        if i[-5:] == "mzXML":
            data.append(mzxml.MzXML(i))
        elif i[-4:] == "mzML":
            data.append(mzml.MzML(i))
        else:
            sys.exit(i+
                   " filename wrong."+
                   " Please, check if file extension is either 'mzML' or 'mzXML'")
    return data

def index_ms1_from_file(files):
    '''Scans the mz(X)ML file and indexes all the MS1 scans, so that you don't have to
    go through the MSn scans as well when doing something only with the MS1.

    Parameters
    ----------
    files : generator
        List with each index containing a generator created by the pyteomics function
        pyteomics.mzxml.MzXML().

    Returns
    -------
    indexes : dict
        Returns a dictionary with each key pointing to a file and each key containing a
        list of indexes of the MS1 spectra.
    '''
    indexes = {}
    for i_i, i in enumerate(files):
        temp_indexes = []
        for j_j, j in enumerate(i):
            if j['msLevel'] == 1:
                temp_indexes.append(j_j)
        indexes[i_i] = temp_indexes
    return indexes

def index_ms2_from_file(files):
    '''Scans the mz(X)ML file and indexes all the MS2 scans, so that you don't have to
    go through the MS1 scans as well when doing something only with the MS2.

    Parameters
    ----------
    files : generator
        List with each index containing a generator created by the pyteomics function
        pyteomics.mzxml.MzXML().

    Returns
    -------
    indexes : dict
        Returns a dictionary with each key pointing to a index of a file in the files
        list and each key containing a list of indexes of the MS2 spectra.
    '''
    indexes = {}
    for i_i, i in enumerate(files):
        temp_indexes = []
        for j_j, j in enumerate(i):
            if j['msLevel'] == 2:
                temp_indexes.append(j_j)
        indexes[i_i] = temp_indexes
    return indexes

def imp_exp_gen_library(full_library,
                        multithreaded_analysis,
                        samples_names,
                        custom_glycans_list,
                        min_max_monos,
                        min_max_hex,
                        min_max_hexnac,
                        min_max_sia,
                        min_max_fuc,
                        min_max_ac,
                        min_max_gc,
                        force_nglycan,
                        max_adducts,
                        max_charges,
                        intensoids_clumping,
                        tag_mass,
                        fast_iso,
                        imp_exp_library): ##need rework
    '''Imports, generates and/or exports a glycans library.

    Parameters
    ----------
    No parameters needed, but must be executed after parameters section of script.

    Uses
    ----
    generate_glycans_library() : list
        A list containing dictionaries of the monosaccharides compositions of all the
        glycans generated.

    full_glycans_library() : dict
        A dictionary with each key containing the glycan formula and each key containing
        a dictionary with monosaccharides composition, atoms composition with tag,
        neutral mass, neutral mass with tag, isotopic distribution and the mzs of the
        glycans with the desired adducts combination.

    Returns
    -------
    full_library : dict
        A dictionary with each key containing the glycan formula and each key containing
        a dictionary with monosaccharides composition, atoms composition with tag,
        neutral mass, neutral mass with tag, isotopic distribution and the mzs of the
        glycans with the desired adducts combination.
    '''
    begin_time = datetime.datetime.now()
    if multithreaded_analysis[0]:
        print('Preparing to split library and execution for multiple threads...')
        if not imp_exp_library[0]:
            monos_library = generate_glycans_library(min_max_monos,
                                                     min_max_hex,
                                                     min_max_hexnac,
                                                     min_max_sia,
                                                     min_max_fuc,
                                                     min_max_ac,
                                                     min_max_gc,
                                                     force_nglycan)
            full_library = full_glycans_library(monos_library,
                                                max_adducts,
                                                max_charges,
                                                intensoids_clumping,
                                                tag_mass,
                                                fast_iso)
        lib_names = []
        split_s = len(full_library)/multithreaded_analysis[1]
        next = 0
        for i in range(multithreaded_analysis[1]):
            lib_names.append('glycans_library_'+str(i)+'.py')
            with open(lib_names[-1], 'w') as f:
                f.write('full_library = {')
                f.close()
            for j_j, j in enumerate(full_library):
                with open(lib_names[-1], 'a') as f:
                    if j_j > next+split_s or j_j == len(full_library)-1:
                        f.write("'"+j+"'"+": "+str(full_library[j])+"}")
                        next = j_j+1
                        f.close()
                        break
                    elif j_j >= next:
                        f.write("'"+j+"'"+": "+str(full_library[j])+", ")
                        f.close()
        with open('GlycoGenius.py', 'r') as f:
            for i in f:
                for j in range(multithreaded_analysis[1]):
                    with open('Multithreaded_'+str(j)+'.py', 'a') as g:
                        if i[:22] == 'custom_glycans_list = ':
                            g.write(i[:22]+"(False, '')\n")
                            g.close()
                            continue
                        if i[:18] == 'imp_exp_library = ':
                            g.write(i[:18]+"(True, False)\n")
                            g.close()
                            continue
                        if i[:25] == 'multithreaded_analysis = ':
                            g.write(i[:25]+"(False, 0)\n")
                            g.close()
                            continue
                        if i[:14] == 'output_file = ':
                            g.write(i[:14]+"'results_"+str(j)+".txt'\n")
                            g.close()
                            continue
                        if i[:26] == 'multithreaded_execution = ':
                            g.write(i[:26]+"True\n")
                            g.close()
                            continue
                        if i[:24] == '    from glycans_library':
                            g.write(i[:24]+"_"+str(j)+" import full_library\n")
                            g.close()
                            continue
                        g.write(i)
                        g.close()
            f.close()
        with open('join_results.py', 'w') as f:
            f.write("with open('results.txt', 'a') as f:\n")
            f.write("   f.write('Glycans\\tAdduct\\tmz\\t")
            for i in samples_names:
                f.write(str(i)+"\\t")
            f.write("\\n')\n   for i in range("+str(multithreaded_analysis[1])+"):\n")
            f.write("      with open('results_'+str(i)+'.txt', 'r') as g:\n")
            f.write("         for i in g:\n")
            f.write("            f.write(i)\n")
            f.write("         g.close()\n")
            f.write("   f.close()\n")
        input("Multithreaded run setup done. Press Enter to exit and run the"+
              "'Multithreaded_n.py' files to execute each part of the script."+
              " After that execution of all the .py files, run 'join_results.p"+
              "y' and access the final results in 'results.txt'")
        sys.exit()
    if custom_glycans_list[0]:
        print('Building custom glycans library...')
        for i in custom_glycans_list[1]:
            custom_glycans_comp.append(form_to_comp(i))
        full_library = full_glycans_library(custom_glycans_comp,
                                            max_adducts,
                                            max_charges,
                                            intensoids_clumping,
                                            tag_mass,
                                            fast_iso)
        print('Done building custom glycans library in '+
              str(datetime.datetime.now()-begin_time))
        return full_library
    if imp_exp_library[0]:
        print('Done importing existing glycans library in '+
              str(datetime.datetime.now()-begin_time))
        return full_library
    print('Building glycans library...')
    monos_library = generate_glycans_library(min_max_monos,
                                             min_max_hex,
                                             min_max_hexnac,
                                             min_max_sia,
                                             min_max_fuc,
                                             min_max_ac,
                                             min_max_gc,
                                             force_nglycan)
    full_library = full_glycans_library(monos_library,
                                        max_adducts,
                                        max_charges,
                                        intensoids_clumping,
                                        tag_mass,
                                        fast_iso)
    if imp_exp_library[1]:
        print('Exporting glycans library...')
        with open('glycans_library.py', 'w') as f:
            f.write('full_library = '+str(full_library))
    print('Done building/exporting glycans library in '+
          str(datetime.datetime.now()-begin_time)+'!')
    return full_library

def sample_names(samples_list):
    '''Extracts the sample names from the file path.
    Parameters
    ----------
    No parameters needed, but must be executed after parameters section of script.

    Returns
    -------
    curated_samples : list
        A list of strings, each string with a sample name.
    '''
    curated_samples = []
    for i in samples_list:
        dot = 0
        backlash = 0
        for j in range(len(i)-1, -1, -1):
            if i[j] == '.':
                dot = j
            if i[j] == '/':
                backlash = j
                break
        curated_samples.append(i[backlash+1:dot])
    return curated_samples

def tolerance(unit,
              value):
    '''A fast, but not super accurate way to convert 'ppm' mass accuracy into 'pw' (peak
    width, aka. mz tolerance).

    Parameters
    ----------
    unit : string
        Can be either "ppm" (particles-per-million) or "pw" (peak width [tolerance]).

    value : float
        Float value of the tolerance, based on the unit inputted.

    Returns
    -------
    tolerance : float
        If unit == "ppm", converts value into "pw", if unit == "pw", outputs value as is.
        ie. 10 ppm gets converted to 0.01 pw tolerance.
    '''
    if unit == "ppm":
        return (1000.0-(-((value*1000)/10**6)+1000))
    elif unit == "pw":
        return value
    else:
        return("Unit for tolerance not 'ppm' or 'pw'.")

def generate_lines_to_print(analyzed_data,
                            samples_names,
                            peak_width,
                            multithreaded_execution): ##Complete
    '''Arrange the raw data to print. Information printed: Glycan formula, Adduct, sum of
    the corresponding peak in all adducts, mz or neutral mass (for total adducts),
    retention time of each peak and AUC of it.

    Parameters
    ----------
    No parameters needed, but must be executed after parameters section of script.

    Returns
    -------
    lines_to_print : list
        A list containing each line that should be printed.
    '''
    begin_time = datetime.datetime.now()
    print('Arranging raw data...')
    data_to_print = [[],[],[]]
    len_adducts = 0
    for i in range(len(samples_names)):
        data_to_print.append([])
    for i_i, i in enumerate(analyzed_data): #i = glycan (key)
        if i_i == 0:
            len_adducts = len(analyzed_data[i]['Adducts_mz'])
        for j in analyzed_data[i]['Adducts_mz_data']: #j = adduct (key)
            data_to_print[0].append(i)
            data_to_print[1].append(j)
            data_to_print[2].append(analyzed_data[i]['Adducts_mz'][j])
            for k in analyzed_data[i]['Adducts_mz_data'][j]: #k = sample (key)
                data_to_print[3+k].append([])
                for l_l, l in enumerate(analyzed_data[i]
                                        ['Adducts_mz_data'][j][k][0]):
                    temp_rt = float("%.2f" % round(l, 2))
                    temp_auc = float("%.2f" % round(analyzed_data[i]
                                                    ['Adducts_mz_data'][j][k][1][l_l], 2))
                    data_to_print[3+k][-1].append((temp_rt, temp_auc))
    lines_to_print = []
    if not multithreaded_execution:
        lines_to_print.append('Glycan\tAdduct\tmz\t')
        for i in samples_names:
            lines_to_print[-1]+=str(i)+'\t'
    for i in range(len(data_to_print[0])):
        found = False
        for j in data_to_print[3:]:
            if len(j[i]) != 0:
                found = True
        if found:
            for j_j, j in enumerate(data_to_print):
                if j_j == 0:
                    lines_to_print.append(str(j[i])+'\t')
                    continue
                if j_j == 1:
                    lines_to_print[-1]+=str(j[i])+'\t'
                if j_j == 2:
                    lines_to_print[-1]+=str("%.4f" % round(j[i], 4))+'\t'
                if j_j >= 3:
                    per_sample = '['
                    for k_k, k in enumerate(j[i]):
                        if len(k) != 0:
                            per_sample+=str(k)
                    per_sample+=']\t'
                    lines_to_print[-1]+= per_sample
        if i != 0 and ((i+1)%(len_adducts) == 0): #make 'total' line
            peaks_sum_sample = []
            for j_j, j in enumerate(data_to_print[3:]): #j = sample
                peaks_sum_sample.append([])
                for k in range(i, i-len_adducts, -1):
                    for l in j[k]: #l = peak in a rt_auc of a sample
                        if len(peaks_sum_sample[j_j]) == 0:
                            peaks_sum_sample[j_j].append([l[0], l[1]])
                        else:
                            found = False
                            for m in peaks_sum_sample[j_j]:
                                if abs(l[0]-m[0]) < peak_width/2:
                                    m[0] = (l[0]+m[0])/2
                                    m[1]+= l[1]
                                    found = True
                            if not found:
                                peaks_sum_sample[j_j].append([l[0], l[1]])
            found = False
            for j in peaks_sum_sample:
                if len(j) != 0:
                    found = True
                    for k in j:
                        k[0] = float("%.2f" % round(k[0], 2))
                        k[1] = float("%.2f" % round(k[1], 2))
            if found:
                lines_to_print.append(str(data_to_print[0][i-1])+'\tTotal\t'+
                                  str("%.4f" % round(analyzed_data[data_to_print[0][i-1]]
                                      ['Neutral_Mass+Tag'], 4))+'\t')
                for j in peaks_sum_sample:
                    lines_to_print[-1]+=str(j)+'\t' #make 'total' line
    print('Data arrangement done in '+str(datetime.datetime.now() - begin_time)+'!')
    return lines_to_print        

def print_to_file(lines_to_print, output_file): ##Complete, just add comment
    '''
    '''
    begin_time = datetime.datetime.now()
    print('Printing to file...')
    with open(output_file, 'w') as f:
        for i in lines_to_print:
            f.write(i+'\n')
        f.close()
    print('Done printing to file in '+str(datetime.datetime.now()-begin_time)+'!')
    print('You may check your data on the file: '+output_file)

def print_sep(): ##Complete
    '''Prints a separator consisting of multiple '-' character.
    '''
    print('----------------------------------------------------------------------------')

def analyze_files(library,
                  lib_size,
                  data,
                  ms1_index,
                  tolerance,
                  threshold,
                  peak_width,
                  ret_time_interval,
                  max_peaks,
                  peaks_relative_intensity,
                  peak_area_threshold): ##Complete
    '''Integrates all the file-accessing associated functions in this script to go
    through the files data, draw eic of hypothetical glycans, check if it is not in
    noise level, based on threshold, does peak-picking, check quality of spectra in
    peaks, deconvolute if checks passed, and calculates AUC of the peaks.

    Parameters
    ----------
    No parameters needed, but must be executed after parameters section of script.

    Uses
    ----
    eic_from_glycan() : dict
        Generates an EIC of the mzs calculated for a glycan.

    peaks_from_eic() : list, tuple
        Does multi-peak-picking in an EIC, based on the most intense peaks in the EIC,
        in the retention time interval of interest. This function works in-line.

    check_monoisotopic_charge() : list, list
        Checks if the peaks identified for a given mz corresponds to monoisotopic,
        correctly charge-assigned peaks.

    deconvoluted_glycan_eic() : list
        Deconvolutes the data around the identified peaks retention time and creates a
        deconvoluted EIC.

    peaks_auc_from_eic() : list
        Calculates the area under curve of the picked peaks based on the EIC given.
        Overlapped portions of peaks are calculated as fractional proportion, based on
        peaks max intensity.

    Returns
    -------
    analyzed_data : dict
        A dictionary similar to the one generated by the full glycans library generating
        function, with the added information of each adducts' peaks' retention time and
        AUC.
    '''
    begin_time = datetime.datetime.now()
    print('Analyzing glycans in samples...')
    analyzed_data = {}
    for i_i, i in enumerate(library):
        glycan_data = library[i]
        glycan_data['Adducts_mz_data'] = {}
        print('Analyzing glycan '+str(i)+': '+str(i_i+1)+'/'+str(lib_size))
        temp_eic = eic_from_glycan(data,
                                   glycan_data,
                                   ms1_index,
                                   tolerance)
        for j in temp_eic:
            glycan_data['Adducts_mz_data'][j] = {}
            for k in temp_eic[j]:
                glycan_data['Adducts_mz_data'][j][k] = [[],[]]
                if max(temp_eic[j][k][1]) < threshold:
                    continue
                temp_peaks = peaks_from_eic(temp_eic[j][k],
                                            peak_width,
                                            threshold,
                                            ret_time_interval,
                                            max_peaks,
                                            peaks_relative_intensity)
                peaks_check = check_monoisotopic_charge(data[k],
                                                        ms1_index[k],
                                                        glycan_data['Adducts_mz'][j],
                                                        glycan_data['Isotopic_Distribution'][1],
                                                        form_to_charge(j),
                                                        threshold,
                                                        temp_peaks[0],
                                                        tolerance)
                to_be_removed = []
                for l_l, l in enumerate(temp_peaks[0]):
                    if not peaks_check[0][l_l] or not peaks_check[1][l_l]:
                        to_be_removed.append(l)
                for l in to_be_removed:
                    temp_peaks[0].remove(l)
                if len(temp_peaks) == 0:
                    continue
                dec_eic = deconvoluted_glycan_eic(data[k],
                                                  ms1_index[k],
                                                  glycan_data['Isotopic_Distribution'],
                                                  glycan_data['Adducts_mz'][j],
                                                  temp_peaks[1],
                                                  threshold,
                                                  form_to_charge(j),
                                                  tolerance)
                temp_peaks_auc = peaks_auc_from_eic(dec_eic,
                                                    temp_peaks[0],
                                                    peak_width)
                for l_l, l in enumerate(temp_peaks[0]):
                    if temp_peaks_auc[l_l] > peak_area_threshold:
                        glycan_data['Adducts_mz_data'][j][k][0].append(l['rt'])
                        glycan_data['Adducts_mz_data'][j][k][1].append(temp_peaks_auc[l_l])
        analyzed_data[i] = glycan_data
    print('Sample analysis done in '+str(datetime.datetime.now() - begin_time)+'!')
    return analyzed_data

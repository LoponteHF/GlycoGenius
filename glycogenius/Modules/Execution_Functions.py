import pathlib
import importlib
import_path = str(pathlib.Path(__file__).parent.resolve())
for i_i, i in enumerate(import_path):
    if i == "\\":
        import_path = import_path[:i_i]+"/"+import_path[i_i+1:]
        
#Absolute import of General_Functions
spec1 = importlib.util.spec_from_file_location("General_Functions", import_path+'/General_Functions.py')
General_Functions = importlib.util.module_from_spec(spec1)
spec1.loader.exec_module(General_Functions)
#Absolute import of Library_Tools
spec2 = importlib.util.spec_from_file_location("Library_Tools", import_path+'/Library_Tools.py')
Library_Tools = importlib.util.module_from_spec(spec2)
spec2.loader.exec_module(Library_Tools)
#Absolute import of File_Accessing
spec3 = importlib.util.spec_from_file_location("File_Accessing", import_path+'/File_Accessing.py')
File_Accessing = importlib.util.module_from_spec(spec3)
spec3.loader.exec_module(File_Accessing)

from pyteomics import mzxml, mzml, mass, auxiliary
from itertools import combinations_with_replacement, islice
from pandas import DataFrame, ExcelWriter
from numpy import percentile
from re import split
from math import inf, isnan
from statistics import mean
from time import sleep
import os
import dill
import sys
import datetime
import traceback

##---------------------------------------------------------------------------------------
##Functions to be used for execution and organizing results data

def print_header(complete = True):
    '''
    '''
    print("\n    GlycoGenius: Glycomics Data Analysis Tool")
    print("    Copyright (C) 2023 by Hector Franco Loponte")
    if not complete:
        print("For more details about the license, run the package stand-alone by")
        print("typing 'python -m glycogenius' in terminal and then type 'license'.")
    if complete:
        print("This program comes with ABSOLUTELY NO WARRANTY; for details type 'warranty'.")
        print("This is free software and can be redistributed under certain conditions.")
        print("If you want to know more details about the licensing, type 'license'.")
    print_sep()
    
def generate_cfg_file(path, comments):
    '''
    '''
    print("Creating settings file...")
    glycogenius_path = str(pathlib.Path(__file__).parent.parent.resolve())
    for i_i, i in enumerate(glycogenius_path):
        if i == "\\":
            glycogenius_path = glycogenius_path[:i_i]+"/"+glycogenius_path[i_i+1:]
    with open(path+'glycogenius_parameters.ini', 'w') as g:
        with open(glycogenius_path+'/Parameters_Template.ini', 'r') as f:
            for line in f:
                if line == "working_path = C:/GlycoGenius/":
                    g.write("working_path = "+path+"\n")
                    continue
                if not comments and line[0] == ';':
                    continue
                g.write(line)
        f.close()
    g.close()
    input("Done! Press Enter to exit.")
    os._exit(1)

def interactive_terminal():
    '''
    '''
    date = datetime.datetime.now()
    begin_time = str(date)[2:4]+str(date)[5:7]+str(date)[8:10]+"_"+str(date)[11:13]+str(date)[14:16]+str(date)[17:19]
    input_order = [None]
    while input_order[0] == None:
        print_header()
        print("1 - Build and output glycans library.\n2 - Analyze sample files in single-threaded mode\n3 - Reanalyze raw results files with new parameters\n4 - Create template parameters file for command-line execution\n")
        var = input("Select your option: ")
        if var == 'warranty':
            print("\nDisclaimer of Warranty.\n")
            print("THERE IS NO WARRANTY FOR THE PROGRAM, TO THE EXTENT PERMITTED BY")
            print("APPLICABLE LAW.  EXCEPT WHEN OTHERWISE STATED IN WRITING THE COPYRIGHT")
            print("HOLDERS AND/OR OTHER PARTIES PROVIDE THE PROGRAM \"AS IS\" WITHOUT WARRANTY")
            print("OF ANY KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO,")
            print("THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR")
            print("PURPOSE.  THE ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE PROGRAM")
            print("IS WITH YOU.  SHOULD THE PROGRAM PROVE DEFECTIVE, YOU ASSUME THE COST OF")
            print("ALL NECESSARY SERVICING, REPAIR OR CORRECTION.\n")
            continue
        if var == 'license':
            license_path = str(pathlib.Path(__file__).parent.parent.resolve())
            for i_i, i in enumerate(import_path):
                if i == "\\":
                    license_path = license_path[:i_i]+"/"+license_path[i_i+1:]
            with open(license_path+"/LICENSE", 'r') as f:
                for line in f:
                    print(line, end = "")
            continue
        try:
            var = int(var)
        except:
            print("Wrong Input")
            continue
        if var < 1 or var > 4:
            print("Wrong Input")
            continue
        if var > 0 and var <= 4:
            input_order[0] = var
    print_sep()
    if input_order[0] == 1 or input_order[0] == 2:
        input_order.append(None)
        while input_order[1] == None:
            print("1 - Use a custom glycans list\n2 - Generate a library based on possible monosaccharides numbers")
            var = input("Select your option: ")
            try:
                var = int(var)
            except:
                continue
            if var < 1 or var > 2:
                continue
            if var > 0 and var <= 2:
                input_order[1] = var
        print_sep()
        if input_order[1] == 1:
            glycans_list = []
            print("Warning: Glycans must be inputted in the format of 'H5N4S1F1G1' where 'H' = Hexose, 'N' = N-Acetylhexosamine, 'S' = N-Acetyl-Neuraminic Acid, 'F' = Deoxyhexose, 'G' = N-Glycolyl-Neuraminic Acid and each number next to it corresponds to the amount of said monosaccharide")
            while True:
                var = input("Insert a glycan, leave blank to end list: ")
                if var == "":
                    print(glycans_list)
                    var2 = input("Proceed with this list? (y/n): ")
                    if var2 == 'y':
                        break
                    if var2 == 'n':
                        print("Emptied glycans list.")
                        glycans_list = []
                        continue
                try:
                    comp = General_Functions.form_to_comp(var)
                except:
                    print('Wrong input')
                    continue
                for i in comp:
                    if i != 'H' and i != 'N' and i != 'S' and i != 'F' and i != 'G':
                        print('Wrong input')
                        continue
                glycans_list.append(var)
                print("Current glycans: ", glycans_list)
        if input_order[1] == 2:
            lib_settings = [None, None, None, None, None, None, None, None, None, None, None, None, None, None, None]
            for i in range(len(lib_settings)):
                if i == 0:
                    while lib_settings[0] == None:
                        var = input("Type the minimum amount of total monosaccharides: ")
                        try:
                            var = int(var)
                        except:
                            continue
                        lib_settings[0] = var
                if i == 1:
                    while lib_settings[1] == None:
                        var = input("Type the maximum amount of total monosaccharides: ")
                        try:
                            var = int(var)
                        except:
                            continue
                        lib_settings[1] = var
                if i == 2:
                    while lib_settings[2] == None:
                        var = input("Type the minimum amount of hexoses: ")
                        try:
                            var = int(var)
                        except:
                            continue
                        lib_settings[2] = var
                if i == 3:
                    while lib_settings[3] == None:
                        var = input("Type the maximum amount of hexoses: ")
                        try:
                            var = int(var)
                        except:
                            continue
                        lib_settings[3] = var
                if i == 4:
                    while lib_settings[4] == None:
                        var = input("Type the minimum amount of N-acetylhexosamines: ")
                        try:
                            var = int(var)
                        except:
                            continue
                        lib_settings[4] = var
                if i == 5:
                    while lib_settings[5] == None:
                        var = input("Type the maximum amount of N-acetylhexosamines: ")
                        try:
                            var = int(var)
                        except:
                            continue
                        lib_settings[5] = var
                if i == 6:
                    while lib_settings[6] == None:
                        var = input("Type the minimum amount of deoxyhexoses (fucoses): ")
                        try:
                            var = int(var)
                        except:
                            continue
                        lib_settings[6] = var
                if i == 7:
                    while lib_settings[7] == None:
                        var = input("Type the maximum amount of deoxyhexoses (fucoses): ")
                        try:
                            var = int(var)
                        except:
                            continue
                        lib_settings[7] = var
                if i == 8:
                    while lib_settings[8] == None:
                        var = input("Type the minimum amount of N-acetyl neuraminic acids: ")
                        try:
                            var = int(var)
                        except:
                            continue
                        lib_settings[8] = var
                if i == 9:
                    while lib_settings[9] == None:
                        var = input("Type the maximum amount of N-acetyl neuraminic acids: ")
                        try:
                            var = int(var)
                        except:
                            continue
                        lib_settings[9] = var
                if i == 10:
                    while lib_settings[10] == None:
                        var = input("Type the minimum amount of N-glycolyl neuraminic acids: ")
                        try:
                            var = int(var)
                        except:
                            continue
                        lib_settings[10] = var
                if i == 11:
                    while lib_settings[11] == None:
                        var = input("Type the maximum amount of N-glycolyl neuraminic acids: ")
                        try:
                            var = int(var)
                        except:
                            continue
                        lib_settings[11] = var
                if i == 12:
                    while lib_settings[12] == None:
                        var = input("Type the minimum amount of total sialic acids: ")
                        try:
                            var = int(var)
                        except:
                            continue
                        lib_settings[12] = var
                if i == 13:
                    while lib_settings[13] == None:
                        var = input("Type the maximum amount of total sialic acids: ")
                        try:
                            var = int(var)
                        except:
                            continue
                        lib_settings[13] = var
                if i == 14:
                    while lib_settings[14] == None:
                        var = input("Force compositions to be closer to N-glycans structure (y/n): ")
                        if var == 'y':
                            lib_settings[14] = True
                        if var == 'n':
                            lib_settings[14] = False
                        else:
                            continue
        print_sep()
        adducts = {}
        while True:
            var = input("Type the first element to calculate as adduct (ie. Na or H). Leave blank to finish: ")
            if var == '':
                print(adducts)
                var2 = input("Proceed with these adducts? (y/n): ")
                if var2 == 'y':
                    break
                if var2 == 'n':
                    print("Emptied adducts list.")
                    adducts = {}
                    continue
            var2 = input("Type the maximum number of such adduct: ")
            try:
                var2 = int(var2)
            except:
                print('Wrong input')
                continue
            adducts[var] = var2
        print_sep()
        max_charges = 0
        while True:
            var = input("Type the maximum amount of charges for the glycans: ")
            try:
                var = int(var)
            except:
                print('Wrong input')
                continue
            max_charges = var
            break
        print_sep()
        tag_mass = 0
        while True:
            var = input("Do you have a reducing end tag attached to your glycans? (y/n): ")
            if var == 'y' or var == 'n':
                break
            else:
                print("Wrong Input")
                continue
            break
        if var == 'y':
            while True:
                var2 = input("Insert the tag added mass (ie. 133.0644 for GirP or 219.1735 for ProA): ")
                try:
                    var2 = float(var2)
                except:
                    print('Wrong input')
                    continue
                tag_mass = var2
                break
        print_sep()
        fast_iso = True
        while True:
            var = input("Do you want to do a quick isotopic distribution calculation? If 'n', then isotopic distribution calculation may take several hours, depending on library size (y/n): ")
            if var == 'y':
                break
            if var == 'n':
                fast_iso = False
                break
            else:
                print('Wrong input')
                continue
        high_res = False
        if not fast_iso:
            while True:
                var = input("Do you need a high resolution isotopic distribution? It may be important for very high accuracy mass spectrometers, such as FT-ICR (y/n): ")
                if var == 'y':
                    high_res = True
                    break
                if var == 'n':
                    break
                else:
                    print('Wrong input')
                    continue
        if input_order[0] == 1: #Outputs of input_order == 1
            path = 'C:\\GlycoGenius_'+begin_time+'\\'
            while True:
                var = input("Insert the path to save the files produced by the script (leave blank for default: C:\\GlycoGenius_[current_date_and_time]): ")
                if var == '':
                    var = path
                print(var)
                var2 = input("Is this path correct? (y/n): ")
                if var2 == 'n':
                    continue
                if var2 == 'y':
                    for i_i, i in enumerate(var):
                        if i == "\\":
                            var = var[:i_i]+"/"+var[i_i+1:]
                    if var[-1] != "/":
                        var = var+"/"
                    path = var
                    break
            if input_order[1] == 1:
                return input_order, glycans_list, adducts, max_charges, tag_mass, fast_iso, high_res, path
            if input_order[1] == 2:
                return input_order, lib_settings, adducts, max_charges, tag_mass, fast_iso, high_res, path
        else:
            print_sep()
            ms2 = [False, False]
            while True:
                var = input("Do you wish to analyze MS2 data? (y/n): ")
                if var == 'y':
                    ms2[0] = True
                    break
                if var == 'n':
                    break
                else:
                    print('Wrong input')
                    continue
            if ms2[0]:
                while True:
                    var = input("Do you want to only output fragments compatible with identified precursor glycan? (y/n): ")
                    if var == 'y':
                        ms2[1] = True
                        break
                    if var == 'n':
                        break
                    else:
                        print('Wrong input')
                        continue
            accuracy_unit = "pw"
            while True:
                var = input("What is the accuracy unit you want to input for mz tolerance? (ppm/pw): ")
                if var == 'ppm':
                    accuracy_unit = var
                    break
                if var == 'pw':
                    break
                else:
                    print('Wrong input')
                    continue
            accuracy_value = 0.0
            while True:
                var = input("Insert the accuracy value for the unit you've chosen (ie. '0.01' or '10'): ")
                try:
                    var = float(var)
                except:
                    print('Wrong input')
                    continue
                accuracy_value = var
                break
            rt_int = [0.0, 999]
            while True:
                var = input("Insert the beggining of the retention time interval at which you want to analyze, in minutes: ")
                try:
                    var = float(var)
                except:
                    print('Wrong input')
                    continue
                rt_int[0] = var
                break
            while True:
                var = input("Insert the end of the retention time interval at which you want to analyze, in minutes: ")
                try:
                    var = float(var)
                except:
                    print('Wrong input')
                    continue
                rt_int[1] = var
                break
            min_isotop = 2
            while True:
                var = input("Insert the minimum amount of detected isotopologue peaks for a mz in a spectrum to be included in the processed EIC: ")
                try:
                    var = int(var)
                except:
                    print('Wrong input')
                    continue
                min_isotop = var
                break
            max_ppm = 10
            while True:
                var = input("Insert the maximum amount of PPM difference that a detected glycan must have in order to show up in results' table: ")
                try:
                    var = int(var)
                except:
                    print('Wrong input')
                    continue
                max_ppm = var
                break
            iso_fit = 0.5
            while True:
                var = input("Insert the minimum isotopic fitting score for a glycan in order for it to show up in the results' table (values between 0.0 and 1.0): ")
                try:
                    var = float(var)
                except:
                    print('Wrong input')
                    continue
                if var < 0.0 or var > 1.0:
                    print('Wrong input')
                    continue
                iso_fit = var
                break
            curve_fit = 0.5
            while True:
                var = input("Insert the minimum curve fitting score for a glycan in order for it to show up in the results' table (values between 0.0 and 1.0): ")
                try:
                    var = float(var)
                except:
                    print('Wrong input')
                    continue
                if var < 0.0 or var > 1.0:
                    print('Wrong input')
                    continue
                curve_fit = var
                break
            sn = 3
            while True:
                var = input("Insert the minimum signal-to-noise ratio that a detected glycan must have in order to show up in results' table: ")
                try:
                    var = int(var)
                except:
                    print('Wrong input')
                    continue
                sn = var
                break
            files = []
            while True:
                var = input("Insert the path to the file to be analyzed (ie. C:\\GlycoGenius\\file.mzxml). Leave blank to finalize: ")
                if var == '':
                    print(files)
                    var2 = input("Proceed with these files? (y/n): ")
                    if var2 == 'y':
                        break
                    if var2 == 'n':
                        print("Emptied files list.")
                        files = []
                        continue
                if var[0] == "'" or var[0] == "\"":
                    var = var[1:-1]
                print(var)
                var2 = input("Is this path correct? (y/n): ")
                if var2 == 'n':
                    continue
                if var2 == 'y':
                    for i_i, i in enumerate(var):
                        if i == "\\":
                            var = var[:i_i]+"/"+var[i_i+1:]
                    if var[-1] != "/":
                        var = var+"/"
                    files.append(var)
                    continue
            path = 'C:\\GlycoGenius_'+begin_time+'\\'
            while True:
                var = input("Insert the path to save the files produced by the script (leave blank for default: C:\\GlycoGenius_[current_date_and_time]): ")
                if var == '':
                    var = path
                print(var)
                var2 = input("Is this path correct? (y/n): ")
                if var2 == 'n':
                    continue
                if var2 == 'y':
                    for i_i, i in enumerate(var):
                        if i == "\\":
                            var = var[:i_i]+"/"+var[i_i+1:]
                    if var[-1] != "/":
                        var = var+"/"
                    path = var
                    break
            if input_order[1] == 1:
                return input_order, glycans_list, adducts, max_charges, tag_mass, fast_iso, high_res, ms2, accuracy_unit, accuracy_value, rt_int, min_isotop, max_ppm, iso_fit, curve_fit, sn, files, path
            if input_order[1] == 2:
                return input_order, lib_settings, adducts, max_charges, tag_mass, fast_iso, high_res, ms2, accuracy_unit, accuracy_value, rt_int, min_isotop, max_ppm, iso_fit, curve_fit, sn, files, path
    if input_order[0] == 3:
        path = 'C:\\Glycogenius\\'
        while True:
            var = input("Insert the working directory (where the 'raw_data' files are, default: C:/Glycogenius/): ")
            print(var)
            var2 = input("Is this path correct? (y/n): ")
            if var2 == 'n':
                continue
            if var2 == 'y':
                for i_i, i in enumerate(var):
                    if i == "\\":
                        var = var[:i_i]+"/"+var[i_i+1:]
                if var[-1] != "/":
                    var = var+"/"
                path = var
                break
        max_ppm = 10
        while True:
            var = input("Insert the maximum amount of PPM difference that a detected glycan must have in order to show up in results' table: ")
            try:
                var = int(var)
            except:
                print('Wrong input')
                continue
            max_ppm = var
            break
        iso_fit = 0.5
        while True:
            var = input("Insert the minimum isotopic fitting score for a glycan in order for it to show up in the results' table (values between 0.0 and 1.0): ")
            try:
                var = float(var)
            except:
                print('Wrong input')
                continue
            if var < 0.0 or var > 1.0:
                print('Wrong input')
                continue
            iso_fit = var
            break
        curve_fit = 0.5
        while True:
            var = input("Insert the minimum curve fitting score for a glycan in order for it to show up in the results' table (values between 0.0 and 1.0): ")
            try:
                var = float(var)
            except:
                print('Wrong input')
                continue
            if var < 0.0 or var > 1.0:
                print('Wrong input')
                continue
            curve_fit = var
            break
        sn = 3
        while True:
            var = input("Insert the minimum signal-to-noise ratio that a detected glycan must have in order to show up in results' table: ")
            try:
                var = int(var)
            except:
                print('Wrong input')
                continue
            sn = var
            break
        return input_order, path, max_ppm, iso_fit, curve_fit, sn
    if input_order[0] == 4:
        commented = False
        while True:
            var = input("Do you want the template file to contain commented information on each parameter? (y/n): ")
            if var == 'y':
                commented = True
                break
            if var == 'n':
                break
            else:
                print('Wrong input')
                continue
        path = 'C:\\GlycoGenius\\'
        while True:
            var = input("Insert the path to the folder to save the template file (Default: C:\\GlycoGenius\\): ")
            if var == "":
                var = path
            print(var)
            var2 = input("Is this path correct? (y/n): ")
            if var2 == 'n':
                continue
            if var2 == 'y':
                for i_i, i in enumerate(var):
                    if i == "\\":
                        var = var[:i_i]+"/"+var[i_i+1:]
                if var[-1] != "/":
                    var = var+"/"
                path = var
                break
        return input_order, commented, path
                
def list_of_data(samples_list): ##complete
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
        if i[-5:] == "mzXML" or i[-5:] == "mzxml" or i[-5:] == "MzXML":
            mzxml_data = mzxml.MzXML(i)
            data.append(mzxml_data)
        elif i[-4:] == "mzML" or i[-4:] == "mzml" or i[-4:] == "MzML":
            mzml_data = File_Accessing.make_mzxml(i)
            data.append(mzml_data)
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
            try:
                if j['msLevel'] == 1:
                    temp_indexes.append(j_j)
            except:
                if j['ms level'] == 1:
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
            try:
                if j['msLevel'] == 2:
                    temp_indexes.append(j_j)
            except:
                if j['ms level'] == 2:
                    temp_indexes.append(j_j)
        indexes[i_i] = temp_indexes
    return indexes

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
            if i[j] == '/' or i[j] == '\\':
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

def imp_exp_gen_library(multithreaded_analysis,
                        multithreaded_execution,
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
                        tag_mass,
                        fast_iso,
                        high_res,
                        imp_exp_library,
                        only_gen_lib,
                        save_path):
    '''Imports, generates and/or exports a glycans library.

    Parameters
    ----------
    No parameters needed, but must be executed after parameters section of script.

    Uses
    ----
    Library_Tools.generate_glycans_library() : list
        A list containing dictionaries of the monosaccharides compositions of all the
        glycans generated.

    Library_Tools.full_glycans_library() : dict
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
    if multithreaded_execution[0]:
        return
    if custom_glycans_list[0] and not imp_exp_library[0]:
        custom_glycans_comp = []
        print('Building custom glycans library...')
        for i in custom_glycans_list[1]:
            custom_glycans_comp.append(General_Functions.sum_monos(General_Functions.form_to_comp(i)))
        full_library = Library_Tools.full_glycans_library(custom_glycans_comp,
                                            max_adducts,
                                            max_charges,
                                            tag_mass,
                                            fast_iso,
                                            high_res)
        if imp_exp_library[1]:
            print('Exporting glycans library...')
            with open(save_path+'glycans_library.py', 'w') as f:
                f.write('full_library = '+str(full_library))
                f.close()
        print('Done building custom glycans library in '+
              str(datetime.datetime.now()-begin_time))
    if multithreaded_analysis[0] and not only_gen_lib and not multithreaded_execution[0]:
        print('Preparing to split library and execution for multiple threads...')
        if not imp_exp_library[0] and not custom_glycans_list[0]:
            monos_library = Library_Tools.generate_glycans_library(min_max_monos,
                                                     min_max_hex,
                                                     min_max_hexnac,
                                                     min_max_sia,
                                                     min_max_fuc,
                                                     min_max_ac,
                                                     min_max_gc,
                                                     force_nglycan)
            full_library = Library_Tools.full_glycans_library(monos_library,
                                                max_adducts,
                                                max_charges,
                                                tag_mass,
                                                fast_iso,
                                                high_res)
            if imp_exp_library[1]:
                print('Exporting glycans library...')
                with open(save_path+'glycans_library.py', 'w') as f:
                    f.write('full_library = '+str(full_library))
                    f.close()
        elif imp_exp_library[0]:
            lib_module = importlib.import_module('glycans_library')
            full_library = lib_module.full_library
        lib_names = []
        full_library_keys_list = list(full_library.keys())
        if len(full_library)%multithreaded_analysis[1] == 0:
            split_s = int(len(full_library)/multithreaded_analysis[1])
        else:
            split_s = int(len(full_library)/multithreaded_analysis[1])+1
        start = 0
        libraries = []
        for i in range(multithreaded_analysis[1]):
            if start >= len(full_library):
                multithreaded_analysis = (True, i)
                break
            temp_library = []
            for j in range(start, start+split_s):
                if j == start+split_s-1 or j == len(full_library)-1:
                    temp_library.append("'"+full_library_keys_list[j]+"'"+": "+str(full_library[full_library_keys_list[j]]))
                    start+=split_s
                    break
                if j < start+split_s-1:
                    temp_library.append("'"+full_library_keys_list[j]+"'"+": "+str(full_library[full_library_keys_list[j]])+", ")
            libraries.append(temp_library)
        mt_path = str(pathlib.Path(__file__).parent.resolve())
        for i_i, i in enumerate(mt_path):
            if i == "\\":
                mt_path = mt_path[:i_i]+"/"+mt_path[i_i+1:]
        with open(mt_path+'/core.py', 'r') as f:
            for i_i, i in enumerate(f):
                for j in range(multithreaded_analysis[1]):
                    with open(save_path+'glycogenius_'+str(j)+'.py', 'a') as g:
                        if i_i == 0:
                            g.write("import importlib\n")
                            g.write("spec1 = importlib.util.spec_from_file_location('Execution_Functions', '"+mt_path+"/Execution_Functions.py')\n")
                            g.write("Execution_Functions = importlib.util.module_from_spec(spec1)\n")
                            g.write("spec1.loader.exec_module(Execution_Functions)\n")
                            continue
                        if i_i == 1:
                            g.write("spec2 = importlib.util.spec_from_file_location('General_Functions', '"+mt_path+"/General_Functions.py')\n")
                            g.write("General_Functions = importlib.util.module_from_spec(spec2)\n")
                            g.write("spec2.loader.exec_module(General_Functions)\n")
                            continue
                        if i[-28:] == "#editted by multithreaded 1\n":
                            g.write("    multithreaded_execution = (True, "+str(j)+", "+str(multithreaded_analysis[1])+")\n")
                            continue
                        if i[-28:] == "#editted by multithreaded 2\n":
                            g.write("        with open('glycogenius_parameters.ini', 'r') as f:\n")
                            continue
                        if i[-28:] == "#editted by multithreaded 3\n":
                            g.write("            for line in f:\n")
                            g.write("                configs+=line\n")
                            continue
                        if i == "#line to add multithreaded library\n":
                            g.write("        full_library = {"+str(libraries[j])[2:-2]+"}\n")
                            continue
                        if i == "#here multithreaded prints execution of main()":
                            g.write("main()")
                            continue
                        g.write(i)
                        g.close()
            f.close()
        print("Multithreaded run setup done. Run the 'glycogenius_n.py' files to execute each part of the script.")
        os._exit(1)
    if imp_exp_library[0]:
        print('Importing existing library...', end = '', flush = True)
        spec = importlib.util.spec_from_file_location("glycans_library", save_path+'glycans_library.py')
        lib_module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(lib_module)
        full_library = lib_module.full_library
        print("Done!")
        return full_library
    if not custom_glycans_list[0]:
        print('Building glycans library...')
        monos_library = Library_Tools.generate_glycans_library(min_max_monos,
                                                 min_max_hex,
                                                 min_max_hexnac,
                                                 min_max_sia,
                                                 min_max_fuc,
                                                 min_max_ac,
                                                 min_max_gc,
                                                 force_nglycan)
        full_library = Library_Tools.full_glycans_library(monos_library,
                                            max_adducts,
                                            max_charges,
                                            tag_mass,
                                            fast_iso,
                                            high_res)
        print('Done building glycans library in '+ 
              str(datetime.datetime.now()-begin_time)+'!')
    if imp_exp_library[1] or only_gen_lib:
        print('Exporting glycans library...')
        with open(save_path+'glycans_library.py', 'w') as f:
            f.write('full_library = '+str(full_library))
            f.close()
        df = {'Glycan' : [], 'Hex' : [], 'HexNAc' : [], 'dHex' : [], 'Neu5Ac' : [], 'Neu5Gc' : [], 'Isotopic Distribution' : [], 'Neutral Mass + Tag' : []}
        for i_i, i in enumerate(full_library):
            df['Glycan'].append(i)
            df['Hex'].append(full_library[i]['Monos_Composition']['H'])
            df['HexNAc'].append(full_library[i]['Monos_Composition']['N'])
            df['dHex'].append(full_library[i]['Monos_Composition']['F'])
            df['Neu5Ac'].append(full_library[i]['Monos_Composition']['S'])
            df['Neu5Gc'].append(full_library[i]['Monos_Composition']['G'])
            temp_isotopic = []
            for j in full_library[i]['Isotopic_Distribution']:
                temp_isotopic.append(float("%.3f" % round(j, 3)))
            df['Isotopic Distribution'].append(str(temp_isotopic)[1:-1])
            df['Neutral Mass + Tag'].append(float("%.4f" % round(full_library[i]['Neutral_Mass+Tag'], 4)))
            for j in full_library[i]['Adducts_mz']:
                if i_i ==0:
                    df[j] = [float("%.4f" % round(full_library[i]['Adducts_mz'][j], 4))]
                else:
                    df[j].append(float("%.4f" % round(full_library[i]['Adducts_mz'][j], 4)))
        df = DataFrame(df)
        with ExcelWriter(save_path+'Glycans_Library.xlsx') as writer:
            df.to_excel(writer, index = False)
    if only_gen_lib:
        print('Library length: '+str(len(full_library)))
        input("Check it in glycans_library.py and Glycans_Library.xlsx, for a readable form. If you wish to analyze files, set 'only_gen_lib' to False and input remaining parameters.\nPress Enter to exit.")
        os._exit(1)
    return full_library
        
def output_filtered_data(curve_fit_score,
                         iso_fit_score,
                         sn,
                         max_ppm,
                         reanalysis,
                         save_path,
                         multithreaded_analysis,
                         analyze_ms2):
    '''
    '''
    date = datetime.datetime.now()
    begin_time = str(date)[2:4]+str(date)[5:7]+str(date)[8:10]+"_"+str(date)[11:13]+str(date)[14:16]+str(date)[17:19]
    if reanalysis[0]:
        print("Reanalyzing raw data with new parameters...")
        results1_list = []
        results2_list = []
        results3_list = []
        results4_list = []
        for i in range(multithreaded_analysis[1]):
            results1_list.append(save_path+'results1_'+str(i))
            results2_list.append(save_path+'results2_'+str(i))
            results3_list.append(save_path+'results3_'+str(i))
            results4_list.append(save_path+'results4_'+str(i))
        try:
            for i_i, i in enumerate(results1_list):
                with open(i, 'rb') as f:
                    file = dill.load(f)
                    df1_import = file[0]
                    df2_import = file[1]
                    if analyze_ms2:
                        fragments_dataframes_import = file[2]
                    f.close()
                with open(results2_list[i_i], 'rb') as f:
                    eic_dataframes_import = dill.load(f)
                    f.close()
                with open(results3_list[i_i], 'rb') as f:
                    smoothed_eic_dataframes_import = dill.load(f)
                    f.close()
                with open(results4_list[i_i], 'rb') as f:
                    curve_fitting_dataframes_import = dill.load(f)
                    f.close()
                if i_i == 0:
                    df1 = df1_import
                    df2 = df2_import
                    eic_dataframes = eic_dataframes_import
                    smoothed_eic_dataframes = smoothed_eic_dataframes_import
                    curve_fitting_dataframes = curve_fitting_dataframes_import
                    if analyze_ms2:
                        fragments_dataframes = fragments_dataframes_import
                else:
                    for j_j, j in enumerate(df1_import):
                        for k in j:
                            df1[j_j][k] = df1[j_j][k] + df1_import[j_j][k]
                    for j_j, j in enumerate(eic_dataframes_import):
                        for k in j:
                            if k[:4] != 'RTs_':
                                eic_dataframes[j_j][k] = eic_dataframes_import[j_j][k]
                                smoothed_eic_dataframes[j_j][k] = smoothed_eic_dataframes_import[j_j][k]
                    for j_j, j in enumerate(curve_fitting_dataframes_import):
                        for k in j:
                            curve_fitting_dataframes[j_j][k] = curve_fitting_dataframes_import[j_j][k]
                    if analyze_ms2:
                        for j_j, j in enumerate(fragments_dataframes_import):
                            for k in j:
                                fragments_dataframes[j_j][k] = fragments_dataframes[j_j][k] + fragments_dataframes_import[j_j][k]
            with open(save_path+'raw_data_1', 'wb') as f:
                if analyze_ms2:
                    dill.dump([df1, df2, fragments_dataframes], f)
                else:
                    dill.dump([df1, df2], f)
                    f.close()
            with open(save_path+'raw_data_2', 'wb') as f:
                dill.dump(eic_dataframes, f)
                f.close()
            with open(save_path+'raw_data_3', 'wb') as f:
                dill.dump(smoothed_eic_dataframes, f)
                f.close()
            with open(save_path+'raw_data_4', 'wb') as f:
                dill.dump(curve_fitting_dataframes, f)
                f.close()
            for i in range(multithreaded_analysis[1]):
                p1 = pathlib.Path(save_path+'Multithreaded_'+str(i)+'.py')
                p2 = pathlib.Path(save_path+'glycans_library_'+str(i)+'.py')
                p3 = pathlib.Path(save_path+'results1_'+str(i))
                p4 = pathlib.Path(save_path+'results2_'+str(i))
                p5 = pathlib.Path(save_path+'results3_'+str(i))
                p6 = pathlib.Path(save_path+'results4_'+str(i))
                p1.unlink(missing_ok=True)
                p2.unlink(missing_ok=True)
                p3.unlink(missing_ok=True)
                p4.unlink(missing_ok=True)
                p5.unlink(missing_ok=True)
                p6.unlink(missing_ok=True)
        except Exception:
            pass
    else:
        print("Analyzing raw data...")
    try:
        with open(save_path+'raw_data_1', 'rb') as f:
            file = dill.load(f)
            df1 = file[0]
            df2 = file[1]
            if analyze_ms2:
                fragments_dataframes = file[2]
            f.close()
    except:
        return
    for i_i, i in enumerate(df2["Sample_Number"]): #QCs cutoff
        for j_j, j in enumerate(df1[i_i]["Adduct"]): 
            temp_rt = df1[i_i]["RT"][j_j]
            temp_auc = df1[i_i]["AUC"][j_j]
            temp_ppm = df1[i_i]["PPM"][j_j]
            temp_sn = df1[i_i]["S/N"][j_j]
            temp_fit = df1[i_i]["Iso_Fitting_Score"][j_j]
            temp_curve = df1[i_i]["Curve_Fitting_Score"][j_j]
            to_remove = []
            to_remove_glycan = []
            to_remove_adduct = []
            for k_k, k in enumerate(temp_sn):
                if float(k) < sn:
                    to_remove.append(k_k)
                    to_remove_glycan.append(df1[i_i]["Glycan"][j_j])
                    to_remove_adduct.append(j)
                    continue
                if float(temp_fit[k_k]) < iso_fit_score:
                    to_remove.append(k_k)
                    to_remove_glycan.append(df1[i_i]["Glycan"][j_j])
                    to_remove_adduct.append(j)
                    continue
                if float(temp_curve[k_k]) < curve_fit_score:
                    to_remove.append(k_k)
                    to_remove_glycan.append(df1[i_i]["Glycan"][j_j])
                    to_remove_adduct.append(j)
                    continue
                if abs(float(temp_ppm[k_k])) > max_ppm:
                    to_remove.append(k_k)
                    to_remove_glycan.append(df1[i_i]["Glycan"][j_j])
                    to_remove_adduct.append(j)
                    continue
            if len(to_remove) != 0:
                to_remove.reverse()
                to_remove_glycan.reverse()
                to_remove_adduct.reverse()
                for k_k in to_remove:
                    del temp_rt[k_k]
                    del temp_auc[k_k]
                    del temp_ppm[k_k]
                    del temp_sn[k_k]
                    del temp_fit[k_k]
                    del temp_curve[k_k]
                    if analyze_ms2:
                        if len(temp_rt) == 0:
                            df1[i_i]["Detected_Fragments"][j_j] = ""
                            for k in range(len(fragments_dataframes[i_i]["Glycan"])-1, -1, -1):
                                if (fragments_dataframes[i_i]["Glycan"][k] == to_remove_glycan[k_k] 
                                    and fragments_dataframes[i_i]["Adduct"][k] == to_remove_adduct[k_k]):
                                    del fragments_dataframes[i_i]["Glycan"][k]
                                    del fragments_dataframes[i_i]["Adduct"][k]
                                    del fragments_dataframes[i_i]["Fragment"][k]
                                    del fragments_dataframes[i_i]["Fragment_mz"][k]
                                    del fragments_dataframes[i_i]["Fragment_Intensity"][k]
                                    del fragments_dataframes[i_i]["RT"][k]
                                    del fragments_dataframes[i_i]["Precursor_mz"][k]
            df1[i_i]["RT"][j_j] = str(temp_rt)[1:-1]
            df1[i_i]["AUC"][j_j] = str(temp_auc)[1:-1]
            df1[i_i]["PPM"][j_j] = str(temp_ppm)[1:-1]
            df1[i_i]["S/N"][j_j] = str(temp_sn)[1:-1]
            df1[i_i]["Iso_Fitting_Score"][j_j] = str(temp_fit)[1:-1]
            df1[i_i]["Curve_Fitting_Score"][j_j] = str(temp_curve)[1:-1]
    to_remove = []
    to_remove_glycan = []
    to_remove_adduct = []
    for i_i, i in enumerate(df1[0]["Adduct"]):
        glycan_good = False
        for j_j, j in enumerate(df2["Sample_Number"]):
            check = df1[j_j]["RT"][i_i]
            if len(check) != 0:
                glycan_good = True
                break
        if not glycan_good:
            to_remove.append(i_i)
            to_remove_glycan.append(df1[j_j]["Glycan"][i_i])
            to_remove_adduct.append(i)
    if len(to_remove) != 0:
        to_remove.reverse()
        to_remove_glycan.reverse()
        to_remove_adduct.reverse()
        for j_j, j in enumerate(df2["Sample_Number"]):
            for i_index, i_i in enumerate(to_remove):
                del df1[j_j]["Glycan"][i_i]
                del df1[j_j]["Adduct"][i_i]
                del df1[j_j]["mz"][i_i]
                del df1[j_j]["RT"][i_i]
                del df1[j_j]["AUC"][i_i]
                del df1[j_j]["PPM"][i_i]
                del df1[j_j]["S/N"][i_i]
                del df1[j_j]["Iso_Fitting_Score"][i_i]
                del df1[j_j]["Curve_Fitting_Score"][i_i]
                if analyze_ms2:
                    del df1[j_j]["Detected_Fragments"][i_i]
                    for k in range(len(fragments_dataframes[j_j]["Glycan"])-1, -1, -1):
                        if (fragments_dataframes[j_j]["Glycan"][k] == to_remove_glycan[i_index] 
                            and fragments_dataframes[j_j]["Adduct"][k] == to_remove_adduct[i_index]):
                            del fragments_dataframes[j_j]["Glycan"][k]
                            del fragments_dataframes[j_j]["Adduct"][k]
                            del fragments_dataframes[j_j]["Fragment"][k]
                            del fragments_dataframes[j_j]["Fragment_mz"][k]
                            del fragments_dataframes[j_j]["Fragment_Intensity"][k]
                            del fragments_dataframes[j_j]["RT"][k]
                            del fragments_dataframes[j_j]["Precursor_mz"][k]
    if analyze_ms2:
        for j_j, j in enumerate(df2["Sample_Number"]):
            for k_k, k in enumerate(df1[j_j]["RT"]):
                splitted = k.split(", ")
                temp_rts = []
                for l in splitted:
                    if len(l) > 0:
                        temp_rts.append(float(l))
                if len(temp_rts) > 0:
                    count = 0
                    for m in temp_rts:
                        for l in range(len(fragments_dataframes[j_j]["Glycan"])-1, -1, -1):
                            if fragments_dataframes[j_j]["Glycan"][l] == df1[j_j]["Glycan"][k_k] and fragments_dataframes[j_j]["Adduct"][l] == df1[j_j]["Adduct"][k_k]:
                                count += 1 
                                if abs(fragments_dataframes[j_j]["RT"][l] - m) > 0.1:
                                    count -= 1
                                    del fragments_dataframes[j_j]["Glycan"][l]
                                    del fragments_dataframes[j_j]["Adduct"][l]
                                    del fragments_dataframes[j_j]["Fragment"][l]
                                    del fragments_dataframes[j_j]["Fragment_mz"][l]
                                    del fragments_dataframes[j_j]["Fragment_Intensity"][l]
                                    del fragments_dataframes[j_j]["RT"][l]
                                    del fragments_dataframes[j_j]["Precursor_mz"][l]
                    if count == 0:
                        df1[j_j]["Detected_Fragments"][k_k] = "No" #QCs cutoff end
    df2 = DataFrame(df2)
    with ExcelWriter(save_path+begin_time+'_Results_'+str(max_ppm)+'_'+str(iso_fit_score)+'_'+str(curve_fit_score)+'_'+str(sn)+'.xlsx') as writer:
        print("Creating results file...", end="", flush=True)
        for i_i, i in enumerate(df1):
            result_df = DataFrame(i)
            result_df.to_excel(writer, sheet_name="Sample_"+str(i_i), index = False)
            if analyze_ms2:
                fragments_df = DataFrame(fragments_dataframes[i_i])
                fragments_df.to_excel(writer, sheet_name="Sample_"+str(i_i)+"_Fragments", index = False)
        df2.to_excel(writer, sheet_name="Index references", index = False)
    print("Done!")
    if (reanalysis[1] and reanalysis[0]) or not reanalysis[0]:
        with open(save_path+'raw_data_2', 'rb') as f:
            eic_dataframes = dill.load(f)
            f.close()
        with open(save_path+'raw_data_3', 'rb') as f:
            smoothed_eic_dataframes = dill.load(f)
            f.close()
        with open(save_path+'raw_data_4', 'rb') as f:
            curve_fitting_dataframes = dill.load(f)
            f.close()
        with ExcelWriter(save_path+begin_time+'_Plot_Data.xlsx') as writer:
            print("Creating data plotting files...", end= "", flush=True)
            for i_i, i in enumerate(eic_dataframes):
                eic_df = DataFrame(i)
                eic_df.to_excel(writer, sheet_name="Sample_"+str(i_i)+"_EICs", index = False)
                smoothed_eic_df = DataFrame(smoothed_eic_dataframes[i_i])
                smoothed_eic_df.to_excel(writer, sheet_name="Sample_"+str(i_i)+"_Smoothed_EICs", index = False)
                if len(curve_fitting_dataframes[i_i]) > 16384:
                    for j in range(int(len(curve_fitting_dataframes[i_i])/16384)+1):
                        if j == 0:
                            curve_df = DataFrame(dict(islice(curve_fitting_dataframes[i_i].items(), 16384)))
                            curve_df.to_excel(writer, sheet_name="Sample_"+str(i_i)+"_Curve_Fits_0", index = False)
                        else:
                            if len(dict(islice(curve_fitting_dataframes[i_i].items(), j*16384, len(curve_fitting_dataframes[i_i])))) <= 16384:
                                curve_df = DataFrame(dict(islice(curve_fitting_dataframes[i_i].items(), j*16384, len(curve_fitting_dataframes[i_i]))))
                                curve_df.to_excel(writer, sheet_name="Sample_"+str(i_i)+"_Curve_Fits_"+str(j), index = False)
                            else:
                                curve_df = DataFrame(dict(islice(curve_fitting_dataframes[i_i].items(), j*16384, (j+1)*16384)))
                                curve_df.to_excel(writer, sheet_name="Sample_"+str(i_i)+"_Curve_Fits_"+str(j), index = False)
                else:
                    curve_df = DataFrame(curve_fitting_dataframes[i_i])
                    curve_df.to_excel(writer, sheet_name="Sample_"+str(i_i)+"_Curve_Fits", index = False)
            df2.to_excel(writer, sheet_name="Index references", index = False)
        print("Done!")

def arrange_raw_data(analyzed_data,
                     samples_names,
                     multithreaded_analysis,
                     multithreaded_execution,
                     analyze_ms2,
                     save_path): ##Complete
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
    print('Arranging raw data...', end='', flush = True)
    df1 = []
    df2 = {"Sample_Number" : [], "File_name" : [], "Noise Level" : []}
    eic_dataframes = []
    smoothed_eic_dataframes = []
    curve_fitting_dataframes = []
    if analyze_ms2:
        fragments_dataframes = []
    for i_i, i in enumerate(samples_names):
        eic_dataframes.append({})
        temp_eic_rt = []
        for j in analyzed_data[1][i_i]:
            temp_eic_rt.append(float("%.4f" % round(j, 4)))
        eic_dataframes[i_i]['RTs_'+str(i_i)] = temp_eic_rt
        smoothed_eic_dataframes.append({})
        smoothed_eic_dataframes[i_i]['RTs_'+str(i_i)] = temp_eic_rt
        curve_fitting_dataframes.append({})
        df2["Sample_Number"].append(i_i)
        df2["File_name"].append(i)
        df2["Noise Level"].append(float("%.1f" % round(analyzed_data[2][i_i],1)))
        if analyze_ms2:
            df1.append({"Glycan" : [], "Adduct" : [], "mz" : [], "RT" : [], "AUC" : [], "PPM" : [], "S/N" : [], "Iso_Fitting_Score" : [], "Curve_Fitting_Score" : [], "Detected_Fragments" : []})
            fragments_dataframes.append({"Glycan" : [], "Adduct" : [], "Fragment" : [], "Fragment_mz" : [], "Fragment_Intensity" : [], "RT" : [], "Precursor_mz" : []})
        else:
            df1.append({"Glycan" : [], "Adduct" : [], "mz" : [], "RT" : [], "AUC" : [], "PPM" : [], "S/N" : [], "Iso_Fitting_Score" : [], "Curve_Fitting_Score" : []})
    for i_i, i in enumerate(analyzed_data[0]): #i = glycan (key)
        for j_j, j in enumerate(analyzed_data[0][i]['Adducts_mz_data']): #j = adduct (key)
            for k_k, k in enumerate(analyzed_data[0][i]['Adducts_mz_data'][j]):
                temp_eic_int = []
                for l in analyzed_data[0][i]['Adducts_mz_data'][j][k][0]:
                    temp_eic_int.append(int(l))
                eic_dataframes[k_k][str(i)+'+'+str(j)+' - '+str(float("%.4f" % round(analyzed_data[0][i]['Adducts_mz'][j], 4)))] = temp_eic_int
                temp_eic_int = []
                for l in analyzed_data[0][i]['Adducts_mz_data'][j][k][2]:
                    temp_eic_int.append(int(l))
                smoothed_eic_dataframes[k_k][str(i)+'+'+str(j)+' - '+str(float("%.4f" % round(analyzed_data[0][i]['Adducts_mz'][j], 4)))] = temp_eic_int
            found = False
            for k_k, k in enumerate(analyzed_data[0][i]['Adducts_mz_data'][j]):
                if len(analyzed_data[0][i]['Adducts_mz_data'][j][k][1]) != 0:
                    found = True
            if not found:
                continue
            for k_k, k in enumerate(analyzed_data[0][i]['Adducts_mz_data'][j]): #k = sample (key)
                df1[k_k]["Glycan"].append(i)
                df1[k_k]["Adduct"].append(j)
                df1[k_k]["mz"].append(float("%.4f" % round(analyzed_data[0][i]['Adducts_mz'][j], 4)))
                temp_rts = []
                temp_aucs = []
                temp_ppm = []
                temp_s_n = []
                temp_iso_score = []
                temp_curve_score = []
                temp_curve_data_total = []
                for l_l, l in enumerate(analyzed_data[0][i]['Adducts_mz_data'][j][k][1]):
                    temp_rts.append(float("%.2f" % round(l['rt'], 2)))
                    temp_aucs.append(float("%.2f" % round(l['AUC'], 2)))
                    temp_ppm.append(float("%.2f" % round(l['Average_PPM'][0], 2)))
                    temp_s_n.append(float("%.1f" % round(l['Signal-to-Noise'], 1)))
                    if isnan(l['Iso_Fit_Score']):
                        temp_iso_score.append(0.0)
                    else:
                        temp_iso_score.append(float("%.4f" % round(l['Iso_Fit_Score'], 4)))
                    if isnan(l['Curve_Fit_Score'][0]):
                        temp_curve_score.append(0.0)
                    else:
                        temp_curve_score.append(float("%.4f" % round(l['Curve_Fit_Score'][0], 4)))
                    temp_curve_data_rt = []
                    temp_curve_data_actual = []
                    temp_curve_data_ideal = []
                    for m_m, m in enumerate(l['Curve_Fit_Score'][1]):
                        temp_curve_data_rt.append(m)
                        temp_curve_data_actual.append(l['Curve_Fit_Score'][2][m_m])
                        temp_curve_data_ideal.append(l['Curve_Fit_Score'][3][m_m])
                    temp_curve_data_total.append((temp_curve_data_rt, temp_curve_data_actual, temp_curve_data_ideal))
                if analyze_ms2:
                    temp_fragments = analyzed_data[3][i][j][k_k]
                if len(temp_rts) == 0:
                    df1[k_k]["RT"].append([0.0])
                    df1[k_k]["AUC"].append([0.0])
                    df1[k_k]["PPM"].append([0.0])
                    df1[k_k]["S/N"].append([0.0])
                    df1[k_k]["Iso_Fitting_Score"].append([0.0])
                    df1[k_k]["Curve_Fitting_Score"].append([0.0])
                    if analyze_ms2:
                        df1[k_k]["Detected_Fragments"].append('Glycan+Adduct not found in sample')
                        temp_fragments = []
                else:
                    df1[k_k]["RT"].append(temp_rts)
                    df1[k_k]["AUC"].append(temp_aucs)
                    df1[k_k]["PPM"].append(temp_ppm)
                    df1[k_k]["S/N"].append(temp_s_n)
                    df1[k_k]["Iso_Fitting_Score"].append(temp_iso_score)
                    df1[k_k]["Curve_Fitting_Score"].append(temp_curve_score)
                    if analyze_ms2:
                        if len(temp_fragments) != 0:
                            for m in temp_fragments:
                                fragments_dataframes[k_k]["Glycan"].append(m[0])
                                fragments_dataframes[k_k]["Adduct"].append(m[1])
                                fragments_dataframes[k_k]["Fragment"].append(m[2])
                                fragments_dataframes[k_k]["Fragment_mz"].append(float("%.4f" % round(m[3], 4)))
                                fragments_dataframes[k_k]["Fragment_Intensity"].append(float("%.2f" % round(m[4], 2)))
                                fragments_dataframes[k_k]["RT"].append(float("%.2f" % round(m[5],2)))
                                fragments_dataframes[k_k]["Precursor_mz"].append(float("%.4f" % round(m[6], 4)))
                            df1[k_k]["Detected_Fragments"].append('Yes')
                        else:
                            df1[k_k]["Detected_Fragments"].append('No')
                    for m_m, m in enumerate(temp_rts):
                        temp_array = []
                        for n in temp_curve_data_total[m_m][0]:
                            temp_array.append(float("%.4f" % round(n, 4)))
                        curve_fitting_dataframes[k_k][str(i)+"+"+str(j)+"_"+str(m)+"_RTs"] = temp_array
                        temp_array = []
                        for n in temp_curve_data_total[m_m][1]:
                            temp_array.append(int(n))
                        curve_fitting_dataframes[k_k][str(i)+"+"+str(j)+"_"+str(m)+"_Found_ints"] = temp_array
                        temp_array = []
                        for n in temp_curve_data_total[m_m][2]:
                            temp_array.append(int(n))
                        curve_fitting_dataframes[k_k][str(i)+"+"+str(j)+"_"+str(m)+"_Ideal_ints"] = temp_array
    biggest_len = 1000
    for i in curve_fitting_dataframes:
        for j in i:
            if len(i[j]) < biggest_len:
                for k in range(biggest_len-len(i[j])):
                    i[j].append(None)
    if multithreaded_execution[0]:
        sleep(multithreaded_execution[2]+multithreaded_execution[1])
        with open(save_path+'results1_'+str(multithreaded_execution[1]), 'wb') as f:
            if analyze_ms2:
                dill.dump([df1, df2, fragments_dataframes], f)
            else:
                dill.dump([df1, df2], f)
            f.close()
        with open(save_path+'results2_'+str(multithreaded_execution[1]), 'wb') as f:
            dill.dump(eic_dataframes, f)
            f.close()
        with open(save_path+'results3_'+str(multithreaded_execution[1]), 'wb') as f:
            dill.dump(smoothed_eic_dataframes, f)
            f.close()
        with open(save_path+'results4_'+str(multithreaded_execution[1]), 'wb') as f:
            dill.dump(curve_fitting_dataframes, f)
            f.close()
        p1 = pathlib.Path(save_path+'glycogenius_'+str(multithreaded_execution[1])+'.py')
        p1.unlink(missing_ok=True)
        results1_list = []
        results2_list = []
        results3_list = []
        results4_list = []
        for i in range(multithreaded_execution[2]):
            results1_list.append(save_path+'results1_'+str(i))
            results2_list.append(save_path+'results2_'+str(i))
            results3_list.append(save_path+'results3_'+str(i))
            results4_list.append(save_path+'results4_'+str(i))
        try:
            for i_i, i in enumerate(results1_list):
                with open(i, 'rb') as f:
                    file = dill.load(f)
                    df1_import = file[0]
                    df2_import = file[1]
                    if analyze_ms2:
                        fragments_dataframes_import = file[2]
                    f.close()
                with open(results2_list[i_i], 'rb') as f:
                    eic_dataframes_import = dill.load(f)
                    f.close()
                with open(results3_list[i_i], 'rb') as f:
                    smoothed_eic_dataframes_import = dill.load(f)
                    f.close()
                with open(results4_list[i_i], 'rb') as f:
                    curve_fitting_dataframes_import = dill.load(f)
                    f.close()
                if i_i == 0:
                    df1 = df1_import
                    df2 = df2_import
                    eic_dataframes = eic_dataframes_import
                    smoothed_eic_dataframes = smoothed_eic_dataframes_import
                    curve_fitting_dataframes = curve_fitting_dataframes_import
                    if analyze_ms2:
                        fragments_dataframes = fragments_dataframes_import
                else:
                    for j_j, j in enumerate(df1_import):
                        for k in j:
                            df1[j_j][k] = df1[j_j][k] + df1_import[j_j][k]
                    for j_j, j in enumerate(eic_dataframes_import):
                        for k in j:
                            if k[:4] != 'RTs_':
                                eic_dataframes[j_j][k] = eic_dataframes_import[j_j][k]
                                smoothed_eic_dataframes[j_j][k] = smoothed_eic_dataframes_import[j_j][k]
                    for j_j, j in enumerate(curve_fitting_dataframes_import):
                        for k in j:
                            curve_fitting_dataframes[j_j][k] = curve_fitting_dataframes_import[j_j][k]
                    if analyze_ms2:
                        for j_j, j in enumerate(fragments_dataframes_import):
                            for k in j:
                                fragments_dataframes[j_j][k] = fragments_dataframes[j_j][k] + fragments_dataframes_import[j_j][k]
            with open(save_path+'raw_data_1', 'wb') as f:
                if analyze_ms2:
                    dill.dump([df1, df2, fragments_dataframes], f)
                else:
                    dill.dump([df1, df2], f)
                    f.close()
            with open(save_path+'raw_data_2', 'wb') as f:
                dill.dump(eic_dataframes, f)
                f.close()
            with open(save_path+'raw_data_3', 'wb') as f:
                dill.dump(smoothed_eic_dataframes, f)
                f.close()
            with open(save_path+'raw_data_4', 'wb') as f:
                dill.dump(curve_fitting_dataframes, f)
                f.close()
            for i in range(multithreaded_execution[2]):
                p3 = pathlib.Path(save_path+'results1_'+str(i))
                p4 = pathlib.Path(save_path+'results2_'+str(i))
                p5 = pathlib.Path(save_path+'results3_'+str(i))
                p6 = pathlib.Path(save_path+'results4_'+str(i))
                p3.unlink(missing_ok=True)
                p4.unlink(missing_ok=True)
                p5.unlink(missing_ok=True)
                p6.unlink(missing_ok=True)
        except Exception:
            pass
    else:
        with open(save_path+'raw_data_1', 'wb') as f:
            if analyze_ms2:
                dill.dump([df1, df2, fragments_dataframes], f)
            else:
                dill.dump([df1, df2], f)
                f.close()
        with open(save_path+'raw_data_2', 'wb') as f:
            dill.dump(eic_dataframes, f)
            f.close()
        with open(save_path+'raw_data_3', 'wb') as f:
            dill.dump(smoothed_eic_dataframes, f)
            f.close()
        with open(save_path+'raw_data_4', 'wb') as f:
            dill.dump(curve_fitting_dataframes, f)
            f.close()
    print("Done!")

def print_sep(): ##Complete
    '''Prints a separator consisting of multiple '-' character.
    '''
    print('-------------------------------------------------')
    
def analyze_files(library,
                  lib_size,
                  data,
                  ms1_index,
                  tolerance,
                  ret_time_interval,
                  min_isotops,
                  min_ppp,
                  max_charges,
                  custom_noise,
                  close_peaks,
                  verbose = False): ##Complete
    '''Integrates all the file-accessing associated functions in this script to go
    through the files data, draw eic of hypothetical glycans, check if it is not in
    noise level, based on threshold, does peak-picking, check quality of spectra in
    peaks, deconvolute if checks passed, and calculates AUC of the peaks.

    Parameters
    ----------
    No parameters needed, but must be executed after parameters section of script.

    Uses
    ----
    File_Accessing.eic_from_glycan() : dict
        Generates an EIC of the mzs calculated for a glycan.

    File_Accessing.peaks_from_eic() : list, tuple
        Does multi-peak-picking in an EIC, based on the most intense peaks in the EIC,
        in the retention time interval of interest. This function works in-line.

    check_monoisotopic_charge() : list, list
        Checks if the peaks identified for a given mz corresponds to monoisotopic,
        correctly charge-assigned peaks.

    deconvoluted_glycan_eic() : list
        Deconvolutes the data around the identified peaks retention time and creates a
        deconvoluted EIC.

    File_Accessing.peaks_auc_from_eic() : list
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
    analyzed_data = {}
    selected_peaks_width = []
    mean_sel_peaks = 0.0
    noise = {}
    rt_array_report = {}
    if not custom_noise[0]:
        print('Analyzing noise level of samples... ', end='', flush = True)
    for i_i, i in enumerate(data):
        rt_array_report[i_i] = []
        temp_noise = []
        if custom_noise[0]:
            temp_noise.append(custom_noise[1][i_i])
            continue
        noise[i_i] = 0
        for j_j, j in enumerate(ms1_index[i_i]):
            rt_array_report[i_i].append(i[j]['retentionTime'])
            if (i[j]['retentionTime'] < ret_time_interval[0] or i[j]['retentionTime'] > ret_time_interval[1]):
                continue
            mz_ints = dict(zip(i[j]['m/z array'], i[j]['intensity array']))
            noise_level = General_Functions.noise_level_calc_mzarray(mz_ints)
            temp_noise.append(noise_level)
        noise[i_i] = percentile(temp_noise, 66.8)
    print('Done!')
    print_sep()
    if verbose:
        print('Noise Level: '+str(noise))
        noise_report = 'Noise Level: '+str(noise)
    print("Analyzing glycans in samples' MS1 spectra...")
    for i_i, i in enumerate(library):
        verbose_info = []
        print('Analyzing glycan '+str(i)+': '+str(i_i+1)+'/'+str(lib_size))
        glycan_data = library[i]
        if verbose:
            verbose_info.append(noise_report)
            verbose_info.append(str(glycan_data))
            verbose_info.append('EIC RT Array: '+str(rt_array_report))
        glycan_data['Adducts_mz_data'] = {}
        temp_eic = File_Accessing.eic_from_glycan(data,
                                   glycan_data,
                                   ms1_index,
                                   ret_time_interval,
                                   tolerance,
                                   min_isotops,
                                   noise,
                                   max_charges,
                                   verbose)
        if verbose:
            with open('log_'+str(i)+'_eic_debug.txt', 'w') as f:
                for j in temp_eic[3]:
                    f.write(j+"\n")
                f.close()
        for j in temp_eic[0]: #moving through adducts
            glycan_data['Adducts_mz_data'][j] = {}
            for k in temp_eic[0][j]: #moving through samples
                if verbose:
                    verbose_info.append('Adduct: '+str(j)+' mz: '+str(glycan_data['Adducts_mz'][j])+' Sample: '+str(k))
                    verbose_info.append('EIC INT Array: '+str(temp_eic[0][j][k][1]))
                temp_eic_smoothed = File_Accessing.eic_smoothing(temp_eic[0][j][k])
                glycan_data['Adducts_mz_data'][j][k] = []
                glycan_data['Adducts_mz_data'][j][k].append(temp_eic[0][j][k][1])
                glycan_data['Adducts_mz_data'][j][k].append([])
                glycan_data['Adducts_mz_data'][j][k].append(temp_eic_smoothed[1])
                if verbose:
                    verbose_info.append('Smoothed EIC INT Array: '+str(temp_eic_smoothed[1]))
                    verbose_info.append('Raw PPM error: '+str(temp_eic[1][j][k])+'\nRaw Isotopic Fitting Score: '+str(temp_eic[2][j][k]))
                if max(temp_eic[0][j][k][1]) < noise[k]:
                    continue
                temp_peaks = File_Accessing.peaks_from_eic(temp_eic[0][j][k],
                                            temp_eic_smoothed,
                                            ret_time_interval,
                                            min_ppp,
                                            close_peaks)
                if len(temp_peaks) == 0:
                    continue
                temp_peaks_auc = File_Accessing.peaks_auc_from_eic(temp_eic[0][j][k],
                                                    ms1_index[k],
                                                    temp_peaks)
                for l_l, l in enumerate(temp_peaks):
                    l['AUC'] = temp_peaks_auc[l_l]
                    l['Average_PPM'] = File_Accessing.average_ppm_calc(temp_eic[1][j][k], tolerance, l)
                    l['Iso_Fit_Score'] = File_Accessing.iso_fit_score_calc(temp_eic[2][j][k], l)
                    l['Signal-to-Noise'] = l['int']/noise[k]
                    l['Curve_Fit_Score'] = File_Accessing.peak_curve_fit(temp_eic_smoothed, l)
                    glycan_data['Adducts_mz_data'][j][k][1].append(l)
                if verbose:
                    verbose_info.append('Adduct: '+str(j)+' mz: '+str(glycan_data['Adducts_mz'][j])+' Sample: '+str(k))
                    verbose_info.append('Peaks found: '+str(temp_peaks))
                    verbose_info.append('Peaks AUC: '+str(temp_peaks_auc))
        analyzed_data[i] = glycan_data
        if verbose:
            with open('log_'+str(i)+'_execution_debug.txt', 'w') as f:
                for i in verbose_info:
                    f.write(i+"\n")
                    f.write("\n")
                f.close()
    print('Sample MS1 analysis done in '+str(datetime.datetime.now() - begin_time)+'!')
    return analyzed_data, rt_array_report, noise
    
def analyze_ms2(ms2_index, 
                data, 
                analyzed_data, 
                rt_interval,
                tolerance,
                min_max_monos,
                min_max_hex,
                min_max_hexnac,
                min_max_sia,
                min_max_fuc,
                min_max_ac,
                min_max_gc,
                max_charges,
                tag_mass,
                filter_output):
    '''
    '''
    begin_time = datetime.datetime.now()
    print('Analyzing MS2 data...')
    fragments = Library_Tools.fragments_library(min_max_monos,
                                  min_max_hex,
                                  min_max_hexnac,
                                  min_max_sia,
                                  min_max_fuc,
                                  min_max_ac,
                                  min_max_gc,
                                  max_charges,
                                  tolerance,
                                  tag_mass)
    fragments_data = {}
    print('Scanning MS2 spectra...')
    for i_i, i in enumerate(analyzed_data[0]):
        fragments_data[i] = {}
        for j_j, j in enumerate(analyzed_data[0][i]['Adducts_mz_data']):
            fragments_data[i][j] = {}
            for k_k, k in enumerate(data):
                fragments_data[i][j][k_k] = []
                if len(analyzed_data[0][i]['Adducts_mz_data'][j][k_k][1]) == 0:
                    continue
                for l in ms2_index[k_k]:
                    if k[l]['retentionTime'] < analyzed_data[0][i]['Adducts_mz_data'][j][k_k][1][0]['peak_interval'][0] or k[l]['retentionTime'] > analyzed_data[0][i]['Adducts_mz_data'][j][k_k][1][-1]['peak_interval'][1]:
                        continue
                    if abs((k[l]['precursorMz'][0]['precursorMz']) - analyzed_data[0][i]['Adducts_mz'][j]) <= (1.0074/General_Functions.form_to_charge(j))+tolerance:
                        for m_m, m in enumerate(k[l]['m/z array']):
                            found = False
                            for n_n, n in enumerate(fragments):
                                if filter_output and ("/" not in n['Formula'] 
                                    and (n['Monos_Composition']['H'] > analyzed_data[0][i]['Monos_Composition']['H']
                                    or n['Monos_Composition']['N'] > analyzed_data[0][i]['Monos_Composition']['N']
                                    or n['Monos_Composition']['S'] > analyzed_data[0][i]['Monos_Composition']['S']
                                    or n['Monos_Composition']['F'] > analyzed_data[0][i]['Monos_Composition']['F']
                                    or n['Monos_Composition']['G'] > analyzed_data[0][i]['Monos_Composition']['G'])):
                                    continue
                                if found:
                                    break
                                for o in n['Adducts_mz']:
                                    if abs(n['Adducts_mz'][o]-m) <= tolerance: #fragments data outputted in the form of (Glycan, Adduct, Fragment, Fragment mz, intensity, retention time, precursor)
                                        if len(n['Formula']) > 3 and n['Formula'][-3] != "_":
                                            fragments_data[i][j][k_k].append((i, j, n['Formula']+'_'+o, n['Adducts_mz'][o], k[l]['intensity array'][m_m], k[l]['retentionTime'], k[l]['precursorMz'][0]['precursorMz']))
                                        else:
                                            fragments_data[i][j][k_k].append((i, j, n['Formula'], n['Adducts_mz'][o], k[l]['intensity array'][m_m], k[l]['retentionTime'], k[l]['precursorMz'][0]['precursorMz']))
                                        found = True
                                        break
    print('Sample MS2 analysis done in '+str(datetime.datetime.now() - begin_time)+'!')
    return analyzed_data[0], analyzed_data[1], analyzed_data[2], fragments_data

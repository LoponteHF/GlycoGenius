# GlycoGenius: Glycomics Data Analysis Tool
# Copyright (C) 2023 by Hector Franco Loponte
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or 
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. It is accessible within the program files
# or by typing 'license' after running it stand-alone in the terminal
# by typing 'glycogenius'. If not, see <https://www.gnu.org/licenses/>.

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
from itertools import combinations_with_replacement, islice, product
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
import pkg_resources
import platform
import copy

try:
    version = pkg_resources.get_distribution("glycogenius").version
except:
    version_path = str(pathlib.Path(__file__).parent.parent.parent.resolve())
    with open(version_path+"/Setup.py", "r") as f: #grabs version from setup.py to add to raw_data files
        for lines in f:
            if lines[:12] == "    version=":
                version = lines[13:-2].strip("'")

##---------------------------------------------------------------------------------------
##Functions to be used for execution and organizing results data

def print_header(complete = True):
    '''Prints a default header to be used in CLI.
    
    Parameters
    ----------
    complete : boolean
        If set to True, produces more complete information used
        when glycogenius is executed stand-alone from the terminal.
    '''
    print("\n    GlycoGenius: Glycomics Data Analysis Tool")
    print("   Copyright (C) 2023 by Hector Franco Loponte")
    if not complete:
        print("\n   For more details about the license, run the")
        print("   package stand-alone by typing 'glycogenius'")
        print("     in terminal and then typing 'license'.")
    if complete:
        print("\n This program comes with ABSOLUTELY NO WARRANTY;")
        print("          for details type 'warranty'.")
        print(" This is free software and can be redistributed")
        print("  under certain conditions. If you want to know")
        print(" more details about the license, type 'license'.")
    print_sep()
    
def generate_cfg_file(path, comments):
    '''Generates a parameter file based on parameters set.
    
    Parameters
    ----------
    path : string
        The path to where the file will be saved.
        
    comments : boolean
        Whether the output file will contain comments or not.
    '''
    print("Creating settings file...")
    glycogenius_path = str(pathlib.Path(__file__).parent.parent.resolve())
    for i_i, i in enumerate(glycogenius_path):
        if i == "\\":
            glycogenius_path = glycogenius_path[:i_i]+"/"+glycogenius_path[i_i+1:]
    with open(path+'glycogenius_parameters.ini', 'w') as g:
        with open(glycogenius_path+'/Parameters_Template.py', 'r') as f:
            for line in f:
                if line[0:14] == "samples_path =":
                    g.write("samples_path = "+path+"Sample Files/\n")
                    continue
                if line[0:14] == "working_path =":
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
    '''This function generates the CLI for user interaction.
    
    Uses
    ----
    pathlib.Path.resolve() : Path object
        Make the path absolute, resolving any symlinks. A new path object is returned
        
    General_Functions.form_to_comp(string) : dict
        Separates a molecular formula or monosaccharides composition of glycans into a
        dictionary with each atom/monosaccharide as a key and its amount as value.
        
    Returns
    -------
    input_order : list
        A list containing all the options chosen in the CLI.
        
    glycans_list : list
        A list of custom glycans to be used by the library generating function.
        
    lib_settings : list
        A list containing all the parameters needed to generate a glycans library.
        
    adducts : dict
        A dictionary with maximum adducts information, to be used for library generation.
        
    max_charges : int
        The maximum amount of charges, to be used for library generation.
        
    tag_mass : float
        The added mass of the tag, to be used for library generation.
        
    fast_iso : boolean
        Whether the isotopic distribution generation should be done quick or not, to be used for library generation.
        
    high_res : boolean
        Whether the isotopic distribution should be of high resolution, to be used for library generation.
        
    ms2 : boolean
        Whether ms2 data should be analyzed.
        
    accuracy_unit : string
        The accuracy unit to be used in the script.
        
    accuracy_value : float, int
        The accuracy value to be used in the script.
    
    rt_int : tuple
        A tuple containing the beggining and end of retention times to analyze.
        
    min_isotop : int
        The minimum amount of isotopologue peaks to consider in the EIC processing.
        
    max_ppm : int
        The maximum amount of PPM difference to consider when outputting results.
        
    iso_fit : float
        The minimum isotopic fitting score to consider when outputting results.
        
    curve_fit : float
        The minimum curve fitting score to consider when outputting results.
        
    sn : int
        The minimum signal-to-noise ration to consider when outputting results.
        
    files : list
        A list of paths to sample files.
        
    path : string
        The working directory of the script.
        
    commented : boolean
        Whether the parameters template file should be commented or not.
        
    number : int
        Number of the execution files to analyze.
    '''
    date = datetime.datetime.now()
    begin_time = str(date)[2:4]+str(date)[5:7]+str(date)[8:10]+"_"+str(date)[11:13]+str(date)[14:16]+str(date)[17:19]
    input_order = [None]
    curr_os = platform.system()
    if curr_os == "Linux":
        default_path = "/home/GlycoGenius/"
    if curr_os == "Windows":
        default_path = "C:/GlycoGenius/"
    if curr_os == "Darwin":
        print("OS not tested for compatibility.")
        default_path = "/home/GlycoGenius/"
    while input_order[0] == None:
        print_header()
        print("1 - Build and output glycans library.\n2 - Analyze sample files in single-threaded mode\n3 - Reanalyze raw results files with new\n    parameters\n4 - Create template parameters file for command-\n    line execution\n5 - Exit")
        var = input("Select your option: ")
        if var == 'warranty':
            print("\nDisclaimer of Warranty.\n")
            print("THERE IS NO WARRANTY FOR THE PROGRAM, TO THE")
            print("EXTENT PERMITTED BY APPLICABLE LAW. EXCEPT WHEN")
            print("OTHERWISE STATED IN WRITING THE COPYRIGHT")
            print("HOLDERS AND/OR OTHER PARTIES PROVIDE THE PROGRAM")
            print("\"AS IS\" WITHOUT WARRANTY OF ANY KIND, EITHER")
            print("EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED")
            print("TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY")
            print("AND FITNESS FOR A PARTICULAR PURPOSE. THE ENTIRE")
            print("RISK AS TO THE QUALITY AND PERFORMANCE OF THE")
            print("PROGRAM IS WITH YOU. SHOULD THE PROGRAM PROVE")
            print("DEFECTIVE, YOU ASSUME THE COST OF ALL NECESSARY")
            print("SERVICING, REPAIR OR CORRECTION.\n")
            continue
        if var == 'version':
            print("\nGlycoGenius version: "+version)
            continue
        if var == 'license':
            license_path = str(pathlib.Path(__file__).parent.parent.parent.resolve())
            for i_i, i in enumerate(license_path):
                if i == "\\":
                    license_path = license_path[:i_i]+"/"+license_path[i_i+1:]
            with open(license_path+"/LICENSE", 'r') as f:
                for line in f:
                    print(line, end = "")
            continue
        if var == 'sneakpeek':
            path = default_path
            while True:
                var = input("Insert the working directory (where the\n'resultsn_' files are, default: "+default_path+"): ")
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
            number = 0
            while True:
                var = input("Type execution number: ")
                try:
                    number = int(var)
                except:
                    print("Wrong Input")
                    continue
                break
            return [69], path, number
        try:
            var = int(var)
        except:
            print("Wrong Input")
            continue
        if var == 5:
            os._exit(1)
        if var < 1 or var > 4:
            print("Wrong Input")
            continue
        if var > 0 and var <= 4:
            input_order[0] = var
    print_sep()
    if input_order[0] == 1 or input_order[0] == 2:
        input_order.append(None)
        while input_order[1] == None:
            print("1 - Use a custom glycans list\n2 - Generate a library based on possible\n    monosaccharides numbers")
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
            print("Warning: Glycans must be inputted in the format\nof 'H5N4S1F1G1' where 'H' = Hexose,'N' =\nN-Acetylhexosamine, 'S' = N-Acetyl-Neuraminic\nAcid, 'F' = Deoxyhexose, 'G' =\nN-Glycolyl-Neuraminic Acid and each number next\nto it corresponds to the amount of said\nmonosaccharide\n")
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
                        var = input("Type the minimum amount of monosaccharides\n (default: 5): ")
                        if var == '':
                            lib_settings[0] = 5
                        try:
                            var = int(var)
                        except:
                            continue
                        lib_settings[0] = var
                if i == 1:
                    while lib_settings[1] == None:
                        var = input("Type the maximum amount of monosaccharides\n (default: 18): ")
                        if var == '':
                            lib_settings[1] = 18
                        try:
                            var = int(var)
                        except:
                            continue
                        lib_settings[1] = var
                if i == 2:
                    while lib_settings[2] == None:
                        var = input("Type the minimum amount of Hex (default: 3): ")
                        if var == '':
                            lib_settings[2] = 3
                        try:
                            var = int(var)
                        except:
                            continue
                        lib_settings[2] = var
                if i == 3:
                    while lib_settings[3] == None:
                        var = input("Type the maximum amount of Hex (default: 10): ")
                        if var == '':
                            lib_settings[3] = 10
                        try:
                            var = int(var)
                        except:
                            continue
                        lib_settings[3] = var
                if i == 4:
                    while lib_settings[4] == None:
                        var = input("Type the minimum amount of HexNAc (default: 2): ")
                        if var == '':
                            lib_settings[4] = 2
                        try:
                            var = int(var)
                        except:
                            continue
                        lib_settings[4] = var
                if i == 5:
                    while lib_settings[5] == None:
                        var = input("Type the maximum amount of HexNAc (default: 8): ")
                        if var == '':
                            lib_settings[5] = 8
                        try:
                            var = int(var)
                        except:
                            continue
                        lib_settings[5] = var
                if i == 6:
                    while lib_settings[6] == None:
                        var = input("Type the minimum amount of dHex (default: 0): ")
                        if var == '':
                            lib_settings[6] = 0
                        try:
                            var = int(var)
                        except:
                            continue
                        lib_settings[6] = var
                if i == 7:
                    while lib_settings[7] == None:
                        var = input("Type the maximum amount of dHex (default: 2): ")
                        if var == '':
                            lib_settings[7] = 2
                        try:
                            var = int(var)
                        except:
                            continue
                        lib_settings[7] = var
                if i == 8:
                    while lib_settings[8] == None:
                        var = input("Type the minimum amount of Neu5Ac (default: 0): ")
                        if var == '':
                            lib_settings[8] = 0
                        try:
                            var = int(var)
                        except:
                            continue
                        lib_settings[8] = var
                if i == 9:
                    while lib_settings[9] == None:
                        var = input("Type the maximum amount of Neu5Ac (default: 4): ")
                        if var == '':
                            lib_settings[9] = 4
                        try:
                            var = int(var)
                        except:
                            continue
                        lib_settings[9] = var
                if i == 10:
                    while lib_settings[10] == None:
                        var = input("Type the minimum amount of Neu5Gc (default: 0): ")
                        if var == '':
                            lib_settings[10] = 0
                        try:
                            var = int(var)
                        except:
                            continue
                        lib_settings[10] = var
                if i == 11:
                    while lib_settings[11] == None:
                        var = input("Type the maximum amount of Neu5Gc (default: 0): ")
                        if var == '':
                            lib_settings[11] = 0
                        try:
                            var = int(var)
                        except:
                            continue
                        lib_settings[11] = var
                if i == 12:
                    while lib_settings[12] == None:
                        var = input("Type the minimum amount of total sialic acids\n(default: 0): ")
                        if var == '':
                            lib_settings[12] = 0
                        try:
                            var = int(var)
                        except:
                            continue
                        lib_settings[12] = var
                if i == 13:
                    while lib_settings[13] == None:
                        var = input("Type the maximum amount of total sialic acids\n(default: 4): ")
                        if var == '':
                            lib_settings[13] = 4
                        try:
                            var = int(var)
                        except:
                            continue
                        lib_settings[13] = var
                if i == 14:
                    while lib_settings[14] == None:
                        var = input("Force compositions to N-glycans structure\n (default: yes) (y/n): ")
                        if var == '':
                            lib_settings[14] = True
                        if var == 'y':
                            lib_settings[14] = True
                        if var == 'n':
                            lib_settings[14] = False
                        else:
                            continue
        print_sep()
        adducts = {}
        while True:
            var = input("Type an element to calculate as adduct\n(ie. Na or H). Leave blank to finish with default (H3): ")
            if var == '':
                if  len(adducts) == 0:
                    adducts = {'H' : 3}
                print(adducts)
                var2 = input("Proceed with these adducts? (default: yes) (y/n): ")
                if var2 == 'y' or var2 == '':
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
        max_charges = 3
        while True:
            var = input("Type the maximum amount of charges (default: 3): ")
            if var == '':
                break
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
            var = input("Do the glycans have a reducing end tag? (y/n): ")
            if var == 'y' or var == 'n':
                break
            else:
                print("Wrong Input")
                continue
        print("")
        if var == 'y':
            while True:
                var2 = input("Insert the tag added mass\nor molecular formula (ie. 133.0644 or C7H7N3\nfor GirP or 219.1735 or C13H21N3 for ProA): ")
                try:
                    var2 = float(var2)
                except:
                    pass
                tag_mass = var2
                break
        permethylated = False
        reduced = False
        if tag_mass == 0:
            while True:
                var = input("Are the glycans permethylated (y/n): ")
                if var == 'y':
                    permethylated = True
                    break
                if var == 'n':
                    break
                else:
                    print('Wrong input')
                    continue
            while True:
                var = input("Are the glycans reduced (y/n): ")
                if var == 'y':
                    reduced = True
                    break
                if var == 'n':
                    break
                else:
                    print('Wrong input')
                    continue
        print_sep()
        fast_iso = True
#        while True:
#            var = input("Do you want to do a quick isotopic distribution\ncalculation? If 'n', then isotopic distribution\ncalculation may take several hours, depending on\nlibrary size (y/n): ")
#            if var == 'y':
#                break
#            if var == 'n':
#                fast_iso = False
#                break
#            else:
#                print('Wrong input')
#                continue
        high_res = False
#        if not fast_iso:
#            print("")
#            while True:
#                var = input("Do you need a high resolution isotopic\ndistribution? It may be important for very high\naccuracy mass spectrometers, such as\nFT-ICR (y/n): ")
#                if var == 'y':
#                    high_res = True
#                    break
#                if var == 'n':
#                    break
#                else:
#                    print('Wrong input')
#                    continue
#            if input_order[0] != 1:
#                print("")
        if input_order[0] == 1: #Outputs of input_order == 1
            path = default_path
            while True:
                var = input("Insert the path to save the files produced by\nthe script (leave blank for default:\n"+default_path+"): ")
                if var == '':
                    var = path
                print(var)
                var2 = input("Is this path correct? (y/n): ")
                if var2 == 'n':
                    continue
                if var2 == 'y' or var2 == '':
                    for i_i, i in enumerate(var):
                        if i == "\\":
                            var = var[:i_i]+"/"+var[i_i+1:]
                    if var[-1] != "/":
                        var = var+"/"
                    path = var
                    break
            if input_order[1] == 1:
                return input_order, glycans_list, adducts, max_charges, tag_mass, fast_iso, high_res, path, permethylated, reduced
            if input_order[1] == 2:
                return input_order, lib_settings, adducts, max_charges, tag_mass, fast_iso, high_res, path, permethylated, reduced
        else:
            print_sep()
            ms2 = [False, False, False]
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
            print("")
            if ms2[0]:
                while True:
                    var = input("Do you want to only output fragments compatible\nwith identified precursor glycan? (y/n): ")
                    if var == 'y':
                        ms2[1] = True
                        break
                    if var == 'n':
                        break
                    else:
                        print('Wrong input')
                        continue
                print("")
            accuracy_unit = "mz"
            while True:
                var = input("What is the accuracy unit you want to input for\nmz tolerance? (default: mz) (ppm/mz): ")
                if var == 'ppm':
                    accuracy_unit = var
                    break
                if var == 'mz' or var == '':
                    break
                else:
                    print('Wrong input')
                    continue
            print("")
            accuracy_value = 0.0
            while True:
                var = input("Insert the accuracy value for the unit you've\nchosen (ie. '0.01' or '10'): ")
                try:
                    var = float(var)
                except:
                    print('Wrong input')
                    continue
                accuracy_value = var
                break
            print("")
            rt_int = [0.0, 999]
            while True:
                var = input("Insert the start of the retention time interval\nat which you want to analyze, in minutes\n (default: 0 mins): ")
                if var == '':
                    break
                try:
                    var = float(var)
                except:
                    print('Wrong input')
                    continue
                rt_int[0] = var
                break
            print("")
            while True:
                var = input("Insert the end of the retention time interval at\nwhich you want to analyze, in minutes\n (default: 999 mins): ")
                if var == '':
                    break
                try:
                    var = float(var)
                except:
                    print('Wrong input')
                    continue
                rt_int[1] = var
                break
            print("")
            min_isotop = 2
            max_ppm = 10
            while True:
                var = input("Insert the maximum PPM error that a detected\nglycan must have in order to show up in\nresults' table (default: 10 ppm): ")
                if var == '':
                    break
                try:
                    var = int(var)
                except:
                    print('Wrong input')
                    continue
                max_ppm = var
                break
            print("")
            iso_fit = 0.9
            while True:
                var = input("Insert the minimum isotopic fitting score for a\nglycan in order for it to show up in the\nresults' table (values between 0.0 and 1.0) (default: 0.9): ")
                if var == '':
                    break
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
            print("")
            curve_fit = 0.9
            while True:
                var = input("Insert the minimum curve fitting score for a\nglycan in order for it to show up in the\nresults' table (values between 0.0 and 1.0) (default: 0.9): ")
                if var == '':
                    break
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
            print("")
            sn = 3
            while True:
                var = input("Insert the minimum signal-to-noise ratio that a\ndetected glycan must have in order to show up in\nresults' table (default: 3): ")
                if var == '':
                    break
                try:
                    var = int(var)
                except:
                    print('Wrong input')
                    continue
                sn = var
                break
            print("")
            files = default_path+"Sample Files/"
            while True:
                var = input("Insert the path to the folder containing the\nsample files to be analyzed ( leave blank for\ndefault: "+default_path+"Sample Files/"+"): ")
                if var == '':
                    var = files
                print(var)
                var2 = input("Is this path correct? (y/n): ")
                if var2 == 'n':
                    continue
                if var2 == 'y' or var2 == '':
                    for i_i, i in enumerate(var):
                        if i == "\\":
                            var = var[:i_i]+"/"+var[i_i+1:]
                    if var[-1] != "/":
                        var = var+"/"
                    path = var
                    break
            print("")
            path = default_path
            while True:
                var = input("Insert the path to save the files produced by\nthe script (leave blank for default:\n"+default_path+"): ")
                if var == '':
                    var = path
                print(var)
                var2 = input("Is this path correct? (y/n): ")
                if var2 == 'n':
                    continue
                if var2 == 'y' or var2 == '':
                    for i_i, i in enumerate(var):
                        if i == "\\":
                            var = var[:i_i]+"/"+var[i_i+1:]
                    if var[-1] != "/":
                        var = var+"/"
                    path = var
                    break
            if input_order[1] == 1:
                return input_order, glycans_list, adducts, max_charges, tag_mass, fast_iso, high_res, ms2, accuracy_unit, accuracy_value, rt_int, min_isotop, max_ppm, iso_fit, curve_fit, sn, files, path, permethylated, reduced
            if input_order[1] == 2:
                return input_order, lib_settings, adducts, max_charges, tag_mass, fast_iso, high_res, ms2, accuracy_unit, accuracy_value, rt_int, min_isotop, max_ppm, iso_fit, curve_fit, sn, files, path, permethylated, reduced
    if input_order[0] == 3:
        path = default_path
        while True:
            var = input("Insert the working directory (where the\n'raw_data' files are, default: "+default_path+"): ")
            if var == "":
                var = path
            print(var)
            var2 = input("Is this path correct? (y/n): ")
            if var2 == 'n':
                continue
            if var2 == 'y' or var2 == '':
                for i_i, i in enumerate(var):
                    if i == "\\":
                        var = var[:i_i]+"/"+var[i_i+1:]
                if var[-1] != "/":
                    var = var+"/"
                path = var
                break
        print("")
        max_ppm = 10
        while True:
            var = input("Insert the maximum amount of PPM difference that\na detected glycan must have in order to show up\nin results' table: ")
            try:
                var = int(var)
            except:
                print('Wrong input')
                continue
            max_ppm = var
            break
        print("")
        iso_fit = 0.9
        while True:
            var = input("Insert the minimum isotopic fitting score for a\nglycan in order for it to show up in the\nresults' table (values between 0.0 and 1.0): ")
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
        print("")
        curve_fit = 0.9
        while True:
            var = input("Insert the minimum curve fitting score for a\nglycan in order for it to show up in the\nresults' table (values between 0.0 and 1.0): ")
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
        print("")
        sn = 3
        while True:
            var = input("Insert the minimum signal-to-noise ratio that\na detected glycan must have in order to show up\nin results' table: ")
            try:
                var = int(var)
            except:
                print('Wrong input')
                continue
            sn = var
            break
        print("")
        return input_order, path, max_ppm, iso_fit, curve_fit, sn
    if input_order[0] == 4:
        commented = False
        while True:
            var = input("Do you want the template file to contain\ncommented information on each parameter? (y/n): ")
            if var == 'y':
                commented = True
                break
            if var == 'n':
                break
            else:
                print('Wrong input')
                continue
        path = default_path
        while True:
            var = input("Insert the path to the folder to save the\ntemplate file (Default: "+default_path+"): ")
            if var == "":
                var = path
            print(var)
            var2 = input("Is this path correct? (y/n): ")
            if var2 == 'n':
                continue
            if var2 == 'y' or var2 == '':
                for i_i, i in enumerate(var):
                    if i == "\\":
                        var = var[:i_i]+"/"+var[i_i+1:]
                if var[-1] != "/":
                    var = var+"/"
                path = var
                break
        return input_order, commented, path
                
def samples_path_to_list(path):
    '''Detects all files in samples path and makes a list of the path to each file for
    internal usage.
    
    Parameters
    ----------
    path : string
        Path to the folder containing mzML and mzXML files to analyze.
        
    Uses
    ----
    itertools.product : 
        Cartesian product of input iterables.
        
    Returns
    -------
    samples_list : list
        A list containing the path to each file to be analyzed.
    '''
    if len(path) == 0:
        return []
    path = path.strip()
    path = path.strip("'")
    path = path.strip("\"")
    mzml_possibilities = list(map(''.join, product(*zip("mzml".upper(), "mzml".lower()))))
    mzxml_possibilities = list(map(''.join, product(*zip("mzxml".upper(), "mzxml".lower()))))
    file_extensions = mzml_possibilities+mzxml_possibilities
    for i_i, i in enumerate(path):
        if i == "\\":
            path = path[:i_i]+"/"+path[i_i+1:]
    if path[-1] != "/":
        path+= "/"
    dir_list = os.listdir(path)
    samples_list = []
    for i_i, i in enumerate(dir_list):
        if i.split('.')[-1] in file_extensions:
            samples_list.append(path+i)
    return samples_list
    
def list_of_data(samples_list): ##complete
    '''Detects if a file is mzXML or mzML and processes it into a generator using
    pyteomics.

    Parameters
    ----------
    samples_list : list
        A list containing the path to each file to be analyzed.

    Uses
    ----
    pyteomics.mzxml.MzXML() : generator
        Indexes the mzXML file into a generator, allowing you to parse the file for
        analysis.

    File_Accessing.make_mzxml : class
        A wrapper that takes the output of pyteomics mzML parser and converts it to
        the mzXML pyteomics parser standard to be used within the script. Allows for 
        full support of mzML.

    Returns
    -------
    data : list
        A list containing the generator of each file name at each index.
    '''
    data = []
    mzml_possibilities = list(map(''.join, product(*zip("mzml".upper(), "mzml".lower()))))
    mzxml_possibilities = list(map(''.join, product(*zip("mzxml".upper(), "mzxml".lower()))))
    for i in samples_list:
        if i.split('.')[-1] in mzxml_possibilities:
            mzxml_data = mzxml.MzXML(i)
            data.append(mzxml_data)
        else:
            mzml_data = File_Accessing.make_mzxml(i)
            data.append(mzml_data)
    return data

def index_ms1_from_file(files):
    '''Scans the mz(X)ML file and indexes all the MS1 scans, so that you don't have to
    go through the MSn scans as well when doing something only with the MS1.

    Parameters
    ----------
    files : list
        List with each index containing a generator created by the pyteomics function
        pyteomics.mzxml.MzXML() or File_Accessing.make_mzxml.

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
    files : list
        List with each index containing a generator created by the pyteomics function
        pyteomics.mzxml.MzXML() or File_Accessing.make_mzxml.

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
    samples_list : list
        A list of strings containing the path to each sample file.

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
                        save_path,
                        internal_standard,
                        permethylated,
                        reduced):
    '''Imports, generates and/or exports a glycans library.

    Parameters
    ----------
    multithreaded_analysis : tuple
        A tuple with two indexes: the first one contains whether or not this execution
        is to produce the tools required for multithreaded analysis and the second one 
        containing the amount of desired threads to split the execution in.
        
    multithread_execution : tuple
        A tuple with three indexes: the first one contains whether or not this execution
        is part of a multithreaded run, the second one contains the current thread number and
        the third one is the total amount of threads.
        
    samples_names : list
        A list containing the extracted sample names from the file names.
        
    custom_glycans_list : tuple
        A tuple with two indexes: The first one is whether or not the library should be built
        off a custom glycans list and the second one is a list containing the glycans you wish
        to analyze.
        
    min_max_monos : tuple
        Minimum and maximum amount of monosaccharides for the hypotethical glycans in the
        library. ie. (5, 20).
        
    min_max_hex : tuple
        Minimum and maximum amount of hexoses for the hypotethical glycans in the library.
        ie. (5, 20).
        
    min_max_hexnac : tuple
        Minimum and maximum amount of N-Acetyl hexosamines for the hypotethical glycans
        in the library. ie. (5, 20).
        
    min_max_sia : tuple
        Minimum and maximum amount of sialic acids for the hypotethical glycans in the
        library. ie. (5, 20).
        
    min_max_fuc : tuple
        Minimum and maximum amount of deoxyhexoses for the hypotethical glycans in the
        library. ie. (5, 20).
        
    min_max_ac : tuple
        Minimum and maximum amount of N-Acetyl Neuraminic acids for the hypotethical
        glycans in the library. ie. (5, 20).
        
    min_max_gc : tuple
        Minimum and maximum amount of N-Glicolyl Neuraminic acids for the hypotethical
        glycans in the library. ie. (5, 20).
        
    force_nglycan : boolean
        Indicates whether the function should force strict conditions based on the
        biological knowledge of glycans in order to avoid possible false positives when
        analysing N-glycans.
        
    max_adducts : dict
        A dictionary with keys containing each possible atomic adducts (ie. 'H', 'Na',
        'K', etc.) and the maximum amount of such adducts as the values.
        
    max_charges : int
        The maximum amount of charges to calculate per glycan.

    tag_mass : float
        The tag's added mass to the glycans, if the glycans are tagged.
        Default = 0 (No Tag).

    fast_iso : boolean
        Makes the isotopic distribution calculation fast (less accurate) or not (more
        accurate).
        Default = True
        
    high_res : boolean
        Decides whether to clump (if set to False) or not (if set to True) the neighbouring
        isotopic envelope peaks. Only works if fast is set to False.
        Default = False

    imp_exp_library : tuple
        A list with two indexes : The first one indicates whether or not the script should try
        to import the library from a glycans_library.py file and the second one indicates
        whether or not the script should export the library to a glycans_libbrary.py file and
        an excel file for visualization.
        
    only_gen_lib : boolean
        A boolean indicating whether or not this execution should stop after generating library.
        
    save_path : string
        A string containing the path to the working directory of the script.
        
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
        
    pathlib.Path.resolve() : Path object
        Make the path absolute, resolving any symlinks. A new path object is returned

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
        print('Building custom glycans library...', end = "", flush = True)
        for i in custom_glycans_list[1]:
            custom_glycans_comp.append(General_Functions.sum_monos(General_Functions.form_to_comp(i)))
        full_library = Library_Tools.full_glycans_library(custom_glycans_comp,
                                            max_adducts,
                                            max_charges,
                                            tag_mass,
                                            fast_iso,
                                            high_res,
                                            internal_standard,
                                            permethylated,
                                            reduced)
        if imp_exp_library[1]:
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
        print('Done!')
    if multithreaded_analysis[0] and not only_gen_lib and not multithreaded_execution[0]:
        print('Splitting execution for multiple threads...')
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
                                                high_res,
                                                internal_standard,
                                                permethylated,
                                                reduced)
            if imp_exp_library[1]:
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
        elif imp_exp_library[0]:
            lib_module = importlib.import_module('glycans_library')
            full_library = lib_module.full_library
            if imp_exp_library[1]:
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
            temp_library = ""
            for j in range(start, start+split_s):
                if j == start+split_s-1 or j == len(full_library)-1:
                    temp_library+= "'"+full_library_keys_list[j]+"'"+": "+str(full_library[full_library_keys_list[j]])
                    start+=split_s
                    break
                if j < start+split_s-1:
                    temp_library+= "'"+full_library_keys_list[j]+"'"+": "+str(full_library[full_library_keys_list[j]])+", "
            libraries.append(temp_library)
        mt_path = str(pathlib.Path(__file__).parent.resolve())
        for i_i, i in enumerate(mt_path):
            if i == "\\":
                mt_path = mt_path[:i_i]+"/"+mt_path[i_i+1:]
        with open(mt_path+'/core.py', 'r') as f:
            for i_i, i in enumerate(f):
                if i_i == 0:
                    mode = 'w'
                else:
                    mode = 'a'
                for j in range(multithreaded_analysis[1]):
                    with open(save_path+'glycogenius_'+str(j)+'.py', mode) as g:
                        if i_i == 18:
                            g.write("import importlib\n")
                            g.write("spec1 = importlib.util.spec_from_file_location('Execution_Functions', '"+mt_path+"/Execution_Functions.py')\n")
                            g.write("Execution_Functions = importlib.util.module_from_spec(spec1)\n")
                            g.write("spec1.loader.exec_module(Execution_Functions)\n")
                            continue
                        if i_i == 19:
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
                            g.write("        full_library = {"+str(libraries[j])+"}\n")
                            continue
                        if i == "#here multithreaded prints execution of main()":
                            g.write("main()")
                            continue
                        g.write(i)
                        g.close()
            f.close()
        print("Setup done. Run 'glycogenius_n.py' files now.")
        os._exit(1)
    if imp_exp_library[0]:
        print('Importing existing library...', end = '', flush = True)
        spec = importlib.util.spec_from_file_location("glycans_library", save_path+'glycans_library.py')
        lib_module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(lib_module)
        full_library = lib_module.full_library
        if imp_exp_library[1]:
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
        print("Done!")
        return full_library
    if not custom_glycans_list[0]:
        print('Building glycans library...', end = "", flush = True)
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
                                            high_res,
                                            internal_standard,
                                            permethylated,
                                            reduced)
        print('Done!')
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
        with open(save_path+'skyline_transitions.csv', 'w') as f:
            f.write('Precursor Name, Precursor Formula, Precursor Adduct, Precursor Charge\n')
            for i_i, i in enumerate(full_library):
                for j_j, j in enumerate(full_library[i]['Adducts_mz']):
                    adduct_comp = General_Functions.form_to_comp(j)
                    if len(adduct_comp) > 1: #can't seem to make skyline work with mixed adducts, so have this in place for now
                        continue
                    adduct = str(adduct_comp[list(adduct_comp.keys())[0]])+str(list(adduct_comp.keys())[0]) #only first adduct
                    del adduct_comp[list(adduct_comp.keys())[0]]
                    formula = General_Functions.comp_to_formula(General_Functions.sum_atoms(full_library[i]['Atoms_Glycan+Tag'], adduct_comp))
                    list_form = [i, str(formula), '[M+'+adduct+']', str(General_Functions.form_to_charge(j))]
                    f.write(",".join(list_form)+'\n')
            f.close()
    if only_gen_lib:
        print('Library length: '+str(len(full_library)))
        print("Check it in Glycans_Library.xlsx.")
        print("If you wish to analyze files,")
        print("set 'only_gen_lib' to False and input")
        print("remaining parameters.")
        try:
            input("\nPress Enter to exit.")
            os._exit(1)
        except:
            os._exit(1)
    return full_library
    
def align_assignments(df, df_type, deltas = None):
    '''
    '''
    dataframe = copy.deepcopy(df)
    if df_type == 'total_glycans': #this is for alignment of total_glycans dataframes... it generates delta values to be used for other alignments
        biggest_df = inf
        biggest_df_size = 0
        for i_i, i in enumerate(dataframe): #this determines the dataframe with most peaks assigned, it will be used as reference
            if len(i['Glycan']) > biggest_df_size:
                biggest_df_size = len(i['Glycan'])
                biggest_df = i_i
        ids_per_sample = []
        deltas_per_sample = [] #this is a list that will contain every delta of aligned glycans per sample, in a dictionary of 'rt : delta'
        skips = 0
        for i_i, i in enumerate(dataframe[biggest_df]['Glycan']): #goes through each glycan in the reference df to align the other samples with
            if skips != 0:
                skips -= 1
                continue
            peaks_number = dataframe[biggest_df]['Glycan'].count(i)
            skips = peaks_number-1
            rts_reference = dataframe[biggest_df]['RT'][i_i:i_i+peaks_number] #this collects the reference rts for the glycan in reference sample
            aucs_reference = []
            for j in range(i_i, i_i+peaks_number):
                aucs_reference.append(dataframe[biggest_df]['AUC'][j]/max(dataframe[biggest_df]['AUC'][i_i:i_i+peaks_number]))
            for j_j, j in enumerate(dataframe): #going through each target sample
                if i_i == 0:
                    ids_per_sample.append([])
                    deltas_per_sample.append({})
                if j_j == biggest_df:
                    continue
                peaks_number_target = j['Glycan'].count(i)
                if peaks_number_target == 0: #glycan not available in target sample
                    continue
                first_glycan_id = 0
                rts = []
                aucs = []
                for k_k, k in enumerate(j['Glycan']): #going through each glycan of target sample
                    if k == i:
                        for l in range(k_k, k_k+peaks_number_target):
                            ids_per_sample[j_j].append(l)
                            rts.append(j['RT'][l])
                            aucs.append(j['AUC'][l]/max(j['AUC'][k_k:k_k+peaks_number_target]))
                        used_peaks = []
                        last_assigned_rt_id = 0
                        for l_l, l in enumerate(aucs):
                            found = False
                            if len(aucs_reference[last_assigned_rt_id:]) != 0: 
                                for m_m in range(last_assigned_rt_id, len(aucs_reference)):
                                    m = aucs_reference[m_m]
                                    if m_m not in used_peaks and abs(l-m) <= 0.5*m:
                                        if m_m > aucs_reference.index(1.0) and l_l < aucs.index(1.0):
                                            break
                                        found = True
                                        last_assigned_rt_id = m_m+1
                                        deltas_per_sample[j_j][rts[l_l]] = [rts_reference[m_m]-rts[l_l], k]
                                        break
                            if not found:
                                deltas_per_sample[j_j][rts[l_l]] = [inf, k]
                        break
        for i_i, i in enumerate(dataframe): #this will apply the alignment based on ids and deltas
            if len(ids_per_sample[i_i]) == 0: #this is the reference sample or blank sample
                continue
            to_fix = []
            for j_j, j in enumerate(i['RT']):
                found = False
                for k_k, k in enumerate(deltas_per_sample[i_i]):
                    if i['RT'][j_j] == k and i['Glycan'][j_j] == deltas_per_sample[i_i][k][1] and deltas_per_sample[i_i][k][0] != inf:
                        i['RT'][j_j] = i['RT'][j_j] + deltas_per_sample[i_i][i['RT'][j_j]][0]
                        found = True
                        break
                if not found:
                    to_fix.append(j_j)
            if len(to_fix) > 0:
                for j_j, j in enumerate(i['RT']):
                    if j_j in to_fix:
                        before_j_j = None
                        after_j_j = None
                        if j_j-1 >= 0 and df[i_i]['Glycan'][j_j-1] == df[i_i]['Glycan'][j_j] and j_j-1 not in to_fix:
                            before_j_j = df[i_i]['RT'][j_j-1]
                        if j_j+1 < len(df[i_i]['Glycan']) and df[i_i]['Glycan'][j_j+1] == df[i_i]['Glycan'][j_j] and j_j+1 not in to_fix:
                            after_j_j = df[i_i]['RT'][j_j+1]
                        x = []
                        y = []
                        for k_k, k in enumerate(deltas_per_sample[i_i]):
                            if deltas_per_sample[i_i][k][0] != inf:
                                x.append(float(k))
                                y.append(deltas_per_sample[i_i][k][0])
                        linear_equation = General_Functions.linear_regression(x, y)
                        fixed_RT = float("%.2f" % round(i['RT'][j_j] + (i['RT'][j_j]*linear_equation[0]) + linear_equation[1], 2))
                        if before_j_j != None and after_j_j != None:
                            scaling_factor = (i['RT'][j_j]-before_j_j)/(after_j_j-before_j_j)
                            fixed_RT = float("%.2f" % round(i['RT'][j_j-1]+(scaling_factor*(i['RT'][j_j+1]-i['RT'][j_j-1])), 2))
                        elif before_j_j == None and after_j_j != None:
                            scaling_factor = i['RT'][j_j]/after_j_j
                            fixed_RT = float("%.2f" % round(i['RT'][j_j+1]*scaling_factor, 2))
                        elif before_j_j != None and after_j_j == None:
                            scaling_factor = i['RT'][j_j]/before_j_j
                            fixed_RT = float("%.2f" % round(i['RT'][j_j-1]*scaling_factor, 2))
                        if j in list(deltas_per_sample[i_i].keys()):
                            deltas_per_sample[i_i][j][0] = i['RT'][j_j] - fixed_RT
                        else:
                            deltas_per_sample[i_i][j] = [i['RT'][j_j] - fixed_RT, i['Glycan'][j_j]]
                        i['RT'][j_j] = fixed_RT
        for i_i, i in enumerate(deltas_per_sample):
            deltas_per_sample[i_i] = dict(sorted(i.items()))
        for i_i, i in enumerate(deltas_per_sample): #cleanup
            positive = []
            negative = []
            for j_j, j in enumerate(i):
                if i[j][0] > 0:
                    positive.append(j)
                if i[j][0] < 0:
                    negative.append(j)
            if len(positive) < len(negative):
                for j in positive:
                    del deltas_per_sample[i_i][j]
            if len(positive) > len(negative):
                for j in negative:
                    del deltas_per_sample[i_i][j]
        return dataframe, deltas_per_sample, biggest_df
        
    if df_type == "chromatograms": #for when you want to align chromatograms based on existing deltas calculated previously
        for i_i, i in enumerate(dataframe): #sample by sample
            if len(deltas[i_i]) > 0:
                chromatogram_length_rt = i['RTs_'+str(i_i)][-1]
                chromatogram_beg_rt = i['RTs_'+str(i_i)][0]
                chromatogram_length = len(i['RTs_'+str(i_i)])
                chromatogram_interval = i['RTs_'+str(i_i)][-1]/len(i['RTs_'+str(i_i)])
                points_per_minute = int(1/chromatogram_interval)
                interval_list = []
                current = 0.0
                minutes_interval = 5
                for j in deltas[i_i]: #rules here still to be discussed
                    if j > current+minutes_interval:
                        current = j
                        closest = 0
                        distance = inf
                        for k_k, k in enumerate(i['RTs_'+str(i_i)]):
                            if abs(k-j) < distance:
                                closest = k_k
                                distance = abs(k-j)
                        interval_list.append(closest-int(points_per_minute*(minutes_interval/5)))
                interval_list.append(interval_list[-1] + points_per_minute*minutes_interval)
                for j_j, j in enumerate(interval_list):
                    x = []
                    y = []
                    if j_j == len(interval_list)-1:
                        last_range = len(i['RTs_'+str(i_i)])-1
                    else:
                        last_range = interval_list[j_j+1]
                    for k_k, k in enumerate(deltas[i_i]):
                        if k > i['RTs_'+str(i_i)][j] and k < i['RTs_'+str(i_i)][last_range]:
                            x.append(float(k))
                            y.append(deltas[i_i][k][0])
                    if len(x) == 0:
                        continue
                    linear_equation = General_Functions.linear_regression(x, y)
                    lower_boundary = i['RTs_'+str(i_i)][j-1]
                    upper_boundary = i['RTs_'+str(i_i)][last_range]
                    lowest = inf
                    lowest_id = 0
                    highest = 0
                    highest_id = inf
                    for k in range(j, last_range):
                        i['RTs_'+str(i_i)][k] = i['RTs_'+str(i_i)][k] + (i['RTs_'+str(i_i)][k]*linear_equation[0])+linear_equation[1]
                        if i['RTs_'+str(i_i)][k] < lowest:
                            lowest = i['RTs_'+str(i_i)][k]
                            lowest_id = k
                        if i['RTs_'+str(i_i)][k] > highest:
                            highest = i['RTs_'+str(i_i)][k]
                            highest_id = k
                    if lowest != inf and lowest < lower_boundary:
                        if j_j > 0:
                            interval = (lowest-i['RTs_'+str(i_i)][interval_list[j_j-1]])/(lowest_id+1-interval_list[j_j-1])
                            for l_l, l in enumerate(range(lowest_id-1, interval_list[j_j-1], -1)):
                                i['RTs_'+str(i_i)][l] = float("%.4f" % round(lowest-(interval*(l_l+1)), 4))
                        else:
                            interval = (lowest)/(lowest_id+1)   
                            for l_l, l in enumerate(range(lowest_id-1, -1, -1)):
                                i['RTs_'+str(i_i)][l] = float("%.4f" % round(lowest-(interval*(l_l+1)), 4))
                    if highest != 0 and highest > upper_boundary:
                        if j_j != len(interval_list)-2:
                            interval = (i['RTs_'+str(i_i)][interval_list[j_j+2]]-highest)/(interval_list[j_j+2]-highest_id)
                            for l_l, l in enumerate(range(highest_id+1, interval_list[j_j+1])):
                                i['RTs_'+str(i_i)][l] = float("%.4f" % round(highest+(interval*(l_l+1)), 4))
                        else:
                            interval = (chromatogram_length_rt-highest)/(chromatogram_length-1-highest_id)
                            for l_l, l in enumerate(range(highest_id+1, chromatogram_length)):
                                i['RTs_'+str(i_i)][l] = float("%.4f" % round(highest+(interval*(l_l+1)), 4))
        return dataframe
        
def output_filtered_data(curve_fit_score,
                         iso_fit_score,
                         sn,
                         max_ppm,
                         reanalysis,
                         save_path,
                         multithreaded_analysis,
                         multithreaded_execution,
                         analyze_ms2,
                         unrestricted_fragments,
                         reporter_ions,
                         plot_metaboanalyst,
                         compositions,
                         align_chromatograms,
                         nglycan,
                         rt_tolerance,
                         rt_tolerance_frag,
                         output_isotopic_fittings,
                         sneakpeek):
    '''This function filters and converts raw results data into human readable
    excel files.
    
    Parameters
    ----------
    curve_fit_score : float
        The minimum curve fitting score to consider when outputting results.
        
    iso_fit_score : float
        The minimum isotopic fitting score to consider when outputting results.
        
    sn : int
        The minimum signal-to-noise ration to consider when outputting results.
        
    max_ppm : int
        The maximum amount of PPM difference to consider when outputting results.
        
    reanalysis : tuple
        Contains two indexes: The first one indicates whether this is only a reanalysis
        execution and the second one indicates whether or not to produce plot data excel
        files. The second one is there because the plot data is more heavy and doesn't 
        change on reanalysis, so should only be set to True if you lost original plot data.
        
    save_path : string
        A string containing the path to the working directory of the script.
        
    multithreaded_analysis : tuple
        A tuple with two indexes: the first one contains whether or not this execution
        is to produce the tools required for multithreaded analysis and the second one 
        containing the amount of desired threads to split the execution in.
        
    multithread_execution : tuple
        A tuple with three indexes: the first one contains whether or not this execution
        is part of a multithreaded run, the second one contains the current thread number and
        the third one is the total amount of threads.
    
    analyze_ms2 : tuple
        A tuple with two indexes: The first one indicates whether to analyze ms2 data and the
        second one indicates whether ms2 data should be forced to fit glycans composition.
        
    reporter_ions : list
        A list containing the reporter ions you wish to filter your data with.
        
    plot_metaboanalyst : tuple
        A tuple with two indexes: The first one indicates whether or not to output a file to be
        used in metaboanalyst and the second one indicates the groups for the samples to be separated
        in.
        
    compositions : boolean
        If set to True, also outputs the compositions analysis.
        
    nglycan : boolean
        Determines whether you're analyzing N-Glycans or not.
    
    rt_tolerance : float
        Tolerance of retention time (in minutes) at which an MS2 feature can be attributed to a 
        specific retention time peak and also for peaks in different samples to be regarded as the
        same peak (and thus be compared with each other).
        
    sneakpeek : boolean
        Allows to peek into incomplete multithreaded execution
        
    Uses
    ----
    datetime.datetime.now : Time object
        Returns the current date and time.
        
    dill.dump : None
        Pickle the current state of __main__ or another module to a file.
    
    dill.load : Module object
        Update the selected module (default is __main__) with the state saved at filename.
        
    pathlib.Path : Path object
        A subclass of PurePath, this class represents concrete paths of the system’s path flavour
        
    os.unlink : None
        Remove (delete) the file path. 
        
    pandas.DataFrame : Dataframe object
        Two-dimensional, size-mutable, potentially heterogeneous tabular data.
        
    pandas.ExcelWriter : ExcelWriter object
        Class for writing DataFrame objects into excel sheets.
    '''
    date = datetime.datetime.now()
    begin_time = str(date)[2:4]+str(date)[5:7]+str(date)[8:10]+"_"+str(date)[11:13]+str(date)[14:16]+str(date)[17:19]
    if reanalysis[0] and not sneakpeek[0]:
        print("Reanalyzing raw data with new parameters...")
        results1_list = []
        results2_list = []
        results3_list = []
        results4_list = []
        results5_list = []
        results6_list = []
        for i in range(multithreaded_analysis[1]):
            results1_list.append(save_path+'results1_'+str(i))
            results2_list.append(save_path+'results2_'+str(i))
            results3_list.append(save_path+'results3_'+str(i))
            results4_list.append(save_path+'results4_'+str(i))
            results5_list.append(save_path+'results5_'+str(i))
            results6_list.append(save_path+'results6_'+str(i))
        try:
            for i_i, i in enumerate(results1_list):
                with open(i, 'rb') as f:
                    file = dill.load(f)
                    df1_import = file[0]
                    df2_import = file[1]
                    if analyze_ms2:
                        fragments_dataframes_import = file[2]
                        total_threads = file[3]
                    else:
                        total_threads = file[2]
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
                with open(results5_list[i_i], 'rb') as f:
                    raw_eic_dataframes_import = dill.load(f)
                    f.close()
                with open(results6_list[i_i], 'rb') as f:
                    isotopic_fits_dataframes_import = dill.load(f)
                    f.close()
                if i_i == 0:
                    df1 = df1_import
                    df2 = df2_import
                    raw_eic_dataframes = raw_eic_dataframes_import
                    eic_dataframes = eic_dataframes_import
                    smoothed_eic_dataframes = smoothed_eic_dataframes_import
                    curve_fitting_dataframes = curve_fitting_dataframes_import
                    isotopic_fits_dataframes = isotopic_fits_dataframes_import
                    if analyze_ms2:
                        fragments_dataframes = fragments_dataframes_import
                else:
                    for j_j, j in enumerate(df1_import):
                        for k in j:
                            df1[j_j][k] = df1[j_j][k] + df1_import[j_j][k]
                    for j_j, j in enumerate(eic_dataframes_import):
                        for k in j:
                            if k[:4] != 'RTs_':
                                raw_eic_dataframes[j_j][k] = raw_eic_dataframes_import[j_j][k]
                                eic_dataframes[j_j][k] = eic_dataframes_import[j_j][k]
                                smoothed_eic_dataframes[j_j][k] = smoothed_eic_dataframes_import[j_j][k]
                    for j_j, j in enumerate(curve_fitting_dataframes_import):
                        for k in j:
                            curve_fitting_dataframes[j_j][k] = curve_fitting_dataframes_import[j_j][k]
                    for j_j, j in enumerate(isotopic_fits_dataframes_import):
                        for k in j:
                            isotopic_fits_dataframes[j_j][k] = j[k]
                    if analyze_ms2:
                        for j_j, j in enumerate(fragments_dataframes_import):
                            for k in j:
                                fragments_dataframes[j_j][k] = fragments_dataframes[j_j][k] + fragments_dataframes_import[j_j][k]
                if i_i == total_threads-1:
                    break
            with open(save_path+'raw_data_1', 'wb') as f:
                if analyze_ms2:
                    dill.dump([df1, df2, fragments_dataframes, version], f)
                else:
                    dill.dump([df1, df2, version], f)
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
            with open(save_path+'raw_data_5', 'wb') as f:
                dill.dump(raw_eic_dataframes, f)
                f.close()
            with open(save_path+'raw_data_6', 'wb') as f:
                dill.dump(isotopic_fits_dataframes, f)
                f.close()
            for i in range(multithreaded_analysis[1]):
                p1 = pathlib.Path(save_path+'results1_'+str(i))
                p2 = pathlib.Path(save_path+'results2_'+str(i))
                p3 = pathlib.Path(save_path+'results3_'+str(i))
                p4 = pathlib.Path(save_path+'results4_'+str(i))
                p5 = pathlib.Path(save_path+'results5_'+str(i))
                p6 = pathlib.Path(save_path+'results6_'+str(i))
                p1.unlink(missing_ok=True)
                p2.unlink(missing_ok=True)
                p3.unlink(missing_ok=True)
                p4.unlink(missing_ok=True)
                p5.unlink(missing_ok=True)
                p6.unlink(missing_ok=True)
        except:
            pass
    if not sneakpeek[0]:
        files_list_in_path = os.listdir(save_path)
        for i in files_list_in_path:
            if "results1" in i:
                return
    try:
        if sneakpeek[0]:
            with open(save_path+'results1_'+str(sneakpeek[1]), 'rb') as f:
                file = dill.load(f)
                df1 = file[0]
                df2 = file[1]
                if analyze_ms2:
                    fragments_dataframes = file[2]
                f.close()
        else:
            with open(save_path+'raw_data_1', 'rb') as f:
                file = dill.load(f)
                df1 = file[0]
                df2 = file[1]
                if type(file[2]) == list:
                    analyze_ms2 = True
                    fragments_dataframes = file[2]
                    if reanalysis[0] and ".".join(version.split('.')[:2]) != ".".join(file[3].split('.')[:2]):
                        print("Raw data files version incompatible with\ncurrent version (Current version: "+version+";\nRaw data version: "+file[3]+")")
                        return
                else:
                    if reanalysis[0] and ".".join(version.split('.')[:2]) != ".".join(file[2].split('.')[:2]):
                        print("Raw data files version incompatible with\ncurrent version (Current version: "+version+";\nRaw data version: "+file[3]+")")
                        return
                f.close()
    except:
        if reanalysis[0]:
            print("\nRaw data files not found. If you're not\nreanalyzing existing raw data, set\n'reanalysis' to 'no' in parameters before\nexecution or choose a different option in the\ncommand-line interface.\n")
        return
    if not reanalysis[0]:
        print("Analyzing raw data...")
    df1_refactor = []
    for i_i, i in enumerate(df2["Sample_Number"]): #QCs cutoff
        if analyze_ms2:
            df1_refactor.append({"Glycan" : [], "Adduct" : [], "mz" : [], "RT" : [], "AUC" : [], "PPM" : [], "S/N" : [], "Iso_Fitting_Score" : [], "Curve_Fitting_Score" : [], "Detected_Fragments" : []})
        else:
            df1_refactor.append({"Glycan" : [], "Adduct" : [], "mz" : [], "RT" : [], "AUC" : [], "PPM" : [], "S/N" : [], "Iso_Fitting_Score" : [], "Curve_Fitting_Score" : []})
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
                if df1[i_i]["Glycan"][j_j] != "Internal Standard":
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
                                    and fragments_dataframes[i_i]["Adduct"][k] == to_remove_adduct[k_k]
                                    and not unrestricted_fragments):
                                    del fragments_dataframes[i_i]["Glycan"][k]
                                    del fragments_dataframes[i_i]["Adduct"][k]
                                    del fragments_dataframes[i_i]["Fragment"][k]
                                    del fragments_dataframes[i_i]["Fragment_mz"][k]
                                    del fragments_dataframes[i_i]["Fragment_Intensity"][k]
                                    del fragments_dataframes[i_i]["RT"][k]
                                    del fragments_dataframes[i_i]["Precursor_mz"][k]
                                    del fragments_dataframes[i_i]["% TIC explained"][k]
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
                            and fragments_dataframes[j_j]["Adduct"][k] == to_remove_adduct[i_index]
                            and not unrestricted_fragments):
                            del fragments_dataframes[j_j]["Glycan"][k]
                            del fragments_dataframes[j_j]["Adduct"][k]
                            del fragments_dataframes[j_j]["Fragment"][k]
                            del fragments_dataframes[j_j]["Fragment_mz"][k]
                            del fragments_dataframes[j_j]["Fragment_Intensity"][k]
                            del fragments_dataframes[j_j]["RT"][k]
                            del fragments_dataframes[j_j]["Precursor_mz"][k]
                            del fragments_dataframes[j_j]["% TIC explained"][k] #QCs cutoff end
                            
    for i_i, i in enumerate(df1): #final arrangement for standard results print
        for j_j, j in enumerate(df1[i_i]["Adduct"]):
            for k_k, k in enumerate(df1[i_i]["RT"][j_j].split(", ")):
                if k != "":
                    df1_refactor[i_i]["Glycan"].append(df1[i_i]["Glycan"][j_j])
                    df1_refactor[i_i]["Adduct"].append(df1[i_i]["Adduct"][j_j])
                    df1_refactor[i_i]["mz"].append(df1[i_i]["mz"][j_j])
                    df1_refactor[i_i]["RT"].append(float(k))
                    df1_refactor[i_i]["AUC"].append(float(df1[i_i]["AUC"][j_j].split(", ")[k_k]))
                    df1_refactor[i_i]["PPM"].append(float(df1[i_i]["PPM"][j_j].split(", ")[k_k]))
                    df1_refactor[i_i]["S/N"].append(float(df1[i_i]["S/N"][j_j].split(", ")[k_k]))
                    df1_refactor[i_i]["Iso_Fitting_Score"].append(float(df1[i_i]["Iso_Fitting_Score"][j_j].split(", ")[k_k]))
                    df1_refactor[i_i]["Curve_Fitting_Score"].append(float(df1[i_i]["Curve_Fitting_Score"][j_j].split(", ")[k_k]))               
                    
    if analyze_ms2:
        to_keep = []
        for j_j, j in enumerate(df1_refactor):
            to_keep.append([])
            for k_k, k in enumerate(df1_refactor[j_j]["RT"]):
                if k_k == 0:
                    df1_refactor[j_j]["%_TIC_explained"] = []
                df1_refactor[j_j]["%_TIC_explained"].append(None)
                found = False
                for l_l, l in enumerate(fragments_dataframes[j_j]["Glycan"]):
                    if l == df1_refactor[j_j]["Glycan"][k_k] and fragments_dataframes[j_j]["Adduct"][l_l] == df1_refactor[j_j]["Adduct"][k_k]:
                        if abs(fragments_dataframes[j_j]["RT"][l_l] - k) <= rt_tolerance_frag:
                            found = True
                            to_keep[j_j].append(l_l)
                if found:
                    df1_refactor[j_j]["Detected_Fragments"].append("Yes")
                else:
                    df1_refactor[j_j]["Detected_Fragments"].append("No")
                    
        if not unrestricted_fragments:  #Filters fragments with RT outside the detected peaks range
            for i_i, i in enumerate(to_keep):
                for k_k in range(len(fragments_dataframes[i_i]["Glycan"])-1, -1, -1):
                    if k_k not in to_keep[i_i]:
                        del fragments_dataframes[i_i]["Glycan"][k_k]
                        del fragments_dataframes[i_i]["Adduct"][k_k]
                        del fragments_dataframes[i_i]["Fragment"][k_k]
                        del fragments_dataframes[i_i]["Fragment_mz"][k_k]
                        del fragments_dataframes[i_i]["Fragment_Intensity"][k_k]
                        del fragments_dataframes[i_i]["RT"][k_k]
                        del fragments_dataframes[i_i]["Precursor_mz"][k_k]
                        del fragments_dataframes[i_i]["% TIC explained"][k_k]
                        
        if len(reporter_ions) != 0: #reporter_ions filtering
            for i_i, i in enumerate(fragments_dataframes):
                to_remove = []
                temp_list_fragments = []
                temp_list_mz = []
                current_checking = ""
                for j_j, j in enumerate(fragments_dataframes[i_i]["Glycan"]):
                    to_check = j+"_"+fragments_dataframes[i_i]["Adduct"][j_j]+"_"+str(fragments_dataframes[i_i]["RT"][j_j])
                    if to_check != current_checking:
                        found = False
                        current_checking = to_check
                        temp_list_fragments = []
                        temp_list_mz = []
                        for k in range(j_j, len(fragments_dataframes[i_i]["Glycan"])):
                            to_check_2 = fragments_dataframes[i_i]["Glycan"][k]+"_"+fragments_dataframes[i_i]["Adduct"][k]+"_"+str(fragments_dataframes[i_i]["RT"][k])
                            if to_check == to_check_2:
                                temp_list_fragments.append(fragments_dataframes[i_i]["Fragment"][k])
                                temp_list_mz.append(fragments_dataframes[i_i]["Fragment_mz"][k])
                            else:
                                break
                        for k_k, k in enumerate(reporter_ions):
                            try:
                                current_reporter = float(k)
                            except:
                                current_reporter = k
                            if type(current_reporter) == float:
                                for l in temp_list_mz:
                                    if abs(current_reporter - l) <= 0.1:
                                        found = True
                                        break
                                if found:
                                    break
                            if type(current_reporter) == str:
                                for l in temp_list_fragments:
                                    if current_reporter == l.split("+")[0].split("-")[0].split("_")[0]:
                                        found = True
                                        break
                                if found:
                                    break
                    if not found and current_checking == to_check:
                        to_remove.append(j_j)
                if len(to_remove) != 0:
                    for k_k in sorted(to_remove, reverse = True):
                        del fragments_dataframes[i_i]["Glycan"][k_k]
                        del fragments_dataframes[i_i]["Adduct"][k_k]
                        del fragments_dataframes[i_i]["Fragment"][k_k]
                        del fragments_dataframes[i_i]["Fragment_mz"][k_k]
                        del fragments_dataframes[i_i]["Fragment_Intensity"][k_k]
                        del fragments_dataframes[i_i]["RT"][k_k]
                        del fragments_dataframes[i_i]["Precursor_mz"][k_k]
                        del fragments_dataframes[i_i]["% TIC explained"][k_k] #end of reporter ions filtering
                        
        for i_i, i in enumerate(fragments_dataframes): #% TIC explained calculation
            fragments_int_sum = 0
            current_checking = ""
            for j_j, j in enumerate(fragments_dataframes[i_i]["Glycan"]):
                to_check = j+"_"+fragments_dataframes[i_i]["Adduct"][j_j]+"_"+str(fragments_dataframes[i_i]["RT"][j_j])
                if to_check != current_checking:
                    current_checking = to_check
                    fragments_int_sum = 0
                    for k in range(j_j, len(fragments_dataframes[i_i]["Glycan"])):
                        to_check_2 = fragments_dataframes[i_i]["Glycan"][k]+"_"+fragments_dataframes[i_i]["Adduct"][k]+"_"+str(fragments_dataframes[i_i]["RT"][k])
                        if to_check == to_check_2:
                            fragments_int_sum += fragments_dataframes[i_i]["Fragment_Intensity"][k]
                        else:
                            break
                if current_checking == to_check and fragments_dataframes[i_i]["% TIC explained"] != 0:
                    fragments_dataframes[i_i]["% TIC explained"][j_j] = float("%.3f" % round(fragments_int_sum/fragments_dataframes[i_i]["% TIC explained"][j_j], 3)) #end of annotated_peaks ratio calculation
                    
        for i_i, i in enumerate(df1_refactor): #start of ms2 score calculation (at the moment its just % TIC explained)
            for j_j, j in enumerate(i['Glycan']):
                if i['Detected_Fragments'][j_j] == 'Yes':
                    glycan = j
                    adduct = i['Adduct'][j_j]
                    rt = i['RT'][j_j]
                    list_tics_explained = []
                    for k_k, k in enumerate(fragments_dataframes):
                        for l_l, l in enumerate(k['Glycan']):
                            if glycan == l and adduct == k['Adduct'][l_l] and abs(rt-k['RT'][l_l]) <= rt_tolerance_frag:
                                list_tics_explained.append(k['% TIC explained'][l_l])
                    if len(list_tics_explained) > 0:
                        i['%_TIC_explained'][j_j] = float('%.2f' % round(max(list_tics_explained), 2))
                else:
                    i['%_TIC_explained'][j_j] = 0.0
                    
        fragments_refactor_dataframes = [] #here it re-structures the fragments data to an updated format
        for i_i, i in enumerate(fragments_dataframes): #moving through samples
            fragments_refactor_dataframes.append({})
            current_scan = ""
            glycan_number = -1
            for j_j, j in enumerate(i['Glycan']):
                scan_name = j+"_"+i['Adduct'][j_j]+"_"+str(i['RT'][j_j])
                if scan_name != current_scan:
                    glycan_number+=1
                    current_scan = scan_name
                    fragments_refactor_dataframes[i_i]['Glycan_'+str(glycan_number)+':'] = [j, 'RT_'+str(glycan_number)+':', i['RT'][j_j], 'Fragment:']
                    fragments_refactor_dataframes[i_i]['Adduct_'+str(glycan_number)+':'] = [i['Adduct'][j_j], '% TIC assigned_'+str(glycan_number)+':', i['% TIC explained'][j_j]*100, 'm/z:']
                    fragments_refactor_dataframes[i_i]['m/z_'+str(glycan_number)+':'] = [i['Precursor_mz'][j_j], None, None, 'Intensity:']
                    fragments_refactor_dataframes[i_i]['Glycan_'+str(glycan_number)+':'].append(i['Fragment'][j_j])
                    fragments_refactor_dataframes[i_i]['Adduct_'+str(glycan_number)+':'].append(i['Fragment_mz'][j_j])
                    fragments_refactor_dataframes[i_i]['m/z_'+str(glycan_number)+':'].append(i['Fragment_Intensity'][j_j])
                else:
                    fragments_refactor_dataframes[i_i]['Glycan_'+str(glycan_number)+':'].append(i['Fragment'][j_j])
                    fragments_refactor_dataframes[i_i]['Adduct_'+str(glycan_number)+':'].append(i['Fragment_mz'][j_j])
                    fragments_refactor_dataframes[i_i]['m/z_'+str(glycan_number)+':'].append(i['Fragment_Intensity'][j_j])
        for i in fragments_refactor_dataframes: #makes all lists in the dataframe equal size so it can be ported to excel
            for j in i:
                while len(i[j]) < 1000:
                    i[j].append(None)
                    
    ambiguity_count = [] #ambiguity indicator
    for i_i, i in enumerate(df1_refactor):
        ambiguity_count.append(0)
        i['Ambiguity'] = []
        for j in i['Glycan']:
            i['Ambiguity'].append([])
        for j_j, j in enumerate(i['Glycan']):
            glycan_j = j+'_'+i['Adduct'][j_j]
            for k_k, k in enumerate(i['Glycan'][j_j+1:]):
                k_k = j_j+k_k+1
                glycan_k = k+'_'+i['Adduct'][k_k]
                if j != k and i['mz'][j_j] == i['mz'][k_k]:
                    ambiguity_count[i_i] += 1
                    i['Ambiguity'][j_j].append(i['Glycan'][k_k]+'_'+i['Adduct'][k_k])
                    i['Ambiguity'][k_k].append(i['Glycan'][j_j]+'_'+i['Adduct'][j_j])
        for j_j, j in enumerate(i['Ambiguity']):
            if len(j) > 0:
                i['Ambiguity'][j_j] = ', '.join(j)
            else:
                i['Ambiguity'][j_j] = 'No'
            
    total_dataframes = [] #total glycans AUC dataframe
    for i_i, i in enumerate(df1_refactor):
        total_dataframes.append({"Glycan": [], "RT": [], "AUC": []})
        current_glycan = ""
        RTs = []
        AUCs = []
        for j_j in range(len(i["Glycan"])):
            j = i["Glycan"][j_j]
            found = False
            if j == current_glycan:
                for k_k, k in enumerate(RTs):
                    if abs(i["RT"][j_j] - k) <= rt_tolerance:
                        RTs[k_k] = (k+i["RT"][j_j])/2
                        AUCs[k_k] = AUCs[k_k]+i["AUC"][j_j]
                        found = True
                        break
                if not found:
                    RTs.append(i["RT"][j_j])
                    AUCs.append(i["AUC"][j_j])
            if j != current_glycan:
                if j_j != 0:
                    for k_k, k in enumerate(RTs):
                            total_dataframes[i_i]["Glycan"].append(current_glycan)
                            total_dataframes[i_i]["RT"].append(k)
                            total_dataframes[i_i]["AUC"].append(AUCs[k_k])
                    RTs = []
                    AUCs = []
                current_glycan = j
                RTs.append(i["RT"][j_j])
                AUCs.append(i["AUC"][j_j])
            if j_j == len(i["Glycan"])-1:
                for k_k, k in enumerate(RTs):
                    total_dataframes[i_i]["Glycan"].append(current_glycan)
                    total_dataframes[i_i]["RT"].append(k)
                    total_dataframes[i_i]["AUC"].append(AUCs[k_k])
                RTs = []
                AUCs = [] #total glycans AUC dataframe
    
    arranged_total_dataframes = []
    for i_i, i in enumerate(total_dataframes):
        arranged_total_dataframes.append({'Glycan' : [], 'RT' : [], 'AUC' : []})
        list_of_glycans = []
        glycans_dict_rt = {}
        glycans_dict_AUC = {}
        current_glycan = ''
        for j_j, j in enumerate(i['Glycan']):
            if j != current_glycan:
                current_glycan = j
                list_of_glycans.append(j)
                glycans_dict_rt[j] = []
                glycans_dict_AUC[j] = []
                glycans_dict_rt[j].append(i['RT'][j_j])
                glycans_dict_AUC[j].append(i['AUC'][j_j])
            else:
                glycans_dict_rt[j].append(i['RT'][j_j])
                glycans_dict_AUC[j].append(i['AUC'][j_j])
        if "Internal Standard" in list_of_glycans:
            highest = glycans_dict_AUC["Internal Standard"].index(max(glycans_dict_AUC["Internal Standard"]))
            to_remove = []
            for j_j, j in enumerate(glycans_dict_AUC["Internal Standard"]):
                if j_j != highest:
                    to_remove.append(j_j)
            for j in sorted(to_remove, reverse = True):
                del glycans_dict_rt["Internal Standard"][j]
                del glycans_dict_AUC["Internal Standard"][j]
        list_of_glycans = sorted(list_of_glycans)
        for j_j, j in enumerate(list_of_glycans):
            for k in range(len(glycans_dict_rt[j])):
                arranged_total_dataframes[i_i]['Glycan'].append(j)
            zipped = zip(glycans_dict_rt[j], glycans_dict_AUC[j])
            zipped = sorted(zipped)
            current_RTs, current_AUCs = zip(*zipped)
            arranged_total_dataframes[i_i]['RT'] += list(current_RTs)
            arranged_total_dataframes[i_i]['AUC'] += list(current_AUCs)
    total_dataframes = arranged_total_dataframes
    
    compositions_dataframes = [] #compositions_dataframes
    for i_i, i in enumerate(total_dataframes):
        compositions_dataframes.append({'Glycan' : [], 'AUC' : []})
        glycans_lib = {}
        for j_j, j in enumerate(i['Glycan']):
            if j not in glycans_lib.keys():
                glycans_lib[j] = i['AUC'][j_j]
            else:
                glycans_lib[j] += i['AUC'][j_j]
        for j_j, j in enumerate(glycans_lib):
            compositions_dataframes[i_i]['Glycan'].append(j)
            compositions_dataframes[i_i]['AUC'].append(glycans_lib[j])
    
    if nglycan: #if N-Glycans, determines its class
        glycan_class = {}
        for i_i, i in enumerate(total_dataframes):
            for j_j, j in enumerate(i['Glycan']):
                if j != 'Internal Standard' and j not in glycan_class.keys():
                    comp = General_Functions.form_to_comp(j)
                    if comp['N'] == 2 and comp['H'] <= 3:
                        glycan_class[j] = 'Paucimannose'
                        continue
                    if comp['H'] > comp['N']+1 and comp['N'] > 2:
                        glycan_class[j] = 'Hybrid'
                        continue
                    if comp['N'] == 2 and comp['H'] > 3:
                        glycan_class[j] = 'High-Mannose'
                        continue
                    else:
                        glycan_class[j] = 'Complex'
                if j == 'Internal Standard':
                    glycan_class[j] = j
        for i_i, i in enumerate(total_dataframes):
            total_dataframes[i_i]['Class'] = []
            for j_j, j in enumerate(i['Glycan']):
                i['Class'].append(glycan_class[j])
        for i_i, i in enumerate(compositions_dataframes):
            compositions_dataframes[i_i]['Class'] = []
            for j_j, j in enumerate(i['Glycan']):
                i['Class'].append(glycan_class[j])
    
    #hook for alignment tool. it'll use the total_dataframes (total_glycans) 
    if align_chromatograms:
        if len(total_dataframes) > 1:
            aligned_total_glycans = align_assignments(total_dataframes, 'total_glycans')
            total_dataframes = aligned_total_glycans[0]
            df2["Average Delta t"] = []
            for i_i, i in enumerate(aligned_total_glycans[1]):
                temp_list = []
                for k_k, k in enumerate(i):
                    temp_list.append(i[k][0])
                if len(temp_list) != 0:
                    average_delta = float("%.2f" % round(sum(temp_list)/len(temp_list), 2))
                    df2["Average Delta t"].append(average_delta)
                else:
                    if i_i == aligned_total_glycans[2]:
                        df2["Average Delta t"].append("Reference Sample")
                    else:
                        df2["Average Delta t"].append("No peaks to align")
    #hook for alignment tool. it'll use the total_dataframes (total_glycans)         
                
    glycans_count = [] #glycans composition counts
    for i_i, i in enumerate(df1_refactor):
        current_glycan = ""
        glycans_count.append(0)
        for j_j, j in enumerate(i["Glycan"]):
            if j == 'Internal Standard':
                continue
            if j != current_glycan:
                glycans_count[i_i]+= 1
                current_glycan = j
    df2["MS1_Glycans_Compositions"] = glycans_count
    if analyze_ms2:
        glycans_count = []
        for i_i, i in enumerate(fragments_dataframes):
            current_glycan = ""
            glycans_count.append(0)
            for j_j, j in enumerate(i["Glycan"]):
                if j != current_glycan:
                    glycans_count[i_i]+= 1
                    current_glycan = j
        df2["MS2_Glycans_Compositions"] = glycans_count #end of glycans counts
        
    df2["Ambiguities"] = ambiguity_count
        
    all_glycans_list = [] #here it makes a list of ALL the glycans found for use in other parts of the data arrangement workflow
    for i in total_dataframes:
        for j_j, j in enumerate(i["RT"]):
            if i["Glycan"][j_j] == "Internal Standard":
                continue
            found = False
            if len(all_glycans_list) > 0:
                for k in all_glycans_list:
                    splitted_glycan = k.split("_")
                    if i["Glycan"][j_j] == splitted_glycan[0] and abs(j-float(splitted_glycan[-1])) <= rt_tolerance:
                        found = True
            if found:
                continue
            if len(all_glycans_list) == 0:
                all_glycans_list.append(i["Glycan"][j_j]+"_"+"%.2f" % j)
                continue
            if not found:
                all_glycans_list.append(i["Glycan"][j_j]+"_"+"%.2f" % j)
                
    if plot_metaboanalyst[0]: #start of metaboanalyst plot
        print("Creating file for metaboanalyst plotting...", end="", flush=True)
        with open(save_path+begin_time+"_metaboanalyst_data.csv", "w") as f:
            samples_line = ["Sample"]
            for i_i, i in enumerate(df2["File_Name"]):
                samples_line.append(i)
            groups_line = ["Group"]
            if len(plot_metaboanalyst[1]) > 0:
                for i in df2["File_Name"]:
                    found = False
                    for j in plot_metaboanalyst[1]:
                        if j in i:
                            found = True
                            groups_line.append(j)
                            break
                    if not found:
                        groups_line.append("Ungrouped")
            else:
                for i in df2["File_Name"]:
                    groups_line.append("Ungrouped")
            found_int_std = False
            for i in total_dataframes:
                if "Internal Standard" in i["Glycan"]:
                    found_int_std = True
                    break
            if found_int_std:
                with open(save_path+begin_time+"_metaboanalyst_data_normalized.csv", "a") as g:
                    g.write(",".join(samples_line)+"\n")
                    g.write(",".join(groups_line)+"\n")
                    g.close()
                is_areas = []
                for i in total_dataframes:
                    if "Internal Standard" in i["Glycan"]:
                        temp_areas = []
                        for j_j, j in enumerate(i["Glycan"]):
                            if j == "Internal Standard":
                                temp_areas.append(i["AUC"][j_j])
                        is_areas.append(max(temp_areas))
                    else:
                        is_areas.append(0.0)
            f.write(",".join(samples_line)+"\n")
            f.write(",".join(groups_line)+"\n")
            for i in all_glycans_list:
                glycan_line = []
                glycan_line_IS = []
                i_splitted = i.split("_")
                glycan_line_IS.append(i)
                glycan_line.append(i)
                for j_j, j in enumerate(total_dataframes): #moving through samples
                    found = False
                    temp_AUC = 0
                    for k_k, k in enumerate(j["Glycan"]):
                        if k == "Internal Standard":
                            continue
                        if k == i_splitted[0] and abs(j["RT"][k_k] - float(i_splitted[-1])) <= rt_tolerance:
                            found = True
                            if "Internal Standard" in j["Glycan"]:
                                if is_areas[j_j] > 0.0:
                                    temp_AUC_IS = j["AUC"][k_k]/is_areas[j_j]
                                else:
                                    temp_AUC_IS = 0.0
                                temp_AUC+= j["AUC"][k_k]
                            else:
                                temp_AUC += j["AUC"][k_k]
                    if found:
                        if "Internal Standard" in j["Glycan"]:
                            glycan_line_IS.append(str(temp_AUC_IS))
                        glycan_line.append(str(temp_AUC))
                        continue
                    if not found:
                        if "Internal Standard" in j["Glycan"]:
                            glycan_line_IS.append("0.0")
                        glycan_line.append("0.0")
                        continue
                if found_int_std:
                    with open(save_path+begin_time+"_metaboanalyst_data_normalized.csv", "a") as g:
                        g.write(",".join(glycan_line_IS)+"\n")
                        g.close()
                f.write(",".join(glycan_line)+"\n")
            f.close()
        if compositions:
            total_glycans_compositions = []
            with open(save_path+begin_time+"_metaboanalyst_data_compositions.csv", "w") as f:
                found_int_std = False
                for i in compositions_dataframes:
                    if "Internal Standard" in i["Glycan"]:
                        found_int_std = True
                        break
                if found_int_std:
                    with open(save_path+begin_time+"_metaboanalyst_data_compositions_normalized.csv", "w") as g:
                        g.write(",".join(samples_line)+"\n")
                        g.write(",".join(groups_line)+"\n")
                        g.close()
                f.write(",".join(samples_line)+"\n")
                f.write(",".join(groups_line)+"\n")
                for i_i, i in enumerate(compositions_dataframes):
                    for j_j, j in enumerate(i['Glycan']):
                        if j not in total_glycans_compositions and j != 'Internal Standard':
                            total_glycans_compositions.append(j)
                for i_i, i in enumerate(sorted(total_glycans_compositions)):
                    glycan_line = [i]
                    glycan_line_IS = [i]
                    for j_j, j in enumerate(compositions_dataframes):
                        if i in j['Glycan']:
                            glycan_line.append(str(j['AUC'][j['Glycan'].index(i)]))
                            if 'Internal Standard' in j['Glycan']:
                                glycan_line_IS.append(str(j['AUC'][j['Glycan'].index(i)]/j['AUC'][j['Glycan'].index('Internal Standard')]))
                        else:
                            glycan_line.append('0.0')
                            glycan_line_IS.append('0.0')
                    f.write(",".join(glycan_line)+"\n")
                    if found_int_std:
                        with open(save_path+begin_time+"_metaboanalyst_data_compositions_normalized.csv", "a") as g:
                            g.write(",".join(glycan_line_IS)+"\n")
                            g.close()
                f.close()
                
        print("Done!") #end of metaboanalyst plot
    
    #start of excel data printing
    df2 = DataFrame(df2) 
    
    with ExcelWriter(save_path+begin_time+'_Results_'+str(max_ppm)+'_'+str(iso_fit_score)+'_'+str(curve_fit_score)+'_'+str(sn)+'.xlsx') as writer:
        print("Creating results file...", end="", flush=True)
        df2.to_excel(writer, sheet_name="Index references", index = False)
        for i_i, i in enumerate(df1_refactor):
            result_df = DataFrame(i)
            result_df.to_excel(writer, sheet_name="Sample_"+str(i_i), index = False)
            total_aucs_df = DataFrame(total_dataframes[i_i])
            total_aucs_df.to_excel(writer, sheet_name="Sample_"+str(i_i)+"_Total_AUCs", index = False)
            if compositions:
                compositions_df = DataFrame(compositions_dataframes[i_i])
                compositions_df.to_excel(writer, sheet_name="Sample_"+str(i_i)+"_Compositions_AUCs", index = False)
            if analyze_ms2:
                if len(fragments_dataframes[i_i]["Glycan"]) > 0:
                    fragments_df = DataFrame(fragments_refactor_dataframes[i_i])
                    fragments_df.to_excel(writer, sheet_name="Sample_"+str(i_i)+"_Fragments", index = False)
    del df1
    del result_df
    del total_dataframes
    del total_aucs_df
    if compositions:
        del compositions_dataframes
        del compositions_df
    
    if analyze_ms2:
        if len(fragments_dataframes[i_i]["Glycan"]) > 0:
            del fragments_df
        del fragments_refactor_dataframes
        del fragments_dataframes
    print("Done!")
    
    found_eic_raw_dataframes = [] #This creates a file with only the found glycan's EIC
    found_eic_processed_dataframes = []
    if sneakpeek[0]:
        with open(save_path+'results5_'+str(sneakpeek[1]), 'rb') as f:
            raw_eic_dataframes = dill.load(f)
            f.close()
    else:
        with open(save_path+'raw_data_5', 'rb') as f:
            raw_eic_dataframes = dill.load(f)
            f.close()
    if sneakpeek[0]:
        with open(save_path+'results3_'+str(sneakpeek[1]), 'rb') as f:
            smoothed_eic_dataframes = dill.load(f)
            f.close()
    else:
        with open(save_path+'raw_data_3', 'rb') as f:
            smoothed_eic_dataframes = dill.load(f)
            f.close()
            
    if align_chromatograms:        
        if len(df2['Sample_Number']) > 1: #aligns the chromatograms, may take some time (around 1 minute per sample, depending on run length)
            print("Aligning chromatograms...", end='', flush=True)
            smoothed_eic_dataframes = align_assignments(smoothed_eic_dataframes, 'chromatograms', aligned_total_glycans[1])
            print("Done!")
        
    for i_i, i in enumerate(df1_refactor): #this selects only the found glycans to draw their EIC
        found_eic_raw_dataframes.append({})
        found_eic_raw_dataframes[i_i]['RTs_'+str(i_i)] = raw_eic_dataframes[i_i]['RTs_'+str(i_i)]
        found_eic_processed_dataframes.append({})
        found_eic_processed_dataframes[i_i]['RTs_'+str(i_i)] = smoothed_eic_dataframes[i_i]['RTs_'+str(i_i)]
        for j_j, j in enumerate(i['Glycan']):
            query = j+"+"+i['Adduct'][j_j]+" - "+str(i['mz'][j_j])
            try:
                found_eic_raw_dataframes[i_i][query] = raw_eic_dataframes[i_i][query]
                found_eic_processed_dataframes[i_i][query] = smoothed_eic_dataframes[i_i][query]
            except:
                pass
    found_eic_processed_dataframes_simplified = [] #combines adducts EICs
    for i_i, i in enumerate(found_eic_processed_dataframes):
        glycans = []
        found_eic_processed_dataframes_simplified.append({})
        for j_j, j in enumerate(i):
            if j_j == 0:
                found_eic_processed_dataframes_simplified[i_i][j] = i[j]
            else:
                glycan = j.split("+")[0]
                if glycan not in glycans:
                    glycans.append(glycan)
                    adduct_no = 1
                    for k_k, k in enumerate(i):
                        if k.split("+")[0] == glycan:
                            if adduct_no == 1:
                                found_eic_processed_dataframes_simplified[i_i][glycan] = i[k]
                                adduct_no += 1
                            else:
                                for l_l, l in enumerate(i[k]):
                                    found_eic_processed_dataframes_simplified[i_i][glycan][l_l] += l
                                adduct_no += 1
    found_eic_processed_dataframes = found_eic_processed_dataframes_simplified
    
    print("Creating data plotting files...", end='', flush=True)
    with ExcelWriter(save_path+begin_time+'_Found_Glycans_EICs.xlsx') as writer:
        for i_i, i in enumerate(found_eic_raw_dataframes):
            found_eic_raw_dataframes_df = DataFrame(i)
            found_eic_processed_dataframes_df = DataFrame(found_eic_processed_dataframes[i_i])
            found_eic_raw_dataframes_df.to_excel(writer, sheet_name="RAW_Sample_"+str(i_i), index = False)
            found_eic_processed_dataframes_df.to_excel(writer, sheet_name="Processed_Sample_"+str(i_i), index = False)
        df2.to_excel(writer, sheet_name="Index references", index = False)
    del found_eic_processed_dataframes
    del found_eic_processed_dataframes_df
    del found_eic_raw_dataframes
    del found_eic_raw_dataframes_df
    
    if output_isotopic_fittings:
        with open(save_path+'raw_data_6', 'rb') as f: #start of isotopic fits output
            isotopic_fits_dataframes = dill.load(f)
            f.close()
        isotopic_fits_dataframes_arranged = []
        for i_i, i in enumerate(isotopic_fits_dataframes): #sample
            temp_fits_dataframes = {}
            for j_j, j in enumerate(i): #glycan
                temp_fits_dataframes[j] = {}
                for k_k, k in enumerate(i[j]): #peaks of glycan
                    temp_fits_dataframes[j]['RT_'+str(k_k)+':'] = []
                    temp_fits_dataframes[j]['Score_'+str(k_k)+':'] = []
                    temp_fits_dataframes[j]['fit_'+str(k_k)] = []
                    temp_fits_dataframes[j]['RT_'+str(k_k)+':'].append(k)
                    temp_fits_dataframes[j]['Score_'+str(k_k)+':'].append(isotopic_fits_dataframes[i_i][j][k][3])
                    temp_fits_dataframes[j]['fit_'+str(k_k)].append(None)
                    temp_fits_dataframes[j]['RT_'+str(k_k)+':'].append('mz:')
                    temp_fits_dataframes[j]['Score_'+str(k_k)+':'].append('Ideal:')
                    temp_fits_dataframes[j]['fit_'+str(k_k)].append('Actual:')
                    temp_fits_dataframes[j]['RT_'+str(k_k)+':'] = temp_fits_dataframes[j]['RT_'+str(k_k)+':']+isotopic_fits_dataframes[i_i][j][k][0]
                    temp_fits_dataframes[j]['Score_'+str(k_k)+':'] = temp_fits_dataframes[j]['Score_'+str(k_k)+':']+isotopic_fits_dataframes[i_i][j][k][1]
                    temp_fits_dataframes[j]['fit_'+str(k_k)] = temp_fits_dataframes[j]['fit_'+str(k_k)]+isotopic_fits_dataframes[i_i][j][k][2]
                    while len(temp_fits_dataframes[j]['RT_'+str(k_k)+':']) < 1000:
                        temp_fits_dataframes[j]['RT_'+str(k_k)+':'].append(None)
                        temp_fits_dataframes[j]['Score_'+str(k_k)+':'].append(None)
                        temp_fits_dataframes[j]['fit_'+str(k_k)].append(None)
            isotopic_fits_dataframes_arranged.append(temp_fits_dataframes)
        for i_i, i in enumerate(isotopic_fits_dataframes_arranged):
            with ExcelWriter(save_path+begin_time+'_Isotopic_Fits_Sample_'+str(i_i)+'.xlsx') as writer:
                for j_j, j in enumerate(i): #navigating glycans
                    isotopic_fits_df = DataFrame(isotopic_fits_dataframes_arranged[i_i][j])
                    isotopic_fits_df.to_excel(writer, sheet_name=j, index = False)
        del isotopic_fits_dataframes
        del isotopic_fits_dataframes_arranged
        del isotopic_fits_df
    
    if (reanalysis[1] and reanalysis[0]) or not reanalysis[0]:
        # if sneakpeek[0]: #commented for now, it's just too much redundant information... smoothed eic will output as processed now
            # with open(save_path+'results2_'+str(sneakpeek[1]), 'rb') as f:
                # eic_dataframes = dill.load(f)
                # f.close()
        # else:
            # with open(save_path+'raw_data_2', 'rb') as f:
                # eic_dataframes = dill.load(f)
                # f.close()
        # with ExcelWriter(save_path+begin_time+'_processed_EIC_Plot_Data.xlsx') as writer:
            # for i_i, i in enumerate(eic_dataframes):
                # eic_df = DataFrame(i)
                # eic_df.to_excel(writer, sheet_name="Sample_"+str(i_i), index = False)
            # df2.to_excel(writer, sheet_name="Index references", index = False)
        # del eic_dataframes
        # del eic_df
        with ExcelWriter(save_path+begin_time+'_processed_EIC_Plot_Data.xlsx') as writer: #smoothed eic, now changed to processed to avoid TMI
            for i_i, i in enumerate(smoothed_eic_dataframes):
                smoothed_eic_df = DataFrame(i)
                smoothed_eic_df.to_excel(writer, sheet_name="Sample_"+str(i_i), index = False)
            df2.to_excel(writer, sheet_name="Index references", index = False)
        del smoothed_eic_dataframes
        del smoothed_eic_df
        
        with ExcelWriter(save_path+begin_time+'_raw_EIC_Plot_Data.xlsx') as writer:
            for i_i, i in enumerate(raw_eic_dataframes):
                raw_eic_df = DataFrame(i)
                raw_eic_df.to_excel(writer, sheet_name="Sample_"+str(i_i), index = False)
            df2.to_excel(writer, sheet_name="Index references", index = False)
        del raw_eic_dataframes
        del raw_eic_df
        
        if sneakpeek[0]:
            with open(save_path+'results4_'+str(sneakpeek[1]), 'rb') as f:
                curve_fitting_dataframes = dill.load(f)
                f.close()
        else:
            with open(save_path+'raw_data_4', 'rb') as f:
                curve_fitting_dataframes = dill.load(f)
                f.close()
        with ExcelWriter(save_path+begin_time+'_curve_fitting_Plot_Data.xlsx') as writer:
            for i_i, i in enumerate(curve_fitting_dataframes):
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
                    curve_df.to_excel(writer, sheet_name="Sample_"+str(i_i), index = False)
            df2.to_excel(writer, sheet_name="Index references", index = False)
        del curve_fitting_dataframes
        del curve_df
        print("Done!")
    elif reanalysis[0] and not reanalysis[1]:
        del smoothed_eic_dataframes
        del raw_eic_dataframes
        print("Done!")

def arrange_raw_data(analyzed_data,
                     samples_names,
                     multithreaded_analysis,
                     multithreaded_execution,
                     analyze_ms2,
                     save_path): ##Complete
    '''Arrange the raw results data into pickled files to be processed by output_filtered_data.

    Parameters
    ----------
    analyzed_data : tuple
        Tuple containing multiple informations about the glycans analysis, as outputted by
        analyze_files or analyze_ms2.
        
    samples_names : list
        List of samples names extracted from file names.
        
    multithreaded_analysis : tuple
        A tuple with two indexes: the first one contains whether or not this execution
        is to produce the tools required for multithreaded analysis and the second one 
        containing the amount of desired threads to split the execution in.
        
    multithread_execution : tuple
        A tuple with three indexes: the first one contains whether or not this execution
        is part of a multithreaded run, the second one contains the current thread number and
        the third one is the total amount of threads.
        
    analyze_ms2 : tuple
        A tuple with two indexes: The first one indicates whether to analyze ms2 data and the
        second one indicates whether ms2 data should be forced to fit glycans composition.
    
    save_path : string
        A string containing the path to the working directory of the script.
        
    Uses
    ----
    dill.dump : None
        Pickle the current state of __main__ or another module to a file.
    
    dill.load : Module object
        Update the selected module (default is __main__) with the state saved at filename.
    '''
    begin_time = datetime.datetime.now()
    print('Arranging raw data...', end='', flush = True)
    df1 = []
    df2 = {"Sample_Number" : [], "File_Name" : [], "Average_Noise_Level" : []}
    raw_eic_dataframes = []
    eic_dataframes = []
    smoothed_eic_dataframes = []
    curve_fitting_dataframes = []
    isotopic_fits_dataframes = []
    if analyze_ms2:
        fragments_dataframes = []
    for i_i, i in enumerate(samples_names):
        isotopic_fits_dataframes.append({})
        raw_eic_dataframes.append({})
        eic_dataframes.append({})
        smoothed_eic_dataframes.append({})
        temp_eic_rt = []
        for j in analyzed_data[1][i_i]:
            temp_eic_rt.append(float("%.4f" % round(j, 4)))
        raw_eic_dataframes[i_i]['RTs_'+str(i_i)] = temp_eic_rt
        eic_dataframes[i_i]['RTs_'+str(i_i)] = temp_eic_rt
        smoothed_eic_dataframes[i_i]['RTs_'+str(i_i)] = temp_eic_rt
        curve_fitting_dataframes.append({})
        df2["Sample_Number"].append(i_i)
        df2["File_Name"].append(i)
        df2["Average_Noise_Level"].append(float("%.1f" % round(analyzed_data[2][i_i],1)))
        if analyze_ms2:
            df1.append({"Glycan" : [], "Adduct" : [], "mz" : [], "RT" : [], "AUC" : [], "PPM" : [], "S/N" : [], "Iso_Fitting_Score" : [], "Curve_Fitting_Score" : [], "Detected_Fragments" : []})
            fragments_dataframes.append({"Glycan" : [], "Adduct" : [], "Precursor_mz" : [], "Fragment" : [], "Fragment_mz" : [], "Fragment_Intensity" : [], "RT" : [], "% TIC explained" : []})
        else:
            df1.append({"Glycan" : [], "Adduct" : [], "mz" : [], "RT" : [], "AUC" : [], "PPM" : [], "S/N" : [], "Iso_Fitting_Score" : [], "Curve_Fitting_Score" : []})
    for i_i, i in enumerate(analyzed_data[0]): #i = glycan (key)
        for j_j, j in enumerate(analyzed_data[0][i]['Adducts_mz_data']): #j = adduct (key)
            for k_k, k in enumerate(analyzed_data[0][i]['Adducts_mz_data'][j]): #k = sample number (key)
                isotopic_fits_dataframes[k_k][i+'_'+j] = analyzed_data[0][i]['Adducts_mz_data'][j][k][4]
                temp_eic_int = []
                for l in analyzed_data[0][i]['Adducts_mz_data'][j][k][3]:
                    temp_eic_int.append(int(l))
                raw_eic_dataframes[k_k][str(i)+'+'+str(j)+' - '+str(float("%.4f" % round(analyzed_data[0][i]['Adducts_mz'][j], 4)))] = temp_eic_int
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
                                fragments_dataframes[k_k]["% TIC explained"].append(float(m[7]))
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
    biggest_len = 10000
    for i in curve_fitting_dataframes:
        for j in i:
            if len(i[j]) < biggest_len:
                for k in range(biggest_len-len(i[j])):
                    i[j].append(None)
    if multithreaded_execution[0]:
#        sleep((multithreaded_execution[1])*(30/multithreaded_execution[2]))
        with open(save_path+'results1_'+str(multithreaded_execution[1]), 'wb') as f:
            if analyze_ms2:
                dill.dump([df1, df2, fragments_dataframes, multithreaded_execution[2]], f)
            else:
                dill.dump([df1, df2, multithreaded_execution[2]], f)
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
        with open(save_path+'results5_'+str(multithreaded_execution[1]), 'wb') as f:
            dill.dump(raw_eic_dataframes, f)
            f.close()
        with open(save_path+'results6_'+str(multithreaded_execution[1]), 'wb') as f:
            dill.dump(isotopic_fits_dataframes, f)
            f.close()
        p1 = pathlib.Path(save_path+'glycogenius_'+str(multithreaded_execution[1])+'.py')
        p1.unlink(missing_ok=True)
        results1_list = []
        results2_list = []
        results3_list = []
        results4_list = []
        results5_list = []
        results6_list = []
        for i in range(multithreaded_execution[2]):
            results1_list.append('results1_'+str(i))
            results2_list.append('results2_'+str(i))
            results3_list.append('results3_'+str(i))
            results4_list.append('results4_'+str(i))
            results5_list.append('results5_'+str(i))
            results6_list.append('results6_'+str(i))
        dir_list = os.listdir(save_path)
        last = True
        for i in results1_list:
            if i not in dir_list:
                last = False
        if last:
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
                with open(results5_list[i_i], 'rb') as f:
                    raw_eic_dataframes_import = dill.load(f)
                    f.close()
                with open(results6_list[i_i], 'rb') as f:
                    isotopic_fits_dataframes_import = dill.load(f)
                    f.close()
                if i_i == 0:
                    df1 = df1_import
                    df2 = df2_import
                    raw_eic_dataframes = raw_eic_dataframes_import
                    eic_dataframes = eic_dataframes_import
                    smoothed_eic_dataframes = smoothed_eic_dataframes_import
                    curve_fitting_dataframes = curve_fitting_dataframes_import
                    isotopic_fits_dataframes = isotopic_fits_dataframes_import
                    if analyze_ms2:
                        fragments_dataframes = fragments_dataframes_import
                else:
                    for j_j, j in enumerate(df1_import):
                        for k in j:
                            df1[j_j][k] = df1[j_j][k] + df1_import[j_j][k]
                    for j_j, j in enumerate(eic_dataframes_import):
                        for k in j:
                            if k[:4] != 'RTs_':
                                raw_eic_dataframes[j_j][k] = raw_eic_dataframes_import[j_j][k]
                                eic_dataframes[j_j][k] = eic_dataframes_import[j_j][k]
                                smoothed_eic_dataframes[j_j][k] = smoothed_eic_dataframes_import[j_j][k]
                    for j_j, j in enumerate(curve_fitting_dataframes_import):
                        for k in j:
                            curve_fitting_dataframes[j_j][k] = curve_fitting_dataframes_import[j_j][k]
                    for j_j, j in enumerate(isotopic_fits_dataframes_import):
                        for k in j:
                            isotopic_fits_dataframes[j_j][k] = j[k]
                    if analyze_ms2:
                        for j_j, j in enumerate(fragments_dataframes_import):
                            for k in j:
                                fragments_dataframes[j_j][k] = fragments_dataframes[j_j][k] + fragments_dataframes_import[j_j][k]
            with open(save_path+'raw_data_1', 'wb') as f:
                if analyze_ms2:
                    dill.dump([df1, df2, fragments_dataframes, version], f)
                else:
                    dill.dump([df1, df2, version], f)
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
            with open(save_path+'raw_data_5', 'wb') as f:
                dill.dump(raw_eic_dataframes, f)
                f.close()
            with open(save_path+'raw_data_6', 'wb') as f:
                dill.dump(isotopic_fits_dataframes, f)
                f.close()
            for i in range(multithreaded_analysis[1]):
                p1 = pathlib.Path(save_path+'results1_'+str(i))
                p2 = pathlib.Path(save_path+'results2_'+str(i))
                p3 = pathlib.Path(save_path+'results3_'+str(i))
                p4 = pathlib.Path(save_path+'results4_'+str(i))
                p5 = pathlib.Path(save_path+'results5_'+str(i))
                p6 = pathlib.Path(save_path+'results6_'+str(i))
                p1.unlink(missing_ok=True)
                p2.unlink(missing_ok=True)
                p3.unlink(missing_ok=True)
                p4.unlink(missing_ok=True)
                p5.unlink(missing_ok=True)
                p6.unlink(missing_ok=True)
    else:
        with open(save_path+'raw_data_1', 'wb') as f:
            if analyze_ms2:
                dill.dump([df1, df2, fragments_dataframes, version], f)
                del df1
                del df2
                del fragments_dataframes
            else:
                dill.dump([df1, df2, version], f)
                del df1
                del df2
            f.close()
        with open(save_path+'raw_data_2', 'wb') as f:
            dill.dump(eic_dataframes, f)
            del eic_dataframes
            f.close()
        with open(save_path+'raw_data_3', 'wb') as f:
            dill.dump(smoothed_eic_dataframes, f)
            del smoothed_eic_dataframes
            f.close()
        with open(save_path+'raw_data_4', 'wb') as f:
            dill.dump(curve_fitting_dataframes, f)
            del curve_fitting_dataframes
            f.close()
        with open(save_path+'raw_data_5', 'wb') as f:
            dill.dump(raw_eic_dataframes, f)
            del raw_eic_dataframes
            f.close()
        with open(save_path+'raw_data_6', 'wb') as f:
            dill.dump(isotopic_fits_dataframes, f)
            del isotopic_fits_dataframes
            f.close()
    print("Done!")

def print_sep(): ##Complete
    '''Prints a separator consisting of multiple '-' character.
    '''
    print('------------------------------------------------')
    
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
                  fast_iso,
                  verbose = False): ##Complete
    '''Integrates all the file-accessing associated functions in this script to go
    through the files data, draw and process eic of hypothetical glycans, does 
    peak-picking and calculates AUC of the peaks.

    Parameters
    ----------
    library : dict
        A glycans library, as generated by Library_Tools.full_glycans_library.
        
    lib_size : int
        The length of the library.
        
    data : list
        A list with each index containing a generator object of the sample file
        to be parsed.
        
    ms1_index : dict
        A dictionary containing the ms1 indexes of each sample file.
        
    tolerance : tuple
        First index contains the unit of the tolerance and the second one is the value of 
        that unit.
        
    ret_time_interval : tuple
        A tuple where the first index contains the beggining time of the retention time
        interval you wish to analyze and the second contains the end time.
        
    min_isotops : int
        The minimum amount of isotopologues required to consider an RT mz peak valid.
    
    min_ppp : tuple
        A tuple where the first index contains a boolean indicating whether or not to
        consider this parameter and the second one containing the minimum amount of
        data points per peak. If min_ppp[0] is set to False, calculates this automatically.
        
    max_charges : int
        The maximum amount of charges the queried mz should have.
        
    custom_noise : tuple
        A tuple containing two indexes: The first one indicates whether or not to use a custom
        noise set by the user and the second one is a list of noise levels, with each index
        indicating the noise level of each set file.
        
    close_peaks : tuple
        A tuple where the first index contains a boolean indicating whether or not to
        consider this parameter and the second one contains the amount of peaks to save.
        If close_peaks[0] is set to True, selects only the most intense peak and its 
        close_peaks[1] surrounding peaks.
        
    verbose : boolean
        Allows to print to file debugging information of each EIC traced so that you
        may find out why some specific pattern or behavior is ocurring.

    Uses
    ----
    General_Functions.noise_level_calc_mzarray : float
        Gathers the mz at the 95th percentile of the mz/int array (which is 
        equivalent to the 3rd SD from the mean.
    
    numpy.percentile : scalar or ndarray
        Compute the q-th percentile of the data along the specified axis.
    
    File_Accessing.eic_from_glycan : tuple
        Generates a very processed EIC for each adduct of each glycan of each sample.
        Removes non-monoisotopic peaks, check charges and calculates multiple quality 
        scoring data.
        
    File_Accessing.eic_smoothing : list
        Smoothes the EIC using the Savitsky Golay algorithm. Smoothed EIC may be
        used for peak-picking and curve-fitting scoring afterwards.

    File_Accessing.peaks_from_eic() : list
        Does multi peak-picking in a given smoothed EIC.

    File_Accessing.average_ppm_calc : tuple
        Calculates the arithmetic mean of the PPM differences of a given peak.
        
    File_Accessing.iso_fit_score_calc : float
        Calculates the mean isotopic fitting score of a given peak.
        
    File_Accessing.peak_curve_fit : tuple
        Calculates the fitting between the actual peak and an ideal peak based
        on a calculated gaussian bell curve.

    Returns
    -------
    analyzed_data : dict
        A dictionary similar to the one generated by the full glycans library generating
        function, with the added information of each adducts' peaks' retention time and
        AUC.
        
    rt_array_report : dict
        A dictionary with lists containing the retention times for each sample.
        
    noise : dict
        A dictionary containing the noise level for each sample.
    '''
    begin_time = datetime.datetime.now()
    analyzed_data = {}
    mean_sel_peaks = 0.0
    noise = {}
    noise_avg = {}
    rt_array_report = {}
    if not custom_noise[0]:
        print('Analyzing noise level of samples...', end='', flush = True)
    for i_i, i in enumerate(data):
        rt_array_report[i_i] = []
        temp_noise = []
        temp_avg_noise = []
        if custom_noise[0]:
            temp_noise.append(custom_noise[1][i_i])
            continue
        for j_j, j in enumerate(ms1_index[i_i]):
            rt_array_report[i_i].append(i[j]['retentionTime'])
            mz_ints = [i[j]['m/z array'], i[j]['intensity array']]
            if len(mz_ints[1]) > 0:
                noise_avg_level = percentile(mz_ints[1], 95)
            else:
                noise_avg_level = 0.0
            if i[j]['retentionTime'] < ret_time_interval[0] or i[j]['retentionTime'] > ret_time_interval[1]:
                noise_level = (0.0, 0.0, 0.0)
            else:
                noise_level = General_Functions.rt_noise_level_parameters_set(mz_ints)
            temp_noise.append(noise_level)
            temp_avg_noise.append(noise_avg_level)
        noise[i_i] = temp_noise
        noise_avg[i_i] = percentile(temp_avg_noise, 66.8)
    print('Done!')
    print_sep()
    if verbose:
        print('Noise Level: '+str(noise_avg))
        noise_report = 'Noise Level: '+str(noise_avg)
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
                                                  i,
                                                  glycan_data,
                                                  ms1_index,
                                                  ret_time_interval,
                                                  tolerance,
                                                  min_isotops,
                                                  noise,
                                                  noise_avg,
                                                  max_charges,
                                                  fast_iso,
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
                if i == "Internal Standard":
                    temp_eic_smoothed = File_Accessing.eic_smoothing(temp_eic[4][j][k])
                else:
                    temp_eic_smoothed = File_Accessing.eic_smoothing(temp_eic[0][j][k])
                glycan_data['Adducts_mz_data'][j][k] = []
                glycan_data['Adducts_mz_data'][j][k].append(temp_eic[0][j][k][1])
                glycan_data['Adducts_mz_data'][j][k].append([])
                glycan_data['Adducts_mz_data'][j][k].append(temp_eic_smoothed[1])
                glycan_data['Adducts_mz_data'][j][k].append(temp_eic[4][j][k][1])
                glycan_data['Adducts_mz_data'][j][k].append(temp_eic[5][j][k]) #all the isotopic fits, have to link them using arrange raw data
                if verbose:
                    verbose_info.append('Smoothed EIC INT Array: '+str(temp_eic_smoothed[1]))
                    verbose_info.append('Raw PPM error: '+str(temp_eic[1][j][k])+'\nRaw Isotopic Fitting Score: '+str(temp_eic[2][j][k]))
                if max(temp_eic[0][j][k][1]) < noise_avg[k] and i != "Internal Standard":
                    continue
                if i == "Internal Standard":
                    temp_peaks = File_Accessing.peaks_from_eic(temp_eic[4][j][k],
                                                               temp_eic_smoothed,
                                                               ret_time_interval,
                                                               min_ppp,
                                                               close_peaks,
                                                               i)
                else:
                    temp_peaks = File_Accessing.peaks_from_eic(temp_eic[0][j][k],
                                                               temp_eic_smoothed,
                                                               ret_time_interval,
                                                               min_ppp,
                                                               close_peaks,
                                                               i)
                if len(temp_peaks) == 0:
                    continue
                if i == "Internal Standard":
                    temp_peaks_auc = File_Accessing.peaks_auc_from_eic(temp_eic[4][j][k],
                                                                       ms1_index[k],
                                                                       temp_peaks)
                else:
                    temp_peaks_auc = File_Accessing.peaks_auc_from_eic(temp_eic[0][j][k],
                                                                       ms1_index[k],
                                                                       temp_peaks)
                for l_l, l in enumerate(temp_peaks):
                    if temp_peaks_auc[l_l] >= noise_avg[k]:
                        l['AUC'] = temp_peaks_auc[l_l]
                        l['Average_PPM'] = File_Accessing.average_ppm_calc(temp_eic[1][j][k], (tolerance[0], tolerance[1], glycan_data['Adducts_mz'][j]), l)
                        l['Iso_Fit_Score'] = File_Accessing.iso_fit_score_calc(temp_eic[2][j][k], l)
                        l['Signal-to-Noise'] = l['int']/(General_Functions.local_noise_calc(noise[k][l['id']], glycan_data['Adducts_mz'][j], noise_avg[k]))
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
    return analyzed_data, rt_array_report, noise_avg
    
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
                nglycan,
                permethylated,
                reduced,
                filter_output,
                unrestricted_fragments,
                rt_tolerance):
    '''Analyzes the MS2 data in the sample files, outputting the found matches.
    
    Parameters
    ----------
    ms2_index : dict
        A dictionary containing the ms2 indexes of each sample file.
        
    data : list
        A list with each index containing a generator object of the sample file
        to be parsed.
        
    analyzed_data : tuple
        Data outputted by the analyze_data function.
        
    ret_time_interval : tuple
        A tuple where the first index contains the beggining time of the retention time
        interval you wish to analyze and the second contains the end time.
        
    tolerance : tuple
        First index contains the unit of the tolerance and the second one is the value of 
        that unit.
        
    min_max_mono : tuple
        Minimum and maximum amount of monosaccharides for the hypotethical glycans in the
        library. ie. (5, 20).

    min_max_hex : tuple
        Minimum and maximum amount of hexoses for the hypotethical glycans in the library.
        ie. (5, 20).

    min_max_hexnac : tuple
        Minimum and maximum amount of N-Acetyl hexosamines for the hypotethical glycans
        in the library. ie. (5, 20).

    min_max_sia : tuple
        Minimum and maximum amount of sialic acids for the hypotethical glycans in the
        library. ie. (5, 20).

    min_max_fuc : tuple
        Minimum and maximum amount of deoxyhexoses for the hypotethical glycans in the
        library. ie. (5, 20).

    min_max_ac : tuple
        Minimum and maximum amount of N-Acetyl Neuraminic acids for the hypotethical
        glycans in the library. ie. (5, 20).

    min_max_gc : tuple
        Minimum and maximum amount of N-Glicolyl Neuraminic acids for the hypotethical
        glycans in the library. ie. (5, 20).
        
    max_charges : int
        The maximum amount of charges to calculate per glycan.
        
    tag_mass : float
        The tag's added mass to the glycans, if the glycans are tagged.
        Default = 0 (No Tag).
        
    filter_output : boolean
        Whether or not to force the output to fit glycans compositions.
        
    unrestricted_fragments : boolean
        Whether or not should take any fragment found
        
    Uses
    ----
    Library_Tools.fragments_library : list
        Generates a list of combinatorial analysis of monosaccharides from the minimum
        amount of monosaccharides to the maximum amount of monosaccharides set, then uses 
        the generated library and increments it with calculations with a series of information,
        clumps the duplicated mzs together as one query and then the script can use this
        library as a fragments library for MS2 analysis.
        
    General_Functions.form_to_charge : int
        Converts adducts formula into raw charge.
        
    Returns
    -------
    analyzed_data : tuple
        Data outputted by the analyze_data function.
        
    fragments_data : dict
        Dictionary containing the fragments data.
    '''
    begin_time = datetime.datetime.now()
    no_ms2 = True
    for i in ms2_index:
        if len(ms2_index[i]) != 0:
            no_ms2 = False
            break
    if no_ms2:
        print('No MS2 data to analyze...')
        dummy_fragment_data = {}
        for i in analyzed_data[0]:
            dummy_fragment_data[i] = {}
            for j in analyzed_data[0][i]['Adducts_mz_data']:
                dummy_fragment_data[i][j] = {}
                for k_k, k in enumerate(data):
                    dummy_fragment_data[i][j][k_k] = []
        return analyzed_data[0], analyzed_data[1], analyzed_data[2], dummy_fragment_data
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
                                  tag_mass,
                                  permethylated,
                                  reduced,
                                  nglycan)
    fragments_data = {}
    print('Scanning MS2 spectra...')
    for i_i, i in enumerate(analyzed_data[0]): #goes through each glycan found in analysis
        fragments_data[i] = {}
        for j_j, j in enumerate(analyzed_data[0][i]['Adducts_mz_data']): #goes through each adduct
            fragments_data[i][j] = {}
            for k_k, k in enumerate(data): #goes through each file
                fragments_data[i][j][k_k] = []
                if len(ms2_index[k_k]) == 0:
                    continue
                if len(analyzed_data[0][i]['Adducts_mz_data'][j][k_k][1]) == 0 and not unrestricted_fragments: #checks if found the adduct
                    continue
                for l in ms2_index[k_k]:
                    if unrestricted_fragments:
                        if k[l]['retentionTime'] < rt_interval[0] or k[l]['retentionTime'] > rt_interval[1]:
                            continue
                    else:
                        if k[l]['retentionTime'] < analyzed_data[0][i]['Adducts_mz_data'][j][k_k][1][0]['peak_interval'][0] - rt_tolerance or k[l]['retentionTime'] > analyzed_data[0][i]['Adducts_mz_data'][j][k_k][1][-1]['peak_interval'][1] + rt_tolerance: #skips spectra outside peak interval of peaks found
                            continue
                    if abs((k[l]['precursorMz'][0]['precursorMz']) - analyzed_data[0][i]['Adducts_mz'][j]) <= (1.0074/General_Functions.form_to_charge(j))+General_Functions.tolerance_calc(tolerance[0], tolerance[1], analyzed_data[0][i]['Adducts_mz'][j]): #checks if precursor matches adduct mz
                        found_count = 0
                        total = sum(k[l]['intensity array'])
                        for m_m, m in enumerate(k[l]['m/z array']):
                            found = False
                            for n_n, n in enumerate(fragments):
                                if 'Monos_Composition' in list(n.keys()):
                                    if (n['Monos_Composition']['H'] == analyzed_data[0][i]['Monos_Composition']['H']
                                        and n['Monos_Composition']['N'] == analyzed_data[0][i]['Monos_Composition']['N']
                                        and n['Monos_Composition']['S'] == analyzed_data[0][i]['Monos_Composition']['S']
                                        and n['Monos_Composition']['F'] == analyzed_data[0][i]['Monos_Composition']['F']
                                        and n['Monos_Composition']['G'] == analyzed_data[0][i]['Monos_Composition']['G']):
                                        continue
                                if found:
                                    break
                                combo = False
                                if filter_output:
                                    if "/" in n['Formula']:
                                        combo = True
                                        fragments_comp = []
                                        formula_splitted = n['Formula'].split("/")
                                        for o in formula_splitted:
                                            for p_p, p in enumerate(o):
                                                if p == "-" or p == "+" or p == "_":
                                                    fragments_comp.append(o[:p_p])
                                                    break
                                        for o_o, o in enumerate(fragments_comp):
                                            fragments_comp[o_o] = General_Functions.form_to_comp(o)
                                        viable = []
                                        for o in fragments_comp:
                                            if 'H' not in o.keys():
                                                o['H'] = 0
                                            if 'N' not in o.keys():
                                                o['N'] = 0
                                            if 'S' not in o.keys():
                                                o['S'] = 0
                                            if 'F' not in o.keys():
                                                o['F'] = 0
                                            if 'G' not in o.keys():
                                                o['G'] = 0
                                            if (o['H'] > analyzed_data[0][i]['Monos_Composition']['H']
                                                or o['N'] > analyzed_data[0][i]['Monos_Composition']['N']
                                                or o['S'] > analyzed_data[0][i]['Monos_Composition']['S']
                                                or o['F'] > analyzed_data[0][i]['Monos_Composition']['F']
                                                or o['G'] > analyzed_data[0][i]['Monos_Composition']['G']):
                                                viable.append(False)
                                                break
                                            else:
                                                viable.append(True)
                                        if True not in viable:
                                            continue
                                        else:
                                            count = 0
                                            new_formula = ""
                                            for o_o, o in enumerate(viable):
                                                if count == 0:
                                                    if o:
                                                        new_formula = formula_splitted[o_o]
                                                        count+= 1
                                                else:
                                                    if o:
                                                        new_formula+= "/"+formula_splitted[o_o]
                                                        count+= 1
                                    elif "/" not in n['Formula'] and not combo:
                                        if (n['Monos_Composition']['H'] > analyzed_data[0][i]['Monos_Composition']['H'] or n['Monos_Composition']['N'] > analyzed_data[0][i]['Monos_Composition']['N'] or n['Monos_Composition']['S'] > analyzed_data[0][i]['Monos_Composition']['S'] or n['Monos_Composition']['F'] > analyzed_data[0][i]['Monos_Composition']['F'] or n['Monos_Composition']['G'] > analyzed_data[0][i]['Monos_Composition']['G']):
                                            continue
                                for o in n['Adducts_mz']:
                                    if abs(n['Adducts_mz'][o]-m) <= General_Functions.tolerance_calc(tolerance[0], tolerance[1], n['Adducts_mz'][o]): #fragments data outputted in the form of (Glycan, Adduct, Fragment, Fragment mz, intensity, retention time, precursor)
                                        if "_" not in n['Formula']:
                                            fragments_data[i][j][k_k].append((i, j, n['Formula']+'_'+o, n['Adducts_mz'][o], k[l]['intensity array'][m_m], k[l]['retentionTime'], k[l]['precursorMz'][0]['precursorMz'], total))
                                        elif "_" in n['Formula'] and combo:
                                            fragments_data[i][j][k_k].append((i, j, new_formula, n['Adducts_mz'][o], k[l]['intensity array'][m_m], k[l]['retentionTime'], k[l]['precursorMz'][0]['precursorMz'], total))
                                        else:
                                            fragments_data[i][j][k_k].append((i, j, n['Formula'], n['Adducts_mz'][o], k[l]['intensity array'][m_m], k[l]['retentionTime'], k[l]['precursorMz'][0]['precursorMz'], total))
                                        found = True
                                        found_count += k[l]['intensity array'][m_m]
                                        break
    print('Sample MS2 analysis done in '+str(datetime.datetime.now() - begin_time)+'!')
    return analyzed_data[0], analyzed_data[1], analyzed_data[2], fragments_data

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

from . import General_Functions
from .Execution_Functions import print_sep
import pkg_resources
import pathlib
import datetime
import platform
import os

version1 = "0.0.0"
version2 = "0.0.0"
try:
    version1 = pkg_resources.get_distribution("glycogenius").version
except:
    pass
try:
    version_path = str(pathlib.Path(__file__).parent.parent.parent.resolve())
    with open(version_path+"/Setup.py", "r") as f: #grabs version from setup.py to add to raw_data files
        for lines in f:
            if lines[:12] == "    version=":
                version2 = lines[13:-2].strip("'")
except:
    pass
if version2 > version1:
    version = version2
else:
    version = version1

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

def interactive_terminal():
    '''This function generates the CLI for user interaction.
    
    Uses
    ----
    pathlib.Path.resolve() : Path object
        Make the path absolute, resolving any symlinks. A new path object is returned
        
    General_Functions.form_to_comp(string) : dict
        Separates a molecular formula or monosaccharides composition of glycans into a
        dictionary with each atom/monosaccharide as a key and its amount as value.
        
    Execution_Functions.print_sep : string
        Prints a fixed size separator.
        
    datetime.datetime.now : time
        Gets current time.
        
    platform.system : string
        Gets current operational system.
        
    os._exit
        Exits the current execution.
        
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
    try:
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
            print("1 - Build and output glycans library.\n2 - Analyze sample files\n3 - Reanalyze raw results files with new\n    parameters\n4 - Create template parameters file for command-\n    line execution\n5 - Exit")
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
                license_path = str(pathlib.Path(__file__).parent.parent.resolve())
                for i_i, i in enumerate(license_path):
                    if i == "\\":
                        license_path = license_path[:i_i]+"/"+license_path[i_i+1:]
                with open(license_path+"/LICENSE.py", 'r') as f:
                    for line in f:
                        print(line, end = "")
                continue
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
                    var = var.strip()
                    if var == "":
                        print(glycans_list)
                        var2 = input("Proceed with this list? (y/n): ")
                        if var2 == 'y' or var2 == "":
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
                n_glycan = True
                while True:
                    var = input("Force compositions to N-glycans structure\n (default: yes) (y/n): ")
                    if var == '':
                        n_glycan = True
                        break
                    if var == 'y':
                        n_glycan = True
                        break
                    if var == 'n':
                        n_glycan = False
                        break
                    else:
                        continue
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
                    var2 = input("Insert the tag added mass\nor molecular formula (ie. 133.0644 or C7H7N3\nfor GirP or 219.1735 or C13H21N3 for ProA) or\ninput 'pep-'+ aminoacids sequence for peptide\nas tag (ie. pep-NK for the peptide NK): ")
                    try:
                        var2 = float(var2)
                    except:
                        if var2.split('-')[0] == 'pep':
                            var2 = dict(mass.Composition(sequence = var2.split('-')[-1]))
                            var2['H'] -= 2
                            var2['O'] -= 1
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
                    var = input("Insert the accuracy value for the unit you've\nchosen (ie. '0.01' for 'mz' or '10' for 'ppm'): ")
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
                        files = var
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
                    return input_order, glycans_list, adducts, max_charges, tag_mass, fast_iso, high_res, ms2, accuracy_unit, accuracy_value, rt_int, min_isotop, max_ppm, iso_fit, curve_fit, sn, files, path, permethylated, reduced, n_glycan
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
    except KeyboardInterrupt:
        print("\n\n----------Execution cancelled by user.----------\n")
        raise SystemExit(1)
        
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

from pyteomics import mzxml, mzml, mass, auxiliary
from itertools import combinations_with_replacement
from numpy import percentile
from re import split
from math import inf, exp, pi
from statistics import stdev, mean
import sys
import datetime

##---------------------------------------------------------------------------------------
##Hard-coded permanent information

monosaccharides = {
    "H": ("Hexose", "C6O6H12", {"C": 6, "O": 5, "N": 0, "H": 10}),
    "N": ("N-Acetyl Hexosamine", "C8O6NH15", {"C": 8, "O": 5, "N": 1, "H": 13}),
    "S": ("Acetyl Neuraminic Acid", "C11O9NH19", {"C": 11, "O": 8, "N": 1, "H": 17}),
    "F": ("Fucose", "C6O5H12", {"C": 6, "O": 4, "N": 0, "H": 10}),
    "G": ("Glycolyl Neuraminic Acid", "C11O10NH19", {"C": 11, "O": 9, "N": 1, "H": 17})
    }
'''A hardcoded dictionary containing each single letter code for monosaccharides as key
and a tuple containing the full monosaccharide name, its full molecular formula and its
residue composition in dict form as value.
'''

h_mass = mass.calculate_mass(composition={'H' : 1})
'''The mass of an hydrogen-1 atom. Pre-calculated here to avoid calculating too many
times during a run.
'''

##---------------------------------------------------------------------------------------
##General functions (these functions use only external libraries, such as itertools and
##pyteomics).

def calculate_ppm_diff(mz, target):
    '''Calculates the PPM difference between a mz and a target mz.
    
    Parameters
    ----------
    mz : float
        A float of the mz you want to check for the difference.
        
    target : float
        A float of the target mz you're comparing your mz to.
        
    Returns
    -------
    float
        A float of the PPM difference between the mz and target mz.
    '''
    return ((target-mz)/target)*(10**6)
    
def tolerance_calc(unit,
                   value, 
                   mz = 1000):
    '''A fast, but not super accurate way to convert 'ppm' mass accuracy into 'pw' (peak
    width, aka. mz tolerance).

    Parameters
    ----------
    unit : string
        Can be either "ppm" (particles-per-million) or "pw" (peak width [tolerance]).

    value : float
        Float value of the tolerance, based on the unit inputted.
        
    mz : float
        mz at which to calculate the PPM difference.

    Returns
    -------
    tolerance : float
        If unit == "ppm", converts value into "pw", if unit == "pw", outputs value as is.
        ie. 10 ppm gets converted to 0.01 pw tolerance.
    '''
    if unit == "ppm":
        return (mz-(-((value*mz)/10**6)+mz))
    elif unit == "mz":
        return value
    else:
        return("Unit for tolerance not 'ppm' or 'mz'.")
    
def rt_noise_level_parameters_set(mz_int):
    '''Gathers the int at the 95th percentile of the mz/int array (which is 
    equivalent to the 3rd SD from the mean.
    
    Parameters
    ----------
    mz_int : list
        A list containing two synchronized lists: the first one contains the
        mzs and the second one the intensity.
        
    Uses
    ----
    numpy.percentile : scalar or ndarray
        Compute the q-th percentile of the data along the specified axis.
        
    Returns
    -------
    scalar
        The intensity of the 95th percentile of the intensity array.
    '''
    if len(mz_int[0]) == 0:
        return 0.0, 0.0, 0.0
    int_list_first_quarter = []
    int_list_last_quarter = []
    first_quarter_end = mz_int[0][-1]/4
    last_quarter_begin = (mz_int[0][-1]/4)*3
    for i_i, i in enumerate(mz_int[0]):
        if i <= first_quarter_end:
            int_list_first_quarter.append(mz_int[1][i_i])
        if i > last_quarter_begin:
            int_list_last_quarter.append(mz_int[1][i_i])
    if len(int_list_first_quarter) == 0 or len(int_list_last_quarter) == 0:
        return 0.0, 0.0, 0.0
    return percentile(int_list_first_quarter, 95), percentile(int_list_last_quarter, 95), mz_int[0][-1]
    
def local_noise_calc(noise_specs, x, avg_noise):
    '''Uses the noise_specs produced by rt_noise_level_parameters_set to 
    calculate the local noise levels. If any of the parameters are considered
    abnormal (ie. absurdly high or no peaks on the array) it defaults to avg_noise.
    
    Parameters
    ----------
    noise_specs : tuple
        A tuple containing scalars, as produced by rt_noise_level_parameters_set
        
    x : float
        The mz at which you want to calculate the noise.
        
    avg_noise : float
        The average noise of the file, used for fallback.
        
    Returns
    -------
    float
        A float containing the local or average noise level.
    '''
    if noise_specs[2] == 0.0 or noise_specs[1] > noise_specs[0]*5 or noise_specs[0] > noise_specs[1]*5:
        return avg_noise
    return noise_specs[0] + (((noise_specs[1]-noise_specs[0])/noise_specs[2])*x)
    
def normpdf(x, mean, sd):
    '''Calculates the intensity of a gaussian bell curve at the x-axis point x
    with the set parameters.
    
    Parameters
    ----------
    x : float
        Given point in the x axis of the gaussian bell curve.
        
    mean : float
        Determines the maximum intensity of the mean point of the gaussian
        bell curve.
        
    sd : float
        Determines the standard deviation of the gaussian.
        
    Uses
    ----
    math.pi : float
        The mathematical constant pi = 3.141592…, to available precision.
        
    math.exp : float
        Return e raised to the power x, where e = 2.718281… is the base of 
        natural logarithms.
        
    Returns
    -------
    float
        A float containing the intensity of the gaussian bell curve at the 
        given x point.
    '''
    var = float(sd)**2
    denom = (2*pi*var)**.5
    num = exp(-(float(x)-float(mean))**2/(2*var))
    return num/denom

def form_to_comp(string):
    '''Separates a molecular formula or monosaccharides formula of glycans into a
    dictionary with each atom/monosaccharide as a key and its amount as value.

    Parameters
    ----------
    string : str
        A string in the form of C6O6N0H12 or H5N4S1F1G1.

    Returns
    -------
    counts : dict
        A dictionary with keys containing the monosaccharides/atoms letter(s) and values
        containing the amounts of each. ie. {"H": 5, "N": 4, "S": 1, "F": 1, "G": 1}.
    '''
    counts = {}
    split_str = split('(\\d+)', string)
    negative = False
    for i_i, i in enumerate(split_str):
        if i != '' and i[-1] == '-':
            split_str[i_i] = i[:-1]
            negative = True
        if i_i%2 != 0 and i != '' and negative:
            split_str[i_i] = '-'+i
    if len(split_str)%2 != 0:
        split_str.append('1')
    for i in range(len(split_str)-1):
        if i%2 == 0:
            counts[split_str[i]] = int(split_str[i+1])
    if '' in counts:
        del counts['']
    return counts

def form_to_charge(string):
    '''Converts adducts formula into raw charge.

    Parameters
    ----------
    string : str
        A string containing adducts formula.

    Uses
    ----
    form_to_comp() : dict
        Separates a molecular formula or monosaccharides formula of glycans into a
        dictionary with each atom/monosaccharide as a key and its amount as value

    Returns
    -------
    charge : int
        The raw charge of the adduct.
    '''
    comp = form_to_comp(string)
    charge = 0
    for i in comp:
        if i != '':
            charge+=comp[i]
    return charge

def glycan_to_atoms(glycan_composition, permethylated):
    '''Calculates the amounts of atoms based on glycan monosaccharides.

    Parameters
    ----------
    glycan_composition : dict
        Accepts as input the glycan monosaccharides formula in the format of {"H": 5,
        "N": 4, "S": 1, "F": 1, "G": 1}.
        
    permethylated : boolean
        Whether the glycans are permethylated or not.

    Uses
    ----
    monosaccharides : dict
        A hardcoded dictionary containing each single letter code for monosaccharides as
        key and a tuple containing the full monosaccharide name, its full molecular
        formula and its residue composition in dict form.

    Returns
    -------
    atoms : dict
        Returns a dictionary with atoms as keys and amounts as values. ie. {"C": 6,
        "O": 6, "N": 0, "H": 12}.
    '''
    atoms = {"C": 0, "O": 0, "N": 0, "H": 0}
    if permethylated:
        for i in monosaccharides:
            if i == 'H':
                monosaccharides[i][2]['C'] += 3
                monosaccharides[i][2]['H'] += 6
            if i == 'N':
                monosaccharides[i][2]['C'] += 3
                monosaccharides[i][2]['H'] += 6
            if i == 'F':
                monosaccharides[i][2]['C'] += 2
                monosaccharides[i][2]['H'] += 4
            if i == 'S' or i == 'G':
                monosaccharides[i][2]['C'] += 5
                monosaccharides[i][2]['H'] += 10
    for i in glycan_composition:
        if i == "T":
            continue
        for j in atoms:
            atoms[j] += monosaccharides[i][2][j]*glycan_composition[i]
    return atoms

def count_seq_letters(string):
    '''If you make anything with itertools for combinatorial analysis, it will produce a
    string that's not very human readable. This converts it into a human readable form.

    Parameters
    ----------
    string : str
        A string of atoms or glycans in the form o CCCCCCOOOONH or HHHHHNNNNFSG.

    Returns
    -------
    friendly_letters : dict
        A dictionary containing the count for each letter in the string. ie. CCCCOONH
        returns {"C": 4, "O": 2, "N": 1, "H": 1}.
    '''
    friendly_letters = {}
    last_letter = ""
    for i in string:
        if i != last_letter:
            friendly_letters[i] = string.count(i)
            last_letter = i
    return friendly_letters

def sum_atoms(*compositions):
    '''Sums the atoms of two compositions.

    Parameters
    ----------
    compositions : dict
        Dictionaries containing the atomic compositions of the items to be summed.
        ie. {"C": 4, "O": 2, "N": 1, "H": 1}.

    Returns
    -------
    summed_comp : dict
        Dictionary containing the sum of each atom of the compositions.
    '''
    summed_comp = {"C": 0, "O": 0, "N": 0, "H": 0}
    for i in compositions:
        for j in i:
            summed_comp[j]+=i[j]
    return summed_comp

def sum_monos(*compositions):
    '''Sums the monosaccharides of two glycan compositions. 'T' stands for TAG, used
    in fragments library calculation.

    Parameters
    ----------
    compositions : dict
        Dictionaries containing the monosaccharrides compositions of the glycans to be
        summed. ie. {"H": 5, "N": 4, "S": 1, "F": 1, "G": 1}.

    Returns
    -------
    summed_comp : dict
        Dictionary containing the sum of each monosaccharides of the compositions.
    '''
    summed_comp = {"H": 0, "N": 0, "S": 0, "F": 0, "G": 0, "T": 0}
    for i in compositions:
        for j in i:
            summed_comp[j]+=i[j]
    return summed_comp

def comp_to_formula(composition):
    '''Transforms a composition dictionary into string formula.

    Parameters
    ----------
    composition : dict
        Dictionary containing the composition of the molecule or glycan.

    Returns
    -------
    formula : string
        Formula of the atomic or monosaccharides composition in string form.
    '''
    formula = ''
    for i in composition:
        if composition[i] != 0:
            formula+=i+str(composition[i])
    return formula

def calculate_comp_from_mass(tag_mass):
    '''Calculates the composition of a molecule based on its mass. Intended to use with
    small tags added to the glycans.

    Parameters
    ----------
    tag_mass : float
        The monoisotopic molecular weight of the tag's added mass.

    Uses
    ----
    pyteomics.mass.calculate_mass(*args, **kwargs) : float
        Calculates the monoisotopic mass of a polypeptide defined by a sequence string,
        parsed sequence, chemical formula or Composition object.

    itertools.combinations_with_replacement : generator
        Return r length subsequences of elements from the input iterable allowing
        individual elements to be repeated more than once.

    Returns
    -------
    closest : tuple
        Returns the proposed composition of the molecule with the calculated mass closest
        to the tag's added mass in dictionary form and the calculated hypothetical mass
        of the molecule as float.
    '''
    closest = ({}, 0)
    test_tag_mass = 0
    atoms_number = int(tag_mass/5)
    while atoms_number > atoms_number/2:
        for i in combinations_with_replacement("CONH", atoms_number):
            seq_readable = count_seq_letters("".join(i))
            test_tag_mass = mass.calculate_mass(composition = seq_readable)
            if abs(test_tag_mass-tag_mass) < abs(closest[1]-tag_mass):
                closest = (seq_readable, test_tag_mass)
        atoms_number -= 1      
    return closest
    
def calculate_isotopic_pattern(glycan_atoms,
                               fast=True,
                               high_res=False):
    '''Calculates up to 5 isotopic pattern peaks relative abundance in relation with the
    monoisotopic one.

    Parameters
    ----------
    glycan_atoms : dict
        A dictionary containing glycan atomic composition in the form of {"C": 4, "O": 2,
        "N": 1, "H": 1}.

    tolerance : float
        Mass tolerance for isotopologues clumping. ie. 0.1 means an isotopologue with
        mass 310.00 and another one with a mass of 310.09 will be clumped as a single
        peak.
    
    fast : boolean
        If True, only calculates the isotopic pattern based on isotopes of carbon and
        nitrogen, thus very inaccurate, but enough for most uses. If fast = False, it
        will take a significantly higher amount of time to produce results (from 10x to
        1000x longer, depending on number of atoms).
        Default = True.

    Uses
    ----
    pyteomics.mass.isotopologues(*args, **kwargs) : iterator
        Iterate over possible isotopic states of a molecule. The molecule can be defined
        by formula, sequence, parsed sequence, or composition.
        
    pyteomics.mass.calculate_mass(*args, **kwargs) : float
        Calculates the monoisotopic mass of a polypeptide defined by a sequence string,
        parsed sequence, chemical formula or Composition object.
    
    Returns
    -------
    relative_isotop_pattern : list
        A list with the isotopic pattern, with the first element being the abundance of
        the monoisotopic peak (1.0, or 100%) and the following ones are the isotopologues
        relative abundance in relation to the monoisotopic peak (around 1 Da apart).
    '''
    if fast:
        isotopologue = mass.isotopologues(glycan_atoms, report_abundance = True,
                                          elements_with_isotopes = ["C"],
                                          overall_threshold = 1e-4)
    else:
        isotopologue = mass.isotopologues(glycan_atoms, report_abundance = True,
                                          elements_with_isotopes = ["C", "N", "O", "H"],
                                          overall_threshold = 1e-4)
    isotop_arranged = []
    relative_isotop_pattern = []
    relative_isotop_mass = []
    for i in isotopologue:
        isotop_arranged.append({'mz' : mass.calculate_mass(i[0]), 'Ab' : i[1]})
    isotop_arranged = sorted(isotop_arranged, key=lambda x: x['mz'])
    for i_i, i in enumerate(isotop_arranged):
        relative_isotop_pattern.append(i['Ab']/isotop_arranged[0]['Ab'])
        relative_isotop_mass.append(i['mz'])
    if not high_res and not fast:
        relative_isotop_pattern_low_res = []
        relative_isotop_mass_low_res = []
        for i_i, i in enumerate(relative_isotop_mass):
            if i_i == 0:
                relative_isotop_pattern_low_res.append(relative_isotop_pattern[i_i])
                relative_isotop_mass_low_res.append(i)
            else:
                if abs(i-relative_isotop_mass[i_i-1]) < h_mass/2:
                    relative_isotop_mass_low_res[-1] = (relative_isotop_mass_low_res[-1]+i)/2
                    relative_isotop_pattern_low_res[-1]+= relative_isotop_pattern[i_i]
                else:
                    relative_isotop_mass_low_res.append(i)
                    relative_isotop_pattern_low_res.append(relative_isotop_pattern[i_i])
        return relative_isotop_pattern_low_res, relative_isotop_mass_low_res
    return relative_isotop_pattern, relative_isotop_mass

def gen_adducts_combo(adducts,
                      max_charge):
    '''Generates a list of dictionaries with compositions of adducts combinations, based
    on parameters set.

    Parameters
    ----------
    adducts : dict
        A dictionary with each key containing a single atom adduct and its value
        containing the maximum amount of such adduct.

    max_charge : int
        The maximum amount of charges for the adducts.

    Uses
    ----
    itertools.combinations_with_replacement : generator
        Return r length subsequences of elements from the input iterable allowing
        individual elements to be repeated more than once.

    Returns
    -------
    adducts_combo_dict : list
        A list of dictionaries containing the composition of each adducts combination.
    '''
    adducts_list = []
    adducts_combo = []
    adducts_combo_dict = []
    for i in adducts:
        adducts_list.append(i)
    for i in range(1, abs(max_charge)+1):
        for j in combinations_with_replacement(adducts_list, i):
            adducts_combo.append(j)
    for i in adducts_combo:
        temp_dict = {}
        for j in i:
            if j not in temp_dict:
                if max_charge > 0:
                    temp_dict[j] = 1
                else:
                    temp_dict[j] = -1
            else:
                if max_charge > 0:
                    temp_dict[j]+= 1
                else:
                    temp_dict[j]-= 1
        adducts_combo_dict.append(temp_dict)
    to_remove = []
    for i in adducts_combo_dict:
        for j in i:
            if abs(i[j]) > abs(adducts[j]):
                to_remove.append(i)
                break
    for i in to_remove:
        adducts_combo_dict.remove(i)
    return adducts_combo_dict

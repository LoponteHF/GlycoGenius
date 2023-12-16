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
    '''
    '''
    return ((target-mz)/target)*(10**6)
    
def noise_level_calc_mzarray(mz_int):
    '''Verifies what is noise by adjusting a slope and the percentile between 1.5 
    standard deviations to 2 standard deviations.
    '''
    int_list = []
    for i in mz_int:
        int_list.append(mz_int[i])
    return percentile(int_list, 95)
    
def normpdf(x, mean, sd):
    '''
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

def glycan_to_atoms(glycan_composition):
    '''Calculates the amounts of atoms based on glycan monosaccharides.

    Parameters
    ----------
    glycan_composition : dict
        Accepts as input the glycan monosaccharides formula in the format of {"H": 5,
        "N": 4, "S": 1, "F": 1, "G": 1}.

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
    '''Sums the monosaccharides of two glycan compositions.

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
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
from pyteomics import mzxml, mzml, mass, auxiliary
from itertools import combinations_with_replacement
from collections import Counter
from re import split
from math import inf
import concurrent.futures
import pathlib
import importlib
import copy
import sys
import datetime
import os

##---------------------------------------------------------------------------------------
##Library generating-associated functions (these functions make vast use of the general
## functions).

def generate_combinations_with_constraints(characters, length, constraints):
    '''
    '''
    def is_valid(counts):
        for char, (min_count, max_count) in constraints.items():
            count = counts[char]
            if count < min_count or count > max_count:
                return False
        return True

    # Generate all possible count combinations within the constraints
    char_indices = list(range(len(characters)))
    valid_count_combinations = []

    for combination in combinations_with_replacement(char_indices, length):
        counts = Counter(combination)
        counts = {characters[idx]: counts[idx] for idx in counts}
        
        # Fill in missing characters with zero counts
        for char in characters:
            if char not in counts:
                counts[char] = 0
        
        if is_valid(counts):
            valid_count_combinations.append(counts)

    # Convert count combinations to strings
    valid_combinations = []
    for counts in valid_count_combinations:
        if 'L' in counts:
            counts['Am'] = counts.pop('L')
        if 'A' in counts:
            counts['AmG'] = counts.pop('A')
        if 'R' in counts:
            counts['EG'] = counts.pop('R')
        if 'M' in counts:
            counts['HN'] = counts.pop('M')
        if 'U' in counts:
            counts['UA'] = counts.pop('U')
        valid_combinations.append(counts)

    return valid_combinations
    
def generate_glycans_library(min_max_mono,
                             min_max_hex,
                             min_max_hexnac,
                             min_max_xyl,
                             min_max_sialics,
                             min_max_fuc,
                             min_max_ac,
                             min_max_gc,
                             min_max_hn,
                             min_max_ua,
                             lactonized_ethyl_esterified,
                             forced):
    '''Generates a list of combinatorial analysis of monosaccharides from the minimum
    amount of monosaccharides to the maximum amount of monosaccharides set, then trims it
    based on the specific monosaccharides criteria inserted.

    Parameters
    ----------
    min_max_mono : tuple
        Minimum and maximum amount of monosaccharides for the hypotethical glycans in the
        library. ie. (5, 20).

    min_max_hex : tuple
        Minimum and maximum amount of hexoses for the hypotethical glycans in the library.
        ie. (5, 20).

    min_max_hexnac : tuple
        Minimum and maximum amount of N-Acetyl hexosamines for the hypotethical glycans
        in the library. ie. (5, 20).

    min_max_xyl : tuple
        Minimum and maximum amount of Xyloses for the hypotethical glycans
        in the library. ie. (5, 20).

    min_max_sialics : tuple
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
        
    min_max_hn : tuple
        Minimum and maximum amount of Hexosamines for the hypotethical glycans in the library. ie. (5, 20).
        
    min_max_ua : tuple
        Minimum and maximum amount of Uronic Acids for the hypotethical glycans in the library. ie. (5, 20).
        
    lactonized_ethyl_esterified : boolean
        Whether the glycans were submitted to lactonization/ethyl-esterification
        derivatization, which differentiates the mass of alpha2-3 and alpha2-6 -bound 
        sialic acids.

    forced : string
        Indicates whether the function should force strict conditions based on the
        biological knowledge of glycans in order to avoid possible false positives when
        analysing N-glycans, O-glycans or GAGs.

    Uses
    ----
    General_Functions.count_seq_letters(string) : string
        If you make anything with itertools for combinatorial analysis, it will produce a
        string that's not very human readable. This converts it into a human readable
        form.

    itertools.combinations_with_replacement : generator
        Return r length subsequences of elements from the input iterable allowing
        individual elements to be repeated more than once.

    General_Functions.sum_monos(*compositions) : dict
        Sums the monosaccharides of two glycan compositions.

    Returns
    -------
    glycans : list
        A list containing dictionaries of the monosaccharides compositions of all the
        glycans generated.
    '''
    glycans = []
    def_glycan_comp = {"H": 0, "N": 0, "X": 0, "S": 0, "Am": 0, "E": 0, "F": 0, "G": 0, "AmG": 0, "EG": 0, "HN": 0, "UA": 0}
    
    if lactonized_ethyl_esterified:
        monos_chars = "HNXLEFARMU" #H = Hexose, N = HexNAc, X = Xylose, L = Amidated Neu5Ac, E = Ethyl-esterified Neu5Ac, F = DeoxyHexose, A = Amidated Neu5Gc, R = Ethyl-esterified Neu5Gc, S = Neu5Ac, G = Neu5Gc, M = Hexosamine, U = Uronic Acids
    else:
        monos_chars = "HNXSFGMU"
        
    constraints = {'H': (min_max_hex[0], min_max_hex[1]),
                   'N': (min_max_hexnac[0], min_max_hexnac[1]),
                   'X': (min_max_xyl[0], min_max_xyl[1]),
                   'F': (min_max_fuc[0], min_max_fuc[1]),
                   'M': (min_max_hn[0], min_max_hn[1]),
                   'U': (min_max_ua[0], min_max_ua[1])}
                   
    if lactonized_ethyl_esterified:
        constraints['L'] = (min_max_ac[0], min_max_ac[1])
        constraints['E'] = (min_max_ac[0], min_max_ac[1])
        constraints['A'] = (min_max_gc[0], min_max_gc[1])
        constraints['R'] = (min_max_gc[0], min_max_gc[1])
    else:
        constraints['S'] = (min_max_ac[0], min_max_ac[1])
        constraints['G'] = (min_max_gc[0], min_max_gc[1])
      
    results = []
    cpu_number = (os.cpu_count())-2 if os.cpu_count() < 60 else 60
    if cpu_number <= 0:
        cpu_number = 1
    with concurrent.futures.ProcessPoolExecutor(max_workers = cpu_number) as executor:  
        for i in range(min_max_mono[0], min_max_mono[1]+1):
            result = executor.submit(generate_combinations_with_constraints,
                                     monos_chars,
                                     i,
                                     constraints)
            results.append(result)
        
        for index, i in enumerate(results):
            result_data = i.result()
            for j in result_data:
                glycans.append(General_Functions.sum_monos(def_glycan_comp, j))
            results[index] = None
            
    to_be_removed = []
    for i_i, i in enumerate(glycans):
        if lactonized_ethyl_esterified:
            if ((i['Am']+i['E']+i['AmG']+i['EG'] < min_max_sialics[0])
                or (i['Am']+i['E']+i['AmG']+i['EG'] > min_max_sialics[1])
                or (i['Am']+i['E'] < min_max_ac[0]) 
                or (i['Am']+i['E'] > min_max_ac[1])
                or (i['AmG']+i['EG'] < min_max_gc[0])
                or (i['AmG']+i['EG'] > min_max_gc[1])):
                to_be_removed.append(i)
        else:
            if ((i['S']+i['G'] < min_max_sialics[0])
                or (i['S']+i['G'] > min_max_sialics[1])):
                to_be_removed.append(i)
    if forced == 'n_glycans':
        if lactonized_ethyl_esterified:
            for i_i, i in enumerate(glycans):
                if ((i['Am']+i['E']+i['AmG']+i['EG'] > (i['N']-2)*2)
                    or (i['F'] >= i['N']) 
                    or (i['Am']+i['E']+i['AmG']+i['EG'] > i['H']-2)
                    or (i['H'] < 2) 
                    or (i['N'] < 2)
                    or (i['X'] > 1)
                    or (i['HN'] > 0)
                    or (i['UA'] > 0)):
                    if i not in to_be_removed:
                        to_be_removed.append(i)
        else:
            for i_i, i in enumerate(glycans):
                if ((i['S']+i['G'] > (i['N']-2)*2)
                    or (i['F'] >= i['N']) 
                    or (i['S']+i['G'] > i['H']-2)
                    or (i['H'] < 2) 
                    or (i['N'] < 2)
                    or (i['X'] > 1)
                    or (i['HN'] > 0)
                    or (i['UA'] > 0)):
                    if i not in to_be_removed:
                        to_be_removed.append(i)
                        
    if forced == 'o_glycans':
        if lactonized_ethyl_esterified:
            for i_i, i in enumerate(glycans):
                if ((i['Am']+i['E']+i['AmG']+i['EG'] > i['N']+i['H'])
                    or (i['H'] > i['N']+1)
                    or (i['N'] > i['H']+3)
                    or (i['X'] > 0)
                    or (i['HN'] > 0)
                    or (i['UA'] > 0)):
                    if i not in to_be_removed:
                        to_be_removed.append(i)
        else:
            for i_i, i in enumerate(glycans):
                if ((i['S']+i['G'] > i['N']+i['H'])
                    or (i['F'] > i['N']+i['H'])
                    or (i['H'] > i['N']+1)
                    or (i['N'] > i['H']+3)
                    or (i['X'] > 0)
                    or (i['HN'] > 0)
                    or (i['UA'] > 0)):
                    if i not in to_be_removed:
                        to_be_removed.append(i)
                        
    if forced == 'gags':
        to_add = []
        for i_i, i in enumerate(glycans):
            monos_count = sum(i.values())
            if ((i['N']+i['HN'] > i['UA']+i['H']+1)
                or (i['H']+i['UA'] > i['N']+i['HN'])
                or (i['H'] > 0 and i['UA'] > 0)
                or (i['S'] + i['G'] + i['AmG']+i['EG']+i['Am']+i['E'] > 0)
                or (i['F'] > 0)):
                if i not in to_be_removed:
                    to_be_removed.append(i)
            # else: # This adds a variant of the gag with its linkage domain; can be useful if analyzing digested proteoglycans peptides containing gags
                # new_glycan = copy.deepcopy(i)
                # new_glycan['X'] += 1
                # new_glycan['H'] += 2
                # new_glycan['UA'] += 1
                # to_add.append(new_glycan)
        # for i_i, i in enumerate(to_add):
            # glycans.append(i)
    
    for i_i, i in enumerate(glycans):
        monos_count = sum(i.values())
        if monos_count == 0:
            if i not in to_be_removed:
                to_be_removed.append(i)
                
    for i in to_be_removed:
        glycans.remove(i)
    return glycans

def full_glycans_library(library,
                         forced,
                         max_adducts,
                         adducts_exclusion,
                         max_charges,
                         tag_mass = 0,
                         fast = True,
                         high_res = False,
                         internal_standard = 0.0,
                         permethylated = False,
                         reduced = False,
                         min_max_sulfation = [0, 0],
                         min_max_phosphorylation = [0, 0]):
    '''Uses the generated glycans library and increments it with calculations of its
    mass, neutral mass with tag mass, its isotopic distribution pattern (intensoids) and
    the mz of the chosen adducts combinations.

    Parameters
    ----------
    library : list
        A list of dictionaries, each dict containing the composition of a hypothetical
        glycan, as generated by generate_glycans_library.

    max_adducts : dict
        A dictionary with keys containing each possible atomic adducts (ie. 'H', 'Na',
        'K', etc.) and the maximum amount of such adducts as the values.
        
    adducts_exclusion : list
        A list containing undesired adducts combinations to be excluded from the function
        result.

    max_charges : int
        The maximum amount of charges to calculate per glycan.
        
    tolerance : tuple
        First index contains the unit of the tolerance and the second one is the value of 
        that unit.

    tag_mass : float, dict or str
        The tag's added mass, molecular formula or peptide sequence, if the glycans 
        are tagged or bound to a peptide.
        Default = 0 (No Tag).

    fast : boolean
        Makes the isotopic distribution calculation fast (less accurate) or not (more
        accurate).
        Default = True
        
    high_res : boolean
        Decides whether to clump (if set to False) or not (if set to True) the neighbouring
        isotopic envelope peaks. Only works if fast is set to False.
        Default = False
        
    internal_standard : float
        If a internal standard is added to the sample, this allows the function
        to calculate its mass based on adducts combination, as well as adding the tag to
        it.
        
    permethylated : boolean
        Whether or not the sample was permethylated.
        
    reduced : boolean
        Whether or not the sample was reduced.
        
    min_max_sulfation : list
        Minimum and maximum amount of sulfations per glycan.
        
    min_max_phosphorylation : list
        Minimum and maximum amount of phosphorylations per glycan.

    Uses
    ----
    General_Functions.calculate_comp_from_mass(tag_mass) : dict
        Calculates the composition of a molecule based on its mass. Intended to use with
        small tags added to the glycans.
        
    General_Functions.gen_adducts_combo(adducts, max_charge) : list
        Generates a list of dictionaries with compositions of adducts combinations,
        based on parameters set.

    General_Functions.comp_to_formula(composition) : string
        Transforms a composition dictionary into string formula.

    General_Functions.sum_atoms(*compositions) : dict
        Sums the atoms of two compositions.

    General_Functions.glycan_to_atoms(glycan_composition) : dict
        Calculates the amounts of atoms based on glycan monosaccharides.

    General_Functions.form_to_comp(string) : dict
        Separates a molecular formula or monosaccharides composition of glycans into a
        dictionary with each atom/monosaccharide as a key and its amount as value.

    pyteomics.mass.calculate_mass(*args, **kwargs) : float
        Calculates the monoisotopic mass of a polypeptide defined by a sequence string,
        parsed sequence, chemical formula or Composition object.

    General_Functions.calculate_isotopic_pattern(glycan_atoms,
                               tolerance,
                               fast=True) : list
        Calculates up to 5 isotopic pattern peaks relative abundance in relation with
        the monoisotopic one.

    Returns
    -------
    full_library : dict
        A dictionary with each key containing the glycan formula and each key containing
        a dictionary with monosaccharides composition, atoms composition with tag,
        neutral mass, neutral mass with tag, isotopic distribution and the mzs of the
        glycans with the desired adducts combination.
    '''
    full_library = {}
    
    try:
        tag_mass = float(tag_mass)
    except:
        if tag_mass.split('-')[0] == 'pep':
            tag_mass = dict(mass.Composition(sequence = tag_mass.split('-')[-1]))
            tag_mass['H'] -= 2
            tag_mass['O'] -= 1
            
    if tag_mass != 0:
        if type(tag_mass) == float:
            tag = General_Functions.calculate_comp_from_mass(tag_mass)
        elif type(tag_mass) != dict:
            comp_tag = General_Functions.form_to_comp(tag_mass)
            tag = (comp_tag, mass.calculate_mass(comp_tag))
        else:
            comp_tag = tag_mass
            tag = (comp_tag, mass.calculate_mass(comp_tag))
        tag_mass = tag[1]
    else:
        tag = ({"C": 0, "O": 0, "N": 0, "H": 0}, 0.0)
        
    adducts_combo = General_Functions.gen_adducts_combo(max_adducts, adducts_exclusion, max_charges)
    
    results = []
    cpu_number = (os.cpu_count())-2 if os.cpu_count() < 60 else 60
    if cpu_number <= 0:
        cpu_number = 1
    with concurrent.futures.ProcessPoolExecutor(max_workers = cpu_number) as executor:
        for i in library:
            result = executor.submit(calculate_one_glycan,
                                     i,
                                     adducts_combo,
                                     tag,
                                     forced,
                                     fast,
                                     high_res,
                                     permethylated,
                                     reduced,
                                     min_max_sulfation,
                                     min_max_phosphorylation)
            results.append(result)
        
        for index, i in enumerate(results):
            result_data = i.result()
            full_library[result_data[1]] = result_data[0]
            results[index] = None
    
    full_library = dict(sorted(full_library.items()))
            
    if internal_standard != '0.0' and internal_standard != 0.0:
        full_library = include_internal_standard(full_library, tag_mass, fast, high_res, internal_standard, permethylated, reduced, adducts_combo)
    return full_library
    
def calculate_one_glycan(i,
                         adducts_combo,
                         tag,
                         forced,
                         fast,
                         high_res,
                         permethylated,
                         reduced,
                         min_max_sulfation,
                         min_max_phosphorylation):
    '''
    '''
    monos_count = sum(i.values())
    number = 1
    if monos_count <= 4 and forced == 'gags':
        number = 2
    for w in range(number):
        for s in range(min_max_sulfation[0], min_max_sulfation[1]+1):
            if s > monos_count*3:
                break
            for p in range(min_max_phosphorylation[0], min_max_phosphorylation[1]+1):
                if p > monos_count*2:
                    break
                sulfations = f"{s}(s)"
                phosphorylations = f"{p}(p)"
                i_formula = General_Functions.comp_to_formula(i)
                if s > 0:
                    i_formula = f"{i_formula}+{sulfations}"
                if p > 0:
                    i_formula = f"{i_formula}+{phosphorylations}"
                if w > 0:
                    i_formula = f"{i_formula}-H2O"
                i_atoms = General_Functions.glycan_to_atoms(i, permethylated)
                if w == 0:
                    i_atoms = General_Functions.sum_atoms(i_atoms, General_Functions.form_to_comp('H2O'))
                if tag[1] == 0.0:
                    if permethylated:
                        i_atoms = General_Functions.sum_atoms(i_atoms, {'C': 2, 'H': 4})
                        if reduced:
                            i_atoms = General_Functions.sum_atoms(i_atoms, {'O': 1})
                    if not permethylated and reduced:
                        i_atoms = General_Functions.sum_atoms(i_atoms, {'H': 2})
                base_mass = mass.calculate_mass(composition=i_atoms)
                i_atoms = General_Functions.sum_atoms(i_atoms, {'S': 1*s, 'O': 3*s}) #sum sulfation atoms
                i_atoms = General_Functions.sum_atoms(i_atoms, {'P': 1*p, 'O': 3*p, 'H': 1*p}) #sum phosphorylation atoms
                i_atoms_tag = General_Functions.sum_atoms(i_atoms, tag[0])
                i_neutral_mass = mass.calculate_mass(composition=i_atoms)
                i_neutral_tag = i_neutral_mass+tag[1]
                i_iso_dist = General_Functions.calculate_isotopic_pattern(i_atoms_tag, fast, high_res)
                iso_corrected = i_iso_dist[0]
                if fast:
                    iso_corrected = []
                    for j_j, j in enumerate(i_iso_dist[0]):
                        if j_j == 1:
                            iso_corrected.append(abs(j*1.02)) #default correction
                            iso_corrected[-1] = iso_corrected[-1]+((((80/base_mass)/2)*iso_corrected[-1])*s)
                            continue
                        if j_j == 2:
                            iso_corrected.append(abs(j*(10.8*(i_neutral_mass**-0.267)))) #default correction
                            iso_corrected[-1] = iso_corrected[-1]+(((((80/base_mass)*6)-0.0733)*iso_corrected[-1])*s)
                            iso_corrected[-1] = iso_corrected[-1]+((((80/base_mass)/2)*iso_corrected[-1])*p)
                            continue
                        if j_j == 3:
                            iso_corrected.append(abs(j*(122.62*(i_neutral_mass**-0.528)))) #default correction
                            iso_corrected[-1] = iso_corrected[-1]+(((((80/base_mass)*6)-0.0733)*iso_corrected[-1])*s)
                            iso_corrected[-1] = iso_corrected[-1]+((((80/base_mass)/2)*iso_corrected[-1])*p)
                            continue
                        if j_j == 4:
                            iso_corrected.append(abs(j*(2192.6*(i_neutral_mass**-0.833)))) #default correction
                            iso_corrected[-1] = iso_corrected[-1]+(((((80/base_mass)*6)-0.0733)*iso_corrected[-1])*s)
                            iso_corrected[-1] = iso_corrected[-1]+((((80/base_mass)/2)*iso_corrected[-1])*p)
                            continue
                        else:
                            iso_corrected.append(j)
                            continue
                glycan_info = {}
                glycan_info['Monos_Composition'] = i
                glycan_info['Atoms_Glycan+Tag'] = i_atoms_tag
                glycan_info['Neutral_Mass'] = i_neutral_mass
                glycan_info['Neutral_Mass+Tag'] = i_neutral_tag
                glycan_info['Isotopic_Distribution'] = iso_corrected
                glycan_info['Isotopic_Distribution_Masses'] = i_iso_dist[1]
                glycan_info['Adducts_mz'] = {}
                for j in adducts_combo:
                    charges = sum(j.values())
                    mz = mass.calculate_mass(i_atoms_tag,
                                             charge = charges,
                                             charge_carrier = j,
                                             carrier_charge = charges)
                    glycan_info['Adducts_mz'][General_Functions.comp_to_formula(j)] = mz
                    
    return glycan_info, i_formula

def include_internal_standard(full_library,
                              tag_mass,
                              fast,
                              high_res,
                              internal_standard,
                              permethylated,
                              reduced,
                              adducts_combo):
    '''Processeses the internal standard inputted into workable library dictionary entry, with isotopic distribution and proper applicable modifications.

    Parameters
    ----------
    full_library : dictionary
        A dictionary containing all the glycans generated for analysis.

    tag_mass : float, dict or str
        The tag's added mass, molecular formula or peptide sequence, if the glycans 
        are tagged or bound to a peptide.
        Default = 0 (No Tag).

    fast : boolean
        Makes the isotopic distribution calculation fast (less accurate) or not (more
        accurate).
        Default = True
        
    high_res : boolean
        Decides whether to clump (if set to False) or not (if set to True) the neighbouring
        isotopic envelope peaks. Only works if fast is set to False.
        Default = False
        
    internal_standard : float
        If a internal standard is added to the sample, this allows the function
        to calculate its mass based on adducts combination, as well as adding the tag to
        it.
        
    permethylated : boolean
        Whether or not the sample was permethylated.
        
    reduced : boolean
        Whether or not the sample was reduced.
        
    adducts_combo : list
        A list containing the combinations fo adducts generated.

    Uses
    ----
    General_Functions.calculate_comp_from_mass(tag_mass) : dict
        Calculates the composition of a molecule based on its mass. Intended to use with
        small tags added to the glycans.

    General_Functions.comp_to_formula(composition) : string
        Transforms a composition dictionary into string formula.

    General_Functions.sum_atoms(*compositions) : dict
        Sums the atoms of two compositions.

    General_Functions.glycan_to_atoms(glycan_composition) : dict
        Calculates the amounts of atoms based on glycan monosaccharides.

    General_Functions.form_to_comp(string) : dict
        Separates a molecular formula or monosaccharides composition of glycans into a
        dictionary with each atom/monosaccharide as a key and its amount as value.

    pyteomics.mass.calculate_mass(*args, **kwargs) : float
        Calculates the monoisotopic mass of a polypeptide defined by a sequence string,
        parsed sequence, chemical formula or Composition object.

    General_Functions.calculate_isotopic_pattern(glycan_atoms,
                               tolerance,
                               fast=True) : list
        Calculates up to 5 isotopic pattern peaks relative abundance in relation with
        the monoisotopic one.

    Returns
    -------
    full_library : dict
        A dictionary with each key containing the glycan formula and each key containing
        a dictionary with monosaccharides composition, atoms composition with tag,
        neutral mass, neutral mass with tag, isotopic distribution and the mzs of the
        glycans with the desired adducts combination.
    '''
    def_glycan_comp = {"H": 0, "N": 0, "X": 0, "S": 0, "Am": 0, "E": 0, "F": 0, "G": 0, "AmG": 0, "EG": 0, "HN": 0, "UA": 0}
    
    i_formula = 'Internal Standard'
    
    if tag_mass != 0:
        if type(tag_mass) == float:
            tag = General_Functions.calculate_comp_from_mass(tag_mass)
        elif type(tag_mass) != dict:
            comp_tag = General_Functions.form_to_comp(tag_mass)
            tag = (comp_tag, mass.calculate_mass(comp_tag))
        else:
            comp_tag = tag_mass
            tag = (comp_tag, mass.calculate_mass(comp_tag))
        tag_mass = tag[1]
    else:
        tag = ({"C": 0, "O": 0, "N": 0, "H": 0}, 0.0)
        
    try:
        i_neutral_mass = float(internal_standard)
        i_atoms = General_Functions.calculate_comp_from_mass(i_neutral_mass)
        i_type = 'direct_mass'
    except:
        i_composition = General_Functions.form_to_comp(internal_standard)
        if "C" in i_composition.keys():
            i_atoms = i_composition
            i_neutral_mass = mass.calculate_mass(composition=i_atoms)
            i_type = 'atomic_comp'
        else:
            i_composition = General_Functions.sum_monos(def_glycan_comp, i_composition)
            i_atoms = General_Functions.glycan_to_atoms(i_composition, permethylated)
            i_atoms = General_Functions.sum_atoms(i_atoms, General_Functions.form_to_comp('H2O'))
            i_neutral_mass = mass.calculate_mass(composition=i_atoms)
            i_type = 'glycan'
            
    if i_type == 'glycan':
        if tag[1] == 0.0:
            if permethylated:
                i_atoms = General_Functions.sum_atoms(i_atoms, {'C': 2, 'H': 4})
                if reduced:
                    i_atoms = General_Functions.sum_atoms(i_atoms, {'O': 1})
            if not permethylated and reduced:
                i_atoms = General_Functions.sum_atoms(i_atoms, {'H': 2})
        i_atoms_tag = General_Functions.sum_atoms(i_atoms, tag[0])
        i_neutral_tag = i_neutral_mass+tag[1]
    else:
        i_atoms_tag = i_atoms
        i_neutral_tag = i_neutral_mass
            
    i_iso_dist = General_Functions.calculate_isotopic_pattern(i_atoms_tag, fast, high_res)
    iso_corrected = i_iso_dist[0]
    if fast:
        iso_corrected = []
        for j_j, j in enumerate(i_iso_dist[0]):
            if j_j == 1:
                iso_corrected.append(abs(j*1.02)) #default correction
                continue
            if j_j == 2:
                iso_corrected.append(abs(j*(10.8*(i_neutral_mass**-0.267)))) #default correction
                continue
            if j_j == 3:
                iso_corrected.append(abs(j*(122.62*(i_neutral_mass**-0.528)))) #default correction
                continue
            if j_j == 4:
                iso_corrected.append(abs(j*(2192.6*(i_neutral_mass**-0.833)))) #default correction
                continue
            else:
                iso_corrected.append(j)
                continue
                
    full_library[i_formula] = {}
    if i_type == 'glycan':
        full_library[i_formula]['Monos_Composition'] = i_composition
    else:
        full_library[i_formula]['Monos_Composition'] = {"H": 0, "N": 0, "X": 0, "S": 0, "Am": 0, "E": 0, "F": 0, "G": 0, "AmG": 0, "EG": 0, "HN": 0, "UA": 0}
    full_library[i_formula]['Neutral_Mass'] = i_neutral_mass
    full_library[i_formula]['Neutral_Mass+Tag'] = i_neutral_tag
    full_library[i_formula]['Isotopic_Distribution'] = iso_corrected
    full_library[i_formula]['Isotopic_Distribution_Masses'] = i_iso_dist[1]
    
    full_library[i_formula]['Adducts_mz'] = {}
    for j in adducts_combo:
        charges = sum(j.values())
        mz = mass.calculate_mass(i_atoms_tag,
                                 charge = charges,
                                 charge_carrier = j,
                                 carrier_charge = charges)
        full_library[i_formula]['Adducts_mz'][General_Functions.comp_to_formula(j)] = mz
        
    return full_library

def fragments_library(min_max_mono,
                      min_max_hex,
                      min_max_hexnac,
                      min_max_xyl,
                      min_max_sialics,
                      min_max_fuc,
                      min_max_ac,
                      min_max_gc,
                      min_max_hn,
                      min_max_ua,
                      max_charges,
                      tolerance,
                      tag_mass,
                      permethylated,
                      reduced,
                      lactonized_ethyl_esterified,
                      nglycan):
    '''Generates a list of combinatorial analysis of monosaccharides from the minimum
    amount of monosaccharides to the maximum amount of monosaccharides set, then uses 
    the generated library and increments it with calculations with a series of information,
    clumps the duplicated mzs together as one query and then the script can use this
    library as a fragments library for MS2 analysis.

    Parameters
    ----------
    min_max_mono : tuple
        Minimum and maximum amount of monosaccharides for the hypotethical glycans in the
        library. ie. (5, 20).

    min_max_hex : tuple
        Minimum and maximum amount of hexoses for the hypotethical glycans in the library.
        ie. (5, 20).

    min_max_hexnac : tuple
        Minimum and maximum amount of N-Acetyl hexosamines for the hypotethical glycans
        in the library. ie. (5, 20).

    min_max_xyl : tuple
        Minimum and maximum amount of Xyloses for the hypotethical glycans
        in the library. ie. (5, 20).

    min_max_sialics : tuple
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
        
    min_max_hn : tuple
        Minimum and maximum amount of Hexosamine for the hypotethical
        glycans in the library. ie. (5, 20).
        
    min_max_ua : tuple
        Minimum and maximum amount of Uronic Acid for the hypotethical
        glycans in the library. ie. (5, 20).

    max_charges : int
        The maximum amount of charges to calculate per glycan.

    tolerance : float
        The mz acceptable tolerance. Used to calculate the intensoids, at this point.

    tag_mass : float
        The tag's added mass to the glycans, if the glycans are tagged.
        Default = 0 (No Tag).
        
    permethylated : boolean
        Whether the glycan is permethylated or not.
        
    reduced : boolean
        Whether the glycan is reduced or not.
        
    lactonized_ethyl_esterified : boolean
        Whether the glycans were submitted to lactonization/ethyl-esterification
        derivatization, which differentiates the mass of alpha2-3 and alpha2-6 -bound 
        sialic acids.
        
    Uses
    ----
    General_Functions.count_seq_letters(string) : string
        If you make anything with itertools for combinatorial analysis, it will produce a
        string that's not very human readable. This converts it into a human readable
        form.

    itertools.combinations_with_replacement : generator
        Return r length subsequences of elements from the input iterable allowing
        individual elements to be repeated more than once.

    General_Functions.sum_monos(*compositions) : dict
        Sums the monosaccharides of two glycan compositions.
        
    General_Functions.calculate_comp_from_mass(tag_mass) : dict
        Calculates the composition of a molecule based on its mass. Intended to use with
        small tags added to the glycans.
        
    General_Functions.gen_adducts_combo(adducts, exclusions, max_charge) : list
        Generates a list of dictionaries with compositions of adducts combinations,
        based on parameters set.

    General_Functions.comp_to_formula(composition) : string
        Transforms a composition dictionary into string formula.

    General_Functions.sum_atoms(*compositions) : dict
        Sums the atoms of two compositions.

    General_Functions.glycan_to_atoms(glycan_composition) : dict
        Calculates the amounts of atoms based on glycan monosaccharides.

    pyteomics.mass.calculate_mass(*args, **kwargs) : float
        Calculates the monoisotopic mass of a polypeptide defined by a sequence string,
        parsed sequence, chemical formula or Composition object.

    Returns
    -------
    glycans : list
        A list containing dictionaries with informations about each fragment generated.
    '''
    time_formatted = str(datetime.datetime.now()).split(" ")[-1].split(".")[0]+" - "
    print(time_formatted+"Building fragments library...", end = "", flush = True)
    glycans = []
    def_glycan_comp = {"H": 0, "N": 0, "X": 0, "S": 0, "Am": 0, "E": 0, "F": 0, "G": 0, "AmG": 0, "EG": 0, "T" : 0, "HN": 0, "UA": 0}
    
    if lactonized_ethyl_esterified:
        monos_chars = "HNXLEFARTMU"
    else:
        monos_chars = "HNXSFGTMU"
        
    constraints = {'H': (0, min_max_hex[1]),
                   'N': (0, min_max_hexnac[1]),
                   'X': (0, min_max_xyl[1]),
                   'F': (0, min_max_fuc[1]),
                   'T': (0, 1),
                   'M': (0, min_max_hn[1]),
                   'U': (0, min_max_ua[1])}
                   
    if lactonized_ethyl_esterified:
        constraints['L'] = (0, min_max_ac[1])
        constraints['E'] = (0, min_max_ac[1])
        constraints['A'] = (0, min_max_gc[1])
        constraints['R'] = (0, min_max_gc[1])
    else:
        constraints['S'] = (0, min_max_ac[1])
        constraints['G'] = (0, min_max_gc[1])
        
    for i in range(1, min_max_mono[1]+2):
        combinations = generate_combinations_with_constraints(monos_chars, i, constraints)
        for j in combinations:
            glycans.append(General_Functions.sum_monos(def_glycan_comp, j))
            
    to_be_removed = []
    for i_i, i in enumerate(glycans):
        if lactonized_ethyl_esterified:
            if ((i['T'] == 1 and i['N'] == 0)
                or (i['Am']+i['E']+i['AmG']+i['EG'] > min_max_sialics[1])
                or (i['Am']+i['E'] > min_max_ac[1])
                or (i['AmG']+i['EG'] > min_max_gc[1])):
                to_be_removed.append(i_i)
        else:
            if ((i['T'] == 1 and i['N'] == 0)
                or (i['S']+i['G'] > min_max_sialics[1])):
                to_be_removed.append(i_i)
        if nglycan and i_i not in to_be_removed: #some rules and hardcoded exceptions for N-Glycans
            if lactonized_ethyl_esterified:
                if ((i['T'] == 1 and sum(i.values()) < 8 and i['Am']+i['E']+i['AmG']+i['EG'] > 0)
                    or (sum(i.values()) < 6 and i['Am']+i['E'] >= 1 and i['N'] > 1)
                    or (i['H'] > 0 and i['T'] == 1 and i['N'] < 2)
                    or (i['H'] == 2 and i['N'] == 1 and i['Am']+i['E'] == 0 and i['F'] == 0 and i['AmG']+i['EG'] == 0 and i['T'] == 1)
                    or (i['H'] == 0 and i['N'] == 1 and (i['Am']+i['E'] == 1 or i['AmG']+i['EG'] == 1) and i['F'] == 0 and i['T'] == 0)
                    or (i['H'] == 1 and i['N'] == 3 and i['Am']+i['E'] == 0 and i['F'] == 0 and i['AmG']+i['EG'] == 0 and i['T'] == 1)
                    or (i['H'] > 1 and i['N'] == 0 and (i['Am']+i['E'] == 1 or i['AmG']+i['EG'] == 1) and i['F'] == 0 and i['T'] == 0)
                    or (i['H'] == 1 and i['N'] == 1 and i['Am']+i['E'] == 0 and i['F'] == 0 and i['AmG']+i['EG'] == 0 and i['T'] == 1)
                    or (i['H'] == 3 and i['N'] == 1 and i['Am']+i['E'] == 0 and i['F'] == 0 and i['AmG']+i['EG'] == 0 and i['T'] == 1)):
                    to_be_removed.append(i_i)
            else:
                if ((i['T'] == 1 and sum(i.values()) < 8 and i['S']+i['G'] > 0)
                    or (sum(i.values()) < 6 and i['S'] >= 1 and i['N'] > 1)
                    or (i['H'] > 0 and i['T'] == 1 and i['N'] < 2)
                    or (i['H'] == 2 and i['N'] == 1 and i['S'] == 0 and i['F'] == 0 and i['G'] == 0 and i['T'] == 1)
                    or (i['H'] == 0 and i['N'] == 1 and (i['S'] == 1 or i['G'] == 1) and i['F'] == 0 and i['T'] == 0)
                    or (i['H'] == 1 and i['N'] == 3 and i['S'] == 0 and i['F'] == 0 and i['G'] == 0 and i['T'] == 1)
                    or (i['H'] > 1 and i['N'] == 0 and (i['S'] == 1 or i['G'] == 1) and i['F'] == 0 and i['T'] == 0)
                    or (i['H'] == 1 and i['N'] == 1 and i['S'] == 0 and i['F'] == 0 and i['G'] == 0 and i['T'] == 1)
                    or (i['H'] == 3 and i['N'] == 1 and i['S'] == 0 and i['F'] == 0 and i['G'] == 0 and i['T'] == 1)):
                    to_be_removed.append(i_i)
    for i in sorted(to_be_removed, reverse = True):
        del glycans[i]
    try:
        tag_mass = float(tag_mass)
    except:
        if tag_mass.split('-')[0] == 'pep':
            tag_mass = dict(mass.Composition(sequence = tag_mass.split('-')[-1]))
            tag_mass['H'] -= 2
            tag_mass['O'] -= 1
            
    if tag_mass != 0:
        if type(tag_mass) == float:
            tag = General_Functions.calculate_comp_from_mass(tag_mass)
        elif type(tag_mass) != dict:
            comp_tag = General_Functions.form_to_comp(tag_mass)
            tag = (comp_tag, mass.calculate_mass(comp_tag))
        else:
            comp_tag = tag_mass
            tag = (comp_tag, mass.calculate_mass(comp_tag))
        tag_mass = tag[1]
    else:
        tag = ({"C": 0, "O": 0, "N": 0, "H": 0}, 0.0)
    adducts_combo = General_Functions.gen_adducts_combo({'H' : 2}, [], max_charges)
    frag_library = []
    combo_frags_lib = []
    adducts_mz = [[], [], []]
    for i_i, i in enumerate(glycans):
        for j_j, j in enumerate(range(-1, 2)):
            if j < 0:
                i_formula = General_Functions.comp_to_formula(i)+str(j)+'H2O'
            elif j > 0:
                i_formula = General_Functions.comp_to_formula(i)+'+'+str(j)+'H2O'
            else:
                i_formula = General_Functions.comp_to_formula(i)
            glycan_atoms = General_Functions.glycan_to_atoms(i, permethylated)
            glycan_atoms['H'] += j*2
            glycan_atoms['O'] += j*1
            i_atoms = glycan_atoms
            if tag[1] == 0.0:
                if permethylated:
                    i_atoms = General_Functions.sum_atoms(i_atoms, {'C': 2, 'H': 4})
                    if reduced:
                        i_atoms = General_Functions.sum_atoms(i_atoms, {'O': 1})
                if not permethylated and reduced:
                    i_atoms = General_Functions.sum_atoms(i_atoms, {'H': 2})
            i_neutral_mass = mass.calculate_mass(composition=i_atoms)
            if i['T'] == 1:
                i_atoms_tag = General_Functions.sum_atoms(i_atoms, tag[0])
                i_neutral_tag = i_neutral_mass+tag[1]
                if tag[1] == 0.0:
                    if permethylated:
                        i_atoms = General_Functions.sum_atoms(i_atoms, {'C': 1, 'H': 2})
                        if reduced:
                            i_atoms = General_Functions.sum_atoms(i_atoms, {'O': 1})
                    if not permethylated and reduced:
                        i_atoms = General_Functions.sum_atoms(i_atoms, {'H': 2})
            else:
                i_atoms_tag = i_atoms
                i_neutral_tag = i_neutral_mass
            frag_library.append({})
            index = (i_i*3)+j_j
            frag_library[index]['Formula'] = i_formula
            frag_library[index]['Monos_Composition'] = i
            frag_library[index]['Adducts_mz'] = {}
            if tag[1] == 0.0: #this skips and marks for removal the reducing end without water
                if j < 0:
                    continue
            for j in adducts_combo:
                charges = sum(j.values())
                mz = mass.calculate_mass(i_atoms_tag,
                                         charge = charges,
                                         charge_carrier = j,
                                         carrier_charge = charges)
                adducts_mz[0].append(mz)
                adducts_mz[1].append(General_Functions.comp_to_formula(j))
                adducts_mz[2].append(index)
                frag_library[index]['Adducts_mz'][General_Functions.comp_to_formula(j)] = {'mz': mz, 'Ambiguities': []}
                
    for i_i, i in enumerate(frag_library):
        for j_j, j in enumerate(i['Adducts_mz']):
            for k_k, k in enumerate(frag_library):
                if k_k == i_i:
                    continue
                for l in k['Adducts_mz']:
                    if abs(k['Adducts_mz'][l]['mz'] - i['Adducts_mz'][j]['mz']) < General_Functions.tolerance_calc(tolerance[0], tolerance[1], k['Adducts_mz'][l]['mz']):
                        i['Adducts_mz'][j]['Ambiguities'].append((k_k, l))
                        
    print("Done!")
    return frag_library
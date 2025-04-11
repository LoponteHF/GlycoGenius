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

def generate_combinations_with_constraints(characters, length, constraints, monosaccharides):
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
        for multi_letter, single_letter in monosaccharides.items():
            single_letter = single_letter[-1]
            if single_letter in counts:
                counts[multi_letter] = counts.pop(single_letter)
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
                             forced,
                             custom_monos = []):
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
    monosaccharides = copy.deepcopy(General_Functions.monosaccharides)
    # Add custom monosaccharides
    if len(custom_monos) > 0:
        for cm in custom_monos:
            monosaccharides[cm['cm_short_code']] = (cm['cm_name'], cm['cm_chem_comp'], General_Functions.sum_atoms({"C": 0, "O": 0, "N": 0, "H": 0}, General_Functions.form_to_comp(cm['cm_chem_comp'])), cm['cm_single_letter_code'])
    
    glycans = []
    def_glycan_comp = {key : 0 for key in monosaccharides}
    sialics = ['S', 'G', 'Am', 'E', 'AmG', 'EG']
    
    # Add sialics from custom monosaccharides
    if len(custom_monos) > 0:
        for cm in custom_monos:
            if cm['sialic']:
                sialics.append(cm['cm_short_code'])
    
    if lactonized_ethyl_esterified:
        monos_chars = "".join([letter[-1] for letter in monosaccharides.values() if (letter[-1] != "S" and letter[-1] != "G")]) #H = Hexose, N = HexNAc, X = Xylose, L = Amidated Neu5Ac, E = Ethyl-esterified Neu5Ac, F = DeoxyHexose, A = Amidated Neu5Gc, R = Ethyl-esterified Neu5Gc, S = Neu5Ac, G = Neu5Gc, M = Hexosamine, U = Uronic Acids
    else:
        monos_chars = "".join([letter[-1] for letter in monosaccharides.values() if (letter[-1] != "R" and letter[-1] != "A" and letter[-1] != "E" and letter[-1] != "L")])
        
    constraints = {'H': (min_max_hex[0], min_max_hex[1]),
                   'N': (min_max_hexnac[0], min_max_hexnac[1]),
                   'X': (min_max_xyl[0], min_max_xyl[1]),
                   'F': (min_max_fuc[0], min_max_fuc[1]),
                   'M': (min_max_hn[0], min_max_hn[1]),
                   'U': (min_max_ua[0], min_max_ua[1])}
    
    # Add constraints from custom monosaccharides
    if len(custom_monos) > 0:
        for cm in custom_monos:
            constraints[cm['cm_single_letter_code']] = (cm['cm_min'], cm['cm_max'])
                   
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
                                     constraints,
                                     monosaccharides)
            results.append(result)
        
        for index, i in enumerate(results):
            result_data = i.result()
            for j in result_data:
                glycans.append(General_Functions.sum_monos(def_glycan_comp, j, monos = monosaccharides))
            results[index] = None
            
    to_be_removed = []
    for i_i, i in enumerate(glycans):
        if lactonized_ethyl_esterified:
            if ((sum([i[sialic] for sialic in sialics]) < min_max_sialics[0])
                or (sum([i[sialic] for sialic in sialics]) > min_max_sialics[1])
                or (i['Am']+i['E'] < min_max_ac[0]) 
                or (i['Am']+i['E'] > min_max_ac[1])
                or (i['AmG']+i['EG'] < min_max_gc[0])
                or (i['AmG']+i['EG'] > min_max_gc[1])):
                to_be_removed.append(i)
        else:
            if ((sum([i[sialic] for sialic in sialics]) < min_max_sialics[0])
                or (sum([i[sialic] for sialic in sialics]) > min_max_sialics[1])):
                to_be_removed.append(i)
    if forced == 'n_glycans':
        if lactonized_ethyl_esterified:
            for i_i, i in enumerate(glycans):
                if ((sum([i[sialic] for sialic in sialics]) > (i['N']-2)*2)
                    or (i['F'] >= i['N']) 
                    or (sum([i[sialic] for sialic in sialics]) > i['H']-2)
                    or (i['H'] < 2) 
                    or (i['N'] < 2)
                    or (i['X'] > 1)
                    or (i['HN'] > 0)
                    or (i['UA'] > 0)):
                    if i not in to_be_removed:
                        to_be_removed.append(i)
        else:
            for i_i, i in enumerate(glycans):
                if ((sum([i[sialic] for sialic in sialics]) > (i['N']-2)*2)
                    or (i['F'] >= i['N']) 
                    or (sum([i[sialic] for sialic in sialics]) > i['H']-2)
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
                if ((sum([i[sialic] for sialic in sialics]) > i['N']+i['H'])
                    or (i['H'] > i['N']+1)
                    or (i['N'] > i['H']+3)
                    or (i['X'] > 0)
                    or (i['HN'] > 0)
                    or (i['UA'] > 0)):
                    if i not in to_be_removed:
                        to_be_removed.append(i)
        else:
            for i_i, i in enumerate(glycans):
                if ((sum([i[sialic] for sialic in sialics]) > i['N']+i['H'])
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
                or (sum([i[sialic] for sialic in sialics]) > 0)
                or (i['F'] > 0)):
                if i not in to_be_removed:
                    to_be_removed.append(i)
            # else: # This adds a variant of the gag with its linkage domain; can be useful if analyzing digested peptides containing gags
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
                         min_max_phosphorylation = [0, 0],
                         lyase_digested = False,
                         custom_monos = []):
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
    monosaccharides = copy.deepcopy(General_Functions.monosaccharides)
    # Add custom monosaccharides
    if len(custom_monos) > 0:
        for cm in custom_monos:
            if cm['cm_short_code'] not in monosaccharides.keys():
                monosaccharides[cm['cm_short_code']] = (cm['cm_name'], cm['cm_chem_comp'], General_Functions.sum_atoms({"C": 0, "O": 0, "N": 0, "H": 0}, General_Functions.form_to_comp(cm['cm_chem_comp'])), cm['cm_single_letter_code'])
                
    full_library = {}
    
    # This codeblock sorts out the tag type and calculates it appropriately
    try:
        tag_mass = float(tag_mass)
    except:
        if tag_mass.split('-')[0] == 'pep':
            sequence = tag_mass.split('-')[-1]
            tag_mass = dict(mass.Composition(sequence = sequence))
            
            if 'C' in sequence: # Alkylation of cysteines by IAA
                cysteines = sequence.count('C')
                tag_mass['C'] += 2*cysteines
                tag_mass['H'] += 3*cysteines
                tag_mass['O'] += cysteines
                tag_mass['N'] += cysteines
                
            tag_mass['H'] -= 2
            tag_mass['O'] -= 1
            
    tag = General_Functions.determine_tag_comp(tag_mass)
        
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
                                     min_max_phosphorylation,
                                     monosaccharides,
                                     lyase_digested)
            results.append(result)
        
        for index, i in enumerate(results):
            result_data = i.result()
            for index_glycan, glycan in enumerate(result_data[1]):
                full_library[glycan] = result_data[0][index_glycan]
            results[index] = None
    
    full_library = dict(sorted(full_library.items()))
            
    if internal_standard != '0.0' and internal_standard != 0.0:
        full_library = include_internal_standard(full_library, tag_mass, fast, high_res, internal_standard, permethylated, reduced, adducts_combo, custom_monos)
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
                         min_max_phosphorylation,
                         monosaccharides,
                         lyase_digested = False):
    '''
    '''
    monos_count = sum(i.values())
    i_formulas = []
    glycan_infos = []
    for s in range(min_max_sulfation[0], min_max_sulfation[1]+1):
        if s > monos_count*3:
            break
        for sodiums in range(s+1):
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
                if lyase_digested:
                    i_formula = f"{i_formula}-H2O"
                i_atoms = General_Functions.glycan_to_atoms(i, permethylated, monosaccharides)
                if not lyase_digested:
                    i_atoms = General_Functions.sum_atoms(i_atoms, {'H': 2, 'O': 1})
                if tag[1] == 0.0:
                    if permethylated:
                        i_atoms = General_Functions.sum_atoms(i_atoms, {'C': 2, 'H': 4})
                        if reduced:
                            i_atoms = General_Functions.sum_atoms(i_atoms, {'C': 1, 'H': 4})
                    if not permethylated and reduced:
                        i_atoms = General_Functions.sum_atoms(i_atoms, {'H': 2})
                base_mass = mass.calculate_mass(composition=i_atoms)
                if permethylated:
                    i_atoms = General_Functions.sum_atoms(i_atoms, {'C': 1*s, 'H': 2*s, 'S': 1*s, 'O': 3*s}) #sum sulfation atoms
                    i_atoms = General_Functions.sum_atoms(i_atoms, {'C': 1*p, 'P': 1*p, 'O': 3*p, 'H': 3*p}) #sum phosphorylation atoms
                else:
                    i_atoms = General_Functions.sum_atoms(i_atoms, {'S': 1*s, 'O': 3*s}) #sum sulfation atoms
                    i_atoms = General_Functions.sum_atoms(i_atoms, {'P': 1*p, 'O': 3*p, 'H': 1*p}) #sum phosphorylation atoms
                
                if forced == "gags" and sodiums > 0:
                    i_formula = f"{i_formula}+{sodiums}Na"
                    i_atoms['H']-= sodiums
                    if "Na" in i_atoms:
                        i_atoms["Na"] += sodiums
                    else:
                        i_atoms["Na"] = sodiums
                    
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
                            if s > 0:
                                iso_corrected[-1] = iso_corrected[-1]*(1+(0.06*s)) #Sulfation correction
                            continue
                        if j_j == 2:
                            iso_corrected.append(abs(j*(10.8*(i_neutral_mass**-0.267)))) #default correction
                            if s > 0:
                                iso_corrected[-1] = iso_corrected[-1]*(1+(0.3*s)) #Sulfation correction
                            continue
                        if j_j == 3:
                            iso_corrected.append(abs(j*(122.62*(i_neutral_mass**-0.528)))) #default correction
                            if s > 0:
                                iso_corrected[-1] = iso_corrected[-1]*s #Sulfation correction
                            continue
                        if j_j == 4:
                            iso_corrected.append(abs(j*(2192.6*(i_neutral_mass**-0.833)))) #default correction
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
                    j, charges = General_Functions.fix_adduct_determine_charge(j)
                    
                    mz = mass.calculate_mass(i_atoms_tag,
                                             charge = charges,
                                             charge_carrier = j,
                                             carrier_charge = charges)
                    glycan_info['Adducts_mz'][General_Functions.comp_to_formula(j)] = mz
                    
                i_formulas.append(i_formula)
                glycan_infos.append(glycan_info)
    
    return glycan_infos, i_formulas

def include_internal_standard(full_library,
                              tag_mass,
                              fast,
                              high_res,
                              internal_standard,
                              permethylated,
                              reduced,
                              adducts_combo,
                              custom_monos = []):
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
    monosaccharides = copy.deepcopy(General_Functions.monosaccharides)
    
    # Add custom monosaccharides
    if len(custom_monos) > 0:
        for cm in custom_monos:
            monosaccharides[cm['cm_short_code']] = (cm['cm_name'], cm['cm_chem_comp'], General_Functions.sum_atoms({"C": 0, "O": 0, "N": 0, "H": 0}, General_Functions.form_to_comp(cm['cm_chem_comp'])), cm['cm_single_letter_code'])
                
    def_glycan_comp = {key : 0 for key in monosaccharides}
    
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
        i_atoms = General_Functions.calculate_comp_from_mass(i_neutral_mass)[0]
        i_type = 'direct_mass'
    except:
        i_composition = General_Functions.form_to_comp(internal_standard)
        
        if "C" in i_composition.keys():
            i_atoms = i_composition
            i_neutral_mass = mass.calculate_mass(composition=i_atoms)
            i_type = 'atomic_comp'
            
        else:
            i_composition = General_Functions.sum_monos(def_glycan_comp, i_composition, monos = monosaccharides)
            i_atoms = General_Functions.glycan_to_atoms(i_composition, permethylated, monosaccharides)
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
        full_library[i_formula]['Monos_Composition'] = {key : 0 for key in monosaccharides}
    full_library[i_formula]['Neutral_Mass'] = i_neutral_mass
    full_library[i_formula]['Neutral_Mass+Tag'] = i_neutral_tag
    full_library[i_formula]['Isotopic_Distribution'] = iso_corrected
    full_library[i_formula]['Isotopic_Distribution_Masses'] = i_iso_dist[1]
    
    full_library[i_formula]['Adducts_mz'] = {}
    for j in adducts_combo:
        j, charges = General_Functions.fix_adduct_determine_charge(j)
        
        mz = mass.calculate_mass(i_atoms_tag,
                                 charge = charges,
                                 charge_carrier = j,
                                 carrier_charge = charges)
        full_library[i_formula]['Adducts_mz'][General_Functions.comp_to_formula(j)] = mz
        
    return full_library

def calculate_glycan_fragments(glycan,
                               glycan_comp,
                               adduct_combos,
                               tolerance,
                               tag,
                               permethylated,
                               reduced,
                               lactonized_ethyl_esterified,
                               forced,
                               fragments_dict = {},
                               fragments_dict_per_glycan = {},
                               custom_monos = []):
    '''
    '''
    monosaccharides = copy.deepcopy(General_Functions.monosaccharides)
    # Add custom monosaccharides
    if len(custom_monos) > 0:
        for cm in custom_monos:
            if cm['cm_short_code'] not in monosaccharides.keys():
                monosaccharides[cm['cm_short_code']] = (cm['cm_name'], cm['cm_chem_comp'], General_Functions.sum_atoms({"C": 0, "O": 0, "N": 0, "H": 0}, General_Functions.form_to_comp(cm['cm_chem_comp'])), cm['cm_single_letter_code'])
                
    fragments_glycan = []
    def_glycan_comp = {key : 0 for key in monosaccharides}
    def_glycan_comp["T"] = 0
    sialics = ['S', 'G', 'Am', 'E', 'AmG', 'EG']
        
    glycolyl_sialics = ['AmG', 'EG', 'G']
    
    # Add sialics from custom monosaccharides
    if len(custom_monos) > 0:
        for cm in custom_monos:
            if cm['sialic']:
                sialics.append(cm['cm_short_code'])
        
    constraints = {data[3]: (0, glycan_comp[mono]) for mono, data in monosaccharides.items() if glycan_comp[mono] != 0}
    constraints["T"] = (0, 1)
    monos_chars = "".join([letter for letter in constraints])
    
    for i in range(1, sum(glycan_comp.values())+2):
        combinations = generate_combinations_with_constraints(monos_chars, i, constraints, monosaccharides)
        for j in combinations:
            fragments_glycan.append(General_Functions.sum_monos(def_glycan_comp, j, monos = monosaccharides))
    
    total_glycolyl_sialics_in_glycan = sum([glycan_comp[sialic] for sialic in sialics if sialic in glycolyl_sialics])
    total_non_glycolyl_sialics_in_glycan = sum([glycan_comp[sialic] for sialic in sialics if sialic not in glycolyl_sialics])
    total_sialics_in_glycan = total_glycolyl_sialics_in_glycan + total_non_glycolyl_sialics_in_glycan
    
    to_be_removed = []
    for i_i, i in enumerate(fragments_glycan):
        
        total_monos = sum(i.values())
        total_glycolyl_sialics = sum([i[sialic] for sialic in sialics if sialic in glycolyl_sialics])
        total_non_glycolyl_sialics = sum([i[sialic] for sialic in sialics if sialic not in glycolyl_sialics])
        total_sialics = total_glycolyl_sialics + total_non_glycolyl_sialics
        
        if ((i['T'] == 1 and total_monos == 1)
            or (total_sialics > total_sialics_in_glycan)
            or (total_non_glycolyl_sialics > total_non_glycolyl_sialics_in_glycan)
            or (total_glycolyl_sialics > total_glycolyl_sialics_in_glycan)
            or ((i['F'] > 0 and total_sialics > 0) and (i['H'] == 0 and i['N'] == 0))):
            to_be_removed.append(i_i)
            
        if forced == 'nglycan' and i_i not in to_be_removed: #some rules and hardcoded exceptions for N-Glycans
            
            if ((i['T'] == 1 and total_monos < 8 and total_sialics > 0)
                or (total_monos < 6 and total_sialics >= 1 and i['N'] > 1)
                or (i['H'] > 0 and i['T'] == 1 and i['N'] < 2)
                or (i['H'] == 2 and i['N'] == 1 and total_sialics == 0 and i['F'] == 0 and i['T'] == 1)
                or (i['H'] == 0 and i['N'] == 1 and total_sialics == 1 and i['F'] == 0 and i['T'] == 0)
                or (i['H'] == 1 and i['N'] == 3 and total_sialics == 0 and i['F'] == 0 and i['T'] == 1)
                or (i['H'] > 1 and i['N'] == 0 and total_sialics == 1 and i['F'] == 0 and i['T'] == 0)
                or (i['H'] == 1 and i['N'] == 1 and total_sialics == 0 and i['F'] == 0 and i['T'] == 1)
                or (i['H'] == 3 and i['N'] == 1 and total_sialics == 0 and i['F'] == 0 and i['T'] == 1)):
                to_be_removed.append(i_i)
                
    for i in sorted(to_be_removed, reverse = True):
        del fragments_glycan[i]
        
    fragments_glycan_formulas = []
    for fragment in fragments_glycan:
        fragment_formula = General_Functions.comp_to_formula(fragment)
        fragments_glycan_formulas.append(fragment_formula)
        
        if fragment_formula not in fragments_dict:
            temp_fragment = {}
            fragment_atoms = General_Functions.glycan_to_atoms(fragment, permethylated, monosaccharides)
            
            if fragment['T'] == 1:
                if tag[1] == 0.0:
                    if permethylated:
                        fragment_atoms = General_Functions.sum_atoms(fragment_atoms, {'C': 1, 'H': 2})
                        if reduced:
                            fragment_atoms = General_Functions.sum_atoms(fragment_atoms, {'O': 1})
                    if not permethylated and reduced:
                        fragment_atoms = General_Functions.sum_atoms(fragment_atoms, {'H': 2})
                fragment_atoms = General_Functions.sum_atoms(fragment_atoms, tag[0])
                
            temp_fragment['Monos_composition'] = fragment
            temp_fragment['Atoms'] = fragment_atoms
            temp_fragment['Adducts_mz'] = {}
            
            for adduct in adduct_combos:
                adduct, charges = General_Functions.fix_adduct_determine_charge(adduct)
                
                mz = mass.calculate_mass(fragment_atoms,
                                         charge = charges,
                                         charge_carrier = adduct,
                                         carrier_charge = charges)
                                         
                temp_fragment['Adducts_mz'][General_Functions.comp_to_formula(adduct)] = mz
                
            fragments_dict[fragment_formula] = temp_fragment
            
    # Generate H2O variants
    fragments_glycan_formulas_h2o = []
    for fragment in fragments_glycan_formulas:
        for water_amount in range(-1, 2, 2):
            new_fragment_name = f"{fragment}{water_amount if water_amount < 0 else '+'+str(water_amount)}H2O"
            fragments_glycan_formulas_h2o.append(new_fragment_name)
            
            if new_fragment_name not in fragments_dict:
                fragment_info = copy.deepcopy(fragments_dict[fragment])
                fragment_info['Atoms'] = General_Functions.sum_atoms(fragment_info['Atoms'], {'O': 1*water_amount, 'H': 2*water_amount})
                
                for adduct in adduct_combos:
                    adduct, charges = General_Functions.fix_adduct_determine_charge(adduct)
                    
                    mz = mass.calculate_mass(fragment_info['Atoms'],
                                             charge = charges,
                                             charge_carrier = adduct,
                                             carrier_charge = charges)
                                             
                    fragment_info['Adducts_mz'][General_Functions.comp_to_formula(adduct)] = mz
                    
                fragments_dict[new_fragment_name] = fragment_info
                
    fragments_glycan_formulas = fragments_glycan_formulas + fragments_glycan_formulas_h2o
    
    # Special fragments for permethylated glycans
    fragments_glycan_formulas_permethylated = []
    if permethylated:
        for fragment in fragments_glycan_formulas:
            new_fragment_name = f"{fragment}+CH2"
            fragments_glycan_formulas_h2o.append(new_fragment_name)
            
            if new_fragment_name not in fragments_dict:
                fragment_info = copy.deepcopy(fragments_dict[fragment])
                fragment_info['Atoms'] = General_Functions.sum_atoms(fragment_info['Atoms'], {'C': 1, 'H': 2})
                
                for adduct in adduct_combos:
                    adduct, charges = General_Functions.fix_adduct_determine_charge(adduct)
                    
                    mz = mass.calculate_mass(fragment_info['Atoms'],
                                             charge = charges,
                                             charge_carrier = adduct,
                                             carrier_charge = charges)
                                             
                    fragment_info['Adducts_mz'][General_Functions.comp_to_formula(adduct)] = mz
                    
                fragments_dict[new_fragment_name] = fragment_info

    fragments_glycan_formulas = fragments_glycan_formulas + fragments_glycan_formulas_permethylated
            
    fragments_dict_per_glycan[glycan] = fragments_glycan_formulas
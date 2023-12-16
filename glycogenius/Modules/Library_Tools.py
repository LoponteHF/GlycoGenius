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

from pyteomics import mzxml, mzml, mass, auxiliary
from itertools import combinations_with_replacement
from re import split
from math import inf
import sys
import datetime

##---------------------------------------------------------------------------------------
##Library generating-associated functions (these functions make vast use of the general
## functions).

def generate_glycans_library(min_max_mono,
                             min_max_hex,
                             min_max_hexnac,
                             min_max_sialics,
                             min_max_fuc,
                             min_max_ac,
                             min_max_gc,
                             n_glycan):
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

    min_max_sialics : tuple
        Minimum and maximum amount of sialic acids for the hypotethical glycans in the
        library. ie. (5, 20).

    min_max_fuc : tuple
        Minimum and maximum amount of deoxyhexoses for the hypotethical glycans in the
        library. ie. (5, 20).

    min_max_ac : tuple
        Minimum and maximum amount of N-Acetyl Neuraminic acids for the hypotethical
        glycans in the library. ie. (5, 20).

    min_max_ac : tuple
        Minimum and maximum amount of N-Glicolyl Neuraminic acids for the hypotethical
        glycans in the library. ie. (5, 20).

    n_glycan : boolean
        Indicates whether the function should force strict conditions based on the
        biological knowledge of glycans in order to avoid possible false positives when
        analysing N-glycans.

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
    def_glycan_comp = {"H": 0, "N": 0, "S": 0, "F": 0, "G": 0}
    for i in range(min_max_mono[0], min_max_mono[1]+1):
        for j in combinations_with_replacement("HNSFG", i):
            glycans.append(General_Functions.sum_monos(def_glycan_comp, General_Functions.count_seq_letters("".join(j))))
    to_be_removed = []
    for i_i, i in enumerate(glycans):
        if ((i['H'] < min_max_hex[0]) or (i['H'] > min_max_hex[1])
            or (i['N'] < min_max_hexnac[0]) or (i['N'] > min_max_hexnac[1])
            or (i['S']+i['G'] < min_max_sialics[0])
            or (i['S']+i['G'] > min_max_sialics[1])
            or (i['F'] < min_max_fuc[0]) or (i['F'] > min_max_fuc[1])
            or (i['S'] < min_max_ac[0]) or (i['S'] > min_max_ac[1])
            or (i['G'] < min_max_gc[0]) or (i['G'] > min_max_gc[1])):
            to_be_removed.append(i)
    if n_glycan:
        for i_i, i in enumerate(glycans):
            if ((i['S']+i['G'] > i['N']-2)
                or (i['F'] >= i['N']) or (i['S']+i['G'] > i['H']-3)
                or (i['H'] < 3) or (i['N'] < 2)):
                if i not in to_be_removed:
                    to_be_removed.append(i)
    for i in to_be_removed:
        glycans.remove(i)
    return glycans

def full_glycans_library(library,
                         max_adducts,
                         max_charges,
                         tag_mass = 0,
                         fast = True,
                         high_res = False):
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

    max_charges : int
        The maximum amount of charges to calculate per glycan.

    tolerance : float
        The mz acceptable tolerance. Used to calculate the intensoids, at this point.

    tag_mass : float
        The tag's added mass to the glycans, if the glycans are tagged.
        Default = 0 (No Tag).

    fast : boolean
        Makes the isotopic distribution calculation fast (less accurate) or not (more
        accurate).
        Default = True

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
    if tag_mass != 0:
        tag = General_Functions.calculate_comp_from_mass(tag_mass)
    else:
        tag = {"C": 0, "O": 0, "N": 0, "H": 0}
    adducts_combo = General_Functions.gen_adducts_combo(max_adducts, max_charges)
    for i in library:
        i_formula = General_Functions.comp_to_formula(i)
        i_atoms = General_Functions.sum_atoms(General_Functions.glycan_to_atoms(i), General_Functions.form_to_comp('H2O'))
        i_atoms_tag = General_Functions.sum_atoms(i_atoms, tag[0])
        i_neutral_mass = mass.calculate_mass(composition=i_atoms)
        i_neutral_tag = i_neutral_mass+tag[1]
        i_iso_dist = General_Functions.calculate_isotopic_pattern(i_atoms_tag, fast, high_res)
        full_library[i_formula] = {}
        full_library[i_formula]['Monos_Composition'] = i
        full_library[i_formula]['Atoms_Glycan+Tag'] = i_atoms_tag
        full_library[i_formula]['Neutral_Mass'] = i_neutral_mass
        full_library[i_formula]['Neutral_Mass+Tag'] = i_neutral_tag
        full_library[i_formula]['Isotopic_Distribution'] = i_iso_dist[0]
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
                      min_max_sialics,
                      min_max_fuc,
                      min_max_ac,
                      min_max_gc,
                      max_charges,
                      tolerance,
                      tag_mass):
    '''
    '''
    print("Building fragments library...", end = "", flush = True)
    glycans = []
    def_glycan_comp = {"H": 0, "N": 0, "S": 0, "F": 0, "G": 0, "T" : 0}
    for i in range(1, 18):
        for j in combinations_with_replacement("HNSFGT", i):
            glycans.append(General_Functions.sum_monos(def_glycan_comp, General_Functions.count_seq_letters("".join(j))))
    to_be_removed = []
    for i_i, i in enumerate(glycans):
        if ((i['T'] > 1) 
            or (i['H'] > min_max_hex[1])
            or (i['N'] > min_max_hexnac[1])
            or (i['S']+i['G'] > min_max_sialics[1])
            or (i['F'] > min_max_fuc[1])
            or (i['S'] > min_max_ac[1])
            or (i['G'] > min_max_gc[1])
            or (i['T'] >= 1 and sum(i.values()) < 8)):
            to_be_removed.append(i_i)
    for i in sorted(to_be_removed, reverse = True):
        del glycans[i]
    if tag_mass != 0:
        tag = General_Functions.calculate_comp_from_mass(tag_mass)
    else:
        tag = {"C": 0, "O": 0, "N": 0, "H": 0}
    adducts_combo = General_Functions.gen_adducts_combo({'H' : 2}, max_charges)
    frag_library = []
    combo_frags_lib = []
    adducts_mz = [[], [], []]
    removed = [[], [], []]
    for i_i, i in enumerate(glycans):
        for j_j, j in enumerate(range(-1, 2)):
            if j < 0:
                i_formula = General_Functions.comp_to_formula(i)+str(j)+'H2O'
            elif j > 0:
                i_formula = General_Functions.comp_to_formula(i)+'+'+str(j)+'H2O'
            else:
                i_formula = General_Functions.comp_to_formula(i)
            glycan_atoms = General_Functions.glycan_to_atoms(i)
            glycan_atoms['H'] += j*2
            glycan_atoms['O'] += j*1
            i_atoms = glycan_atoms
            i_neutral_mass = mass.calculate_mass(composition=i_atoms)
            if i['T'] == 1:
                i_atoms_tag = General_Functions.sum_atoms(i_atoms, tag[0])
                i_neutral_tag = i_neutral_mass+tag[1]
            else:
                i_atoms_tag = i_atoms
                i_neutral_tag = i_neutral_mass
            frag_library.append({})
            index = (i_i*3)+j_j
            frag_library[index]['Formula'] = i_formula
            frag_library[index]['Monos_Composition'] = i
            frag_library[index]['Adducts_mz'] = {}
            for j in adducts_combo:
                charges = sum(j.values())
                mz = mass.calculate_mass(i_atoms_tag,
                                         charge = charges,
                                         charge_carrier = j,
                                         carrier_charge = charges)
                found = False
                for k_k, k in enumerate(adducts_mz[0]):
                    if abs(mz - k) <= tolerance*2:
                        combo_frags_lib.append({})
                        combo_frags_lib[-1]['Formula'] = frag_library[adducts_mz[2][k_k]]['Formula']+'_'+adducts_mz[1][k_k]+'/'+frag_library[index]['Formula']+'_'+General_Functions.comp_to_formula(j)
                        combo_frags_lib[-1]['Adducts_mz'] = {}
                        combo_frags_lib[-1]['Adducts_mz'][General_Functions.comp_to_formula(j)] = mz
                        combo_frags_lib[-1]['Adducts_mz'][adducts_mz[1][k_k]] = mz
                        try:
                            del frag_library[adducts_mz[2][k_k]]['Adducts_mz'][adducts_mz[1][k_k]]
                        except:
                            pass
                        adducts_mz[0].append(mz)
                        adducts_mz[1].append(General_Functions.comp_to_formula(j))
                        adducts_mz[2].append(index)
                        found = True
                        break
                if not found:
                    adducts_mz[0].append(mz)
                    adducts_mz[1].append(General_Functions.comp_to_formula(j))
                    adducts_mz[2].append(index)
                    frag_library[index]['Adducts_mz'][General_Functions.comp_to_formula(j)] = mz
                    
    frag_library = frag_library+combo_frags_lib
    print("Done!")
    return frag_library
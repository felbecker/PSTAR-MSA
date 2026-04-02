#PSTAR-MSA
#Author: Marcel Gabor
#Institute for Mathematics and Computer Science, University of Greifswald
#11/04/2024

import argparse
#from codecs import ignore_errors
import sys
from datetime import datetime
import re
import os
import time
import pandas
import math
import random
import requests
import subprocess
import gzip
import shutil
import json
#import svg_stack
from itertools import tee
#import rsync
#import lxml
import traceback
import contextlib
try:
    from rcsbsearchapi.search import Query, SequenceQuery
    _rcsbsearchapi_available = True
except Exception:
    _rcsbsearchapi_available = False
from Bio import AlignIO
from Bio import SeqIO
from Bio import Seq
from Bio import SeqRecord
import matplotlib
import multiprocessing

#defining path to script
def set_script_path():
    
    global SCRIPT_DIR
    SCRIPT_DIR = os.path.abspath(os.path.dirname(__file__))
    os.chdir(SCRIPT_DIR)


#writing error messages to log
#void (str)
def errorFct(errorMsg):
    try:
        print(errorMsg)
        log_path = os.path.join(SCRIPT_DIR, "errorLog.txt")
        with open(log_path, 'a') as errFile:
            logEntry = "\n [%s] %s" % ((datetime.now()).strftime("%d/%m/%Y %H:%M:%S"), errorMsg)
            errFile.write(logEntry)
    except:
        open_file_err_msg = "Error log file could not be opened."
        if errorMsg == 'mp':
            return open_file_err_msg
        print(open_file_err_msg)
        exit(1)

#safely closing files, writing to error log if it fails
#void (file handle, str, str)
def close_file_safely(file, file_path, errMsgForward):
    try:
        file.close()
    except:
        errorMsg = "%s\nCould not close %s." % (errMsgForward, file_path)
        if errMsgForward == "mp":
            return errorMsg
        errorFct(errorMsg)
        exit(1)
        
#create dir in specified path
#void str        
def create_dir_safely(path): 
    
    try:
        os.mkdir(path)
    except FileExistsError:
        pass
    except:
        errMsg = "Failed to create \"%s\" directory. Insufficient writing priviliges in specified output path?" % os.path.split(path)[1]
        errorFct(errMsg)
        exit(1)

#writing list of dictionaries to custom CSV file importable by pandas into dataframe
#void ([{},{},...], str, str)
def write_list_of_dicts_to_csv(mydict_list, output_full_path):
    
    try:
        csv_file = open(output_full_path, 'w')
    except:
        errMsg = "CSV file %s could not be opened for writing." % output_full_path
        errorFct(errMsg)
        return
    #writing the csv, likely to be optimizable
    #remove trailing delimiters, consider using more robust delimiter
    try:
        #usually all dicts should have same size. if function is to be generalized, a check needs to be implemented
        for key in mydict_list[0].keys():
            csv_file.write(key + ",")
        csv_file.write("\n")
        for mydict in mydict_list:
            for value in mydict.values():
                csv_file.write(str(value)+",")
            csv_file.write("\n")
    except:
        errMsg = "An error occured while writing to %s." % output_full_path
        close_file_safely(csv_file, output_full_path, errMsg)
        errorFct(errMsg)

    close_file_safely(csv_file, output_full_path, '')

#function to select alignment sequence alphabet (e.g. protein, RNA, DNA)
#currently without choice    
#(str) -> [str,str,...]
def sel_alphabet(choice):
 
    if choice == "AA":
        #depending on the length of the protein sequence, different orders might be more efficient
        #see "Amino acid composition and protein dimension", Protein Sci. 17(12): 2187-2191 (2008)
        #maybe make into set. then rework of build_regex is required
        #this list of letters might change depending on what the database or the given structure aligner requires
        #for example, DALI ignores all "X" in sequences, that's why it's left out here for now
        #this results in cleaned sequences that have no "X" and will produce a match, albeit not a perfect match with sequences in the diamond
        #database that contains gaps (formerly "X"), excluding those sequences by applying our exact match definition will result in DALI never seeing sequences,
        #which contain "X" somewhere in the middle of the sequence (edges are ok)
        valid_symbols = ["A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","O","S","U","T","W","Y","V","B","Z","J"]
        #IUPAC-IUB Joint Commission on Biochemical Nomenclature.Nomenclature and Symbolism for Amino Acids and Peptides. Eur. J. Biochem. 138: 9-37 (1984)
    elif choice == "DNA/RNA":
        valid_symbols = ["A","G","C","T","U"]
    
    return valid_symbols

#takes string of invalid symbols and returns a list
#this function mostly exists only for factoring consistency
#str -> [str,str,...]
def set_non_alphabet(non_alphabet):
    
    return list(non_alphabet)

#build regex from list of valid symbols for clean_seq() function
#([str,str,...]) -> str
def build_regex_for_seq_cleaning_whitelist(valid_symbols):

    regex_str_list = []
    regex_str_list = valid_symbols.copy()
    valid_symbols_lower = list(map(str.lower, valid_symbols.copy()))
    regex_str_list.extend(valid_symbols_lower)
    regex_str_list.append("]")
    regex_str_list.insert(0,"^")
    regex_str_list.insert(0,"[")
    regex_str = ''.join(regex_str_list)
    
    return regex_str

def build_regex_for_seq_cleaning_blacklist(invalid_symbols):
    
    regex_str_list = []
    regex_str_list = invalid_symbols.copy()
    invalid_symbols_lower = list(map(str.lower, invalid_symbols.copy()))
    regex_str_list.append("]")
    regex_str_list.insert(0,"[")
    regex_str = ''.join(regex_str_list)
    
    return regex_str

#removing gaps and other symbols not in "valid_symbols" from given alignment sequence data
#preserves case
#(str, str) -> str
def clean_seq(seq):

    cleaned_seq = re.sub(CLEAN_SEQ_REGEX_STR,'',seq)
    
    return cleaned_seq

#function for custom selection of relevant metadata contained in the response body to be exported into CSV for further evaluation
#([str,str,...], {}) -> [{},{},...]
def parse_relevant_search_metadata(sample_seq_list, header_pdb_dict):
   
    relevant_item_dict = {}
    relevant_item_dict_list = []
    header_id_list = [row[1:3] for row in sample_seq_list]
    
    for header,internal_id in header_id_list:
        for entry in header_pdb_dict[header]:
            relevant_item_dict.clear()
            relevant_item_dict['identifier'] = entry['identifier']
            relevant_item_dict['original_fasta_header'] = header
            relevant_item_dict['score'] = entry['score']
            relevant_item_dict['sequence_identity'] = entry['services'][0]['nodes'][0]['match_context'][0]['sequence_identity']
            relevant_item_dict['evalue'] = entry['services'][0]['nodes'][0]['match_context'][0]['evalue']
            relevant_item_dict['bitscore'] = entry['services'][0]['nodes'][0]['match_context'][0]['bitscore']
            relevant_item_dict['alignment_length'] = entry['services'][0]['nodes'][0]['match_context'][0]['alignment_length']
            relevant_item_dict['mismatches'] = entry['services'][0]['nodes'][0]['match_context'][0]['mismatches']
            relevant_item_dict['gaps_opened'] = entry['services'][0]['nodes'][0]['match_context'][0]['gaps_opened']
            relevant_item_dict['query_beg'] = entry['services'][0]['nodes'][0]['match_context'][0]['query_beg']
            relevant_item_dict['query_end'] = entry['services'][0]['nodes'][0]['match_context'][0]['query_end']
            relevant_item_dict['subject_beg'] = entry['services'][0]['nodes'][0]['match_context'][0]['subject_beg']
            relevant_item_dict['subject_end'] = entry['services'][0]['nodes'][0]['match_context'][0]['subject_end']
            relevant_item_dict['query_length'] = entry['services'][0]['nodes'][0]['match_context'][0]['query_length']
            relevant_item_dict['subject_length'] = entry['services'][0]['nodes'][0]['match_context'][0]['subject_length']
            relevant_item_dict['internal_id'] = internal_id
            relevant_item_dict_list.append(relevant_item_dict.copy())

    return relevant_item_dict_list
 
#building sequence search query with set parameters and begin search. cast results to list.
#([str,str,int], str) -> []
def search_pdb_for_sequences(seq, aln_file_path, identity_cutoff):

    if not _rcsbsearchapi_available:
        errMsg = "rcsbsearchapi is not available (import failed). Cannot perform sequence search."
        errorFct(errMsg)
        return []

    raw_response_list = []
    result_list = []
    e_val_cutoff = 1
    try:
        results = SequenceQuery(seq[0],e_val_cutoff,identity_cutoff)
        #request_options is currently not accessible via exec(), e.g. request_options["sort_by"]="sequence_identity"...
        result_list = list(results(results_verbosity="verbose", return_type="polymer_entity"))
    except:
        errMsg = "Failed to get results for sequence with header %s in file %s with %s and %s in sequence search." %(seq[1], aln_file_path, e_val_cutoff, identity_cutoff)
        errorFct(errMsg)

    for entry in result_list:
        raw_response_list.append(entry)
        
    if len(raw_response_list) == 0:
        return []

    return raw_response_list

#extracting sequences and headers from alignment file in FASTA format and adding entry indices for good measure
# [sequence, header, number, cleaned_sequence]
#(str) -> [[str,str,int,str],[str,str,int,str],...]
def import_seq_list_from_fasta_aln(aln_file_path):

    seq_list = []

    try:
        aln_file = open(aln_file_path, "r")
    except:
        errMsg = "Could not open alignment file %s." % aln_file_path
        errorFct(errMsg)
        raise
    try:
        alignment = AlignIO.read(aln_file, "fasta")
    except ValueError:
        close_file_safely(aln_file, aln_file_path, '')
        try:
            alignment = SeqIO.parse(aln_file_path, "fasta")
        except:
            errMsg = "Could not read fasta file %s." %aln_file_path
            errorFct(errMsg)
            raise
        finally:
            #implement this in a way that warnings are not alwaays written when a cleaned fasta file is written
            warn_msg = "Warning: Sequences in %s might not all have the same length." % aln_file_path
            close_file_safely(aln_file, aln_file_path, warn_msg)
            errorFct(warn_msg)


    internal_id = 0

    #dict for duplicate check
    seen = {}

    for record in alignment:

        #check if key has been seen before
        try:
            if seen[str(record.id)]:
                warn_msg = "Warning: Duplicate header %s in file %s." % (str(record.id),aln_file_path)
                errorFct(warn_msg)
        except KeyError:
            pass
        internal_id += 1
        cleaned_seq = clean_seq(str(record.seq))
        seq_list.append([str(record.seq), str(record.id), internal_id, cleaned_seq])

        seen[str(record.id)] = True

    return seq_list

#reading file names in specified path and returning list of full paths
#(str) -> [str,str,...]
def get_file_paths_in_dir(path):

    path_file_list = []
 
    for element in os.listdir(path):
        full_path = os.path.join(path, element)
        if os.path.isfile(full_path):
            path_file_list.append(full_path)

    return path_file_list

#function that returns a list of PDB-IDs, one for each fasta header in CSV. default is the one with the highest SEQUENCE IDENTITY
#(str) -> ([str,str,...],{})
def get_ids_from_csv(csv_file_path, sort_by):
    
    id_list = []
    internal_id_pdb_name_dict = {}
    df = pandas.read_csv(csv_file_path)
    #remove next line, once trailing delimiters have been removed
    df = df.drop(['Unnamed: 16'], axis=1)
    #split df into one dataframe per fasta_header
    df_list = [x for _,x in df.groupby('original_fasta_header')]
    for split_df in df_list:
        if sort_by == "identity":
            split_df = split_df.sort_values(by=['sequence_identity'], ascending=False)
        elif sort_by == "score":
            split_df = split_df.sort_values(by=['score'], ascending=False)
        else:
            errMsg = "Can't sort by %s. Exiting." % sort_by
            errorFct(errMsg)
            exit(1)
        pdb_identifier = (split_df['identifier'].iloc[0])[0:4]
        id_list.append(pdb_identifier)
        internal_id_pdb_name_dict[pdb_identifier] = split_df['internal_id'].iloc[0]
    
    return id_list, internal_id_pdb_name_dict

#function to automatically send download requests to PDB according to a list of PDB-IDs,
#gathered from the generated CSV files
#inspired by https://www.rcsb.org/docs/programmatic-access/batch-downloads-with-shell-script
#(str,str,int) -> {}
def download_pdb_files_by_CSV(csv_file_path, out_path, sort_by, verbosity):

    counter = 1
    id_list, internal_id_pdb_name_dict = get_ids_from_csv(csv_file_path, sort_by)
    id_list_len = len(id_list)
    
    for identifier in id_list:
        if verbosity > 0:
            msg = "Downloading \"%s\"... [%s/%s]" % (identifier, counter, id_list_len)
            print(msg)
            counter += 1
            
        ext_id = identifier + ".pdb"
        url = "https://files.rcsb.org/download/%s" % ext_id
        try:
            r = requests.get(url)
        except:
            errMsg = "Failed to get download request for ID %s." % identifier
            errorFct(errMsg)
            continue
        try:
            with open(os.path.join(out_path, ext_id), 'wb') as f:
                f.write(r.content)
        except:
            errMsg = "Failed to write to %s." % ext_id
            errorFct(errMsg)
    
    #maybe refactor so that id_list and internal_id_name_dict is passed to download_pdb_files_by_CSV directly
    return internal_id_pdb_name_dict

#calculates the reference "sum-of-pairs-set" for structural alignments (lowercase chars mean not equivalent)
#which is a set containing double-nested 2-tuples ( i.e. ((x,y),(a,b)) )
#where x is the row index of the query structure, a is the row index of the sbjct structure
#y is the column index of the query structure, b is the column index of the sbjct structure
#if ((x,y),(a,b)) is in ref_ind_set, then this means that according to the structure alignment
#the sequence with index x at residue number y is equivalent to the sequence with index a at residue number b
#((int,int),(str,str)) -> set( ((),()) )
def calc_ref_SP_set(sequence_ids, aligned_sequences):

    query = aligned_sequences[0]
    sbjct = aligned_sequences[1]
    query_gaps = 0
    sbjct_gaps = 0
    ref_ind_set = set()

    #raw query and sbjct strings should have the same length in alignment (incl. gaps)
    if len(query) == len(sbjct):
        max_len = len(query)
    elif len(query) < len(sbjct):
        errMsg = "Warning: Aligned subject sequence is longer than query sequence."
        errorFct(errMsg)
        max_len = len(sbjct)
    else:
        errMsg = "Warning: Aligned query sequence is longer than subject sequence."
        errorFct(errMsg)
        max_len = len(query)
    for i in range(0, max_len):
        
        try:
            equiv_sbjct = True
            equiv_query = True
            #need to define a set of expected symbols or set of gap/unknown symbols
            if sbjct[i] == "-" or sbjct[i] == ".":
                sbjct_gaps += 1
                equiv_sbjct = False
            if query[i] == "-"  or sbjct[i] == ".":
                query_gaps += 1
                equiv_query = False
            #if one of the sequences has a gap, it's a mismatch
            if equiv_sbjct is False or equiv_query is False:
                continue
            if sbjct[i].islower():
                continue
            if query[i].islower():
                continue
            #is the following is an inefficient approach, when we know what symbols to expect?
            #if sbjct[i] not in valid_symbols:
             #   continue
            #if query[i] not in valid_symbols:
             #   continue
            ref_ind_set.add(((sequence_ids[0],i-query_gaps),(sequence_ids[1],i-sbjct_gaps)))
        except IndexError:
            break
    
    return ref_ind_set

#function to calculate MSA sum-of-pairs-set (different equivalence relation)
#its partly duplicate to "calc_ref_SP_set" and will most likely be refactored together with latter function
#((int,int),(str,str)) -> set( ((),()) )
def calc_MSA_SP_set(sequence_ids, aligned_sequences):
    
    query = aligned_sequences[0]
    sbjct = aligned_sequences[1]
    query_gaps = 0
    sbjct_gaps = 0
    MSA_ind_set = set()
    #raw query and sbjct strings should have the same length in alignment (incl. gaps)
    if len(query) == len(sbjct):
        max_len = len(query)
    elif len(query) < len(sbjct):
        errMsg = "Warning: Aligned subject sequence is longer than query sequence."
        errorFct(errMsg)
        max_len = len(sbjct)
    else:
        errMsg = "Warning: Aligned query sequence is longer than subject sequence."
        errorFct(errMsg)
        max_len = len(query)
    for i in range(0, max_len):
        try:
            #define a set of expected symbols or set of gap/unknown symbols
            sbjct_gap = False
            query_gap = False
            if sbjct[i] == "-" or sbjct[i] == ".":
                sbjct_gaps += 1
                sbjct_gap = True
            if query[i] == "-"  or query[i] == ".":
                query_gaps += 1
                query_gap = True
            if sbjct_gap or query_gap:
                continue
            MSA_ind_set.add(((sequence_ids[0],i-query_gaps),(sequence_ids[1],i-sbjct_gaps)))
        except IndexError:
            break
    
    return MSA_ind_set
        

#function for naive approach to maximizing intersection between two sets by shifting all element values in set_2 by some constant
#sets need to contain mutable elements
# (set(),set()) -> (set(),int)
def argmax_intersection(set_1, set_2):
    
    #make copies of original sets
    orig_set_1 = set_1.copy()
    orig_set_2 = set_2.copy()

    #initial shift bounds +-1000
    shift_bounds = range(1000,-1001,-1)
    max_intersection = orig_set_1.intersection(orig_set_2)
    cur_max = len(max_intersection)
    cur_argmax = orig_set_2
    max_shift = 0
    for shift in shift_bounds:  
        new_set_2 = set()
        set_intersection = set()
        for element in set_2:
            y = element[0][1] + shift
            b = element[1][1] + shift
            x = element[0][0]
            a = element[1][0]
            new_set_2.add(((x,y),(a,b)))
        set_intersection = set_1.intersection(new_set_2)
        if len(set_intersection) > cur_max:
            max_intersection = set_intersection.copy()
            cur_max = len(max_intersection)
            cur_argmax = new_set_2.copy()
            max_shift = shift
        set_1 = orig_set_1
        set_2 = orig_set_2

    return cur_argmax, max_shift

#writes SP sets in job and resulting SPS to CSV for subsequent statistical evaluation
# 
def write_job_SPS_CSV(ref_set_list, aln_set_list, job_SPS):
    
    pass

#calculates SP scores for job and returns list of dictionaries containing information about each alignment
#highly specific signature
#([[{},{},...],[{},{},...],...], {}, {}) -> ([{},{},...], float)
def calc_job_SPS(job_list_of_dict_lists, internal_id_pdb_name_dict, internal_id_raw_seq_dict):
    current_max_ref_SP_set = set()
    save_dict = {}
    SPS_dict_list = []
    max_ref_SP_set_list = []
    argmax_aln_SP_set_list = []
    job_SPS = 0
    #for each alignment folder in job, i.e. for each sub job
    #for-else statement to continue outer loop in case of breaking of inner loop
    for dict_list in job_list_of_dict_lists:
        aln_SP_set = set()
        pair_SPS = 0
        pdb_name_1_const_check = dict_list[0]["query"][0:4]
        pdb_name_2_const_check = dict_list[0]["sbjct"][0:4]
        pdb_internal_id_1_const = internal_id_pdb_name_dict[pdb_name_1_const_check]
        pdb_internal_id_2_const = internal_id_pdb_name_dict[pdb_name_2_const_check]
        sub_job_counter = 0
        #for each chain alignment in PDB alignment 
        #calculating reference SP set for DALI alignment
        for dic in dict_list:
            ref_SP_set = set()
            pdb_name_1 = dic["query"][0:4]
            pdb_name_2 = dic["sbjct"][0:4]
            if (pdb_name_1 != pdb_name_1_const_check) or (pdb_name_2 != pdb_name_2_const_check):
                errMsg = "Differing PDB names detected within sub job. Skipping sub job with data:\n %s" % str(dict_list)
                errorFct(errMsg)
                break
            ref_SP_set = calc_ref_SP_set((pdb_internal_id_1_const, pdb_internal_id_2_const), (dic["DALI_query_seq"], dic["DALI_sbjct_seq"]))
            #DALI generates an unknown number of TXT files or alignments per TXT file depending on the number of chains per PDB file and what file acts as query or subject
            #for example if PDB1 contains 2 chains and PDB2 contains 3 chains, if PDB1 is query and PDB2 is subject, DALI will generate
            #2 TXT files, namely "PDB1A.txt" and "PDB1B.txt". each will contain 3 alignments. Thus we are left with 6 alignments.
            #namely: PDB1A vs. PDB2A, PDB1A vs. PDB2B, PDB1A vs. PDB2C, PDB1B vs. PDB2A, ... and so forth
            #the following check chooses the set with maximum length out of all possible chain alignments, so that we dont accidentally align
            #dissimilar chains with each other 
            if (len(ref_SP_set) > len(current_max_ref_SP_set)) or sub_job_counter == 0:
                current_max_ref_SP_set.clear()
                save_dict.clear()
                current_max_ref_SP_set = ref_SP_set.copy()
                save_dict = dic.copy()
            sub_job_counter += 1
        #calculate SP set for original alignment
        #outer loop    
        else:
            max_ref_SP_set_list.append(current_max_ref_SP_set.copy())
            aln_seq_1 = internal_id_raw_seq_dict[pdb_internal_id_1_const]
            aln_seq_2 = internal_id_raw_seq_dict[pdb_internal_id_2_const]
            aln_SP_set = calc_MSA_SP_set((pdb_internal_id_1_const, pdb_internal_id_2_const), (aln_seq_1, aln_seq_2))
            #extend saved aln dict with additional data
            save_dict["SP_DALI"] = len(current_max_ref_SP_set)
            save_dict["SP_orig_aln"] = len(aln_SP_set)
            save_dict["query_internal_id"] = pdb_internal_id_1_const
            save_dict["sbjct_internal_id"] = pdb_internal_id_2_const
            save_dict["orig_query_seq"] = aln_seq_1
            save_dict["orig_sbjct_seq"] = aln_seq_2
            #calculate SPS
            
            #try to maximize intersection size by shifting y and b in ((x,y),(a,b)) in aln_SP_set by some constants
            #try to find shift pattern to minimize maximization effort
            #this could actually be viable. but since we know from the search meta data where in the PDB the particular sequence match starts
            #we should use that information first
            #in any way, maximizing intersection increases values
            argmax_aln_SP_set, max_shift = argmax_intersection(current_max_ref_SP_set, aln_SP_set)
            argmax_aln_SP_set_list.append(argmax_aln_SP_set.copy())
            save_dict["max_shift"] = max_shift
            try:
                pair_SPS = len(argmax_aln_SP_set.intersection(current_max_ref_SP_set))/len(current_max_ref_SP_set)
            except ZeroDivisionError:
                pair_SPS = 0
            save_dict["pair_SPS"] = pair_SPS
            SPS_dict_list.append(save_dict.copy())
            
    #first unite all current_max_ref_SP_sets. unite aln_SP_sets. then intersect unions and divide by union of cur_max_ref_SP_sets.
    overall_ref_SP_set = set()
    overall_MSA_SP_set = set()
    for s in max_ref_SP_set_list:
        overall_ref_SP_set = overall_ref_SP_set.union(s)
    for s in argmax_aln_SP_set_list:
        overall_MSA_SP_set = overall_MSA_SP_set.union(s)
    try:
        job_SPS = len(overall_ref_SP_set.intersection(overall_MSA_SP_set))/len(overall_ref_SP_set)
    except ZeroDivisionError:
        #if length of denominator set is 0, it means that the intersection must be empty as well
        job_SPS = 0
        
    return SPS_dict_list, job_SPS

#import function for DALI outputs. one pass, line per line. might need rework,
#because it's hard to read and looks ugly. hopefully it's fast(-ish), at least.
#(str) -> ([{},{},...])
def import_DALI_aln(aln_file_path):

    try:
        alignment = open(aln_file_path, 'r')
    except:
        errorMsg = "Could not open %s for reading." % aln_file_path
        errorFct(errorMsg)
        exit(1)

    #read file line by line
    lineCount = 0
    queryFound = False
    sbjctFound = False
    query_list = []
    sbjct_list = []
    aln_dict_list = []
    aln_dict = {}
    aln_counter = 0
    while True:
        lineCount += 1
        line = alignment.readline()
        if not line:
            break

        #import jobname
        job_match = re.search(r"^#[ ]*Job:[ ]*(.*)", line)
        if job_match is not None:
            job_name = job_match.group(1)
            job_match = None

        #import header lines
        header_match = re.search(r"^[ ]*[0-9]+:[ ]*[a-zA-Z0-9-]*[ ]+([0-9.]*)[ ]+([0-9.]*)[ ]+([0-9.]*)[ ]+([0-9.]*)[ ]+([0-9.]*)[ ]+([A-Za-z0-9.: -]*)", line)
        if header_match is not None:
            aln_dict.clear()
            aln_dict['job_name'] = job_name
            aln_dict['Z-score'] = header_match.group(1)
            aln_dict['rmsd'] = header_match.group(2)
            aln_dict['lali'] = header_match.group(3)
            aln_dict['nres'] = header_match.group(4)
            aln_dict['perc_id'] = header_match.group(5)
            aln_dict['pdb_descr'] = header_match.group(6)
            aln_dict_list.append(aln_dict.copy())
            header_match = None
            continue
        
        
        #we're not using structural equivalences or transrot at the moment,
        #thus we can end parsing once we reached this point
        if aln_counter > 0 and (line.startswith("# S") or line.startswith("# T")):
            overall_query = ''.join(query_list)
            overall_sbjct = ''.join(sbjct_list)
            aln_dict_list[aln_counter-1]['DALI_query_seq'] = overall_query
            aln_dict_list[aln_counter-1]['DALI_sbjct_seq'] = overall_sbjct
            close_file_safely(alignment, aln_file_path, "")
            return aln_dict_list


        #set aln_counter
        aln_count_match = re.search(r"^No[ ]+([0-9]+):[ ]*Query=([a-zA-Z0-9]+)[ ]*Sbjct=([a-zA-Z0-9]+)", line)
        if aln_count_match is not None:
            aln_counter = int(aln_count_match.group(1))
            aln_dict_list[aln_counter-1]['query'] = aln_count_match.group(2)
            aln_dict_list[aln_counter-1]['sbjct'] = aln_count_match.group(3)
            #check if there's previous alignment. if we're here, it means that we're at the start of a new alignment
            if aln_counter > 1:
                overall_query = ''.join(query_list)
                overall_sbjct = ''.join(sbjct_list)
                aln_dict_list[aln_counter-2]['DALI_query_seq'] = overall_query
                aln_dict_list[aln_counter-2]['DALI_sbjct_seq'] = overall_sbjct
                query_list = []
                sbjct_list = []
            aln_count_match = None
            continue

        #search for Query sequence
        if queryFound is False:
            query_match = re.search(r"^Query ([a-zA-Z-]+)", line)
            if query_match is None:
                continue
            else:
                queryFound = True
                continue
        #search for Sbjct sequence
        if sbjctFound is False:
            sbjct_match = re.search(r"^Sbjct ([a-zA-Z-]+)", line)
            if sbjct_match is None:
                continue
            else:
                sbjctFound = True
                continue
            
        if (queryFound and sbjctFound) is True:
            query_list.append(query_match.group(1))
            queryFound = False
            sbjct_list.append(sbjct_match.group(1))
            sbjctFound = False
            continue

    #still assemble sequences, if end of file is reached
    if aln_counter <= 0:
        aln_dict.clear()
        #get query and sbjct name from dir name
        try:
            head_tail = os.path.split(aln_file_path)
            head_tail2 = os.path.split(head_tail[0])
            #if chain identifiers are used for naming of dirs
            if len(head_tail2[1]) == 11:
                query_name = head_tail2[1][0:5]
                sbjct_name = head_tail2[1][6:11]
            elif len(head_tail2[1]) == 9:
                query_name = head_tail2[1][0:4]
                sbjct_name = head_tail2[1][5:9]
            else:
                errMsg = "Sub job dir name %s is not of expected format." % head_tail[0]
                errorFct(errMsg)
                raise NameError
        except:
            query_name = "NA"
            sbjct_name = "NA"
            errMsg = "Path \"%s\" does not contain any information about query and subject name at expected location." % aln_file_path
            errorFct(errMsg)
        aln_dict['query'] = query_name
        aln_dict['sbjct'] = sbjct_name
        #set remaining dict entries
        aln_dict['job_name'] = job_name
        aln_dict['Z-score'] = "NA"
        aln_dict['rmsd'] = "NA"
        aln_dict['lali'] = "NA"
        aln_dict['nres'] = "NA"
        aln_dict['perc_id'] = "NA"
        aln_dict['pdb_descr'] = "NA"
        aln_dict['DALI_query_seq'] = "0"
        aln_dict['DALI_sbjct_seq'] = "0"
        aln_dict_list.append(aln_dict.copy())
        errMsg = "File \'%s\' contains no alignments." % aln_file_path
        errorFct(errMsg)
        close_file_safely(alignment, aln_file_path, errMsg)
        return aln_dict_list
    
    if (len(query_list) == len(sbjct_list)) and (queryFound == sbjctFound):
        overall_query = ''.join(query_list)
        overall_sbjct = ''.join(sbjct_list)
        aln_dict_list[aln_counter-1]['DALI_query_seq'] = overall_query
        aln_dict_list[aln_counter-1]['DALI_sbjct_seq'] = overall_sbjct           
    else:
        errorMsg = "Differing amounts of regex matches for Query and Sbjct sequences in %s" % aln_file_path
        close_file_safely(alignment, aln_file_path, errorMsg)
        errorFct(errorMsg)
        exit(1)

    close_file_safely(alignment, aln_file_path, '')
    return aln_dict_list

#import all alignment files in a dir, appends job_list_of_dict_list exactly ONCE
#(str) -> [{},{},...]
def import_aln_files_in_dir(aln_path):
    
    #dict_list for alignment dir
    aln_dict_list = []

    #get all TXT files
    txt_path_list = []
    for file in os.listdir(aln_path):
        if file.endswith(".txt"):
            txt_path_list.append(os.path.join(aln_path, file))
    #import all TXT files in dir and extend list
    for txt_path in txt_path_list:
        aln_dict_list.extend(import_DALI_aln(txt_path))

    return aln_dict_list

#import all alignments in one job
#(str, str) -> [[{},{},...], [{},{},...],...]
def import_aln_files_in_job(job_title, out_path):

    job_list_of_dict_lists = []
    job_file_list = []
    job_aln_path = os.path.join(out_path, job_title, "ALN")
    #list all subdirs in job dir
    job_file_list = os.listdir(job_aln_path)
    for file in job_file_list:
        if not os.path.isdir(os.path.join(job_aln_path, file)):
            job_file_list.remove(file)
    for aln_dir in job_file_list:
        aln_path = os.path.join(job_aln_path, aln_dir)
        job_list_of_dict_lists.append(import_aln_files_in_dir(aln_path))
        
    return job_list_of_dict_lists

#build regex from list of valid symbols for clean_seq() function
#([str,str,...]) -> str
def build_regex_for_seq_cleaning_whitelist(valid_symbols):

    regex_str_list = []
    #all upper case alphabet letters
    regex_str_list = valid_symbols.copy()
    valid_symbols_lower = list(map(str.lower, valid_symbols.copy()))
    #and all lower case alphabet letters
    regex_str_list.extend(valid_symbols_lower)
    regex_str_list.append("]")
    regex_str_list.insert(0,"^")
    regex_str_list.insert(0,"[")
    regex_str = ''.join(regex_str_list)
    
    return regex_str

def build_regex_for_seq_cleaning_blacklist(invalid_symbols):
    
    regex_str_list = []
    regex_str_list = invalid_symbols.copy()
    invalid_symbols_lower = list(map(str.lower, invalid_symbols.copy()))
    regex_str_list.append("]")
    regex_str_list.insert(0,"[")
    regex_str = ''.join(regex_str_list)
    
    return regex_str

#takes a list of sequence-header pairs and returns a random sample of given size
#([[str,str,int],[str,str,int],...],int) -> [[str,str,int],[str,str,int],...]
def choose_random_sample_from_aln_file(seq_list, sample_size):

    sample_seq_list = []
    random_int_list = []
    counter = 0

    try:
        sample_size = int(sample_size)
    except:
        errMsg = "Sample size \"%s\" must be given as an integer." % sample_size
        errorFct(errMsg)
        exit(0)

    if sample_size >= len(seq_list):
        sample_seq_list = seq_list
        return sample_seq_list

    while counter < sample_size:

        random_int = random.randint(0, len(seq_list)-1)
        if random_int in random_int_list:
            continue
        random_int_list.append(random_int)
        sample_seq_list.append(seq_list[random_int].copy())
        counter += 1

    return sample_seq_list

#counts number of "real" proteins in alignment
#irrelevant function
#(str,str) -> (str,int,int,float) 
def count_homologues(aln_file_path, hom_aln_file_path):
    
    seq_list = import_seq_list_from_fasta_aln(aln_file_path)
    hom_seq_list = import_seq_list_from_fasta_aln(hom_aln_file_path)
    
    no_of_proteins = len(seq_list)
    no_of_proteins_and_homologues = len(hom_seq_list)
    prot_frac = no_of_proteins/no_of_proteins_and_homologues
    head, tail = os.path.split(hom_aln_file_path)
    
    #might want to write this to CSV
    msg = "There are a total of %s sequences in file %s. %s of them are proteins. Fraction: %s" % (no_of_proteins_and_homologues, tail, no_of_proteins, prot_frac)
    print(msg)
    
    return (tail, no_of_proteins_and_homologues, no_of_proteins, prot_frac)

#wrapper function to parse fasta alignment sequences and headers, query a PDB search and write relevant metadata to CSV in one function call
#verbosity option allows progress tracking in terminal
#([str,str,...],str,str,int,int) -> {}
def write_aln_file_search_hits_to_csv(valid_symbols, aln_file_path, output_path, verbosity, sample_size):

    seq_list = []
    regex_str = ""
    header_pdb_dict = {}

    if(verbosity>0):
        msg = "Parsing fasta sequences from %s." % aln_file_path
        print(msg)

    #rename seq_list, because it's not only a list of sequences 
    seq_list = import_seq_list_from_fasta_aln(aln_file_path)

    if sample_size != 0:
        sample_seq_list = choose_random_sample_from_aln_file(seq_list, sample_size)
    else:
        sample_seq_list = seq_list
    
    #delete next line, if everything runs fine
    regex_str = build_regex_for_seq_cleaning_whitelist(valid_symbols)
    
    #make a copy of seq list before sequences are cleaned. returned by function.
    raw_sample_seq_list = sample_seq_list.copy()
    internal_id_raw_seq_dict = {}
    #fill dictionary to get sequence by internal sequence id or header later, returned by function. not elegant.
    for raw_sample in raw_sample_seq_list:
        internal_id_raw_seq_dict[raw_sample[2]] = raw_sample[0]
    
    if(verbosity>0):
        msg = "Searching for sequences in alignment file %s." % (os.path.split(aln_file_path)[1])
        seq_no = len(sample_seq_list)
        counter = 0
        print(msg)

    for seq in sample_seq_list:
        
        identity_cutoff = 1.00
        
        #pass regex_str to clean_seq() again, if something goes wrong
        seq[0] = clean_seq(str(seq[0]), valid_symbols)
        if seq[1] in header_pdb_dict.keys():
            errMsg = "Duplicate header \"%s\" in alignment file %s." %(seq[1], os.path.split(aln_file_path)[1])
            errorFct(errMsg)
        
        header_pdb_dict[seq[1]] = search_pdb_for_sequences(seq, aln_file_path, identity_cutoff)

        #this loop keeps lowering the identity cutoff until proteins have been found
        #currently theres no good and easy way to sort search hits via the python interface
        while not header_pdb_dict[seq[1]]:
            if (verbosity>0):
                notification = "No matches for header \'%s\' with identity cutoff %s." %(seq[1], identity_cutoff)
                print(notification)

            header_pdb_dict[seq[1]] = search_pdb_for_sequences(seq, aln_file_path, identity_cutoff)
            if identity_cutoff < 0 or math.isclose(identity_cutoff, 0.00,abs_tol=1e-03):
                header_pdb_dict[seq[1]] = []
                break
            if identity_cutoff < 0.30:
                identity_cutoff = 0.00
            if math.isclose(identity_cutoff,1.00,abs_tol=1e-03): 
                identity_cutoff -= 0.02
            if identity_cutoff > 0.90:
                identity_cutoff -= 0.05
            if identity_cutoff < 0.90:
                identity_cutoff -= 0.2
        
        if(verbosity>0):
            counter += 1
            msg = "%s/%s. Header: \"%s\". Done!" % (counter, seq_no, seq[1])
            print(msg)
    
    #extract list of headers from list of sample and headers
    #header_list = [row[1] for row in sample_seq_list]

    if(verbosity>0):
            print("Parsing search metadata.")

    relevant_search_results_list = parse_relevant_search_metadata(sample_seq_list, header_pdb_dict)

    if(verbosity>0):
            msg = "Writing CSV to %s" % (output_path)
            print(msg)
    
        
    aln_file_name = os.path.basename(aln_file_path)
    csv_output_file_name = "%s.csv" % aln_file_name
    csv_output_full_path = os.path.join(output_path, csv_output_file_name)
    write_list_of_dicts_to_csv(relevant_search_results_list, csv_output_full_path)

    return internal_id_raw_seq_dict

#function to write error log. for errors encountered during subprocesses involving the shell
#void str        
def shell_err_fct(errMsg):
    try:
        log_path = os.path.join(SCRIPT_DIR, "shell_err_log.txt")
        with open(log_path, 'a') as errFile:
            logEntry = "\n [%s] %s" % ((datetime.now()).strftime("%d/%m/%Y %H:%M:%S"), errMsg)
            errFile.write(logEntry)
    except:
        open_err_file_msg = "Error log file could not be opened."
        if errMsg == "mp":
            return open_err_file_msg
        print(open_err_file_msg)
        exit(1)

#returns list of PDB names from given dir of PDB files
#str -> [str,str,...]        
def get_pdb_name_list(pdb_path):
    
    ent_gz = False
    pdb_path_list = get_file_paths_in_dir(pdb_path)
    pdb_name_list = []
    for pdb_path in pdb_path_list:
        head, file_name = os.path.split(pdb_path)
        #this loop works for .pdb and for .ent.gz formats
        while "." in file_name:
            file_name = os.path.splitext(file_name)[0]
        #standard rsync pdb database has file names of format "pdbxxxx.ent.gz"
        if len(file_name) == 7:
            file_name = file_name[3:8]
            ent_gz = True
        pdb_name_list.append(file_name)
        
    return pdb_name_list, ent_gz

#http://ekhidna2.biocenter.helsinki.fi/dali/README.v5.html
#creates subprocess to run import.pl of DaliLite.v5 to import PDBs in a given dir to DAT file format used by DALI
#void (str, str, str, int)
def DALI_import_PDBs(pl_bin_path, pdb_path, DAT_path, verbosity):

    pdb_name_list, ent_gz = get_pdb_name_list(pdb_path)
        
    import_pl_full_path = os.path.join(pl_bin_path, "import.pl")

    #create dir for output and err
    create_dir_safely(os.path.join(DAT_path, "err_logs"))
    
    create_dir_safely(os.path.join(DAT_path, "out_logs"))
    
    #this loop is necessary for clean error logging, since DaliLite.v5 seems to have bad logging
    counter = 0
    no_of_files = len(pdb_name_list)
    for pdb_name in pdb_name_list:
        
        shell_input = []
        if ent_gz:
            pdb_full_path = os.path.join(pdb_path, "pdb"+pdb_name+".ent.gz")
        else:
            pdb_full_path = os.path.join(pdb_path, pdb_name+".pdb")
        err_file_path = os.path.join(DAT_path, "err_logs/", pdb_name)
        out_file_path = os.path.join(DAT_path, "out_logs/", pdb_name)

        counter += 1
        if verbosity>0:
            msg = "Importing PDBs... [%s/%s] " % (counter, no_of_files)
            print(msg, end="\r")
        
        #/.../<dali_dir>/import.pl --pdbfile <path> --pdbid <id> --dat <path> --verbose --clean
        shell_input = [import_pl_full_path, '--pdbfile', pdb_full_path, '--pdbid', pdb_name.upper(), '--dat', DAT_path, '--clean', '1', '>', out_file_path, '2', '>', err_file_path]
        
        #create file handles for stdout and stderr
        try:
            err_file = open(err_file_path, 'w+')
        except:
            errMsg = "Failed to open file %s." % err_file_path
            errorFct(errMsg)
            exit(1)
        try:
            out_file = open(out_file_path, 'w+')
        except:
            errMsg = "Failed to open file %s." % out_file_path
            errorFct(errMsg)
            exit(1)

        #subprocess for importing PDBs
        #shell input already sends stdout to out_file_path and stderr to err_file_path
        #probably just pipe output in the next line    
        process = subprocess.Popen(shell_input, stdout=out_file, stderr=err_file)
        
        try:
            out, err = process.communicate()
            #its really just a lot to print...
            #if out is not None and verbosity>1:
                #print(out)
        except:
            process.kill()
            out, err = process.communicate()
            errMsg = "Communication with import.pl subprocess timed out. Killed it. Forwarding error to shell_err_log.txt."
            close_file_safely(err_file, err_file_path, errMsg)
            close_file_safely(out_file, out_file_path, errMsg)
            shell_err_fct(err)
            errorFct(errMsg)
            if out is not None:
                print(out)
            exit(1)
            
        #close files
        close_file_safely(err_file, err_file_path, "")
        close_file_safely(out_file, out_file_path, "")

        if os.path.exists(os.path.join(SCRIPT_DIR,'dali.lock')):
            err_msg = "Error while importing pdb file %s with import.pl. Skipping." % pdb_name.upper()
            errorFct(err_msg)
            os.remove(os.path.join(SCRIPT_DIR,'dali.lock'))
    
    #rewrite a sufficient check of whether import has been successful
    for file_name in os.listdir(DAT_path):
        if file_name.endswith(".d"):
            errMsg = "Failed to convert PDB %s to DAT format. Exiting." % file_name
            errorFct(errMsg)
            exit(1)

#http://ekhidna2.biocenter.helsinki.fi/dali/README.v5.html
#creates subprocess to run dali.pl of DaliLite.v5. Uses "pairwise" functionality of DALI for better control over what alignments are being made.
#writes stderr and stdout to files for each job and alignment. improved output dir structure, since DaliLite.v5 doesn't seem to support specifying an output dir
#dir structure: <jobtitle>/ALN/<queryPDBname_sbjctPDBname>/
#if chain_list is passed, only the specified chains of the PDBs will be used in alignment       
#void (str,str,str,str,str,int,[str,str,...])
def DALI_all_vs_all_query(pl_bin_path, pdb_path, out_path, DAT_path, job_title, verbosity, chain_list=[]):

    #for some weird implementation reasons (probably Fortran though) DALI doesn't generate alignment files for job titles longer than 12 chars
    #and generates only empty alignments for job titles exactly 12 chars long
    #this is a DALI problem and can't be changed at the moment
    #i could, however probably find a workaround, but this isn't really worth the effort right now
    #maybe it's also just because job_title needs to be passed with "" around it.
    if len(job_title) > 11:
        errMsg = "Note: Job title is longer than 11 characters. Will be cut off in DALI output files due to title length limitations."
        job_title = job_title[0:11]
    
    #saving cwd for later file system navigation
    wd = os.getcwd()
    
    #changing relative paths to absolute
    if os.path.isabs(pl_bin_path):
        dali_pl_full_path = os.path.join(pl_bin_path, "dali.pl")
    else:
        dali_pl_full_path = os.path.join(wd, pl_bin_path, "dali.pl")
        
    if not os.path.isabs(pdb_path):
        pdb_path = os.path.join(wd, pdb_path)
        
    if not os.path.isabs(out_path):
        out_path = os.path.join(wd, out_path)
    
    if not os.path.isabs(DAT_path):
        DAT_path = os.path.join(wd, DAT_path)
    
    chains_used = False
    #if optional argument chain_list is passed, use that instead    
    if len(chain_list) != 0:
        pdb_name_list = list(chain_list).copy()
        chains_used = True
    else:
        pdb_name_list = get_pdb_name_list(pdb_path)
        
    
    #creating job dir
    job_path = os.path.join(out_path, job_title)
    create_dir_safely(job_path)
    
    #Dali's all-against-all feature behaves in unexpected ways. 
    #That's why manual looping through the PDB list for generation of pairwise alignments.
    #this is also beneficial, since it allows full control over what alignments are to be generated
    #and individual sub job error logging
    n = len(pdb_name_list)
    no_of_alns = sum(range(1,n))
    if n<2:
        errMsg = "Please provide more than one PDB file for alignment."
        errorFct(errMsg)
        exit(0)
    
    #we dont want redundant or trivial alignments:
    # e.g. (1,1) is trivial, (1,2) = (2,1) are (effectively?) redundant
    
    counter = 0    
    for i in range(0, n):
        k = i+1
        for j in range(k, n):
            
            if chains_used:
                try:
                    chain_id1 = str(pdb_name_list[i][0:4])+str(pdb_name_list[i][5])
                    chain_id2 = str(pdb_name_list[j][0:4])+str(pdb_name_list[j][5])
                except IndexError:
                    errMsg = "Chain identifier %s or %s is shorter than 6 chars." % (pdb_name_list[i], pdb_name_list[j])
                    errorFct(errMsg)
                    raise

            shell_input = []
            sub_job = "%s_%s" % (chain_id1, chain_id2)
            
            counter += 1 
            if verbosity>0:
                msg = "Calculating alignment %s. [%s/%s] " % (sub_job, counter, no_of_alns)
                print(msg, end="\r")
            
            #creating ALN dir
            aln_path = os.path.join(job_path, "ALN")
            create_dir_safely(aln_path)
            
            #creating sub job dir
            sub_job_path = os.path.join(aln_path, sub_job)    
            create_dir_safely(sub_job_path)
                          
            err_file_path = os.path.join(sub_job_path, "err_log")
            out_file_path = os.path.join(sub_job_path, "output_log")
            if chains_used:
                #/.../<dali_dir>/dali.pl --cd1 <chain_id1> --cd2 <chain_id2> --dat1 <path> --dat2 <path> --title <string> \
                # --outfmt "summary,alignments,equivalences,transrot" --clean 1> <out_log_file> 2> <err_log_file>
                shell_input = [dali_pl_full_path, '--cd1', chain_id1, '--cd2', chain_id2, 
                           '--dat1', DAT_path, '--dat2', DAT_path, '--title', job_title,
                           '--outfmt', "summary,alignments,equivalences,transrot",
                           '--clean', '1', '>', out_file_path, '2', '>', err_file_path]
            else:
                #/.../<dali_dir>/dali.pl --pdbfile1 <path> --pdbid1 <id> --pdbfile2 <path> --pdbid2 <id> --dat1 <path> --dat2 <path> --title <string> \
                # --outfmt "summary,alignments,equivalences,transrot" --clean 1> <out_log_file> 2> <err_log_file>
                #technically --dat argument makes no sense here
                pdb_full_path_query = os.path.join(pdb_path, pdb_name_list[i]+".pdb")
                pdb_full_path_sbjct = os.path.join(pdb_path, pdb_name_list[j]+".pdb")    
                shell_input = [dali_pl_full_path, '--pdbfile1', pdb_full_path_query, '--pdbid1', pdb_name_list[i], 
                               '--pdbfile2', pdb_full_path_sbjct, '--pdbid2', pdb_name_list[j], '--dat1', DAT_path,
                               '--dat2', DAT_path, '--title', job_title, '--outfmt', "summary,alignments,equivalences,transrot",
                               '--clean', '1', '>', out_file_path, '2', '>', err_file_path]

            #navigate to output dir in loop
            os.chdir(sub_job_path)
            
            #create file handles for stdout and stderr, only needed for Popen. rework maybe
            try:
                err_file = open(err_file_path, 'w+')
            except:
                errMsg = "Failed to open file %s." % err_file_path
                errorFct(errMsg)
                exit(1)
            try:
                out_file = open(out_file_path, 'w+')
            except:
                errMsg = "Failed to open file %s." % out_file_path
                errorFct(errMsg)
                exit(1)

            #subprocess for alignment generation
            #need to PIPE stdout and stderr for output forwarding    
            process = subprocess.Popen(shell_input, stdout=out_file, stderr=err_file)
    
            try:
                #timeout may need to be longer (or shorter), depending on size of alignments
                #for the calculation of pairwise alignments of up to a few chains
                #500 seconds seems to be reasonable though
                out, err = process.communicate()
                #if out is not None and verbosity>1:
                    #print(out)
            except subprocess.TimeoutExpired:
                process.kill()
                out, err = process.communicate()
                errMsg = "Communication with dali.pl subprocess timed out. Killed it. Forwarding error to shell_err_log.txt"
                close_file_safely(err_file, err_file_path, errMsg)
                close_file_safely(out_file, out_file_path, errMsg)
                shell_err_fct(err)
                errorFct(errMsg)
                if out is not None:
                    print(out)
                exit(1)
            
            #close files
            close_file_safely(err_file, err_file_path, "")
            close_file_safely(out_file, out_file_path, "")
            
            #change dir back to previous wd
            os.chdir(wd)
            
#wrapper function to calculate SPS of a (sub-)set of aligned sequences 
#(str, str, str, str, int, [str,str,...], int) -> (float, [{},{},...])
def calc_SPS_from_aln_sample(output_path, aln_file_path, job_title, dali_path, sample_size, valid_symbols, verbosity):
    
    job_out_path = os.path.join(output_path, job_title)
    create_dir_safely(job_out_path)
    
    job_data_out_path = os.path.join(job_out_path, "DATA")
    create_dir_safely(job_data_out_path)
    
    CSV_out_path = os.path.join(job_data_out_path, "CSV")
    create_dir_safely(CSV_out_path)
    
    #refactor write_aln_file_search_hits_to_csv
    internal_id_raw_seq_dict = write_aln_file_search_hits_to_csv(valid_symbols, aln_file_path, CSV_out_path, verbosity, sample_size)
    
    pdb_out_path = os.path.join(job_data_out_path, "PDB_lib")
    create_dir_safely(pdb_out_path)
    
    #refactor download_pdb_files_by_CSV
    aln_file_name = os.path.split(aln_file_path)[1]
    CSV_file_path = os.path.join(CSV_out_path, aln_file_name+".csv")
    internal_id_pdb_name_dict = download_pdb_files_by_CSV(CSV_file_path, pdb_out_path, "identity", verbosity)
    
    dat_out_path = os.path.join(job_data_out_path, "DAT_lib")
    create_dir_safely(dat_out_path)

    DALI_import_PDBs(dali_path, pdb_out_path, dat_out_path, verbosity)
    
    DALI_all_vs_all_query(dali_path, pdb_out_path, output_path, dat_out_path, job_title, verbosity)
    
    job_list_of_dict_lists = import_aln_files_in_job(job_title, output_path)
    
    SPS_dict_list, job_SPS = calc_job_SPS(job_list_of_dict_lists, internal_id_pdb_name_dict, internal_id_raw_seq_dict)
    
    print(SPS_dict_list)
    print(job_SPS)

    SPS_out_path = os.path.join(job_data_out_path, "SPS")
    create_dir_safely(SPS_out_path)
    SPS_csv_full_out_path = os.path.join(SPS_out_path, job_title+"_SPS.csv")
    #write SPS_dict_list to csv
    write_list_of_dicts_to_csv(SPS_dict_list, SPS_csv_full_out_path)
    
    return job_SPS, SPS_dict_list

#function to take a list of svg images and concatenates them to a large svg image
#void ([str,str,...],str)    
def append_svg_images(svg_path, out_img_full_path):
    
    svg_full_path_list = []
    svg_file_name_list = os.listdir(svg_path)
    for svg_file_name in svg_file_name_list:
        svg_full_path = os.path.join(svg_path, svg_file_name)
        if os.path.isfile(svg_full_path) is False:
            continue
        svg_full_path_list.append(svg_full_path)

    doc = svg_stack.Document()
    
    layout1 = svg_stack.HBoxLayout()
    for svg_full_path in svg_full_path_list:
        layout1.addSVG(svg_full_path, alignment=svg_stack.AlignCenter)
    doc.setLayout(layout1)
    doc.save(out_img_full_path)

#evaluate and create overview of search data
#void (str, str, str)    
def evaluate_rcsb_search_data(csv_path, output_path, sort_by):
    
    ex_match_counter = 0
    total_counter = 0
    csv_file_list = get_file_paths_in_dir(csv_path)
    
    for csv_file_path in csv_file_list:
        df = pandas.read_csv(csv_file_path)
        new_df = pandas.DataFrame(columns=['original_fasta_header','sequence_identity','score','mismatches','gaps_opened'])
        row_list = []
        #remove next line, once trailing delimiters have been removed
        df = df.drop(['Unnamed: 16'], axis=1)
        #split df into one dataframe per fasta_header
        df_list = [x for _,x in df.groupby('original_fasta_header')]
        for split_df in df_list:
            if sort_by == "identity":
                split_df = split_df.sort_values(by=['sequence_identity'], ascending=False)
            elif sort_by == "score":
                split_df = split_df.sort_values(by=['score'], ascending=False)
            else:
                errMsg = "Can't sort by %s. Exiting." % sort_by
                errorFct(errMsg)
                exit(1)
            total_counter += 1
            #append highest score/identity row to df list
            if math.isclose(split_df.iloc[[0]]['sequence_identity'],1,abs_tol=10e-5):
                ex_match_counter += 1
            row_list.append(split_df.iloc[[0]].copy())
        #combine list into one df
        new_df = pandas.concat(row_list, ignore_index=True)
        plot_df = pandas.DataFrame(new_df[['sequence_identity','score']])
        plot_df['mean_identity'] = new_df['sequence_identity'].mean()
        csv_file_name = os.path.splitext(os.path.split(csv_file_path)[1])[0]
        #pandas.plot returns matplotlib.axes.Axes object
        axes_obj = plot_df.plot(title=csv_file_name)
        #we can get figure object from axes_obj using get_figure() method and save it to vector image
        axes_obj.get_figure().savefig(os.path.join(output_path,csv_file_name+".svg"))
    print("Exact matches: %s.Total number of sequences: %s." % (ex_match_counter, total_counter))

#evaluate diamond blastp TSV files and create overview plots for structure availability
#optimizable (?) iloc is massively slowing things down.
#try with df.set_value(), almost 3 times faster probably, but its deprecated. annoying  
# void (str,str,str,int)    
def evaluate_diamond_search_data(tsv_path, output_path, sort_by, verbosity):
    
    tsv_file_list = get_file_paths_in_dir(tsv_path)
    exact_match_dict = {}
    file_counter = 1
    no_of_files = len(tsv_file_list)

    for tsv_file_path in tsv_file_list:
        if verbosity>0:
            msg = "[%s/%s] Plotting %s..." % (file_counter, no_of_files, tsv_file_path)
            print(msg)
        tsv_ex_match_counter = 0
        tsv_total_counter = 0
        df = pandas.read_csv(tsv_file_path, sep='\t')
        new_df = pandas.DataFrame(columns=['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore','qlen'])
        row_list = []
        df['is_exact_match'] = 0
        ex_match_col_num = df.columns.get_loc('is_exact_match')
        mismatch_col_num = df.columns.get_loc('mismatch')
        len_col_num = df.columns.get_loc('length')
        qlen_col_num = df.columns.get_loc('qlen')
        #split df into one dataframe per fasta_header
        df_list = [x for _,x in df.groupby('qseqid')]
        for split_df in df_list:
            if sort_by == "identity":
                split_df = split_df.sort_values(by=['pident'], ascending=False)
            elif sort_by == "score":
                split_df = split_df.sort_values(by=['bitscore'], ascending=False)
            else:
                errMsg = "Can't sort by %s. Exiting." % sort_by
                errorFct(errMsg)
                exit(1)
            tsv_total_counter += 1
            #append highest score/identity row to df list
            #if math.isclose(split_df.iloc[[0]]['pident'],1,abs_tol=10e-1):
             #   ex_match_counter += 1
            if int(split_df.iloc[0,mismatch_col_num]) == 0 and int(split_df.iloc[0,len_col_num]) == int(split_df.iloc[0,qlen_col_num]):
                tsv_ex_match_counter += 1
                split_df.iloc[0,ex_match_col_num] = 100
            row_list.append(split_df.iloc[[0]].copy())
        tsv_file_name = os.path.splitext(os.path.split(tsv_file_path)[1])[0]
        exact_match_dict[tsv_file_name] = (tsv_ex_match_counter/tsv_total_counter)*100
        #combine list into one df
        new_df = pandas.concat(row_list, ignore_index=True)
        plot_df = pandas.DataFrame(new_df[['pident','is_exact_match']])
        plot_df['mean_identity'] = new_df['pident'].mean()
        plot_df['exact_match_perc'] = (tsv_ex_match_counter/tsv_total_counter)*100

        plot_df.to_csv(os.path.join(output_path,"data_avail.tsv"), sep='\t')

        #pandas.plot returns matplotlib.axes.Axes object
        plot_df_len = len(plot_df)
        xtick_list = list(range(0,plot_df_len,math.floor(plot_df_len/5)))
        axes_obj = plot_df[['pident','is_exact_match',]].plot(title=tsv_file_name, kind="bar", yticks=[0,2.5,5,7.5,10,15,20,25,30,40,50,60,70,80,90,100],
                                                             xticks=xtick_list)
        plot_df[['mean_identity','exact_match_perc']].plot(kind="line", ax=axes_obj)
        #we can get figure object from axes_obj using get_figure() method and save it to vector image
        axes_obj.get_figure().savefig(os.path.join(output_path, tsv_file_name+".svg"))
        #matplot.pyplot.figure objects need to be explicitly closed. otherwise they are all open at the same time
        matplotlib.pyplot.close()
        file_counter += 1
    exact_perc_df = pandas.DataFrame.from_dict(data=exact_match_dict, orient='index', columns=['ex_match_perc'])
    exact_perc_df['mean_exact_perc'] = exact_perc_df['ex_match_perc'].mean()
    axes_obj = exact_perc_df[['ex_match_perc']].plot(title="Exact match distribution among test alignments", kind="bar")
    exact_perc_df[['mean_exact_perc']].plot(kind="line", ax=axes_obj, color='orange')
    axes_obj.set_xticklabels(exact_perc_df.index.values, rotation = 45)
    axes_obj.tick_params(axis='x', which='major', labelsize=2)
    axes_obj.get_figure().savefig(os.path.join(output_path, "mean_ex_match_perc.svg"))
    matplotlib.pyplot.close()

#recursively get all pdb file paths in a specified path and return them as list
# str -> [str, str, ...]        
def get_pdb_path_list_recursively(loc_pdb_db_path):
    
    pdb_path_list = []
    
    for root, dirs, files in os.walk(loc_pdb_db_path):
        for file in files:
            if ".pdb" in file.lower() or ".ent" in file.lower():
                pdb_path_list.append(os.path.join(root, file))
                
    return pdb_path_list

#extract pdb ID from pdb file path and write it to index file as line
# (str,file_handle) -> (str)
def write_pdb_index_entry_from_file_path(raw_line, index_file_handle):

    #this might need to be generalized more for other structure formats downloaded by rsync
    pdb_id_match = re.search(r'(.+\/)*pdb(....)\..*', raw_line)
    pdb_id = pdb_id_match.group(2)
    index_line = pdb_id+"\n"
    index_file_handle.write(index_line)
    
    return pdb_id

#wrapper function to write pdb index entries according to list of paths to pdb files
# void (str, [str,str,...])
def write_index_file_from_file_path_list(old_index_file_full_path, pdb_path_list):
    
    try:
        old_index_file_handle = open(old_index_file_full_path, "w")
    except:
        errMsg = "Cannot open file %s." % old_index_file_full_path
        errorFct(errMsg)
        exit(1)
    try:
        old_index_file_handle.truncate()
        for pdb_path in pdb_path_list:
            write_pdb_index_entry_from_file_path(pdb_path, old_index_file_handle)
    except:
        errMsg = "Error while writing PDB IDs to index file %s." % old_index_file_full_path
        close_file_safely(old_index_file_handle, old_index_file_full_path, errMsg)
        exit(1)
            
    close_file_safely(old_index_file_handle, old_index_file_full_path, "")

#writes an index file from the contents of a set
# void (str, set())   
def write_index_file_from_set(index_file_full_path, pdb_id_set):
    
    try:
        index_file_handle = open(index_file_full_path, "w")
    except:
        errMsg = "Cannot open file %s." % index_file_full_path
        errorFct(errMsg)
        exit(1)
    try:
        index_file_handle.truncate()
        for pdb_id in pdb_id_set:
            index_line = str(pdb_id).lower()+"\n"
            index_file_handle.write(index_line)
    except:
        errMsg = "Error while writing PDB IDs to index file %s." % index_file_full_path
        close_file_safely(index_file_handle, index_file_full_path, errMsg)
        exit(1)
            
    close_file_safely(index_file_handle, index_file_full_path, "")
        
    

#read index file and return set with pdb IDs
# (file_handle) -> set()    
def read_index_file(old_index_file_handle):
    
    old_index_set = set()

    pdb_ids = old_index_file_handle.readlines()
    for pdb_id in pdb_ids:
        pdb_id = pdb_id.rstrip().lower()
        old_index_set.add(pdb_id)
        
    return old_index_set

#function to be called before editing the database fasta file
#creates a copy of the fasta file in the same location
# void (str,str)
def backup_fasta_file(fasta_file_full_path):
     
    dest_backup = os.path.join(os.path.split(fasta_file_full_path)[0], "pdb_db_backup.fasta")
    cp_process = subprocess.Popen(['cp', fasta_file_full_path, dest_backup], stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    try:
        out, err = cp_process.communicate()
    except:
        cp_process.kill()
        out, err = cp_process.communicate()
        shell_err_fct(err)
        raise
    
    return dest_backup

#called when editing fasta file has been completed successfully
#double checks to see, if original fasta file is present in the same location
# void (str,str)    
def remove_fasta_backup(dest_backup, fasta_file_full_path):
    
    #just to be extra safe
    if os.path.exists(fasta_file_full_path):
        os.remove(dest_backup)

#called when exception is raised during editing of fasta file
#simply overwrites the (potentially) broken, edited fasta file with the backup and renames it
# void (str,str)        
def restore_fasta_file_from_backup(dest_backup, fasta_file_full_path):
    
    mv_process = subprocess.Popen(['mv', dest_backup, fasta_file_full_path], stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    try:
        out, err = mv_process.communicate()
    except:
        mv_process.kill()
        out, err = mv_process.communicate()
        shell_err_fct(err)
        raise

#function for removal of specific entries in fasta file
#unfortunately theres no perfect way to do it at the moment
#lines need to be rewritten when we want to delete some
#but maybe the biopython parser is what makes it slow
#just try it later and see how performance does
# void (set(),str)    
def remove_fasta_db_entries_by_header(removal_set, diamond_db_path):

    fasta_file_full_path = os.path.join(diamond_db_path, "pdb_db.fasta")
    with open(fasta_file_full_path, 'r+') as fasta_file:
        for line in fasta_file:
            if line[0] == ">":
                if line[1:5].lower() in removal_set:
                    pass

#takes a SeqRecord iterator, if chain.id does not have expected format and rewrites the whole iterator
#fixing those entries that don't adhere to XXXX:X format in the header by getting pdb id from file name
#if not even the chain identifier can be read (accessing first char of chain.id throws IndexError)
#the header will simply contain the pdb id without chain information
# (iterator, str) -> (iterator)                
def update_chains_iterator(iterator, pdb_file_path):
         
    new_chains_list = []

    for chain in iterator:
        try:
            chain_id_str = str(chain.id)
            if chain_id_str[4] != ":":
                pdb_id_match = re.search(r'(.+\/)*pdb(....)\..*', pdb_file_path)
                pdb_id = pdb_id_match.group(2)
                new_chain_id = pdb_id+":"+chain_id_str[0].upper()
                new_chain_seq = str(chain.seq)
                new_chain = SeqIO.SeqRecord(Seq.Seq(new_chain_seq), id=new_chain_id)
                new_chains_list.append(new_chain)
            else:
                #new_chain_id = str(chain.id)
                #new_chain_seq = Seq.Seq(str(chain.seq))
                #new_chain = SeqIO.SeqRecord(new_chain_seq, id=new_chain_id)
                new_chains_list.append(chain)
        except IndexError:
            pdb_id_match = re.search(r'(.+\/)*pdb(....)\..*', pdb_file_path)
            pdb_id = pdb_id_match.group(2)
            try:
                new_chain_id = (pdb_id+":"+chain_id_str[0]).upper()
                new_chain_seq = str(chain.seq)
                new_chain = SeqIO.SeqRecord(Seq.Seq(new_chain_seq), id=new_chain_id)
                new_chains_list.append(new_chain)
            except IndexError:
                new_chain_id = pdb_id
                new_chain_seq = str(chain.seq)
                new_chain = SeqIO.SeqRecord(Seq.Seq(new_chain_seq), id=new_chain_id)
                new_chains_list.append(new_chain)
                
    return iter(new_chains_list)
                
#extracts chain numbers, pdbid and sequences from pdb file and writes them to FASTA file
#also writes biopython output to dedicated log file
#consider parallelization of pdb file parser, it's very slow like this
#i suspect the biopython pdb parser     
#void (str, file_handle, file_handle)
def write_loc_pdb_to_fasta_file(pdb_file_path, fasta_file_handle, biopython_log_file_handle):
    
    #make sure that this also works with uncompressed files
    try:
        pdb_file_handle = gzip.open(pdb_file_path, 'rt')
    except:
        errMsg = "Failed to open %s." % pdb_file_path
        errorFct(errMsg)
        exit(1)
    #supressing warnings for this line, specifically Bio.PDB.PDBExceptions.PDBConstructionWarning does NOT work in ANY way
    #not with BiopythonWarnings nor with BiopythonExperimentalWarnings
    #redirecting all stdout and stderr streams to a different file instead
    with contextlib.redirect_stdout(biopython_log_file_handle):
        with contextlib.redirect_stderr(sys.stdout):
            try:
                #version with SeqResIterator
                #chains = SeqIO.PdbIO.PdbSeqresIterator(pdb_file_handle)

                #version with AtomIterator
                #atom iterator is speed killer.
                #IMPORTANT REMARK: AtomIterator ends up producing FASTA files down the line that
                #confuse DIAMOND. The current FTP path of rsync also downloads DNA files.
                #something about the handling of HETATOMS in AtomIterator causes problems (wrong Biopython parser).
                #stick to using SeqResIterator for now
                #UPDATE: what makes the use of AtomIterator problematic is that it automatically fills the first n-1
                #AAs with "X", if the AA count starts at values n > 1. We'd also need to cover the edge case, where the first
                #listed AA is unrecognized, thus gets marked as "X"
                #we'll need to remove those leading "X" for it to work properly with DALI, since DALIs importer
                #does not use them, and additionally NOT remove "X", if its not recognized and the first listed AA
                chains = SeqIO.PdbIO.PdbAtomIterator(pdb_file_handle)  
            except:
                errMsg = "Failed to parse chain data from %s." % pdb_file_path
                close_file_safely(pdb_file_handle, pdb_file_path, errMsg)
                errorFct(errMsg)
                return
    
    update = False
    #copies of iterators are required for a) counting length, b) looping to check for incomplete headers
    #c) passing to update_chains_iterator() function 
    chains, chains_copy, chains_copy_2, chains_copy_3 = tee(chains, 4)
    chains_choice = chains

    #if parser fails to get valid chain ID, try to get PDB ID from file name
    #this is not a particularly elegant implementation, but at least it requires only
    #a minimal amount of rewriting iterators  
    for chain in chains_copy_2:
        try:
            chain_id_str = str(chain.id)
            if chain_id_str[4] != ":":
                new_chains = update_chains_iterator(chains_copy_3, pdb_file_path)
                update = True
                break
        except IndexError:
            new_chains = update_chains_iterator(chains_copy_3, pdb_file_path)
            update = True
            break
        
    if update:
        #check len of old iterator
        old_iter_len = 0
        for chain in chains_copy:
            old_iter_len +=1
        new_chains, new_chains_copy = tee(new_chains)
        #check len of new iterator
        new_iter_len = 0
        for chain in new_chains_copy:
            new_iter_len += 1
        if old_iter_len == new_iter_len:
            chains_choice = new_chains
        else:
            errMsg = "Warning: Failed to update chains for %s." % (pdb_file_path)
            errorFct(errMsg)

    try:
        SeqIO.write(chains_choice, fasta_file_handle, "fasta")
    except:
        #try:
            #   repair_recent_fasta_entries(fasta_file_handle, chains)
        #except:
            #   errMsg = "Error while attempting to repair most recent FASTA entries."
            #  errorFct(errMsg)
            # raise Exception
        print(traceback.format_exc())
        errMsg = "Error while writing %s to FASTA file. Skipping." % pdb_file_path
        close_file_safely(pdb_file_handle, pdb_file_path, errMsg)
        errorFct(errMsg)
    close_file_safely(pdb_file_handle, pdb_file_path, "")

#takes a set of pdb IDs and appends existing fasta file with the elements of this set
#throws FileNotFoundException, if indexed PDB ID is not in local PDB database
# void (set(),str,str,str,int)                
def add_fasta_db_entries(added_set, fasta_file_full_path, loc_pdb_db_path, biopython_log_full_path, verbosity, backup_freq):
    
    if verbosity>0:
        print("Appending FASTA file...")
    try:
        dest_backup = backup_fasta_file(fasta_file_full_path)
    except:
        head, tail = os.path.split(fasta_file_full_path)
        errMsg = "Failed to back up %s. Insufficient writing priviliges in %s?" % (tail, head)
        errorFct(errMsg)
        exit(1)
    
    try:
        biopython_log_file_handle = open(biopython_log_full_path)
    except:
        errMsg = "File %s could not be opened for writing." % biopython_log_full_path
        errorFct(errMsg)
        exit(1)

    try:
        counter = 1
        number_of_pdbs = len(added_set)
        with open(fasta_file_full_path, 'a+') as fasta_file_handle:
            for pdb_id in added_set:
                if verbosity>0:
                    msg = "[%s/%s]              " % (str(counter), str(number_of_pdbs))
                    print(msg, end="\r")
                pdb_file_path = os.path.join(loc_pdb_db_path, pdb_id[1:3], "pdb"+pdb_id+".ent.gz")
                write_loc_pdb_to_fasta_file(pdb_file_path, fasta_file_handle, biopython_log_file_handle)
                counter += 1
                #this is a bit shaky. maybe change it in the future
                if backup_freq != 0 and counter % backup_freq == 0:
                    fasta_file_handle.seek(0)
                    dest_backup = backup_fasta_file(fasta_file_full_path)
                    fasta_file_handle.read()
                    fasta_file_handle.write('\n')
    except:
        try:
            restore_fasta_file_from_backup(dest_backup, fasta_file_full_path)
        except:
            print(traceback.format_exc())
            errMsg = "Error while restoring FASTA file from backup %s." % dest_backup
            errorFct(errMsg)
        finally:
            print(traceback.format_exc())
            errMsg = "Error while writing to %s." % fasta_file_full_path
            close_file_safely(fasta_file_handle, fasta_file_full_path, errMsg)
            close_file_safely(biopython_log_file_handle, biopython_log_full_path, errMsg)
            errorFct(errMsg)
            exit(1)

    remove_fasta_backup(dest_backup, fasta_file_full_path)
    close_file_safely(biopython_log_file_handle, biopython_log_full_path, "")

#updates index file by reading the old index file and a new index file.
#the new index file is generated by recursively listing all pdb files that are in the local database.
#this ensures that the up-to-date index really only contains IDs that exist locally.    
#the index entries are then added to two different sets. the resulting difference set is subsequently used
#to write only the new entries to the fasta file (s. function "sync_pdb_copy()")    
#this funciton gets called for every rsync call
# (str,str,str,str) -> set()    
def update_pdb_index(loc_pdb_db_path, diamond_db_path, tmp_index_file_full_path, old_index_file_full_path):

    tmp_index_set = set()
    old_index_set = set()
    pdb_id = ""

    try:
        tmp_index_file_handle = open(tmp_index_file_full_path, 'w')
    except:
        errMsg = "Cannot open file %s." % tmp_index_file_full_path
        errorFct(errMsg)
        exit(1)
    try:      
        pdb_path_list = get_pdb_path_list_recursively(loc_pdb_db_path)
        for pdb_path in pdb_path_list:
            pdb_id = write_pdb_index_entry_from_file_path(pdb_path, tmp_index_file_handle)
            tmp_index_set.add(pdb_id.lower())
    except:
        errMsg = "Error while writing index file %s." % tmp_index_file_full_path
        close_file_safely(tmp_index_file_handle, tmp_index_file_full_path, errMsg)
    
    #if old index file does not exist, then write it  
    if os.path.exists(old_index_file_full_path) is False:
        write_index_file_from_file_path_list(old_index_file_full_path, pdb_path_list)
     
    try:
        old_index_file_handle = open(old_index_file_full_path, "r")
    except:
        errMsg = "Cannot open file %s." % old_index_file_full_path
        close_file_safely(tmp_index_file_handle, tmp_index_file_full_path, errMsg)
        errorFct(errMsg)
        exit(1)
            
    try:
        old_index_set = read_index_file(old_index_file_handle)
    except:
        errMsg = "Error while reading index file %s." % old_index_file_full_path
        close_file_safely(tmp_index_file_handle, tmp_index_file_full_path, errMsg)
        close_file_safely(old_index_file_handle, old_index_file_full_path, errMsg)
        exit(1)
        
    close_file_safely(tmp_index_file_handle, tmp_index_file_full_path, "")
    close_file_safely(old_index_file_handle, old_index_file_full_path, "")
  
    #for all entries that have been removed in recent sync
    #this will probably be a rare occurence
    #removal_set = old_index_set.difference(tmp_index_set)
    #if len(removal_set) != 0:
     #   remove_fasta_db_entries_by_header(removal_set)
        
    #set of all entries that got added in recent sync
    added_set = tmp_index_set.difference(old_index_set)
    
    return added_set


#replaces old index file with new one
#called after appending of fasta file finished successfully
# void (str,str)
def replace_old_index_file(tmp_index_file_full_path, old_index_file_full_path):
    
    mv_process = subprocess.Popen(['mv', tmp_index_file_full_path, old_index_file_full_path], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    try:
        out, err = mv_process.communicate()
    except:
        mv_process.kill()
        out, err = mv_process.communicate()
        shell_err_fct(err)
        errMsg = "Communcation with mv_process to replace old index file failed. Killed it."
        errorFct(errMsg)
        exit(1)

#extra. biopython should handle fasta write errors just fine.
#the idea would be to traverse the file backwards with the cursor and delete everything
#between the first occurence of the PDB ID contained within "chains" and the end of the file
# void (file_handle,Bio.SeqRecord)
def repair_recent_fasta_entries(fasta_file_handle, chains):
    
    offset = 0
    while True:
        fasta_file_handle.seek(offset, 2)
        if fasta_file_handle.read(1) == ">":
            pass

#using this function probably does not have a big advantage over using "get_pdb_path_list_recursively()"
#though tying everything to the index file instead of to the file system structure might be useful later
# (str,str) -> [str,str,...]                      
def get_pdb_path_list_by_index_file(index_file_full_path, loc_pdb_db_path):
    
    pdb_path = ""
    pdb_path_list = []

    try:
        index_file_handle = open(index_file_full_path, 'r')
    except:
        errMsg = "Could not open %s." % index_file_full_path
        errorFct(errMsg)
        exit(1)
    try:
        for pdb_id in index_file_handle:
            pdb_id_stripped = pdb_id.rstrip()
            pdb_path = os.path.join(loc_pdb_db_path, pdb_id_stripped[1:3], "pdb"+pdb_id_stripped+".ent.gz")
            pdb_path_list.append(pdb_path)
    except:
        errMsg = "Error while reading %s." % index_file_full_path
        close_file_safely(index_file_handle, index_file_full_path, errMsg)
        errorFct(errMsg)
        exit(1)

    close_file_safely(index_file_handle, index_file_full_path, "")
    
    return pdb_path_list

#if there is an exception during the writing of the FASTA file in create_fasta_file_from_index_file(),
#then this function is called to remove the most recently added entries belonging to the PDB that was read
#while the exception occured.
#starts with the file handle at the end of file, works its way up to the last occurence of most_recent_pdb_id
#and deletes everything from there to end of file
# void (str, str)
def rm_most_recent_pdb_from_fasta_file(fasta_file_full_path, most_recent_pdb_path):
    
    #get PDB ID from file path
    pdb_id_match = re.search(r'(.+\/)*pdb(....)\..*', most_recent_pdb_path)
    pdb_id = pdb_id_match.group(2).upper()

    try:
        fasta_file_handle = open(fasta_file_full_path, 'rb+')
    except:
        errMsg = "Failed to open %s." % fasta_file_full_path
        errorFct(errMsg)
        raise
    try:
        offset = 0
        first_occurence_offset = 0
        while True:
            fasta_file_handle.seek(offset, 2)
            char_byte = fasta_file_handle.read(1)
            #every char in UTF-8 that is only decoded by exactly one byte has a 0 as a first bit
            if int.from_bytes(char_byte, byteorder='little') & (1 << 7):
                errMsg = "FASTA file contains char that is not encoded by exactly 1 byte in UTF-8. Exiting."
                errorFct(errMsg)
                raise UnicodeDecodeError('utf-8',char_byte,offset,offset+1,'Char is encoded using more than 1 byte.')
            if char_byte == ">".encode('utf-8'):
                header = fasta_file_handle.read(4)
                if header == str(pdb_id).encode('utf-8'):
                    first_occurence_offset = offset
                    offset -= 1
                    continue
                else:
                    break
            else:
                offset -= 1
                continue
        #sets the file handle end of line before first occurence of matching pdb header
        fasta_file_handle.seek(first_occurence_offset-1, 2)
        #deletes everything from there to eof
        fasta_file_handle.truncate()
            
    except:
        errMsg = "Error during repairing of %s." % fasta_file_full_path
        close_file_safely(fasta_file_handle, fasta_file_full_path, errMsg)
        errorFct(errMsg)
        raise
        
    close_file_safely(fasta_file_handle, fasta_file_full_path, "")

#reads a fasta file and returns all pdb ids (not chain names) in the headers
# (str) -> set() 
def read_fasta_pdb_ids(fasta_file_full_path):
    
    try:
        fasta_file_handle = open(fasta_file_full_path, 'r')
    except:
        errMsg = "Failed to open %s." % fasta_file_full_path
        errorFct(errMsg)
        raise

    pdb_id_set = set()

    for line in fasta_file_handle:
        if line.startswith(">"):
            pdb_id_set.add(line[1:5])
    
    return pdb_id_set     

#creates fasta file from index file
#void (str, str, str, str, int)
def create_fasta_file_from_index_file(loc_pdb_db_path, fasta_file_full_path, index_file_full_path, biopython_log_full_path, verbosity):
    
    pdb_path_list = get_pdb_path_list_by_index_file(index_file_full_path, loc_pdb_db_path)

    counter = 1
    number_of_pdbs = len(pdb_path_list)

    if verbosity>0:
        print("Writing Fasta file...")
    try:
        fasta_file_handle = open(fasta_file_full_path, 'w')
    except:
        errMsg = "File %s could not be opened for writing." % fasta_file_full_path
        errorFct(errMsg)
        exit(1)
    try:
        biopython_log_file_handle = open(biopython_log_full_path, 'a+')
    except:
        errMsg = "File %s could not be opened for writing." % biopython_log_full_path
        errorFct(errMsg)
        exit(1)
    try:
        for pdb_path in pdb_path_list:
            most_recent_pdb_path = pdb_path
            write_loc_pdb_to_fasta_file(pdb_path, fasta_file_handle, biopython_log_file_handle)
            if verbosity>0:
                msg = "[%s/%s]" % (str(counter), str(number_of_pdbs))
                print(msg, end="\r")
            counter += 1
    except:
        print(traceback.format_exc())
        errMsg1 = "Error while writing to %s." % fasta_file_full_path
        close_file_safely(biopython_log_file_handle, biopython_log_full_path, errMsg1)
        #if an exception is thrown during writing of fasta file this will flush fasta file to memory
        close_file_safely(fasta_file_handle, fasta_file_full_path, errMsg1)
        try:
            dest_backup = backup_fasta_file(fasta_file_full_path)
        except:
            errMsg = "Error while backing up FASTA file. Exiting."
            errorFct(errMsg)
            exit(1)
        #removes most recent pdb chains from fasta file
        try:    
            rm_most_recent_pdb_from_fasta_file(fasta_file_full_path, most_recent_pdb_path)
            fasta_pdb_id_set = read_fasta_pdb_ids(fasta_file_full_path)
            #overwrite index file
            write_index_file_from_set(index_file_full_path, fasta_pdb_id_set)
        except:
            try:
                restore_fasta_file_from_backup(dest_backup, fasta_file_full_path)
            except:
                print(traceback.format_exc())
                errMsg3 = "Error while restoring FASTA file backup %s." % dest_backup
                errorFct(errMsg3)
            finally:
                print(traceback.format_exc())
                errMsg2 = "Error while attempting to repair FASTA file %s." % fasta_file_full_path
                errorFct(errMsg2)
        else:
            remove_fasta_backup(dest_backup, fasta_file_full_path)
        finally:
            errorFct(errMsg1)
            exit(1)
          
    close_file_safely(biopython_log_file_handle, biopython_log_full_path, '')
    close_file_safely(fasta_file_handle, fasta_file_full_path, '')


#creates diamond database from fasta file
#void (str, str, str, int)        
def create_diamond_database(diamond_file_path, diamond_db_path, fasta_file_full_path, verbosity):
    
    diamond_db_full_path = os.path.join(diamond_db_path, "diamond_db")

    #<path_to_diamond_file> makedb --in <path_to_FASTA> --db <db_out_path>
    shell_input = []
    shell_input = [diamond_file_path, "makedb", "--in", fasta_file_full_path, "--db", diamond_db_full_path]
    
    if verbosity>0:
        print("Creating local diamond database...")
    if verbosity>1:
        process = subprocess.Popen(shell_input, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        for line in process.stdout:
            print(line, end='\n')
    else:
        process = subprocess.Popen(shell_input, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
    try:
        out, err = process.communicate()
        shell_err_fct(err)
    except:
        process.kill()
        out, err = process.communicate()
        errMsg = "Communication with \"diamond makedb\" subprocess failed. Killed it."
        shell_err_fct(err)
        errorFct(errMsg)
        exit(1)

#think about when to call this.
#when the rsync process fails, e.g. internet connection is lost,
#broken pdb files might remain that might not be cleaned up by calling rsync again, since
#they're already indexed as present in the directory tree.
#the assumption here is that the pdb files are being written from top to bottom.
#when the writing process fails halfway though the file will be truncated, but might still seem intact.
#thus, in truncated files the trailing "END" statement will be missing. check if its there.
# void (str)        
def rm_incomplete_pdbs(loc_pdb_db_path):
    
    pdb_path_list = get_pdb_path_list_recursively(loc_pdb_db_path)
        
#synchronizes copy of remote PDB archive with local copy
#void (str,str,str,int)  
def sync_pdb_copy(diamond_file_path, loc_pdb_db_path, diamond_db_path, verbosity, backup_freq):
    
    tmp_index_file_full_path = os.path.join(diamond_db_path, "tmp_pdb_db.index")
    old_index_file_full_path = os.path.join(diamond_db_path, "pdb_db.index")
    fasta_file_full_path = os.path.join(diamond_db_path, "pdb_db.fasta")
    biopython_log_full_path = os.path.join(diamond_db_path, "biopython_log.txt")

    create_dir_safely(loc_pdb_db_path)
    create_dir_safely(diamond_db_path)
    #https://www.wwpdb.org/ftp/pdb-ftp-sites
    #rsync -rlpt -v -z --delete --port=33444 rsync.rcsb.org::ftp_data/structures/divided/pdb/ ./pdb
    #there's also DNA/RNA in there, look into it to maybe only sync protein structures
    shell_input = []
    shell_input = ["rsync", "-rlpt", "-v", "-z", "--delete", "--port=33444", "rsync.rcsb.org::ftp_data/structures/divided/pdb/", loc_pdb_db_path]
    
    if verbosity>0:
        print("Synchronizing local PDB database with remote. This may take a while.")
    
    process = subprocess.Popen(shell_input, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)

    if verbosity>0:
        #counter = 1
        for line in process.stdout:
            #REMOVE THIS, testing
            #if counter>1000:
             #   process.kill()
              #  break
            line = line.rstrip()+"                   "
            print(line, end="\r")
            #counter += 1

    try:
        out, err = process.communicate()
    except:
        process.kill()
        #rm_incomplete_pdbs(loc_pdb_db_path)
        out, err = process.communicate()
        errMsg = "Communication with rsync subprocess failed. Killed it."
        shell_err_fct(err)
        errorFct(errMsg)
        exit(1)
                    
    if verbosity>0:
        print("Updating PDB index file...")
    
    added_set = update_pdb_index(loc_pdb_db_path, diamond_db_path, tmp_index_file_full_path, old_index_file_full_path)
    if len(added_set)>0 and os.path.exists(fasta_file_full_path):
        add_fasta_db_entries(added_set, fasta_file_full_path, loc_pdb_db_path, biopython_log_full_path, verbosity, backup_freq)
        replace_old_index_file(tmp_index_file_full_path, old_index_file_full_path)
    elif os.path.exists(fasta_file_full_path) is False:
        create_fasta_file_from_index_file(loc_pdb_db_path, fasta_file_full_path, old_index_file_full_path, biopython_log_full_path, verbosity)
        replace_old_index_file(tmp_index_file_full_path, old_index_file_full_path)
    else:
        if verbosity>0:
            print("No new PDB entries for update.")
       
    create_diamond_database(diamond_file_path, diamond_db_path, fasta_file_full_path, verbosity)

#removes gaps and other symbols not in valid_symbols and writes new fasta file with cleaned sequences
# void (str,str,[str,str,...], int)    
def clean_fasta_file(aln_file_full_path, out_file_full_path, verbosity):
    
    if verbosity>0:
        msg = "Cleaning %s..." % os.path.split(aln_file_full_path)[1]
        print(msg)

    try:
        aln_file_handle = open(aln_file_full_path)
    except:
        errMsg = "Could not open alignment file %s." % aln_file_full_path
        errorFct(errMsg)
        exit(1)
    try:
        records = list(SeqIO.parse(aln_file_handle, "fasta"))
    except:
        errMsg = "Could not read fasta file %s." % aln_file_full_path
        close_file_safely(aln_file_handle, aln_file_full_path, errMsg)
        errorFct(errMsg)
        exit(1)

    close_file_safely(aln_file_handle, aln_file_full_path, "")
    
    counter = 1
    no_of_records = len(records)
    for record in records:
        if verbosity>0:
            msg = "[%s/%s]" % (str(counter), str(no_of_records))
            print(msg, end='\r')
        cleaned_seq = clean_seq(str(record.seq))
        #if this fails use seqrecord constructor
        record.seq = Seq.Seq(cleaned_seq)
        counter += 1

    try:
        SeqIO.write(records, out_file_full_path, "fasta")
    except:
        errMsg = "Error while writing %s." % out_file_full_path
        errorFct(errMsg)
        exit(1)

#clean sequences of all fasta files in given dir. writes new ones to out_path with the same name
# void (str,str,[str,str,...],int)        
def clean_fasta_files_in_dir(path_to_fasta_files, out_path, valid_symbols, verbosity):
    
    for fasta_file_name in os.listdir(path_to_fasta_files):
        fasta_file_full_path = os.path.join(path_to_fasta_files, fasta_file_name)
        out_file_full_path = os.path.join(out_path, os.path.split(fasta_file_full_path)[1])
        clean_fasta_file(fasta_file_full_path, out_file_full_path, verbosity)

def init_parallel_diamond_run(diamond_exe_full_path, diamond_db_file_full_path, query_fasta_file_full_path, out_file_full_path,
                              diamond_block_size, diamond_tmp_dir, diamond_mp_tmp_dir, diamond_threads):

    shell_input = [str(diamond_exe_full_path), 'blastp', '--db', str(diamond_db_file_full_path), '--query', str(query_fasta_file_full_path), '-o', str(out_file_full_path),
                   '--multiprocessing', '--mp-init', '--tmpdir', str(diamond_tmp_dir), '--parallel-tmpdir', str(diamond_mp_tmp_dir)]
    
    blastp_process = subprocess.Popen(shell_input, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
    try:
        out, err = blastp_process.communicate()
    except:
        blastp_process.kill()
        out, err = blastp_process.communicate()
        shell_err_fct(err)
        errMsg = "Failed to communicate with diamond blastp subprocess. Killed it."
        errorFct(errMsg)
        raise

# ./diamond blastp -d <diamond_db_file_full_path> -q <query_fasta_file_full_path> -o <out_path> --iterate
#launches a subprocess for a diamond blastp query without restriction of sequence identity
# void (str,str,str,str,int)        
def diamond_blastp_query(diamond_exe_full_path, diamond_db_file_full_path, query_fasta_file_full_path, out_file_full_path,
                         verbosity, diamond_block_size, diamond_threads, diamond_tmp_dir, diamond_mp_tmp_dir, diamond_mp):

    shell_input = []
    shell_input = [str(diamond_exe_full_path), 'blastp', '--db', str(diamond_db_file_full_path), '--query', str(query_fasta_file_full_path), '-o', str(out_file_full_path),
                   '--log', '--header', 'simple', '--outfmt', '6', 'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
                   'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qlen', '--id', '99.9']

    if int(diamond_threads) > 0:
        shell_input.append('--threads')
        shell_input.append(str(diamond_threads))
    if diamond_tmp_dir != "":
        shell_input.append('--tmpdir')
        shell_input.append(str(diamond_tmp_dir))
    if diamond_mp_tmp_dir != "":
        shell_input.append('--parallel-tmpdir')
        shell_input.append(str(diamond_mp_tmp_dir))
    if diamond_mp > 0:
        shell_input.append('--multiprocessing')
        shell_input.append('--mid-sensitive')
        if int(diamond_threads) <= 0:
            shell_input.append('--threads')
            try:
                diamond_threads = os.environ['SLURM_CPUS_PER_TASK']
            except:
                warn_msg = "Warning: Could not get number of CPUs per task from environment variable SLURM_CPUS_PER_TASK."
                errorFct(warn_msg)
            print("SLURM_CPUS_PER_TASK="+str(diamond_threads))
            shell_input.append(str(diamond_threads))
        if verbosity>0:
            msg = "Initializing parallel run for file %s..." % query_fasta_file_full_path
            print(msg)
        try:
            init_parallel_diamond_run(diamond_exe_full_path, diamond_db_file_full_path, query_fasta_file_full_path, out_file_full_path,
                                      diamond_block_size, diamond_tmp_dir, diamond_mp_tmp_dir, diamond_threads)
        except:
            err_msg = "Error while initializing parallel diamond run for file %s!" % query_fasta_file_full_path
            errorFct(err_msg)
            exit(1)
    else:
        shell_input.append('--iterate')
        if int(diamond_threads) <= 0:
            shell_input.append('--threads')
            try:
                diamond_threads = os.environ['SLURM_CPUS_PER_TASK']
            except:
                warn_msg = "Warning: Could not get number of CPUs per task from environment variable SLURM_CPUS_PER_TASK."
                errorFct(warn_msg)
            print("SLURM_CPUS_PER_TASK="+str(diamond_threads))
            shell_input.append(str(diamond_threads))
    if diamond_block_size != None:
        shell_input.append('--block-size')
        shell_input.append(str(diamond_block_size))
    if verbosity>0:
        msg = "Performing blastp query for file %s." % os.path.split(query_fasta_file_full_path)[1]
        print(msg)
    
    blastp_process = subprocess.Popen(shell_input, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
    try:
        out, err = blastp_process.communicate()
    except:
        blastp_process.kill()
        out, err = blastp_process.communicate()
        shell_err_fct(err)
        errMsg = "Failed to communicate with diamond blastp subprocess. Killed it."
        errorFct(errMsg)
        exit(1)

#starts a diamond blastp query, but with an entire dir of fasta files, writes them with the same file name to out_path
#void (str,str,str,str,int)        
def batch_blastp_query(diamond_exe_full_path, diamond_db_file_full_path, query_fasta_file_dir_path, out_path, verbosity):
    
    create_dir_safely(out_path)

    for query_fasta_file_name in os.listdir(query_fasta_file_dir_path):
        query_fasta_file_full_path = os.path.join(query_fasta_file_dir_path, query_fasta_file_name)
        if os.path.isfile(query_fasta_file_full_path) is False:
            continue
        out_file_full_path = os.path.join(out_path, os.path.splitext(query_fasta_file_name)[0]+".tsv")
        diamond_blastp_query(diamond_exe_full_path, diamond_db_file_full_path, query_fasta_file_full_path, out_file_full_path, verbosity)

#checks if TSV file generated by diamond is properly tab-separated
# (str) -> ([int,int,...])        
def check_TSV_integrity(TSV_file_full_path):

    broken_line_list = []
    with open(TSV_file_full_path, 'r') as TSV_file_handle:
        header_match_number = 0
        line_count = 1
        for header in TSV_file_handle:
            header_match = re.findall(r'([^\t\n]+)(\t|\n)', header)
            header_match_number = len(header_match)
            break
        for line in TSV_file_handle:
            line_match = re.findall(r'([^\t\n]+)(\t|\n)', line)
            line_match_number = len(line_match)
            if line_match_number != header_match_number:
                errMsg = "Warning: Inconsistent data in file %s, line %s:\n%s" % (TSV_file_full_path, str(line_count), line_match)
                errorFct(errMsg)
            broken_line_list.append(line_count)
            line_count += 1
            
    return broken_line_list

#this function defines what an exact match is
#currently it's 0 mismatches and alignment over full length of query
# (str) -> (pandas.DataFrame)
def filter_exact_matches(TSV_file_full_path):
    
    df = pandas.read_csv(TSV_file_full_path, sep='\t')
    #new_df = pandas.DataFrame(columns=['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore','qlen'])
    row_list = []
    for index, row in df.iterrows():
        #exact query match definition. no mismatches and no additional gaps and qlen=send-sstart+1 =>(?) pident = 100
        if row.mismatch == 0 and row.gapopen == 0:
            aln_len = row.send - row.sstart + 1
            if row.qlen == aln_len:
                row_list.append(row.values)
    new_df = pandas.DataFrame(row_list, columns=df.columns)

    return new_df

#definition of exact match needs to be applied here
# (str) -> { {},{},... },{ {},{},... }
def get_exact_matches(TSV_file_full_path, aln_entry_list):
    
    #optional check, if I trust diamond to generate good data
    #broken_line_list = check_TSV_integrity(TSV_file_full_path)
    #print(broken_line_list)
    #if len(broken_line_list) != 0:
     #   attempt_TSV_repair()
    #apply filter according to exact match definition to get at most one chain id for each query header
    #we also need to make sure that if the input sequences are the same that we always choose the SAME 
    #chain as exact query match
    exact_matches_df = filter_exact_matches(TSV_file_full_path)
    #write to file
    TSV_path, TSV_name = os.path.split(TSV_file_full_path)
    name, ext = os.path.splitext(TSV_name)
    exact_matches_df.to_csv(os.path.join(TSV_path, name+"_exact_matches.tsv"), sep="\t")
    
    #keys = FASTA_header; values = [row,row,...] 
    exact_match_dict = {}
    #initializing exact_match_dict with empty lists
    for aln_entry in aln_entry_list:
        exact_match_dict[str(aln_entry[1])] = []
    #filling lists
    for index, row in exact_matches_df.iterrows():
        exact_match_dict[str(row.qseqid)].append(row)    

    return exact_match_dict
 
def remove_pdb_entry_from_fasta(pdb_id):
    pass


#copies PDB files from local PDB database to <job_dir>/DALI/PDB_lib
#mainly implemented for consistency, no real advantage
# void ([str,str,...],str,str)
def cp_PDB_files_to_job_dir(chain_list, pdb_db_path, pdb_out_path):
    
    for chain in chain_list:
        try:
            pdb_id = chain[0:4].lower()
        except IndexError:
            errMsg = "Chain name has less than 4 chars."
            errorFct(errMsg)
            raise
        pdb_file_name = "pdb"+str(pdb_id)+".ent.gz"
        src_pdb_full_path = os.path.join(pdb_db_path, pdb_id[1:3], pdb_file_name)
        dest_pdb_full_path = os.path.join(pdb_out_path, pdb_file_name)
        try:
            shutil.copy(src_pdb_full_path, dest_pdb_full_path)
        except:
            errMsg = "Error while copying file %s to %s." % (src_pdb_full_path, dest_pdb_full_path)
            errorFct(errMsg)
            remove_pdb_entry_from_fasta(pdb_id)

#assigns numbers (temporary ids) to headers and sequences in alignment file, return dict
#this function is very similar to import_seq_list_from_fasta_aln(), but i dont want to change the latter at the moment
# (str) -> { (str,int),(str,int),... }        
def enumerate_aln_headers(aln_file_full_path):
    
    aln_dict = {}
    
    with open(aln_file_full_path, 'r') as aln_file_handle:
        tmp_id = 0
        try:
            alignment = AlignIO.read(aln_file_handle, "fasta")
        except:
            errMsg = "Failed to read alignment file %s." % aln_file_handle
            errorFct(errMsg)
            raise
        for record in alignment:
            aln_dict[record.id] = (str(record.seq), tmp_id)
            tmp_id += 1
            
    return aln_dict
                
def diamond_calc_job_SPS(job_list_of_dict_lists, exact_match_dict, aln_dict, job_title, verbosity):
    
    if verbosity>0:
        msg = "Calculating SPS for job %s." % job_title
        print(msg)

    ref_DALI_set = set()
    MSA_set = set()
    
    #each list in job_list_of_dict_lists represents a sub job, i.e. a dir of format xxxxX_xxxxX
    #for each sub job in job
    for sub_job_dict_list in job_list_of_dict_lists:
        #for each alignment in sub job
        for dali_aln in sub_job_dict_list:
            DALI_aln_index_set = set()
            query_range = set()
            sbjct_range = set()
            DALI_aln_index_set.clear()
            query_range.clear()
            sbjct_range.clear()
            try:
                query_chain_id = dali_aln["query"][0:4]+":"+dali_aln["query"][4]
                sbjct_chain_id = dali_aln["sbjct"][0:4]+":"+dali_aln["sbjct"][4]
            except IndexError:
                errMsg = "Query or subject name in DALI alignment %s contain less than 5 characters."
                errorFct(errMsg)
                exit(1)
            DALI_query_seq = dali_aln["DALI_query_seq"]
            DALI_sbjct_seq = dali_aln["DALI_sbjct_seq"]
            orig_query_header = exact_match_dict[query_chain_id]["qseqid"]
            orig_sbjct_header = exact_match_dict[sbjct_chain_id]["qseqid"]
            tmp_query_id = aln_dict[orig_query_header][1]
            tmp_sbjct_id = aln_dict[orig_sbjct_header][1]
            orig_query_seq = aln_dict[orig_query_header][0]
            orig_sbjct_seq = aln_dict[orig_sbjct_header][0]
            #calculate interval intersection
            query_sstart = exact_match_dict[query_chain_id]["sstart"]
            query_send = exact_match_dict[query_chain_id]["send"]
            query_range = set(range(query_sstart,query_send+1))
            sbjct_sstart = exact_match_dict[sbjct_chain_id]["sstart"]
            sbjct_send = exact_match_dict[sbjct_chain_id]["send"]
            sbjct_range = set(range(sbjct_sstart,sbjct_send+1))
            overlap_interval = query_range.intersection(sbjct_range)
            DALI_tmp_aln_index_set = calc_ref_SP_set((tmp_query_id, tmp_sbjct_id), (DALI_query_seq, DALI_sbjct_seq))
            #add only those elements to ref_DALI_set that overlap in the original chain, since non-overlapping regions
            #in the original MSA cannot be measured
            for element in DALI_tmp_aln_index_set:
                if element[0][1] in overlap_interval and element[1][1] in overlap_interval:
                    DALI_aln_index_set.add(element)
            ref_DALI_set = ref_DALI_set.union(DALI_aln_index_set.copy())
            orig_aln_index_set = calc_MSA_SP_set((tmp_query_id, tmp_sbjct_id), (orig_query_seq, orig_sbjct_seq))
            MSA_set = MSA_set.union(orig_aln_index_set.copy())
    print("DALI_set: ")
    print(ref_DALI_set)
    print("MSA_set: ")
    print(MSA_set)
    intersection_set = ref_DALI_set.intersection(MSA_set)
    if len(ref_DALI_set) == 0:
        job_SPS = 0
    else:
        job_SPS = len(intersection_set)/len(ref_DALI_set)
    
    return job_SPS

#takes alignment and calculates SPS score from pairwise structural alignments
#        
def calc_aln_score_with_diamond(aln_file_full_path, diamond_exe_full_path, diamond_db_file_full_path, pdb_db_path, out_path,
                                dali_path, sample_size, job_title, valid_symbols, verbosity):
    
    job_out_path = os.path.join(out_path, job_title)
    create_dir_safely(job_out_path)
    
    job_data_out_path = os.path.join(job_out_path, "DATA")
    create_dir_safely(job_data_out_path)
    
    aln_file_name = os.path.splitext(os.path.split(aln_file_full_path)[1])[0]
    TSV_out_path = os.path.join(job_data_out_path, "TSV")
    TSV_file_full_path = os.path.join(TSV_out_path, aln_file_name+".tsv")
    create_dir_safely(TSV_out_path)

    clean_fasta_out_path = os.path.join(job_data_out_path, "CLEAN_ALN")
    create_dir_safely(clean_fasta_out_path)
    clean_fasta_full_path = os.path.join(clean_fasta_out_path, aln_file_name)
    clean_fasta_file(aln_file_full_path, clean_fasta_full_path, valid_symbols, verbosity)

    #it's a bit of a waste to read the alignment file yet again and not import this information while the
    #file is read the first time, but it's definitely clearer this way
    #aln_dict has original fasta headers as keys, check for duplicate keys
    aln_dict = enumerate_aln_headers(aln_file_full_path)

    diamond_blastp_query(diamond_exe_full_path, diamond_db_file_full_path, clean_fasta_full_path, TSV_file_full_path, verbosity)
    #exact_match_dict has chains as keys, i.e. XXXX:X
    #be extra careful and watch out for leading spaces and wrong or missing \t separations
    #check for duplicate keys
    
    #this function was changed!
    exact_match_dict = get_exact_matches(TSV_file_full_path)
    
    if sample_size != 0:
        #randomize key popping and write function
        key_list = list(exact_match_dict.keys())
        key_list_len = len(key_list)
        no_to_remove = key_list_len - sample_size
        if no_to_remove < 0:
            no_to_remove = 0
        i = 1
        while i <= no_to_remove:
            #generate random number between 0 and len(key_list)
            random_int = random.randint(0, len(key_list)-1)
            key = key_list[random_int]
            try:
                exact_match_dict.pop(key)
            except KeyError:
                continue
            i += 1
    
    pdb_out_path = os.path.join(job_data_out_path, "PDB_lib")
    create_dir_safely(pdb_out_path)

    cp_PDB_files_to_job_dir(exact_match_dict.keys(), pdb_db_path, pdb_out_path)

    dat_out_path = os.path.join(job_data_out_path, "DAT_lib")
    create_dir_safely(dat_out_path)

    DALI_import_PDBs(dali_path, pdb_out_path, dat_out_path, verbosity)
    
    #use xxxxX format for all dali related things, instead of xxxx:X
    DALI_all_vs_all_query(dali_path, pdb_out_path, out_path, dat_out_path, job_title, verbosity, chain_list=exact_match_dict.keys())

    job_list_of_dict_lists = import_aln_files_in_job(job_title, out_path)

    job_SPS = diamond_calc_job_SPS(job_list_of_dict_lists, exact_match_dict, aln_dict, job_title, verbosity)
    
    print(job_SPS)

#calculation rework
#__________________________________________________________________________________________________________________________________________________
def get_mantissa_and_exponent_from_number(scientific_number):
    scientific_regex_str = "([0-9.,-]+)([eE]|(\*[0-9]+\^))([0-9-]+)"
    match = re.search(scientific_regex_str, str(scientific_number))
    try:
        mant = float(match.group(1))
    except:
        if match.group(1) is not None:
            err_msg = "Could not get mantissa from match group 1 %s." % match.group(1)
            errorFct(err_msg)
        raise
    try:
        exp = int(match.group(4))
    except:
        if match.group(4) is not None:
            err_msg = "Could not get exponent from match group 4 %s." % match.group(4)
            errorFct(err_msg)
        raise

    return mant, exp

#this function serves to detect and resolve duplicate values in dictionary for any one key
#right now it takes FASTA headers as keys and checks if the list in values is of length <= 1
#its important to note that ambiguous rows are rows that are objectively equivalent in respect to
#exact match selection. this means that there is no direct advantage in choosing one hit over the other
#BUT we want to make sure that the proteins, which the hits are embedded in are the same among all measurements
#because this impacts structural alignment generation and thus the reference set and in consequence also the score
def resolve_ambivalences(exact_match_dict, aln_file_full_path):

    for qseqid in exact_match_dict.keys():
        if len(exact_match_dict[qseqid]) > 1:
            cur_best_entry = exact_match_dict[qseqid][0]
            for row in exact_match_dict[qseqid][1:]:
                #try to select entry by evalue
                try:
                    mant_best, exp_best = get_mantissa_and_exponent_from_number(cur_best_entry.evalue)
                except:
                    mant_best = float('inf')
                    exp_best = float('inf')
                try:
                    mant_row, exp_row = get_mantissa_and_exponent_from_number(row.evalue)
                except:
                    mant_row = float('inf')
                    exp_row =  float('inf')
                if exp_row < exp_best:
                    cur_best_entry = row
                    continue
                elif exp_row == exp_best and mant_row < mant_best:
                    cur_best_entry = row
                    continue
                #try to select entry by bitscore
                elif int(row.bitscore) > int(cur_best_entry.bitscore):
                    cur_best_entry = row
                    continue
                elif int(row.bitscore) == int(cur_best_entry.bitscore):
                    #if bitscore is the same, try lower subject match start
                    if int(row.sstart) < int(cur_best_entry.sstart):
                        cur_best_entry = row
                        continue
                    elif int(row.sstart) == int(cur_best_entry.sstart):
                        #if lower match start is the same, we must have different chain IDs
                        for i in range(0,6):
                            if ord(str(row.sseqid)[i]) < ord(str(cur_best_entry.sseqid)[i]):
                                cur_best_entry = row
                                break
                            elif ord(str(row.sseqid)[i]) > ord(str(cur_best_entry.sseqid)[i]):
                                break
                        #if we have the same chain IDs, we'll have the same query match restriction
                        #and the choice is irrelevant. This SHOULD NOT happen usually.
                        else:
                            err_msg = """Warning: Identical exact matches suspected for headers %s and 
                                       %s in %s.""" % (str(cur_best_entry.qseqid), str(row.qseqid), aln_file_full_path)
                            errorFct(err_msg)
            #we'll use a list for consistency
            exact_match_dict[qseqid] = []
            exact_match_dict[qseqid].append(cur_best_entry)

    return exact_match_dict

def create_MSA_ex_match_df(aln_entry_list, unique_exact_match_dict):

    MSA_df = pandas.DataFrame(columns=['fasta_entry_index','fasta_entry_header','ungapped_seq','gapped_seq',
                                       'chain_id','sstart','send','evalue','bitscore','qlen'])
    row_list = []
    #aln_entry = [sequence, header, entry_number, ungapped_sequence]
    for aln_entry in aln_entry_list:
        row = []
        qseqid = aln_entry[1]
        if len(unique_exact_match_dict[qseqid]) == 1:
            row.append(aln_entry[2])
            row.append(qseqid)
            row.append(aln_entry[3])
            row.append(aln_entry[0])
            exact_match_tuple = unique_exact_match_dict[qseqid][0]
            row.append(exact_match_tuple.sseqid)
            row.append(exact_match_tuple.sstart)
            row.append(exact_match_tuple.send)
            row.append(exact_match_tuple.evalue)
            row.append(exact_match_tuple.bitscore)
            row.append(exact_match_tuple.qlen)

            row_list.append(row.copy())
    
    MSA_df = pandas.DataFrame(row_list, columns=['fasta_entry_index','fasta_entry_header','ungapped_seq',
                              'gapped_seq','chain_id','sstart','send','evalue','bitscore','qlen'])
    return MSA_df

#writes index files for each input alignment that tie together diamond blastp query data
# with the headers and the aligned sequences in the original input files
# void ([str,str,...], str, str, str, int)
def create_aln_index_files(aln_file_full_path_list, job_path, diamond_exe_full_path, diamond_db_file_full_path,
                           verbosity, diamond_block_size, diamond_threads, diamond_tmp_dir, diamond_mp_tmp_dir, diamond_mp):

    DATA_path = os.path.join(job_path, "DATA")
    create_dir_safely(DATA_path)

    for aln_file_full_path in aln_file_full_path_list:

        path, aln_file_name_ext = os.path.split(aln_file_full_path)
        aln_file_name, ext = os.path.splitext(aln_file_name_ext)

        #sub-subdir of jobdir for each alignment, contains all temporary or permanent files important for measurment
        #of each individual MSA file
        aln_misc_path = os.path.join(DATA_path, aln_file_name)
        create_dir_safely(aln_misc_path)
        #aln_entry_list is a list of lists each of length 4, containing [sequence, header, entry_number, ungapped_sequence]
        aln_entry_list = import_seq_list_from_fasta_aln(aln_file_full_path)
        
        #VERY CAREFUL HERE, DOUBLE CHECK, TRIPLE CHECK, IF IT MAKES SENSE TO DO IT LIKE THIS
        #IF WE FORGOT A SINGLE LETTER, THIS RUINS ALL RESULTS
        #HOW DO WE HANDLE UNCOMMON AAs OR 'X' AND THE LIKE?
        #___________________________________________________________________________________________
        #path for fasta output
        clean_fasta_file_full_path = os.path.join(aln_misc_path, aln_file_name+"_cleaned.fasta")
        _ = clean_fasta_file(aln_file_full_path, clean_fasta_file_full_path, verbosity)
        #-------------------------------------------------------------------------------------------

        #path for TSV output, return file for db query
        TSV_file_full_path = os.path.join(aln_misc_path, aln_file_name+"_diamond.tsv")
        diamond_blastp_query(diamond_exe_full_path, diamond_db_file_full_path, clean_fasta_file_full_path, TSV_file_full_path,
                             verbosity, diamond_block_size, diamond_threads, diamond_tmp_dir, diamond_mp_tmp_dir, diamond_mp)

        #get exact query matches. this is shaky, because we dont know how diamond processes large
        # FASTA files internally and if this has an impact on what results we are presented with
        #however, diamond excels at processing large FASTA files and it would be wasteful to feed in 
        #the individual FASTA entries one by one (it might be necessary though, if exact query match
        #output turns out to differ between different alignments with the same underlying sequences)
        exact_match_dict = get_exact_matches(TSV_file_full_path, aln_entry_list)
        unique_exact_match_dict = resolve_ambivalences(exact_match_dict, aln_file_full_path)

        #create dataframe for summary of MSA data
        MSA_df = create_MSA_ex_match_df(aln_entry_list, unique_exact_match_dict)
        #write MSA_df to file
        MSA_df.to_csv(os.path.join(aln_misc_path, aln_file_name+"_exact_match.index"), sep="\t")

def read_MSA_index_files_in_job(aln_file_full_path_list, job_path):

    #extract aln_file_names from paths and build list of full paths to index files
    DATA_path = os.path.join(job_path, "DATA")
    aln_index_full_path_list = []
    for aln_file_full_path in aln_file_full_path_list:
        aln_path, aln_name_ext = os.path.split(aln_file_full_path)
        aln_file_name, ext = os.path.splitext(aln_name_ext)
        aln_index_full_path = os.path.join(DATA_path, aln_file_name, aln_file_name+"_exact_match.index")
        aln_index_full_path_list.append(aln_index_full_path)

    #load index files into list of DataFrames
    MSA_df_dict = {}
    for aln_index_full_path in aln_index_full_path_list:
        MSA_df = pandas.read_csv(aln_index_full_path, sep="\t")
        MSA_df_dict[aln_index_full_path] = MSA_df.copy()

    return MSA_df_dict

#this function is given a set of ungapped sequences that are shared among all input MSAs
#it is also given a dict of dictionaries. latter contain ungapped_sequences of ALL queries that have 
#exact query matches (the ungapped sequences) as keys and the corresponding rows in the aln_index as values as a list.
# each key of the higher level dict represents a different input MSA
# (set(), [{},{},...]) -> set()
def verify_common_seq_data(common_seq_set, seq_prim_key_dict_dict):

    exclude_set = set()
    #following three sets should sufficiently define the exact query match restriction
    chain_id_set = set()
    sstart_set = set()
    send_set = set()
    #this is just an extra check, it should be the same, but it needn't be
    bitscore_set = set()
    for common_seq in common_seq_set:
        chain_id_set.clear()
        sstart_set.clear()
        send_set.clear()
        bitscore_set.clear()
        #the next loop adds the respective column values to sets for all aln_index files
        #if a set of the relevalent data has more than one element, then it means that
        #there had to be a difference in relevant data (i.e. chain_id, sstart, send)
        #somewhere among all the aln_index files for rows that share the same unngapped sequence
        #this means that there would be NOT at most one P_{ref,C} for a given C of ungapped sequences
        #contradicting Lemma 2 
        prev_aln_index_full_path = ""
        for aln_index_full_path in seq_prim_key_dict_dict.keys():
            for i in range(0, len(seq_prim_key_dict_dict[aln_index_full_path][common_seq])):
                chain_id_set.add(str(seq_prim_key_dict_dict[aln_index_full_path][common_seq][i].chain_id))
                sstart_set.add(str(seq_prim_key_dict_dict[aln_index_full_path][common_seq][i].sstart))
                send_set.add(str(seq_prim_key_dict_dict[aln_index_full_path][common_seq][i].send))
                bitscore_set.add(str(float(seq_prim_key_dict_dict[aln_index_full_path][common_seq][i].bitscore)))
                if len(chain_id_set)>1 or len(sstart_set)>1 or len(send_set)>1 or len(bitscore_set)>1:
                    exclude_set.add(common_seq)
                    warn_msg = """Warning: Inconsistency in common sequence set between %s and %s. Excluding sequence.\n
                                  Chain_ids: %s\n Sstart: %s\n Send: %s\n Bitscore: %s""" % (prev_aln_index_full_path,aln_index_full_path,
                                                                                             str(chain_id_set),str(sstart_set),str(send_set),str(bitscore_set))
                    errorFct(warn_msg)
            prev_aln_index_full_path = aln_index_full_path
    red_common_seq_set = common_seq_set.difference(exclude_set)
    
    return red_common_seq_set


def generate_common_seq_set(MSA_df_dict, verbosity):

    #load cleaned exact query match sequences into sets
    seq_set_list = []
    seq_prim_key_dict_dict = {}
    for aln_index_full_path in MSA_df_dict.keys():
        seq_set = set()
        seq_prim_key_dict = {}
        seq_prim_key_dict.clear()
        for index, row in MSA_df_dict[aln_index_full_path].iterrows():
            upper_ungapped_seq = str(row.ungapped_seq).upper()
            seq_set.add(upper_ungapped_seq)
            #initialize seq_prim_key_dict with empty list.
            #these lists will contain every row, where the ungapped sequences
            #are the same (yes, this actually happens for whatever reason)
            if upper_ungapped_seq not in seq_prim_key_dict.keys():
                seq_prim_key_dict[upper_ungapped_seq] = []
            seq_prim_key_dict[upper_ungapped_seq].append(row)
        seq_set_list.append(seq_set.copy())
        seq_prim_key_dict_dict[aln_index_full_path] = seq_prim_key_dict.copy()

    #splat seq_set_list and intersect every set in it
    common_seq_set = set.intersection(*seq_set_list)

    red_common_seq_set = verify_common_seq_data(common_seq_set, seq_prim_key_dict_dict)

    if verbosity>0:
        msg = "[%s/%s] common sequences remaining after verification." % (str(len(common_seq_set)),str(len(red_common_seq_set)))
        print(msg)

    return red_common_seq_set, seq_prim_key_dict_dict

#create reduced MSA dataframe and write to file. return list of full paths to reduced MSA index files.
#i'm personally dissatisfied with this function, but I had to work with duplicate underlying sequences somehow
#thus this function is a bit of a frankenstein's creation
#rewrite at appropriate moment in time
def write_red_aln_index_files(common_seq_set, seq_prim_key_dict_dict, unique_seq_only):

    red_MSA_index_full_path_list = []
    #assign ID to each common seq, e.g. convert set to list.
    common_seq_list = list(common_seq_set)
    #initialize total row dict of dicts for each alignment file
    total_row_dict = {}
    total_row_dict.clear()
    for aln_index_full_path in seq_prim_key_dict_dict.keys():
        total_row_dict[aln_index_full_path] = {}
        total_row_dict[aln_index_full_path].clear()
    
    while_com_seq_list_len = len(common_seq_list)
    common_seq_ID = 0
    canon_ind = 0
    while common_seq_ID < while_com_seq_list_len:
        rows_per_com_seq_dict = {}
        rows_per_com_seq_dict.clear()
        for aln_index_full_path in seq_prim_key_dict_dict.keys():
            #list of rows where ungapped sequences are the same per alignment file
            rows_per_com_seq_dict[aln_index_full_path] = seq_prim_key_dict_dict[aln_index_full_path][common_seq_list[common_seq_ID]]
        #making sure that for each alignment to be compared we have the same amount
        # of aligned sequences for each common sequence
        len_dict = {}
        len_dict.clear()
        for aln_index_full_path in seq_prim_key_dict_dict.keys():
            len_dict[aln_index_full_path] = len(rows_per_com_seq_dict[aln_index_full_path])
        #determine minimum of values in len_dict to restrict number of aligned sequences to this int
        #if unique_seq_only is true, then take only one aligned sequence per common sequence
        if unique_seq_only:
            min_value = 1
        else:
            min_value = min(list(len_dict.values()))
        #RANDOM ELEMENTS! This is controversial in my mind. Ultimately, this only affects the MSA index sets
        # NOT the reference sets, because the underlying sequences are EXACTLY the same. This means, that if 
        # DIAMOND and the subsequent processing of data is deterministic, we should only see a difference in the 
        # reduced index files in the columns "fasta_header", "fasta_index" and "gapped_sequence". 
        # Columns "sstart", "send" and "ungapped_seq" should be the same, which means that our
        # exact query match restriction should remain constant, no matter the choice of the next rows
        # implement check to be extra sure
        new_rows_per_com_seq_dict = {}
        new_rows_per_com_seq_dict.clear()
        #choosing subset for each list of rows with gapped_seq, where ungapped sequences are the same 
        for aln_index_full_path in rows_per_com_seq_dict.keys():
            rand_ind_list = random.sample(range(0,len(rows_per_com_seq_dict[aln_index_full_path])), min_value)
            new_rows_per_com_seq_dict[aln_index_full_path] = [rows_per_com_seq_dict[aln_index_full_path][index] for index in rand_ind_list]

        #assign canonical index
        for aln_index_full_path in new_rows_per_com_seq_dict.keys():
            for row in new_rows_per_com_seq_dict[aln_index_full_path]:
                total_row_dict[aln_index_full_path][canon_ind] = row
                canon_ind += 1
            #new_rows_per_com_seq_dict[key] for all aln_index_full_paths in keys contains exactly min_value rows
            #thus we can subtract min_value from canon_ind to keep the canonical index the same between all red_index_files
            canon_ind -= min_value
        canon_ind += min_value
        common_seq_ID += 1

    #create dfs
    for aln_index_full_path in total_row_dict.keys():
        red_MSA_df = pandas.DataFrame.from_dict(total_row_dict[aln_index_full_path], orient='index')
        #construct path
        path, name_ext = os.path.split(aln_index_full_path)
        red_name_ext = "red_"+name_ext
        red_MSA_index_full_path = os.path.join(path, red_name_ext)
        #write to file
        red_MSA_df.to_csv(red_MSA_index_full_path, sep='\t')
        red_MSA_index_full_path_list.append(red_MSA_index_full_path)

    return red_MSA_index_full_path_list

#elaborate test function to test integrity of reduced index files
def test_red_index_file(aln_dir, out_path, job_name, verbosity, differing_sequences, unique_seq_only):

    #get all alignment files in dir
    aln_dir_list = os.listdir(aln_dir)
    aln_file_full_path_list = []
    for thing in aln_dir_list:
        if os.path.isfile(os.path.join(aln_dir, thing)):
            aln_file_full_path_list.append(os.path.join(aln_dir, thing))

    job_path = os.path.join(out_path, job_name)
    create_dir_safely(job_path)

    #use index files to determine set of common underlying sequences
    #MSA_df_dict will serve as a central information hub
    MSA_df_dict = read_MSA_index_files_in_job(aln_file_full_path_list, job_path)
    #common_seq_set is our C. seq_prim_key_dict_dict is utility structure helping to generate
    #the reduced alignment index files
    common_seq_set, seq_prim_key_dict_dict = generate_common_seq_set(MSA_df_dict, verbosity)
    if differing_sequences.issubset(common_seq_set) and len(differing_sequences) != 0:
        warn_msg = ("Warning: Integrity test failed for reduced index files."
                    "Common sequences %s are not present, although they should be.") % (str(differing_sequences))
        errorFct(warn_msg)
    red_MSA_index_full_path_list = write_red_aln_index_files(common_seq_set, seq_prim_key_dict_dict, unique_seq_only)
    #test contents of reduced index files
    red_MSA_df_dict = import_red_aln_index_files(red_MSA_index_full_path_list)

    test_dict = {}
    test_dict.clear()
    i = 0
    while True:
        try:
            test_dict[i] = []
            for aln_index_full_path in red_MSA_df_dict.keys():
                test_dict[i].append(red_MSA_df_dict[aln_index_full_path].iloc[i])
            i += 1
        except IndexError:
            max_i = i-1
            break
    #are the rows almost the same for each canonical index?
    seq_set_red = set()
    for i in range(0,max_i+1):
        last_row = test_dict[i][0]
        for row in test_dict[i]:
            if row.ungapped_seq.upper() != last_row.ungapped_seq.upper():
                msg = "Warning: Rows with index %s have differing ungapped sequences." % (str(i))
                print(row.ungapped_seq)
                print(last_row.ungapped_seq)
                errorFct(msg)
            if row.sstart != last_row.sstart:
                msg = "Warning: Rows with index %s have differing sstart." % (str(i))
                print(row.sstart)
                print(last_row.sstart)
                errorFct(msg)
            if row.send != last_row.send:
                msg = "Warning: Rows with index %s have differing send." % (str(i))
                print(row.send)
                print(last_row.send)
                errorFct(msg)
            if row.evalue != last_row.evalue:
                msg = "Warning: Rows with index %s have differing evalue." % (str(i))
                print(row.evalue)
                print(last_row.evalue)
                errorFct(msg)
            if row.bitscore != last_row.bitscore:
                msg = "Warning: Rows with index %s have differing bitscore." % (str(i))
                print(row.bitscore)
                print(last_row.bitscore)
                errorFct(msg)

            seq_set_red.add(row.ungapped_seq.upper())

    seq_set_index_file = set()
    test_dict_2 = {}
    test_dict_2.clear()
    i = 0
    while True:
        try:
            test_dict_2[i] = []
            for aln_index_full_path in MSA_df_dict.keys():
                test_dict_2[i].append(MSA_df_dict[aln_index_full_path].iloc[i])
            i += 1
        except IndexError:
            max_i = i-1
            break

    for i in range(0, max_i+1):
        for row in test_dict_2[i]:
            seq_set_index_file.add(row.ungapped_seq.upper())

    differing_sequences = seq_set_red.symmetric_difference(seq_set_index_file)

    return differing_sequences

def random_aln_seq(clean_seq, base_aln_len):
    clean_seq_len = len(clean_seq)
    cur_seq_len = clean_seq_len
    gaps_to_insert = base_aln_len - clean_seq_len
    if gaps_to_insert <= 0:
        return '0'
    while gaps_to_insert > 0:
        #cur_seq_len+1 places to insert gap 
        rand_int = random.randint(0,cur_seq_len)
        if rand_int != 0 and rand_int != cur_seq_len:
            clean_seq = clean_seq[:rand_int]+'-'+clean_seq[rand_int:]
            cur_seq_len += 1
            gaps_to_insert -= 1
            continue
        elif rand_int == 0:
            clean_seq = '-'+clean_seq
            cur_seq_len += 1
            gaps_to_insert -= 1
            continue
        elif rand_int == cur_seq_len:
            clean_seq = clean_seq+'-'
            cur_seq_len += 1
            gaps_to_insert -= 1
            continue

    return clean_seq

def write_baseline_aln_file(aln_file_full_path_list, job_path, verbosity):

    cleaned_seq_set = set()
    cleaned_seq_set.clear()
    aln_len_set = set()
    aln_len_set.clear()
    for aln_file_full_path in aln_file_full_path_list:
        if verbosity>0:
            msg = "Cleaning file %s for baseline." % (aln_file_full_path)
            print(msg, end="\n")
        #why is this here?
        #entry is [sequence, header, entry_number, ungapped_sequence]
        #aln_entry_list = import_seq_list_from_fasta_aln(aln_file_full_path)
        #cleaning sequences
        try:
            aln_file_handle = open(aln_file_full_path)
        except:
            errMsg = "Could not open alignment file %s." % aln_file_full_path
            errorFct(errMsg)
            exit(1)
        try:
            records = AlignIO.read(aln_file_handle, "fasta")
        except:
            errMsg = "Could not read fasta file %s." % aln_file_full_path
            close_file_safely(aln_file_handle, aln_file_full_path, errMsg)
            errorFct(errMsg)
            exit(1)

        close_file_safely(aln_file_handle, aln_file_full_path, "")

        counter = 0
        for record in records:
            counter += 1
            if verbosity>0:
                msg = "[%s]" % (counter)
                print(msg, end="\r")
            if counter<2:
                aln_len_set.add(len(str(record.seq)))
            cleaned_seq_set.add(clean_seq(str(record.seq).upper()))
    
    base_aln_len_sorted_list = iter(sorted(list(aln_len_set)))
    cur_base_aln_len = next(base_aln_len_sorted_list)
    #create SeqRecord list
    seq_record_list = []
    counter = 1
    if verbosity>0:
        msg = "Generating random alignment file."
        print(msg, end="\n")
    cleaned_seq_list = list(cleaned_seq_set)
    cleaned_seq_list_len = len(cleaned_seq_list)
    i = 0
    while i < cleaned_seq_list_len:
        if verbosity>0:
            msg = "[%s/%s]" % (counter, cleaned_seq_list_len)
            print(msg, end="\r")
        aln_seq = random_aln_seq(cleaned_seq_list[i], cur_base_aln_len)
        #if the chosen alignment length is too short for at least one sequence, then use
        #a bigger alignment length and restart the loop 
        if aln_seq == '0':
            seq_record_list.clear()
            cur_base_aln_len = next(base_aln_len_sorted_list)
            i = 0
            counter = 1
            continue
        else:
            seq_record = SeqRecord.SeqRecord(Seq.Seq(aln_seq), id=str(counter))
            seq_record_list.append(seq_record)
            i += 1
            counter += 1

    baseline_aln_full_path = os.path.join(job_path, "baseline_aln.fasta")
    SeqIO.write(seq_record_list, baseline_aln_full_path, "fasta")

    aln_file_full_path_list.insert(0, baseline_aln_full_path)

    return aln_file_full_path_list

#creates job dir from a list of alignments
def create_job(aln_dir, out_path, job_name, diamond_exe_full_path, diamond_db_file_full_path, verbosity, unique_seq_only,
               diamond_block_size, diamond_threads, diamond_tmp_dir, baseline, diamond_mp_tmp_dir, diamond_mp):

    #get all alignment files in dir
    aln_dir_list = os.listdir(aln_dir)
    aln_file_full_path_list = []
    for thing in aln_dir_list:
        if os.path.isfile(os.path.join(aln_dir, thing)):
            aln_file_full_path_list.append(os.path.join(aln_dir, thing))

    job_path = os.path.join(out_path, job_name)
    create_dir_safely(job_path)

     #generate baseline alignment file
    if baseline > 0:
        if unique_seq_only < 1:
            err_msg = "Using --baseline currently requires use of --uniqueonly. Please provide the latter to use baseline."
            errorFct(err_msg)
            exit(1)
        aln_file_full_path_list = write_baseline_aln_file(aln_file_full_path_list, job_path, verbosity)

    create_aln_index_files(aln_file_full_path_list, job_path, diamond_exe_full_path, diamond_db_file_full_path, verbosity,
                           diamond_block_size, diamond_threads, diamond_tmp_dir, diamond_mp_tmp_dir, diamond_mp)
    
    if verbosity>0:
        msg = "Determining maximal common sequence set..."
        print(msg)

    #use index files to determine set of common underlying sequences
    #MSA_df_dict will serve as a central information hub
    MSA_df_dict = read_MSA_index_files_in_job(aln_file_full_path_list, job_path)
    #common_seq_set is our C. seq_prim_key_dict_dict is utility structure helping to generate
    #the reduced alignment index files
    common_seq_set, seq_prim_key_dict_dict = generate_common_seq_set(MSA_df_dict, verbosity)
    red_MSA_index_full_path_list = write_red_aln_index_files(common_seq_set, seq_prim_key_dict_dict, unique_seq_only)

    #testing reduced index files, especially the canonical index for correctness
    #this can be left out, once sufficient trust is built
    differing_sequences = set()
    differing_sequences = test_red_index_file(aln_dir, out_path, job_name, 0, differing_sequences, unique_seq_only)
    _ = test_red_index_file(aln_dir, out_path, job_name, 0, differing_sequences, unique_seq_only)

def import_red_aln_index_files(red_aln_index_full_path_list):
    red_MSA_df_dict = {}
    for red_aln_index_full_path in red_aln_index_full_path_list:
        red_MSA_df_dict[red_aln_index_full_path] = pandas.read_csv(red_aln_index_full_path, sep="\t")
    return red_MSA_df_dict

def get_aln_path_list(DATA_path):
    aln_path_list = []
    aln_dir_list = os.listdir(DATA_path)
    for thing in aln_dir_list:
        if not os.path.isfile(thing):
            aln_path_list.append(os.path.join(DATA_path,thing))
    return aln_path_list

def get_red_aln_index_full_path_list(aln_path_list):
    red_aln_index_full_path_list = []
    for aln_path in aln_path_list:
        head, aln_name = os.path.split(aln_path)
        red_aln_index_full_path = os.path.join(aln_path, "red_"+aln_name+"_exact_match.index")
        red_aln_index_full_path_list.append(red_aln_index_full_path)
    return red_aln_index_full_path_list

def get_max_seq_num_for_job(red_MSA_df_dict):
    len_set = set()
    for red_aln_index_full_path in red_MSA_df_dict.keys():
        len_set.add(len(red_MSA_df_dict[red_aln_index_full_path].index))
    if len(len_set) > 1:
        err_msg = "Reduced alignment index files are of unequal length."
        errorFct(err_msg)
        try:
            #fix_red_aln_index_files()
            raise
        except:
            err_msg = "Fixing alignment index failed. Try to create new job."
            errorFct(err_msg)
    else:
        max_seq_num = len_set.pop()
        return max_seq_num

#by design the comparison strategy is a set of 2-tuples
#to be precise: comp_strat \subseteq sample_set x sample_set
#lower triangle means that, if we had a square matrix A=(a_{ik}) with
#i \in sample_set and k \in sample_set, we exclude the diagonal elements
#and also assume symmetry a_{ik}=a_{ki}. Of course this leads to a large number of
#pairwise reference alignments of (n^2-n)/2, where n=len(sample_set) and also is not true
#for arbitrary structure aligners. i.e. alignments a_{ik}=a_{ki} are not identical in general.
#i'm not sure what exactly this means for measurement though.
def generate_low_triangle(max_seq_num):
    comp_strat = set()
    range_list = range(0, max_seq_num)
    for i in range_list:
        for k in range_list:
            if k >= i:
                break
            comp_strat.add((i,k))
    return comp_strat

def load_comp_strat_from_file(comp_strat_file_full_path):
    with open(comp_strat_file_full_path, "r") as file:
        comp_strat_list_list = json.load(file)
    comp_strat = set(tuple(x) for x in comp_strat_list_list)

    return comp_strat

def generate_comp_strat(max_seq_num, sample_size, comp_strat_name="low_triangle", comp_strat_file_path=""):

    if comp_strat_file_path == "":
        int_sample_size = int(sample_size)
        #select from list of different comparison strategies
        if comp_strat_name == "low_triangle":
            comp_strat = generate_low_triangle(max_seq_num)
        elif comp_strat_name == "random":
            pass
        elif comp_strat_name == "complete_linear_coverage":
            pass

        if int_sample_size > 0:
            #random sample on sets is deprecated since python 3.9, that's why i need to waste a bit of resources here
            #to keep using this function
            max_sample_size = len(comp_strat)
            if int_sample_size > max_sample_size:
                int_sample_size = max_sample_size
            sampled_comp_strat_list = random.sample(list(comp_strat), int_sample_size)
            sampled_comp_strat = set(sampled_comp_strat_list)

            return sampled_comp_strat

        return comp_strat
    #else if path to comparison strategy file is provided, then load it up
    else:
        comp_strat = load_comp_strat_from_file(comp_strat_file_path)
        return comp_strat

def get_chain_IDs_from_red_MSA_df(red_MSA_df_dict, comp_strat):
    index_set = set()
    for i,k in comp_strat:
        index_set.add(i)
        index_set.add(k)
    chain_ID_dict = {}
    #index is the common sequence ID
    for index in index_set:
        small_chain_id_set = set()
        small_chain_id_set.clear()
        for key in red_MSA_df_dict.keys():
            chain_id_col_num = red_MSA_df_dict[key].columns.get_loc('chain_id')
            chain_id = red_MSA_df_dict[key].iloc[index, chain_id_col_num]
            small_chain_id_set.add(str(chain_id))
        if len(small_chain_id_set)>1:
            err_msg = "Differing chain IDs for the same common sequence index. Attempting to fix..."
            errorFct(err_msg)
            try:
                #fix_red_aln_index_files()
                raise
            except:
                err_msg = "Fixing alignment index failed. Try to create new job."
                errorFct(err_msg)
        else:
            #maybe decolonize chain ID (?)
            chain_id = small_chain_id_set.pop()
            chain_ID_dict[index] = chain_id
    
    return chain_ID_dict

def write_set_to_file(comp_strat, dest_full_path):
    comp_strat_list = list(comp_strat)
    comp_strat_list_list = [list(element) for element in comp_strat_list]
    with open(dest_full_path, "w") as file:
        json.dump(comp_strat_list_list, file)

#listener function that writes errors to err_log given a queue of error messages
def err_queue_writer(err_queue):
    while True:
        queued_err = err_queue.get()
        if queued_err == "kill":
            break
        open_file_err = errorFct(queued_err)
        if open_file_err is not None:
            raise

#listener function that writes errors to shell_err_log given a queue of error messages
def shell_err_queue_writer(shell_err_queue):
    while True:
        queued_err = shell_err_queue.get()
        if queued_err == "kill":
            break
        open_file_err = shell_err_fct(queued_err)
        if open_file_err is not None:
            raise

def sub_job_queue_handler(no_of_alns, verbosity, backup_freq, json_exclude_path, sub_job_queue):
    counter = 0 
    comp_strat_exclude_set = set()
    while True:
        pair = sub_job_queue.get()
        counter += 1
        if pair == 'kill':
            return comp_strat_exclude_set
        if verbosity>0:
            msg = "Calculated alignment %s_%s. [%s/%s] " % (str(pair[0]), str(pair[1]), str(counter), str(no_of_alns))
            print(msg, end="\r") 

        #add (i,k) tuple to exclude set, since alignment was successful
        comp_strat_exclude_set.add((pair[0], pair[1]))
        
        if backup_freq == 0: 
            continue
        if counter % backup_freq == 0:
            write_set_to_file(comp_strat_exclude_set, json_exclude_path)

#pre-initialize constant data, e.g. not (i,k)
#a single DALI process to be launched in parallel based on i,k
def DALI_comp_strat_child(arg_tuple):
    try:
        #unpack arg_tuple
        err_queue = arg_tuple[9]
        i = arg_tuple[0]
        k = arg_tuple[1]
        ALN_path = arg_tuple[2]
        dali_pl_full_path = arg_tuple[3]
        dali_dat_lib_path = arg_tuple[4]
        wd = arg_tuple[5]
        job_name = arg_tuple[6]
        chain_id_dict = arg_tuple[7]
        shell_err_queue = arg_tuple[8]
    except:
        msg = "Error while unpacking argument tuple in DALI_comp_strat_child process."
        print(msg)
        err_queue.put(msg)
        raise

    try:
        chain_id1 = str(chain_id_dict[i][0:4])+str(chain_id_dict[i][5])
        chain_id2 = str(chain_id_dict[k][0:4])+str(chain_id_dict[k][5])
    except IndexError:
        errMsg = "Chain identifier %s or %s is shorter than 6 chars." % (chain_id_dict[i], chain_id_dict[k])
        err_queue.put(errMsg)
        raise

    try:
        shell_input = []
        sub_job = "%s_%s" % (str(i), str(k))
            
        #creating sub job dir
        sub_job_path = os.path.join(ALN_path, sub_job)    
        create_dir_safely(sub_job_path)
                          
        err_file_path = os.path.join(sub_job_path, "err_log")
        out_file_path = os.path.join(sub_job_path, "output_log")

        #/.../<dali_dir>/dali.pl --cd1 <chain_id1> --cd2 <chain_id2> --dat1 <path> --dat2 <path> --title <string> \
        # --outfmt "summary,alignments,equivalences,transrot" --clean 1> <out_log_file> 2> <err_log_file>
        shell_input = [dali_pl_full_path, '--cd1', chain_id1, '--cd2', chain_id2, 
                    '--dat1', dali_dat_lib_path, '--dat2', dali_dat_lib_path, '--title', job_name,
                    '--outfmt', "summary,alignments,equivalences,transrot",
                    '--clean', '1', '>', out_file_path, '2', '>', err_file_path]

        #navigate to output dir
        os.chdir(sub_job_path)
    except:
        raise
            
    #create file handles for stdout and stderr, only needed for Popen. rework maybe
    #or remove
    try:
        err_file = open(err_file_path, 'w+')
    except:
        errMsg = "Failed to open file %s." % err_file_path
        err_queue.put(errMsg)
        raise
    try:
        out_file = open(out_file_path, 'w+')
    except:
        errMsg = "Failed to open file %s." % out_file_path
        err_queue.put(errMsg)
        raise

    #subprocess for alignment generation
    #need to PIPE stdout and stderr for output forwarding
    try:  
        process = subprocess.Popen(shell_input, stdout=out_file, stderr=err_file)
    except:
        raise
    
    try:
        #timeout may need to be longer (or shorter), depending on size of alignments
        #for the calculation of pairwise alignments of up to a few chains
        #500 seconds seems to be reasonable though
        #remove killing of the process, if it causes problems. Instead, add a warning.
        out, err = process.communicate(timeout=1000)
        #if out is not None and verbosity>1:
            #print(out)
    except subprocess.TimeoutExpired:
        process.kill()
        out, err = process.communicate()
        errMsg = "Communication with dali.pl subprocess timed out during %s. Killed it. Forwarding error to shell_err_log.txt" % (sub_job)
        shell_err_queue.put(err)
        err_queue.put(errMsg)
        raise
    #this exception should not occur when using Popen().communicate()
    #instead errors are logged in specific err_file via the shell_err_watcher process
    #except subprocess.CalledProcessError as e:
     #   process.kill()
      #  out, err = process.communicate()
       # errMsg = "Called process exited with non-zero return code (Code:%s) during dali.pl subprocess of alignment %s." % (str(e.returncode), sub_job)
        #shell_err_queue.put(err)
        #err_queue.put(errMsg)
        #raise
    except:
        process.kill()
        errMsg = "Exception raised during dali.pl subprocess of alignment %s." % (sub_job)
        shell_err_queue.put(err)
        err_queue.put(errMsg)
        raise
    finally:
        close_file_err = close_file_safely(err_file, err_file_path, 'mp')
        close_file_err2 = close_file_safely(out_file, out_file_path, 'mp')
        if close_file_err is not None:
            err_queue.put(close_file_err)
        if close_file_err2 is not None:
            err_queue.put(close_file_err2)
        #if out is not None:
            #return out

    return i,k

def DALI_comp_strat_query_mp(pl_bin_path, small_pdb_lib_path, out_path, dali_dat_lib_path, job_name, comp_strat, chain_id_dict, verbosity, backup_freq, threads):
    #for some weird implementation reasons (probably Fortran though) DALI doesn't generate alignment files for job titles longer than 12 chars
    #and generates only empty alignments for job titles exactly 12 chars long
    #this is a DALI problem and can't be changed at the moment
    #i could, however probably find a workaround, but this isn't really worth the effort right now
    #maybe it's also just because job_title needs to be passed with "" around it. (?) test it!
    comp_strat_exclude_set = set()
    if len(job_name) > 11:
        errMsg = "Job title is longer than 11 characters. Will be cut off in DALI output due to length limitations."
        job_name = job_name[0:11]
    
    #saving cwd for later file system navigation
    wd = os.getcwd()
    
    #changing relative paths to absolute
    if os.path.isabs(pl_bin_path):
        dali_pl_full_path = os.path.join(pl_bin_path, "dali.pl")
    else:
        dali_pl_full_path = os.path.join(wd, pl_bin_path, "dali.pl")
        
    if not os.path.isabs(small_pdb_lib_path):
        small_pdb_lib_path = os.path.join(wd, small_pdb_lib_path)
        
    if not os.path.isabs(out_path):
        out_path = os.path.join(wd, out_path)
    
    if not os.path.isabs(dali_dat_lib_path):
        dali_dat_lib_path = os.path.join(wd, dali_dat_lib_path)

    #creating ALN subdir
    ALN_path = os.path.join(out_path, job_name, "DALI", "ALN")
    create_dir_safely(ALN_path)

    json_exclude_path = os.path.join(ALN_path, "continue_exclude.json")
    try:
        #set up queue manager and process pool
        manager = multiprocessing.Manager()
        shell_err_queue = manager.Queue()
        err_queue = manager.Queue()
        sub_job_queue = manager.Queue()
        if threads == 0:
            threads = multiprocessing.cpu_count()
            pool = multiprocessing.Pool(threads)
        else:
            pool = multiprocessing.Pool(threads)
    
        #launch sub-jobs as child processes
        no_of_alns = len(comp_strat)
        #intentionally underestimating chunk_size to account for variance in alignment generation time
        chunk_size = int(round(no_of_alns/(threads+2)))
        #launch error-handling processes
        shell_err_watcher = pool.apply_async(shell_err_queue_writer, (shell_err_queue,))
        err_watcher = pool.apply_async(err_queue_writer, (err_queue,))
        sub_job_watcher = pool.apply_async(sub_job_queue_handler, (no_of_alns, verbosity, backup_freq, json_exclude_path, sub_job_queue,))

        arg_tuple_list = []
        arg_tuple_list.clear()
        for i,k in comp_strat:
            arg_tuple_list.append(tuple([i, k, ALN_path, dali_pl_full_path, dali_dat_lib_path,
                                  wd, job_name, chain_id_dict, shell_err_queue, err_queue]))
        #fire off workers with task list
        res_iter = pool.imap_unordered(DALI_comp_strat_child, arg_tuple_list, chunksize=chunk_size)
        #loop through results, until every result is returned
        #look into what blocking means exactly and if this loop makes sense
        for j in range(0, len(comp_strat)):
            #next() call SHOULD block until next element is available, unless timeout is given
            #in which case we should raise a TimeoutException
            res_pair = next(res_iter)
            #put result in a queue to be handled by the sub_job_watcher process
            sub_job_queue.put(res_pair)

    except:
        raise
    finally:
        #kill watcher processes
        shell_err_queue.put('kill')
        err_queue.put('kill')
        comp_strat_exclude_set = sub_job_queue.put('kill')
        pool.close()
        pool.join()

def DALI_comp_strat_query(pl_bin_path, small_pdb_lib_path, out_path, dali_dat_lib_path, job_name, comp_strat, chain_id_dict, verbosity, backup_freq, threads, mpi_path):
    #for some weird implementation reasons (probably Fortran though) DALI doesn't generate alignment files for job titles longer than 12 chars
    #and generates only empty alignments for job titles exactly 12 chars long
    #this is a DALI problem and can't be changed at the moment
    #no matter, if the job title is shortened (splicing job_name[0:11]) beforehand and then passed on to DALI, DALI will refuse to generate alignments, if the
    # the name that was INITIALLY specified at the beginning of the script under --jobname exceeds 11 or 12 characters
    #thus, to fix this issue it'd be necessary to look at everything that is being done with the job_name in this script to find out exactly why DALI
    #does not work the expected way with it. The only plausible explanation for this could be that DALI somehow references generated files by looking at the paths
    # in the background, which DO indeed contain the full job name and thus is unable to find intermediate files, if a directory name contained within this path
    # exceeds the character limit. This could be tested by executing DALI manually in a directory which name is longer than 12 chars.
    comp_strat_exclude_set = set()

    #saving cwd for later file system navigation
    wd = os.getcwd()

    #changing relative paths to absolute
    if os.path.isabs(pl_bin_path):
        dali_pl_full_path = os.path.join(pl_bin_path, "dali.pl")
    else:
        dali_pl_full_path = os.path.join(wd, pl_bin_path, "dali.pl")

    if not os.path.isabs(small_pdb_lib_path):
        small_pdb_lib_path = os.path.join(wd, small_pdb_lib_path)

    if not os.path.isabs(out_path):
        out_path = os.path.join(wd, out_path)

    if not os.path.isabs(dali_dat_lib_path):
        dali_dat_lib_path = os.path.join(wd, dali_dat_lib_path)

    #creating ALN subdir
    ALN_path = os.path.join(out_path, job_name, "DALI", "ALN")
    create_dir_safely(ALN_path)

    json_exclude_path = os.path.join(ALN_path, "continue_exclude.json")

    counter = 0
    no_of_alns = len(comp_strat)
    for i,k in comp_strat:

        try:
            chain_id1 = str(chain_id_dict[i][0:4])+str(chain_id_dict[i][5])
            chain_id2 = str(chain_id_dict[k][0:4])+str(chain_id_dict[k][5])
        except IndexError:
            errMsg = "Chain identifier %s or %s is shorter than 6 chars." % (chain_id_dict[i], chain_id_dict[k])
            errorFct(errMsg)
            raise

        shell_input = []
        sub_job = "%s_%s" % (str(i), str(k))

        counter += 1 
        if verbosity>0:
            msg = "Calculating alignment %s. [%s/%s] " % (sub_job, counter, no_of_alns)
            print(msg, end="\r")

        #creating sub job dir
        sub_job_path = os.path.join(ALN_path, sub_job)    
        create_dir_safely(sub_job_path)

        err_file_path = os.path.join(sub_job_path, "err_log")
        out_file_path = os.path.join(sub_job_path, "output_log")


        #/.../<dali_dir>/dali.pl --cd1 <chain_id1> --cd2 <chain_id2> --dat1 <path> --dat2 <path> --title <string> \
        # --outfmt "summary,alignments,equivalences,transrot" --clean 1> <out_log_file> 2> <err_log_file>
        shell_input = [dali_pl_full_path, '--cd1', chain_id1, '--cd2', chain_id2, 
                    '--dat1', dali_dat_lib_path, '--dat2', dali_dat_lib_path, '--title', job_name,
                    '--outfmt', "summary,alignments,equivalences,transrot",
                    '--clean']
        if mpi_path != "":
            shell_input.append('--MPIRUN_EXE')
            shell_input.append('--np')
            if threads > 0:
                shell_input.append(str(threads))
            else:
                shell_input.append(str(multiprocessing.cpu_count()))

        shell_input.extend(['1', '>', out_file_path, '2', '>', err_file_path])

        #navigate to output dir in loop
        os.chdir(sub_job_path)

        #create file handles for stdout and stderr, only needed for Popen. rework maybe
        #or remove
        try:
            err_file = open(err_file_path, 'w+')
        except:
            errMsg = "Failed to open file %s." % err_file_path
            errorFct(errMsg)
            raise
        try:
            out_file = open(out_file_path, 'w+')
        except:
            errMsg = "Failed to open file %s." % out_file_path
            errorFct(errMsg)
            raise

        #subprocess for alignment generation
        #need to PIPE stdout and stderr for output forwarding    
        process = subprocess.Popen(shell_input, stdout=out_file, stderr=err_file)

        try:
            #timeout may need to be longer (or shorter), depending on size of alignments
            #for the calculation of pairwise alignments of up to a few chains
            #500 seconds seems to be reasonable though
            #remove killing of the process, if it causes problems. Instead, add a warning.
            out, err = process.communicate()
            #if out is not None and verbosity>1:
                #print(out)
        except subprocess.TimeoutExpired:
            process.kill()
            out, err = process.communicate()
            errMsg = "Communication with dali.pl subprocess timed out. Killed it. Forwarding error to shell_err_log.txt"
            close_file_safely(err_file, err_file_path, errMsg)
            close_file_safely(out_file, out_file_path, errMsg)
            shell_err_fct(err)
            errorFct(errMsg)
            if out is not None:
                print(out)
            raise

        #close files
        close_file_safely(err_file, err_file_path, "")
        close_file_safely(out_file, out_file_path, "")

        #change dir back to previous wd
        os.chdir(wd)

        #add (i,k) tuple to exclude set, since alignment was successful
        comp_strat_exclude_set.add((i,k))

        if backup_freq == 0:
            continue
        if counter % backup_freq == 0:
            write_set_to_file(comp_strat_exclude_set, json_exclude_path)

def restore_comp_strat_backup(comp_strat, out_path, job_name):
    exclude_set_file_full_path = os.path.join(out_path, job_name, "DALI", "ALN", "continue_exclude.json")
    with open(exclude_set_file_full_path, "r") as file:
        exclude_comp_strat_list_list = json.load(file)
    exclude_comp_strat_set = set(tuple(x) for x in exclude_comp_strat_list_list)
    comp_strat = comp_strat.difference(exclude_comp_strat_set)

    return comp_strat

#after job has been set up, this function uses the data to launch structural alignment generation
# void (str, str, str, str, int, int, str, int, str, str, int, int)
def run_job(out_path, job_name, pl_bin_path, pdb_db_path, verbosity, threads, mpi_path,
            sample_size, comp_strat_name="low_triangle", comp_strat_file_path="", cont=0, backup_freq=100):
    #set and get paths
    DATA_path = os.path.join(out_path, job_name, "DATA")

    auto_gen_comp_strat_file_full_path = os.path.join(out_path, job_name, "comp_strat.json")

    aln_path_list = get_aln_path_list(DATA_path)
    red_aln_index_full_path_list = get_red_aln_index_full_path_list(aln_path_list)
    #import reduced alignment index files as dict of dataframes (keys: red_aln_index_full_path; values: dataframe)
    red_MSA_df_dict = import_red_aln_index_files(red_aln_index_full_path_list)
    #get max sample size
    max_seq_num = get_max_seq_num_for_job(red_MSA_df_dict)
    #if --continuejob flag is passed (cont>0), load comp_strat from auto-generated file and
    #exclude elements that have already been used for pairwise structure alignment
    if cont > 0:
        try:
            comp_strat = load_comp_strat_from_file(auto_gen_comp_strat_file_full_path)
        except:
            err_msg = ("Error loading auto-generated comparison strategy file. Please make sure that the comparison strategy"
                       "specified under --compstratname or --compstratfile matches the one that was initially used to run the job"
                       "when proceeding.")
            errorFct(err_msg)
            comp_strat = generate_comp_strat(max_seq_num, sample_size, comp_strat_name=comp_strat_name, comp_strat_file_path=comp_strat_file_path)
            try:
                write_set_to_file(comp_strat, auto_gen_comp_strat_file_full_path)
            except:
                err_msg = "Warning: Could not write comparison strategy to file for path %s." % auto_gen_comp_strat_file_full_path
                errorFct(err_msg)
        try:
            comp_strat = restore_comp_strat_backup(comp_strat, out_path, job_name)
        except:
            err_msg = ("Unable to read backup file for comparison strategy. By default, it is being generated every 100 structure alignments."
                       "You can specify backup frequency with the --backupfreq flag.")
            errorFct(err_msg)
            raise
    #else generate comparison strategy the normal way and write to file
    else:
        comp_strat = generate_comp_strat(max_seq_num, sample_size, comp_strat_name=comp_strat_name, comp_strat_file_path=comp_strat_file_path)
        try:
            write_set_to_file(comp_strat, auto_gen_comp_strat_file_full_path)
        except:
            err_msg = "Warning: Could not write comparison strategy to file for path %s." % auto_gen_comp_strat_file_full_path
            errorFct(err_msg)

    #everything below this line can be continued after aborting with a reduced comparison strategy
    #-----------------------------------------------------------------------------------------------------------------------------------------
    #generating dir structure for DALI
    dali_job_dir = os.path.join(out_path, job_name, "DALI")
    create_dir_safely(dali_job_dir)

    small_pdb_lib_path = os.path.join(dali_job_dir, "PDB_lib")
    create_dir_safely(small_pdb_lib_path)
    #get chain IDs from reduced MSA dataframes and associate them with common seq ID
    chain_id_dict = get_chain_IDs_from_red_MSA_df(red_MSA_df_dict, comp_strat)
    #create a small copy of PDB files that are relevant to the specific job
    cp_PDB_files_to_job_dir(list(chain_id_dict.values()), pdb_db_path, small_pdb_lib_path)

    dali_dat_lib_path = os.path.join(dali_job_dir, "DAT_lib")
    create_dir_safely(dali_dat_lib_path)

    #retrying import might be necessary sometimes, if DALI behaves unexpectedly
    #for sbatch it is advisable to create a copy of the DALI program for each
    #sbatch call to prevent them from interfering with each other, if they are doing so at all
    retry_import_counter = 0
    small_pdb_lib_list = os.listdir(small_pdb_lib_path)
    small_pdb_id_set = set()
    small_pdb_id_set.clear()
    for file_name in small_pdb_lib_list:
        small_pdb_id_set.add(file_name[3:7].upper())
    
    #remove dali.lock, if it was accidentally written to SCRIPT_DIR
    if os.path.exists(os.path.join(SCRIPT_DIR,'dali.lock')):
        os.remove(os.path.join(SCRIPT_DIR,'dali.lock'))

    #import data for DALI
    DALI_import_PDBs(pl_bin_path, small_pdb_lib_path, dali_dat_lib_path, verbosity)

    #start of a quite elaborate check, if there is a mismatch between PDB files to do alignments with and
    # ones that have successfully been imported into DAT format
    dat_file_set = set()
    dat_file_set.clear()
    for file_name in os.listdir(dali_dat_lib_path):
        if file_name.endswith('.dat'):
            dat_file_set.add(file_name[0:4].upper())
    if len(dat_file_set.difference(small_pdb_id_set)) != 0 or len(small_pdb_id_set.difference(dat_file_set)) != 0:
        if len(dat_file_set.difference(small_pdb_id_set)) != 0:
            mis_msg_1 = "Too many DAT file(s) present. %s" % str(dat_file_set.difference(small_pdb_id_set))
            errorFct(mis_msg_1)
        if len(small_pdb_id_set.difference(dat_file_set)) != 0:
            mis_msg_2 = "Skipped DAT files for PDB entries %s." % str(small_pdb_id_set.difference(dat_file_set))
            errorFct(mis_msg_2)
        err_msg = "Could not import all PDB files into DALI DAT format. [%s/%s] DAT files present." % (len(dat_file_set), len(small_pdb_id_set))
        errorFct(err_msg)
    
    #start pairwise structural alignment
    try:
        DALI_comp_strat_query(pl_bin_path, small_pdb_lib_path, out_path, dali_dat_lib_path, job_name, comp_strat, chain_id_dict, verbosity, backup_freq, threads, mpi_path)
    except:
        err_msg = ("Exception raised during alignment generation. Use flag --continuejob to resume"
                   "alignment from where it failed, after error has been identified.")
        errorFct(err_msg)
        raise

#not implemented.
def guess_comp_strat_from_DALI_dir(DALI_ALN_path):
    err_msg = "Failed to guess comparison strategy from DALI alignments."
    errorFct(err_msg)
    raise

#parses DALI structural alignment information as requested by each entry in the comparison strategy and
#writes it to a dict
# (set(), str, int) -> dict()
def import_DALI_alns_by_comp_strat(comp_strat, DALI_ALN_dir, verbosity):
    #keys are all tuples in comp_strat and values are
    #lists of dictionaries returned by the function import_DALI_aln()
    #each dict represents an individual alignment within a TXT file
    #with keys being 'job_name','Z-score','rmsd','lali','nres','perc_id','pdb_descr','query'(chain id),'sbjct'(chain id),'DALI_query_seq','DALI_sbjct_seq'
    #query is the first index in comp_strat index pair (i), sbjct is the second (k)  
    comp_strat_struct_aln_dict = {}
    no_of_alns = len(comp_strat)
    counter = 0
    for i,k in comp_strat:
        counter += 1
        sub_job_name = "%s_%s" % (str(i),str(k))
        sub_job_path = os.path.join(DALI_ALN_dir, sub_job_name)
        txt_file_counter = 0
        for item in os.listdir(sub_job_path):
            #there should exactly be one .txt file in the dir. it can be empty.
            if item.endswith('.txt'):
                txt_file_counter += 1
                struct_aln_full_path = os.path.join(sub_job_path, item)
        if txt_file_counter != 1:
            warn_msg = "Warning: Not exactly one structural alignment in sub job dir %s. Counted %s." % (sub_job_path, str(txt_file_counter))
            errorFct(warn_msg)
        #import_DALI_aln() function returns a list of dictionaries. each dict corresponds to one alignment
        #here we should have at most one alignment per file, i.e. at most lists of length one
        if txt_file_counter > 0:
            struct_aln_dict_list = import_DALI_aln(struct_aln_full_path)
        else:
            comp_strat_struct_aln_dict[(i,k)] = ["NA"]
        no_of_alns = len(struct_aln_dict_list)
        if no_of_alns != 1:
            warn_msg = "Warning: Not exactly one structural alignment in file %s. Counted %s." % (struct_aln_full_path, str(no_of_alns))
            errorFct(warn_msg)
        comp_strat_struct_aln_dict[(i,k)] = struct_aln_dict_list.copy()

        if verbosity>0:
            msg = "Importing DALI alignments... [%s/%s] " % (str(counter), str(no_of_alns))
            print(msg, end="\r")


    return comp_strat_struct_aln_dict

#acquires information about the restriction bounds of each pair of aligned sequences by reading
#the respective sstart and send values. this function is inefficient because it is practically not necessary to
#iterate through all elements in the comparison strategy. All N sequence would suffice.
#It's just easier this way.
# (set(), dict()) -> dict()
def get_restr_bounds(comp_strat, red_MSA_df_dict):
    restr_bounds_dict = {}
    for i,k in comp_strat:
        sstart_i_set = set()
        send_i_set = set()
        sstart_k_set = set()
        send_k_set = set()
        sstart_i_set.clear()
        send_i_set.clear()
        sstart_k_set.clear()
        send_k_set.clear()
        last_key = ""
        for key in red_MSA_df_dict.keys():
            sstart_col_num = red_MSA_df_dict[key].columns.get_loc('sstart')
            send_col_num = red_MSA_df_dict[key].columns.get_loc('send')
            sstart_i = red_MSA_df_dict[key].iloc[i, sstart_col_num]
            send_i = red_MSA_df_dict[key].iloc[i, send_col_num]
            sstart_k = red_MSA_df_dict[key].iloc[k, sstart_col_num]
            send_k = red_MSA_df_dict[key].iloc[k, send_col_num]
            #watch out for zero-based indexing
            sstart_i_set.add(int(sstart_i)-1)
            send_i_set.add(int(send_i)-1)
            sstart_k_set.add(int(sstart_k)-1)
            send_k_set.add(int(send_k)-1)
            if len(sstart_i_set) != 1 or len(send_i_set) != 1:
                err_msg = ("Differing subject-embedded query index sets ('sstart','send') among common sequence matches"
                           "at index %s among files %s and %s.") % (str(i),last_key,key)
                errorFct(err_msg)
                exit(1)
            if len(sstart_k_set) != 1 or len(send_k_set) != 1:
                err_msg = ("Differing subject-embedded query index sets ('sstart','send') among common sequence matches"
                           "at index %s in files %s and %s.") % (str(k),last_key,key)
                errorFct(err_msg)
                exit(1)
            last_key = key
        restr_bounds_dict[i] = (sstart_i_set.pop(),send_i_set.pop())
        restr_bounds_dict[k] = (sstart_k_set.pop(),send_k_set.pop())

    return restr_bounds_dict

#applies current definition of reference set restriction to individual reference index sets
#as for right now, we chose to only work with the intersection of the two intervals bounds_i and bounds_k
#for reasons that are explained in detail in the adjacent work
#in short: it can be that lack of specificity in MSA algorithm is either punished too harshly or not at all,
#depending on whether or not any given pair of indices is in the union of bounds_i and bounds_k 
# (set(),(int,int),(int,int)) -> set()
def apply_ref_set_restr(small_ref_ind_set, bounds_i, bounds_k):
    restr_small_ref_ind_set = set()
    new_min = max(bounds_i[0], bounds_k[0])
    new_max = min(bounds_i[1], bounds_k[1])
    if new_max < new_min:
        return set()
    for elem in small_ref_ind_set:
        if elem[0][1] >= new_min-1 and elem[0][1] <= new_max and elem[1][1] >= new_min-1 and elem[1][1] <= new_max:
            new_elem = ((elem[0][0],elem[0][1]-bounds_i[0]),(elem[1][0],elem[1][1]-bounds_k[0]))
            restr_small_ref_ind_set.add(new_elem)

    return restr_small_ref_ind_set

#calculates reference sets based on comparison strategy
# (set(), str, dict(), int) -> ( dict(), dict() )
def generate_ref_sets(comp_strat, aln_dir, red_MSA_df_dict, verbosity):
    
    #import structural alignments for all entries in comp_strat
    comp_strat_struct_aln_dict = import_DALI_alns_by_comp_strat(comp_strat, aln_dir, verbosity)

    #get subject-embedded query index set for all (i,k) in comp_strat
    #for that note that each common sequence C_i is associated with exactly
    #one subject_embedded query index set, which can be found in files red_<aln_name>_exact_match.index
    #column 'sstart' represents sms_i and 'send' sms_i+m_i-1 with m_i='qlen'
    #indices 'sstart'-1 and 'send'-1 are both included in the subject-embedded query index set
    #these files are already readily imported into red_MSA_df_dict
    #restr_bounds_dict has common sequences indices i as keys and tuple ('sstart'_i-1, 'send'_i-1) as values
    #property 'qlen' emerges from that as 'qlen'=('send'-1)-('sstart'-1)+1 with min('sstart')=1
    restr_bounds_dict = get_restr_bounds(comp_strat, red_MSA_df_dict)

    #calculate index sets for DALI alignments
    #keys are (i,k) and values are restr_small_ref_ind_sets (reference index sets for structural alignments between i and k)
    restr_small_ref_ind_set_dict = {}
    for i,k in comp_strat:
        aln_query_seq = comp_strat_struct_aln_dict[(i,k)][0]['DALI_query_seq']
        aln_sbjct_seq = comp_strat_struct_aln_dict[(i,k)][0]['DALI_sbjct_seq']
        small_ref_ind_set = calc_ref_SP_set((i,k), (aln_query_seq, aln_sbjct_seq))

        #this set should be smaller or equal to 'qlen'
        restr_small_ref_ind_set = apply_ref_set_restr(small_ref_ind_set, restr_bounds_dict[i], restr_bounds_dict[k])

        #add to dict
        restr_small_ref_ind_set_dict[(i,k)] = restr_small_ref_ind_set.copy()

    return restr_small_ref_ind_set_dict, comp_strat_struct_aln_dict

#function to calculate MSA index sets for comparison strategy
# (set(), dict(), int) -> dict()
def generate_MSA_sets(comp_strat, red_MSA_df_dict, verbosity):

    comp_strat_MSA_aln_dict = {}
    comp_strat_MSA_aln_dict.clear()
    for red_index_file_full_path in red_MSA_df_dict.keys():
        if verbosity>0:
            msg = "Generating MSA set for file %s ..." % red_index_file_full_path
            print(msg)

        orig_fasta_index_col_num = red_MSA_df_dict[red_index_file_full_path].columns.get_loc('fasta_entry_index')
        fasta_header_col_num = red_MSA_df_dict[red_index_file_full_path].columns.get_loc('fasta_entry_header')
        gapped_seq_col_num = red_MSA_df_dict[red_index_file_full_path].columns.get_loc('gapped_seq')
        chain_id_col_num = red_MSA_df_dict[red_index_file_full_path].columns.get_loc('chain_id')
        evalue_col_num = red_MSA_df_dict[red_index_file_full_path].columns.get_loc('evalue')
        bitscore_col_num = red_MSA_df_dict[red_index_file_full_path].columns.get_loc('bitscore')
        sstart_col_num = red_MSA_df_dict[red_index_file_full_path].columns.get_loc('sstart')
        send_col_num = red_MSA_df_dict[red_index_file_full_path].columns.get_loc('send')
        MSA_ind_set_dict = {}
        comp_strat_red_ind_dict = {}
        MSA_ind_set_dict.clear()
        comp_strat_red_ind_dict.clear()
        no_of_alns = len(comp_strat)
        counter = 0
        for i,k in comp_strat:
            counter += 1
            if verbosity>0:
                msg = "[%s/%s]" % (str(counter), str(no_of_alns))
                print(msg, end="\r")
            #generate MSA index sets
            gapped_seq_i = red_MSA_df_dict[red_index_file_full_path].iloc[i, gapped_seq_col_num]
            gapped_seq_k = red_MSA_df_dict[red_index_file_full_path].iloc[k, gapped_seq_col_num]
            MSA_ind_set = calc_MSA_SP_set((i,k), (gapped_seq_i, gapped_seq_k))
            MSA_ind_set_dict[(i,k)] = MSA_ind_set.copy()
            #generate additional info
            row_dict = {}
            row_dict.clear()
            orig_fasta_index_i = red_MSA_df_dict[red_index_file_full_path].iloc[i, orig_fasta_index_col_num]
            orig_fasta_index_k = red_MSA_df_dict[red_index_file_full_path].iloc[k, orig_fasta_index_col_num]
            fasta_header_i = red_MSA_df_dict[red_index_file_full_path].iloc[i, fasta_header_col_num]
            fasta_header_k = red_MSA_df_dict[red_index_file_full_path].iloc[k, fasta_header_col_num]
            chain_id_i = red_MSA_df_dict[red_index_file_full_path].iloc[i, chain_id_col_num]
            chain_id_k = red_MSA_df_dict[red_index_file_full_path].iloc[k, chain_id_col_num]
            evalue_i = red_MSA_df_dict[red_index_file_full_path].iloc[i, evalue_col_num]
            evalue_k = red_MSA_df_dict[red_index_file_full_path].iloc[k, evalue_col_num]
            bitscore_i = red_MSA_df_dict[red_index_file_full_path].iloc[i, bitscore_col_num]
            bitscore_k = red_MSA_df_dict[red_index_file_full_path].iloc[k, bitscore_col_num]
            sstart_i = red_MSA_df_dict[red_index_file_full_path].iloc[i, sstart_col_num]
            sstart_k = red_MSA_df_dict[red_index_file_full_path].iloc[k, sstart_col_num]
            send_i = red_MSA_df_dict[red_index_file_full_path].iloc[i, send_col_num]
            send_k = red_MSA_df_dict[red_index_file_full_path].iloc[k, send_col_num]
            row_dict['fasta_entry_index'] = (orig_fasta_index_i, orig_fasta_index_k)
            row_dict['fasta_entry_header'] = (fasta_header_i, fasta_header_k)
            row_dict['chain_id'] = (chain_id_i, chain_id_k)
            row_dict['evalue'] = (evalue_i, evalue_k)
            row_dict['bitscore'] = (bitscore_i, bitscore_k)
            row_dict['sstart'] = (sstart_i, sstart_k)
            row_dict['send'] = (send_i, send_k)
            comp_strat_red_ind_dict[(i,k)] = row_dict.copy()
        comp_strat_MSA_aln_dict[red_index_file_full_path] = (MSA_ind_set_dict.copy(), comp_strat_red_ind_dict.copy())
    
    return comp_strat_MSA_aln_dict

#this function calculates the pstar-score for each pair of sequences given by the comparison strategy
#this function can be exchanged for any function that calculates values for each pair of sequences
# (set(), dict(), dict()) -> dict()
def calc_pair_SPS_by_comp_strat(comp_strat, ref_set_dict, comp_strat_MSA_aln_dict):
    pair_SPS_dict = {}
    pair_SPS_dict.clear()
    for red_index_file_full_path in comp_strat_MSA_aln_dict.keys():
        MSA_ind_set_dict = comp_strat_MSA_aln_dict[red_index_file_full_path][0]
        small_pair_SPS_dict = {}
        small_pair_SPS_dict.clear()
        for i,k in comp_strat:
            intersect_set = MSA_ind_set_dict[(i,k)].intersection(ref_set_dict[(i,k)])
            if len(ref_set_dict[(i,k)]) != 0:
                pair_SPS = len(intersect_set)/len(ref_set_dict[(i,k)])
            else:
                pair_SPS = 0
            small_pair_SPS_dict[(i,k)] = pair_SPS
        pair_SPS_dict[red_index_file_full_path] = small_pair_SPS_dict.copy()

    return pair_SPS_dict

#calculate SPS for each pair of sequences indexed by comparison strategy and write a detailed output file
# void (str, set(), dict(), dict(), dict(), int)
def write_output_tsv(eval_path, comp_strat, ref_set_dict, comp_strat_struct_aln_dict, comp_strat_MSA_aln_dict, verbosity):
    #pair_SPS_dict: keys are paths to red_index_file; values are small_pair_SPS_dict
    #small_pair_SPS_dict: keys are each element in comp_strat (i,k); values are the SPS between sequences i and k ("pair_SPS")
    pair_SPS_dict = calc_pair_SPS_by_comp_strat(comp_strat, ref_set_dict, comp_strat_MSA_aln_dict)
    
    #merge comp_strat_MSA_aln_dict with pair_SPS_dict and comp_strat_struct_aln_dict
    output_dict = {}
    output_dict.clear()
    no_of_alns = len(comp_strat)
    for red_ind_file_full_path in comp_strat_MSA_aln_dict.keys():
        counter = 0
        small_output_dict = {}
        small_output_dict.clear()
        if verbosity>0:
            msg = "Preparing output for file %s..." % (red_ind_file_full_path)
            print(msg)
        for i,k in comp_strat:
            counter += 1
            if verbosity>0:
                msg = "[%s/%s]" % (str(counter), str(no_of_alns))
                print(msg, end="\r")
            row_dict = {}
            row_dict.clear()
            row_dict['fasta_entry_index'] = comp_strat_MSA_aln_dict[red_ind_file_full_path][1][(i,k)]['fasta_entry_index']
            row_dict['fasta_entry_header'] = comp_strat_MSA_aln_dict[red_ind_file_full_path][1][(i,k)]['fasta_entry_header']
            row_dict['chain_id'] = comp_strat_MSA_aln_dict[red_ind_file_full_path][1][(i,k)]['chain_id']
            row_dict['evalue'] = comp_strat_MSA_aln_dict[red_ind_file_full_path][1][(i,k)]['evalue']
            row_dict['bitscore'] = comp_strat_MSA_aln_dict[red_ind_file_full_path][1][(i,k)]['bitscore']
            row_dict['sstart'] = comp_strat_MSA_aln_dict[red_ind_file_full_path][1][(i,k)]['sstart']
            row_dict['send'] = comp_strat_MSA_aln_dict[red_ind_file_full_path][1][(i,k)]['send']
            row_dict['Z-score'] = comp_strat_struct_aln_dict[(i,k)][0]['Z-score']
            row_dict['rmsd'] = comp_strat_struct_aln_dict[(i,k)][0]['rmsd']
            row_dict['lali'] = comp_strat_struct_aln_dict[(i,k)][0]['lali']
            row_dict['nres'] = comp_strat_struct_aln_dict[(i,k)][0]['nres']
            row_dict['perc_id'] = comp_strat_struct_aln_dict[(i,k)][0]['perc_id']
            row_dict['pdb_descr'] = comp_strat_struct_aln_dict[(i,k)][0]['pdb_descr']
            #next two are only for verification
            row_dict['query_id'] = comp_strat_struct_aln_dict[(i,k)][0]['query']
            row_dict['sbjct_id'] = comp_strat_struct_aln_dict[(i,k)][0]['sbjct']
            #include pair_SPS
            row_dict['pair_SPS'] = pair_SPS_dict[red_ind_file_full_path][(i,k)]
            small_output_dict[(i,k)] = row_dict.copy()

        output_dict[red_ind_file_full_path] = small_output_dict.copy()

    #use output dict to write dataframe to csv for each key
    for red_ind_file_full_path in output_dict.keys():
        if verbosity>0:
            msg = "Writing output file for %s..." % (red_ind_file_full_path)
            print(msg)
        path, file = os.path.split(red_ind_file_full_path)
        path2, aln_name = os.path.split(path)
        tsv_output_file_full_path = os.path.join(eval_path, aln_name+".tsv")
        output_df = pandas.DataFrame.from_dict(output_dict[red_ind_file_full_path], orient='index')
        output_df.to_csv(tsv_output_file_full_path, sep='\t')

#aggregates all small reference and MSA index sets into one large set each
# (dict, dict) -> (set, set) 
def pstar_union(ref_set_dict, comp_strat_MSA_aln_dict):
    
    union_ref_set = set.union(*list(ref_set_dict.values()))

    union_test_set_dict = {}
    union_test_set_dict.clear()
    for red_aln_index_full_path in comp_strat_MSA_aln_dict.keys():
        union_test_set_dict[red_aln_index_full_path] = set.union(*list(comp_strat_MSA_aln_dict[red_aln_index_full_path][0].values())).copy()

    return union_ref_set, union_test_set_dict

def write_scores_to_file(union_test_set_dict, eval_path, union_ref_set, verbosity):

    scores_file_full_path = os.path.join(eval_path, 'scores_per_alignment.tsv')
    jacc_scores_file_full_path = os.path.join(eval_path, 'jacc_scores_per_alignment.tsv')

    scores_dict = {}
    jacc_scores_dict = {}
    for key in union_test_set_dict.keys():
        #set paths and write individual sets to files
        head, file_name = os.path.split(key)
        head, aln_name = os.path.split(head)
        union_test_set_full_path = os.path.join(eval_path, "union_"+aln_name+"_MSA_set.json")
        write_set_to_file(union_test_set_dict[key], union_test_set_full_path)
        intersect_union = union_test_set_dict[key].intersection(union_ref_set)
        union_A_C = union_test_set_dict[key].union(union_ref_set)
        if len(union_ref_set) > 0:
            SPS = len(intersect_union)/len(union_ref_set)
        else:
            SPS = 0
        if len(union_A_C) > 0:
            jacc_score = len(intersect_union)/len(union_A_C)
        else:
            jacc_score = 0
        scores_dict[os.path.split(key)[1]] = SPS
        jacc_scores_dict[os.path.split(key)[1]] = jacc_score
        if verbosity > 0:
            sps_msg = "Score for alignment file %s: %s" % (key, str(SPS))
            jacc_msg = "Jaccard index for alignment file %s: %s" % (key, str(jacc_score))
            print(sps_msg)
            print(jacc_msg)
    jacc_scores_df = pandas.DataFrame.from_dict(jacc_scores_dict, orient='index', columns=['Jaccard index'])
    scores_df = pandas.DataFrame.from_dict(scores_dict, orient='index', columns=['score'])
    scores_df.to_csv(scores_file_full_path, sep='\t')
    jacc_scores_df.to_csv(jacc_scores_file_full_path, sep='\t')

def write_PSM_and_jacc_DM_to_file(union_test_set_dict, eval_path, union_ref_set):
    jacc_dm_file_full_path = os.path.join(eval_path, 'jacc_dm.tsv')
    PSM_dm_file_full_path = os.path.join(eval_path, 'psm_dm.tsv')

    jacc_dict = {}
    PSM_dict = {}
    for key_A in union_test_set_dict.keys():
        for key_B in union_test_set_dict.keys():

            intersect_A_C = union_test_set_dict[key_A].intersection(union_ref_set)
            intersect_B_C = union_test_set_dict[key_B].intersection(union_ref_set)
            intersect_A_B = union_test_set_dict[key_A].intersection(union_test_set_dict[key_B])
            union_A_C = union_test_set_dict[key_A].union(union_ref_set)
            union_B_C = union_test_set_dict[key_B].union(union_ref_set)
            union_A_B = union_test_set_dict[key_A].union(union_test_set_dict[key_B])
            if len(union_test_set_dict[key_A]) == 0 and len(union_test_set_dict[key_B]) == 0:
                PSM = float(0)
                JACC = float(0)
            elif (len(union_test_set_dict[key_A]) == 0 and len(union_ref_set) == 0) or (len(union_test_set_dict[key_B]) == 0 and len(union_ref_set) == 0):
                PSM = float(0.5)
                JACC = 1-(len(intersect_A_B)/len(union_A_B))
            else:
                JACC = 1-(len(intersect_A_B)/len(union_A_B))
                PSM = 0.5*(abs(((len(intersect_A_C))/(len(union_A_C)))-((len(intersect_B_C))/(len(union_B_C))))+1-(len(intersect_A_B)/len(union_A_B)))

            jacc_dict[(os.path.split(key_A)[1],os.path.split(key_B)[1])] = JACC
            PSM_dict[(os.path.split(key_A)[1],os.path.split(key_B)[1])] = PSM

    jacc_df = pandas.DataFrame.from_dict(jacc_dict, orient='index', columns=['Jaccard-distance'])
    PSM_df = pandas.DataFrame.from_dict(PSM_dict, orient='index', columns=['PSTAR-distance'])
    jacc_df.to_csv(jacc_dm_file_full_path, sep='\t')
    PSM_df.to_csv(PSM_dm_file_full_path, sep='\t')


#calculates numbers based off precalculated structural alignments and corresponding comparison strategy
# void (str, str, int)
def evaluate_job(out_path, job_name, verbosity):

    auto_gen_comp_strat_file_full_path = os.path.join(out_path, job_name, "comp_strat.json")
    job_path = os.path.join(out_path, job_name)
    DATA_path = os.path.join(job_path, "DATA")
    DALI_path = os.path.join(job_path, "DALI")
    DALI_ALN_path = os.path.join(DALI_path, "ALN")

    #import reduced alignment index files as dict of dataframes (keys: red_aln_index_full_path; values: dataframe)
    aln_path_list = get_aln_path_list(DATA_path)
    red_aln_index_full_path_list = get_red_aln_index_full_path_list(aln_path_list)
    red_MSA_df_dict = import_red_aln_index_files(red_aln_index_full_path_list)

    try:
        comp_strat = load_comp_strat_from_file(auto_gen_comp_strat_file_full_path)
    except:
        err_msg = ("Error loading auto-generated comparison strategy file. Trying to guess comparison strategy from generated alignments.")
        errorFct(err_msg)
        comp_strat = guess_comp_strat_from_DALI_dir(DALI_ALN_path)
        exit(1)

    #generate reference sets
    #ref_set_dict: keys are (i,k); values are reference index sets
    #comp_strat_struct_aln_dict: keys are (i,k); values are (one-element) lists of dict with struct alignment data (s. function "import_DALI_alns_by_comp_strat()")
    ref_set_dict, comp_strat_struct_aln_dict = generate_ref_sets(comp_strat, DALI_ALN_path, red_MSA_df_dict, verbosity)


    #generate MSA index sets
    #comp_strat_MSA_aln_dict: keys are paths to red_index_file (<path>/DATA/<aln_name>/red_<aln_name>_exact_match.index)
    #values are 2-tuples. at pos [0] they have MSA_index_set_dict, at pos [1] they have comp_strat_red_ind_dict

        #MSA_index_set_dict: keys are each element in comp_strat (i,k); values are MSA_index_sets
        # for gapped sequences i and k

        #comp_strat_red_ind_dict: keys are each element in comp_strat (i,k); values are 2-tuples with
        # pos [0] representing row i and pos [1] representing row k in red_index_file as dictionaries
   
    comp_strat_MSA_aln_dict = generate_MSA_sets(comp_strat, red_MSA_df_dict, verbosity)

    #create output folder and 
    eval_path = os.path.join(job_path, "EVAL")
    create_dir_safely(eval_path)

    write_output_tsv(eval_path, comp_strat, ref_set_dict, comp_strat_struct_aln_dict, comp_strat_MSA_aln_dict, verbosity)

    union_ref_set, union_test_set_dict = pstar_union(ref_set_dict, comp_strat_MSA_aln_dict)
    #write sets to files
    union_ref_set_full_path = os.path.join(eval_path, "union_ref_set.json")
    write_set_to_file(union_ref_set, union_ref_set_full_path)
    #write scores to file
    write_scores_to_file(union_test_set_dict, eval_path, union_ref_set, verbosity)
    #write PSTAR metric distance matrix to file
    write_PSM_and_jacc_DM_to_file(union_test_set_dict, eval_path, union_ref_set)

def test_SPS(aln_file_full_path):

    #[sequence, header, number, cleaned_sequence]
    entry_list = import_seq_list_from_fasta_aln(aln_file_full_path)

    #((int,int),(str,str)) -> set( ((),()) )
    sequence_ids_1 = (1,2)
    aligned_sequences_1 = ((entry_list[0][0]),(entry_list[1][0]))
    sequence_ids_2 = (1,2)
    aligned_sequences_2 = ((entry_list[2][0]),(entry_list[3][0]))
    MSA_set_1 = calc_MSA_SP_set(sequence_ids_1, aligned_sequences_1)
    MSA_ref_set_2 = calc_MSA_SP_set(sequence_ids_2, aligned_sequences_2)

    print("SET 1:\n")
    print(MSA_set_1)
    print("\nSET 2:\n")
    print(MSA_ref_set_2)
    print("\nINTERSECTION:\n")
    print(MSA_set_1.intersection(MSA_ref_set_2))
    print("\nSPS:\n")

    SPS = len(MSA_set_1.intersection(MSA_ref_set_2))/len(MSA_ref_set_2)

    print(SPS)
    print("\n")

def test_duplicate_clean_seq(fasta_file_full_path):
    #[sequence, header, number, cleaned_sequence]
    entry_list = import_seq_list_from_fasta_aln(fasta_file_full_path)

    seq_set = set()
    header_set = set()
    seq_set.clear()
    header_set.clear()
    seq_header_dict = {}
    entry_no = len(entry_list)
    last_seq_set_size = len(seq_set)
    for entry in entry_list:
        upper_seq = str(entry[0]).upper()
        seq_set.add(upper_seq)
        header_set.add(entry[1])
        if last_seq_set_size == len(seq_set):
            msg = "Seq\n %s \n with header %s is the same as seq with header %s." %(str(entry[0]),str(entry[1]),str(seq_header_dict[upper_seq]))
            print(msg)
        last_seq_set_size = len(seq_set)
        seq_header_dict[upper_seq] = entry[1]
    seq_msg = "Sequences: [%s/%s]" % (str(len(seq_set)), str(entry_no))
    header_msg = "Headers: [%s/%s]" % (str(len(header_set)), str(entry_no))
    print(seq_msg)
    print(header_msg)

def test_PDB_AtomIterator(pdb_file_full_path_1, pdb_file_full_path_2, output_path, dali_pl_full_path, import_pl_full_path, chain_letter_1, chain_letter_2):

    try:
        pdb_file_handle_1 = gzip.open(pdb_file_full_path_1, 'rt')
        pdb_file_handle_2 = gzip.open(pdb_file_full_path_2, 'rt')
    except:
        errMsg = "Failed to open file."
        errorFct(errMsg)
        exit(1)
    try:
        chains_1 = SeqIO.PdbIO.PdbAtomIterator(pdb_file_handle_1)
        chains_2 = SeqIO.PdbIO.PdbAtomIterator(pdb_file_handle_2)
    
        for record in chains_1:
            msg = "Record id %s, sequence %s, length %s" % (str(record.id), str(record.seq), str(len(str(record.seq))))
            print(msg)
        for record in chains_2:
            msg = "Record id %s, sequence %s, length %s" % (str(record.id), str(record.seq), str(len(str(record.seq))))
            print(msg)
    except:
        raise
    finally:
        close_file_safely(pdb_file_handle_1,pdb_file_full_path_1,'')
        close_file_safely(pdb_file_handle_2,pdb_file_full_path_2,'')

    wd = os.getcwd()
    output_path = os.path.join(wd, output_path)
    dali_pl_full_path = os.path.join(wd, dali_pl_full_path)
    import_pl_full_path = os.path.join(wd, import_pl_full_path)
    DAT_path = os.path.join(output_path, 'DAT_path')
    pdb_file_full_path_1 = os.path.join(wd, pdb_file_full_path_1)
    pdb_file_full_path_2 = os.path.join(wd, pdb_file_full_path_2)
    create_dir_safely(DAT_path)

    head, ext = os.path.splitext(pdb_file_full_path_1)
    head2, ext = os.path.splitext(head)
    path, file_name = os.path.split(head2)
    pdb_name_1 = file_name[3:7]
    head, ext = os.path.splitext(pdb_file_full_path_2)
    head2, ext = os.path.splitext(head)
    path, file_name = os.path.split(head2)
    pdb_name_2 = file_name[3:7]
    pdb_name_dict = {}
    pdb_name_dict[pdb_name_1] = pdb_file_full_path_1
    pdb_name_dict[pdb_name_2] = pdb_file_full_path_2

    for key in pdb_name_dict.keys():
        #try with and without importing first
        shell_input_import = [import_pl_full_path, '--pdbfile', pdb_name_dict[key], '--pdbid', key.upper(),
                       '--dat', DAT_path,'1>','/dev/null','2>','/dev/null']

        #subprocess for importing PDBs
        #shell input already sends stdout to out_file_path and stderr to err_file_path
        #probably just pipe output in the next line    
        process = subprocess.Popen(shell_input_import)
        
        try:
            out, err = process.communicate()
            #its really just a lot to print...
            #if out is not None and verbosity>1:
                #print(out)
        except:
            process.kill()
            out, err = process.communicate()
            errMsg = "Communication with import.pl subprocess timed out. Killed it."
            print(errMsg)
            exit(1)

    chain_id1 = pdb_name_1+chain_letter_1
    chain_id2 = pdb_name_2+chain_letter_2
    job_title = 'TestAtIt'

    print(chain_id1)
    print(chain_id2)

    #/.../<dali_dir>/dali.pl --cd1 <chain_id1> --cd2 <chain_id2> --dat1 <path> --dat2 <path> --title <string> \
    # --outfmt "summary,alignments,equivalences,transrot" --clean 1> <out_log_file> 2> <err_log_file>
    shell_input_align = [dali_pl_full_path, '--cd1', chain_id1.upper(), '--cd2', chain_id2.upper(), 
                '--dat1', DAT_path, '--dat2', DAT_path, '--title', job_title,
                '--outfmt', "summary,alignments,equivalences,transrot",
                '--clean','1>','/dev/null','2>','/dev/null']

    os.chdir(output_path)

    process = subprocess.Popen(shell_input_align)
        
    try:
        out, err = process.communicate()
        #its really just a lot to print...
        #if out is not None and verbosity>1:
            #print(out)
    except:
        process.kill()
        out, err = process.communicate()
        errMsg = "Communication with dali.pl subprocess timed out. Killed it."
        print(errMsg)
        exit(1)

    os.chdir(wd)

def transform_alignment(fasta_file_full_path, output_file_full_path):
    with open(fasta_file_full_path, 'r') as fasta_file_handle:
        alignment = AlignIO.read(fasta_file_handle, "fasta")
    for record in alignment:
        seq = str(record.seq)
        seq = re.sub('\.', '-', seq)
        seq = seq.upper()
        record.seq = Seq.Seq(seq)
    
    SeqIO.write(alignment, output_file_full_path, "fasta")

#maps a single input FASTA to PDB files based on exact DIAMOND matches against the PDB-derived database.
#writes per-sequence filter statistics to <out_path>/<fasta_name>_filter_stats.tsv and copies matching
#PDB files (as .ent.gz) to <out_path>/PDB/.
#failure reasons:
#  no_diamond_hit  - no DIAMOND hit at all for this sequence
#  no_exact_match  - DIAMOND hit exists but fails exact-match filter (0 mismatches, 0 gaps, full length)
# void (str, str, str, str, str, int, float, int, str, str, int)
def map_fasta_to_pdb(fasta_file_full_path, out_path, diamond_exe_full_path, diamond_db_file_full_path,
                     loc_pdb_db_path, verbosity, diamond_block_size, diamond_threads,
                     diamond_tmp_dir, diamond_mp_tmp_dir, diamond_mp):

    create_dir_safely(out_path)
    fasta_name = os.path.splitext(os.path.basename(fasta_file_full_path))[0]

    # clean input FASTA (strip gaps and non-alphabet characters)
    clean_fasta_full_path = os.path.join(out_path, fasta_name + "_cleaned.fasta")
    clean_fasta_file(fasta_file_full_path, clean_fasta_full_path, verbosity)

    # import sequence list: [gapped_seq, header, entry_number, ungapped_seq]
    aln_entry_list = import_seq_list_from_fasta_aln(fasta_file_full_path)
    total_count = len(aln_entry_list)
    if verbosity > 0:
        print("Total sequences in input: %d" % total_count)

    # run DIAMOND blastp
    TSV_full_path = os.path.join(out_path, fasta_name + "_diamond.tsv")
    diamond_blastp_query(diamond_exe_full_path, diamond_db_file_full_path, clean_fasta_full_path,
                         TSV_full_path, verbosity, diamond_block_size, diamond_threads,
                         diamond_tmp_dir, diamond_mp_tmp_dir, diamond_mp)

    # collect best raw DIAMOND hit per header (for reporting on sequences that fail the exact-match filter)
    best_raw_hit = {}
    try:
        raw_df = pandas.read_csv(TSV_full_path, sep='\t')
        for index, row in raw_df.iterrows():
            qseqid = str(row.qseqid)
            if qseqid not in best_raw_hit:
                best_raw_hit[qseqid] = row
            else:
                if float(row.bitscore) > float(best_raw_hit[qseqid].bitscore):
                    best_raw_hit[qseqid] = row
    except Exception:
        pass

    # apply exact-match filter and resolve multi-hit ambivalences
    exact_match_dict = get_exact_matches(TSV_full_path, aln_entry_list)
    unique_exact_match_dict = resolve_ambivalences(exact_match_dict, fasta_file_full_path)

    # build per-sequence statistics and collect chain IDs for PDB copying
    stat_rows = []
    passed_chain_list = []
    for entry in aln_entry_list:
        header = entry[1]
        row_stat = {"fasta_header": header}
        if len(unique_exact_match_dict[header]) == 1:
            match = unique_exact_match_dict[header][0]
            row_stat["status"] = "passed"
            row_stat["fail_reason"] = ""
            row_stat["chain_id"] = match.sseqid
            row_stat["sstart"] = match.sstart
            row_stat["send"] = match.send
            row_stat["pident"] = match.pident
            row_stat["evalue"] = match.evalue
            row_stat["bitscore"] = match.bitscore
            passed_chain_list.append(str(match.sseqid))
        else:
            if header in best_raw_hit:
                fail_reason = "no_exact_match"
                best = best_raw_hit[header]
                row_stat["pident"] = best.pident
                row_stat["evalue"] = best.evalue
                row_stat["bitscore"] = best.bitscore
            else:
                fail_reason = "no_diamond_hit"
                row_stat["pident"] = ""
                row_stat["evalue"] = ""
                row_stat["bitscore"] = ""
            row_stat["status"] = "failed"
            row_stat["fail_reason"] = fail_reason
            row_stat["chain_id"] = ""
            row_stat["sstart"] = ""
            row_stat["send"] = ""
        stat_rows.append(row_stat)

    # write filter statistics TSV
    stats_full_path = os.path.join(out_path, fasta_name + "_filter_stats.tsv")
    stats_df = pandas.DataFrame(stat_rows, columns=["fasta_header", "status", "fail_reason",
                                                     "chain_id", "sstart", "send",
                                                     "pident", "evalue", "bitscore"])
    stats_df.to_csv(stats_full_path, sep='\t', index=False)

    # print summary
    passed_count = sum(1 for r in stat_rows if r["status"] == "passed")
    no_hit_count = sum(1 for r in stat_rows if r["fail_reason"] == "no_diamond_hit")
    no_exact_count = sum(1 for r in stat_rows if r["fail_reason"] == "no_exact_match")
    if verbosity > 0:
        print("Passed (exact PDB match found):  %d / %d" % (passed_count, total_count))
        print("Failed - no DIAMOND hit:         %d / %d" % (no_hit_count, total_count))
        print("Failed - hit but not exact match: %d / %d" % (no_exact_count, total_count))
        print("Filter statistics written to: %s" % stats_full_path)

    # copy matching PDB files to <out_path>/PDB/
    pdb_out_path = os.path.join(out_path, "PDB")
    create_dir_safely(pdb_out_path)
    cp_PDB_files_to_job_dir(passed_chain_list, loc_pdb_db_path, pdb_out_path)
    if verbosity > 0:
        print("PDB files written to: %s" % pdb_out_path)

def main():

    tic = time.perf_counter()
    
    set_script_path()

    #multiprocessing.set_start_method('spawn')

    argParser = argparse.ArgumentParser()
    argParser.add_argument("-p", "--createsearchplots", default=0, action="count", help="Creates search plots from CSV search metadata. Not a feature at the moment.", required=False)
    argParser.add_argument("-all", "--wrapjob", default=0, action="count", help="Sets up job, runs structural alignments and evaluates results.", required=False)
    argParser.add_argument("-sync", "--syncdb", default=0, action="count", help="Sync local PDB database with remote.", required=False)
    argParser.add_argument("-db", "--locpdbdb", help="Path to local PDB database.", required=False)
    argParser.add_argument("-af", "--alignmentfile", help="MSA file in FASTA format. No purpose at the moment.", required=False)
    argParser.add_argument("-bat", "--msadir", help="Path to MSA files to be compared; in FASTA format", required=False)
    argParser.add_argument("-ss", "--samplesize", default=0, type = int,
                           help="Give sample size of random comparisons taken from comparison strategy as integer. Default:0, means full sample.", required=False)
    argParser.add_argument("-bf", "--backupfreq", default=100, type = int, help="Specify how often you'd like to save progress of alignment generation.", required=False)
    argParser.add_argument("-csn", "--compstratname", choices=["low_triangle"], default="low_triangle",
                           help="Select comparison strategy by name. Default: low_triangle", required=False)
    argParser.add_argument("-csf", "--compstratfile", default="", help="Path to comparison strategy set in file format.", required=False)
    argParser.add_argument("-sn", "--samplenumber",type = int, default = 1, help="Provide integer for how many times you want to sample. Default:1", required=False)
    argParser.add_argument("-dali", "--dalibindir", help="Path to DaliLite.v5/bin directory", required=False)
    argParser.add_argument("-ti", "--jobname", help="Job title. At most 11 characters.", required=False)
    argParser.add_argument("-o", "--outputdir", help="Path to output dir", required=False)
    argParser.add_argument("-a", "--alphabet", choices=["AA", "DNA/RNA"], default="AA", help="Select alphabet: (AA, DNA/RNA). Default: AA", required=False)
    argParser.add_argument("-na", "--nonalphabet", default="-.", help="Select non-alphabet. Default: -.", required=False)
    argParser.add_argument("-v", "--verbose", default=0, action="count",
                           help="Print progress to terminal. No effect on error logging.", required=False)
    argParser.add_argument("-cnt", "--continuejob", default=0, action="count", help="Continue job from last checkpoint", required=False)
    argParser.add_argument("-dia", "--diamondfile", help="Path to diamond program file.", required=False)
    argParser.add_argument("-clf", "--cleanfasta", action="count", default=0, help="testing fasta cleaning", required=False)
    argParser.add_argument("-ddb", "--diamonddbfile", help="Path to diamond database file.", required=False)
    argParser.add_argument("-cj", "--createjob", default=0, action="count", help="Set up job folder and required data.", required=False)
    argParser.add_argument("-rj", "--runjob", default=0, action="count", help="Running job by specifying output folder and job title.", required=False)
    argParser.add_argument("-ej", "--evaljob", default=0, action="count", help="Evaluating finished job.", required=False)
    argParser.add_argument("-uo", "--uniqueonly", default=1, action="count", help="Provide this flag to only use each common sequence at most once.", required=False)
    argParser.add_argument("-bz", "--diamondblocksize", type = float, default = None, help="Provide float for diamond block size. Default is 2.0.", required=False)
    argParser.add_argument("-dtd", "--diamondtmpdir", default = "", help="Provide temporary storage directory for diamond. Default is output dir.", required=False)
    argParser.add_argument("-mptd", "--diamondmptmpdir", default = "", help="Provide temporary shared storage directory for diamond multiprocessing.", required=False)
    argParser.add_argument("-cpu", "--threads", type = int, default = 0, help=("Provide int for number of CPU threads to use."
                           "If not provided, diamond will try to auto-detect number of available cores and DALI will run on one core only."), required=False)
    argParser.add_argument("-dmp", "--diamondmp", default=0, action="count", help="Provide this flag to use diamond in a multiprocessing cluster environment across multiple nodes.", required=False)
    argParser.add_argument("-mpi", "--mpibin", default = "", help="Provide path to \"openmpi/bin\" path.", required=False)
    argParser.add_argument("-bl", "--baseline", default=1, action="count", help="Provide this flag to use a random alignment as baseline score.", required=False)

    argParser.add_argument("-gp", "--getpdb", default=0, action="count",
                           help="Map input FASTA (--alignmentfile) to PDB files. Writes per-sequence filter statistics "
                                "to <outputdir>/<name>_filter_stats.tsv and copies matched PDB files to <outputdir>/PDB/. "
                                "Requires --alignmentfile, --outputdir, --diamondfile, --diamonddbfile, --locpdbdb.",
                           required=False)
    argParser.add_argument("-t", "--test", default=0, action="count", help="Testing ", required=False)

    args = argParser.parse_args()

    if args.jobname is not None and len(args.jobname) > 11:
        errMsg = "Job title must be shorter than 12 characters. Due to DALI longer names are not possible at the moment."
        errorFct(errMsg)
        exit(1)
    
    #setting list of valid symbols and constructing regex for sequence cleaning as a global variable
    valid_symbols = sel_alphabet(args.alphabet)
    global CLEAN_SEQ_REGEX_STR
    CLEAN_SEQ_REGEX_STR = build_regex_for_seq_cleaning_whitelist(valid_symbols)
    #invalid_symbols = set_non_alphabet(args.nonalphabet)

    #not a feature
    if args.cleanfasta>0:
        clean_fasta_files_in_dir(args.msadir, args.outputdir, valid_symbols, args.verbose)

    if args.syncdb>0:
        sync_pdb_copy(args.diamondfile, args.locpdbdb, args.outputdir, args.verbose, args.backupfreq)

    #not a feature
    if args.createsearchplots > 0:
        evaluate_diamond_search_data(args.csvdir, args.outputdir, "identity", args.verbose)
    
    if args.wrapjob > 0:
        create_job(args.msadir, args.outputdir, args.jobname, args.diamondfile, args.diamonddbfile, args.verbose, args.uniqueonly,
                    args.diamondblocksize, args.threads, args.diamondtmpdir, args.baseline, args.diamondmptmpdir, args.diamondmp)

        run_job(args.outputdir, args.jobname, args.dalibindir, args.locpdbdb, args.verbose, args.threads, args.mpibin, args.samplesize,
                args.compstratname, args.compstratfile, args.continuejob, args.backupfreq)

        evaluate_job(args.outputdir, args.jobname, args.verbose)

    if args.createjob > 0:
        create_job(args.msadir, args.outputdir, args.jobname, args.diamondfile, args.diamonddbfile, args.verbose, args.uniqueonly,
                    args.diamondblocksize, args.threads, args.diamondtmpdir, args.baseline, args.diamondmptmpdir, args.diamondmp)
    if args.runjob > 0:
        run_job(args.outputdir, args.jobname, args.dalibindir, args.locpdbdb, args.verbose, args.threads, args.mpibin, args.samplesize,
                args.compstratname, args.compstratfile, args.continuejob, args.backupfreq)
    if args.evaljob > 0:
        evaluate_job(args.outputdir, args.jobname, args.verbose)
    if args.getpdb > 0:
        map_fasta_to_pdb(args.alignmentfile, args.outputdir, args.diamondfile, args.diamonddbfile,
                         args.locpdbdb, args.verbose, args.diamondblocksize, args.threads,
                         args.diamondtmpdir, args.diamondmptmpdir, args.diamondmp)
    if args.test > 0:
        #test_SPS(args.alignmentfile)
        #test_duplicate_clean_seq(args.alignmentfile)
        #differing_sequences = set()
        #differing_sequences = test_red_index_file(args.batchsearch, args.output, args.title, args.verbose, differing_sequences)
        #_ = test_red_index_file(args.batchsearch, args.output, args.title, args.verbose, differing_sequences)
        transform_alignment(args.alignmentfile, args.outputdir)

    toc = time.perf_counter()
    print(f"Done in {toc - tic:0.4f} seconds.")
    exit(0)

if __name__ == "__main__":
    main()
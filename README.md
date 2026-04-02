Feel free to read about the theory behind PSTAR and some results in the file "PSTAR-theory-results.pdf". Disclaimer: Use with caution. This is unpublished and not peer-reviewed. As such it might contain errors.

# Introduction
*PSTAR-MSA* (Pairwise STructural Alignment Reference for Multiple Sequence Alignment) is an MSA benchmarking tool that aims to generalize the comparative scoring of multiple sequence alignment (MSA) algorithms. This is tried to be achieved by using pairwise structural alignments as reference data for each pair of aligned sequences in a given MSA. Tertiary protein structure data is sourced from the *"Protein Data Bank"* (PDB) and is planned to be expanded by high-certainty protein structure predictors like *AlphaFold*. Novel measure indices like the **PSTAR-score**, which is the cumulative pairwise structure alignment equivalent of the **sum-of-pairs-score** (SPS), have been itroduced in order to evaluate the similarity between MSA algorithms.
Usually MSA algorithm quality is measured through the use of manually curated reference alignment databases such as *HOMSTRAD*, or *HomFam* for that matter. However, all of these databases suffer from limited data availability and it is not clear, if their status as references is legitimate. Because of this limited data availability it is likely that algorithms perform worse in general than they perform on the small select data set.
An in-depth explanation of the scoring and some first results can be found [<ins>here</ins>](link to pdf not available yet)

# Installation and setup
**Linux** environment required.
List of dependencies:
```
sudo apt-get install rsync
sudo apt-get install python3.6
sudo apt install python3-pip
pip install argparse DateTime pandas requests rcsbsearchapi biopython matplotlib
```
-install [DaliLite.v5](http://ekhidna2.biocenter.helsinki.fi/dali/) from http://ekhidna2.biocenter.helsinki.fi/dali/

-install [DIAMOND](https://github.com/bbuchfink/diamond) from https://github.com/bbuchfink/diamond

-download **PSTAR-MSA.py** from the master branch.

In case you're not in a hurry, you are done now. Downloading the PDB database and setting up the DIAMOND database can take up to a few hours. You don't need to create the specified directories manually, but you can. You can start the process via
```
python3 -m PSTAR-MSA --syncdb --diamondfile <path_to_DIAMOND_exe_file> --locpdbdb <dir_path_to_create_local_pdb_copy_in> --outputdir <dir_path_to_place_DIAMOND_output_db_files_in> --verbose
```
| Option   | Description | Required |
| --------- | ------- | ------- |
| --diamondfile  |  Specify the path to the DIAMOND program file. | Yes |
| --locpdbdb  | Specify path to directory where local PDB copy should be placed. | Yes |
| --outputdir | Specify a directory path to write DIAMOND-database-adjacent files to.  | Yes |
| --backupfreq | Specify an integer to control "checkpoint frequency". Default is 100. Which means that every 100 PDB files a copy of the FASTA file is made. Increase this value, if hard-drive speed is a problem. | No |
| --verbose    | Enable terminal output. | No |

##### **<ins>Waiting is recommended to ensure flawless database generation.</ins>**

However, you can download the file "/PSTAR-MSA/data/atom_diamond_db.7z", manually create a directory you can later refer to under **--outputdir** of the **--syncdb** option, then extract the two files in the zip archive into this directory and launch the above command with these two files already in place and specifying the directory name under the **--outputdir** flag. The generation of the FASTA file based on the PDB copy is then skipped. But it can happen, if entries were removed from the PDB that manual deletions in the FASTA file and manual subsequent DIAMOND database generation have to be carried out.

# Usage
There are three main parts to alignment benchmarking with PSTAR-MSA, which can be done independently of each other. But we also provide the user with a wrapper function to do everything in one call. First: **Setting up the job**. This can be done via
```
python3 -m PSTAR-MSA --createjob --msadir <dir_path_containing_alignments_to_be_compared> --outputdir <dir_path_where_to_create_job_dirs_in> --jobname <name_of_job> --diamondfile <path_to_DIAMOND_exe_file> --diamonddbfile <path_to_DIAMOND_seq_db_file> --verbose
```
| Option   | Description | Required |
| --------- | ------- | ------- |
| --msadir  | Specify a directory path containing an arbitrary number of FASTA files to be scored against each other.  | Yes |
| --outputdir | Specify a directory path to write the outputs to. For each new job there will be a subdirectory created in this directory, according to the job name.  | Yes |
| --jobname    | Specify a job name. Needs to be shorter than 12 character. At least in our environment DALI does not like job titles that exceed 11 characters for some reason.  | Yes |
| --diamondfile  |  Specify the path to the DIAMOND program file. | Yes |
| --diamonddbfile  | Specify the path to the DIAMOND database file (extension DMND). By default this can be found in the directory specified under --outputdir during the --syncdb call. See section "Installation and setup". | Yes |
| --diamondblocksize | Specify DIAMOND block size as float. Refer to the DIAMOND GitHub page for more information on the effect of this option. | No |
| --threads | Specify the number of threads used. By default this will use all available threads by reading from a environmal variable. However, if you are in a cluster environment it might be better to specify this value manually. | No |
| --diamondtmpdir | Specify the directory to use for intermediate files to be placed during DIAMOND database searches. Useful to avoid bottlenecks in a cluster environment. Refer to DIAMOND GitHub page. | No |
| --verbose | Enable terminal output. | No |

Next part: **Launching the job**
```
python3 -m PSTAR-MSA --runjob --outputdir <dir_path_where_to_create_job_dirs_in> --jobname <name_of_job> --dalibindir <path_to_DALI_bin_dir> --locpdbdb <dir_path_to_local_PDB_copy> --verbose
```
| Option   | Description | Required |
| --------- | ------- | ------- |
| --outputdir  | Specify a directory path to write the outputs to. Same directory as in --createjob.  | Yes |
| --jobname    | Specify the job name you want to run. Job has to be set up before running it.  | Yes |
| --dalibindir  |  Specify the path to the DALI "bin" directory, which contains "import.pl" and "dali.pl". | Yes |
| --locpdbdb  | Specify path to directory containing your local PDB copy that was created during the setup process. | Yes |
| --samplesize | Specify an integer to give an upper bound of pairwise alignments to be calculated. By default PSTAR-MSA calculates all possible pairwise alignments. This number can get large. | No |
| --compstratfile | Specify the path to a JSON file containing a custom comparison strategy in a python-readable format, such as lists of lists. Experimental feature, might not work. | No |
| --continuejob | Provide this flag to continue a job that had been running before, but was interrupted. | No |
| --backupfreq | Specify an integer to control "checkpoint frequency". Default is 100. Which means that a job's progress is saved every 100 pairwise structure alignments. | No |
| --verbose | Enable terminal output. | No |

Last part: **Evaluating the job**
```
python3 -m PSTAR-MSA --evaljob --outputdir <dir_path_where_finished_jobs_reside> --jobname <name_of_job_to_evaluate> --verbose 
```
There is also the flag "**--wrapjob**", which wraps all 3 of the above steps into one call. You still need to specify each required flag as before. The optional arguments are also supported. The order does not matter.

# Output
After running through the 3 steps in the last section (or by using the --wrapjob flag) you will acquire the results for each job. If you specified --verbose during the runs, you will have the results printed to stdout. Additionally PSTAR-MSA generates a lot of different outputs that might be more or less useful. For each job the outputs can be found in
```
<outputdir>/<jobname>/EVAL/
```
For each input alignment there will be a file called "<input_alignment_name>.tsv", which can be imported into a pandas dataframe and it contains a plethora of information on each pairwise alignment that was carried out including their pairwise PSTAR-score, making the contribution of each pair of sequences to the overall PSTAR-score transparent. The user is also provided with set representations for all of the alignments in the form of JSON files and also the set representation for the reference. This allows calculation of different measure indices at a later point without needing to rerun the job entirely. The two main output files however are "scores_per_alignment.tsv" and "jacc_scores_per_alignment.tsv", which simply contain a list of scores per input file. The former lists the individual PSTAR-scores, the latter the corresponding Jaccard indices. The output directory also contains two files "jacc_dm.tsv" and "psm_dm.tsv", which contain lists of pairs of set representations between the individual alignments and their respective Jaccard distances in the first file and PSTAR-metric distances in the second file. These two files are just additional information and it's not clear, whether the space of all MSAs together with either of these distance measures is a metric space at all. 











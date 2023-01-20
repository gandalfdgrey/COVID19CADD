The algorithms require the creation of Conda environment, with an installation of the latest 
RDKit package - including the installation of the Postgresql cartridge. For the purpose of
these experiments, the RDKit version 2022.03.1 was used.

The following modules are necessary for the correct running of the scripts:
psycopg2 - PostgreSQL database adapter
os - Module that allows the use of operating system dependent functionalities
utils - Module for the execution of scripts
collections - Module that implements specialized container datatypes
numpy - Python library used for working with arrays
datetime - Module that supplies classes for manipulating dates and times
pandas - Library that facilitates data analysis and data manipulation
statistics - Module that allows the calculation of mathematical statistics of numeric data

## Unzip the file "202210_MMB5010_TristanCamilleri_codenddata.zip" into the selected folder

The Disco_db reational database needs to be restored using the following command from the command line:
"psql disco_db < disco_db.sql"

The extraction of the molecules of interest from the database and filtering is performed by running the 
python script 1_file_importer.py

The extracted molecules, in the respective file, are then filtered using the python script 2_filters.py

The molecular clustering of the filtered molecules is performed by running the python script 
3_excl_clustering_tc.py on the lists of filtered molecules.

The molecular fingerprints for the molecules and respective similarity searches (using Tanimoto and 
Dice coefficients) are obtained by executing the python scripts 5_lbvs_tanimoto.py and 5a_lbvs_tanimoto.py 
on each filtered list of molecules.

The USRCAT similarity search for the molecules, including conformer generation are obtained by executing 
the python scripts 4_usrcat_conf_gen.py and 4a_usrcat_conf_gen_act.py on each filtered list of molecules.
The process for the generation of the lowest energy conformer is time intensive, with long waiting time for
the processing of an entire molecule file set.

The enrichment factors are calculated for all of the generated lists by running the python script 
6_stat_enr_factor.py on each of the generated lists.

The results are generated into the respective sub-folder in the /data sub-folder.

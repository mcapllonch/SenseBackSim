MANUSCRIPT: "Modelling the Effects of Ephaptic Coupling on Selectivity and Response Patterns during Artificial Stimulation of Peripheral Nerves", by Miguel Capllonch-Juan and Francisco Sepulveda (2020). PLOS Computational Biology.

SECTION: "Results: Effects of ephaptic coupling on axon recruitment and selectivity"

MODEL: Nerve 1

DATASET:
This dataset contains the code and the data for the subsection "Results: Effects of ephaptic coupling on axon recruitment and selectivity". In particular, this folder contains the results for the simulations with no EC.

 - Data: recruitment_data.txt and *.json
 - Code: ./code/

INSTRUCTIONS:

Running the simulations:
1. Prepare the simulations by running "create_folders.py".
2. If necessary, modify any files in the folders with "modify_files_all.sh".
3. Run all the simulations with "run_all_simulations.sh". Alternatively, a selected list of simulations can be run with "run_all_simulations_specific.sh". Note that there are two steps to run simulations without EC:
	3.1. First, replace "sim_launcher.py" by "get_extstim.py" in "run_all_simulations[...].sh" in order to get the stimulating fields from the electrodes.
	3.2. Second, restore the string "sim_launcher.py" in the bash script of choice and run it.

Data processing:
4. To collect the recruitment results from the simulation outputs, run "get_stim_results.py". Note: the json files also contain information about AP delays from stimulation onset.
5. Once this has been done both for EC and noEC, "cd .." and build "recruitment_data.csv" manually from their corresponding "recruitment_data.txt" files (i.e., "./ec/recruitment_data.txt", and "./noec/recruitment_data.txt"). Although the process can be programmed, this was done manually while processing the data for the manuscript.

Data visualization:
6. Run "fig3.py" or "fig4.py" for representations of figures 3 and 4 in the manuscript, respectively.

TROUBLESHOOTING
 - Any problems? Please search for your issue or open a new one on the software's GitHub repository: https://github.com/mcapllonch/SenseBackSim/issues
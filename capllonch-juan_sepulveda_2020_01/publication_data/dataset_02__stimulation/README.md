MANUSCRIPT: "Modelling the Effects of Ephaptic Coupling on Selectivity and Response Patterns during Artificial Stimulation of Peripheral Nerves", by Miguel Capllonch-Juan and Francisco Sepulveda (2020). PLOS Computational Biology.

SECTION: "Results: Effects of ephaptic coupling on propagation"

MODEL: Nerve 1

DATASET:
This dataset contains the code and the data for the subsection "Results: Effects of ephaptic coupling on propagation". Here are separate folders for EC and no-EC simulations.


 - Data (combined for both EC and no EC): recruitment_data.csv
 - Specific code and data for EC and no EC: "./ec/" and "./noec/", respectively.

INSTRUCTIONS:

Running the simulations:
1. Enter "./ec" and "./noec" and follow the instructions in the "README.md" in each folder to run the simulations.

Data processing:
2. Once this has been done both for EC and noEC, "cd .." and build "recruitment_data.csv" manually from their corresponding "recruitment_data.txt" files (i.e., "./ec/recruitment_data.txt", and "./noec/recruitment_data.txt"). Although the process can be programmed, this was done manually while processing the data for the manuscript.

Visualizing data:
3. Run "figs_3_5.py", or "fig4.py" for representations of figures 3 and 5, or 4 in the manuscript, respectively. "figs_3_5.py" generates Figs 3 and 5 in the manuscript, apart from other figures.

TROUBLESHOOTING
 - Any problems? Please search for your issue or open a new one on the software's GitHub repository: https://github.com/mcapllonch/SenseBackSim/issues
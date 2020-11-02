This folder contains specific data and code for the results in the article <b>Modelling the effects of ephaptic coupling on selectivity and response patterns during artificial stimulation of peripheral nerves</b>, by Miguel Capllonch-Juan and Francisco Sepulveda (2020). PLOS Computational Biology.

<b>Installation:</b>

Install the required Python dependencies:
`pip install -r requirements.txt`

If there are problems due to the use of a different Python version, try installing the dependencies without version numbers.

NEURON:
For instructions on how to download and install NEURON, visit:
https://www.neuron.yale.edu/neuron/download/
or:
https://www.neuron.yale.edu/neuron/getstd

<b>Notes:</b>

Before running any simulation as indicated in the instructions of each dataset, the MRG and Gaines & al. (2018) models need to be compiled. For this, go to the `code` folder of each dataset, and run:

`nrnivmodl`

This will create a folder, whose name will depend on your computer's architecture, containing the compiled models for NEURON to use them. This only needs to be done once on each folder and computer.
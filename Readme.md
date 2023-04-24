# AlphaBSP

In this research project a standardised AlphaFold 2-based molecular replacement strategy will be developed and implemented in an existing biomolecule structure solution pipeline at MAX IV Laboratory

## How to run

### 1. Run prediction
The first script creates a SLURM job and submits it to the MAX IV cluster. A fasta file for the protein of interest should be placed in the same directory as the 1_AF_prediction python script. The prediction is then automatically started for this protein fasta. The user will be reported when the job is done and the output file with all created predictions is completed.

To start the prediction run:
```
python 1_AF_prediction.py
```
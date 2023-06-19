# mrALFO

In this research project a standardised AlphaFold 2-based molecular replacement strategy will be developed and implemented in an existing biomolecule structure solution pipeline at MAX IV Laboratory

## Requirements
To run the pipeline following resources and requirements are necessary to produce confident prediction:

## How to run

### 1. Run prediction & molecular replacement based on AlphaFold 2

To start the automated mrALFOtask:

```
python ALFOtask <fasta_file>
```

#### 1.1 Run prediction
The first step creates a SLURM job and submits it to the MAX IV cluster. A fasta file for the protein of interest should be placed in the same directory as the 1_AF_prediction python script. The prediction is then automatically started for this protein fasta. The user will be reported when the job is done and the output file with all created predictions is completed.

#### 1.2 Make molecular replacement
From the previous AlphaFold 2 predicted structures, the model with the highest pLDDT (confidence for per-residue estimate) is used as a search model for molecular replacement.



Please refer to the following sources:

[Jumper, J et al. Highly accurate protein structure prediction with AlphaFold. Nature (2021).](https://www.nature.com/articles/s41586-021-03819-2)

[Varadi, M et al. AlphaFold Protein Structure Database: massively expanding the structural coverage of protein-sequence space with high-accuracy models. Nucleic Acids Research (2022).](https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gkab1061/6430488)

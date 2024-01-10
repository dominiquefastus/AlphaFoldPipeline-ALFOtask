#!/usr/bin/env python

""" 
Small program to quickly split a PDB file into multiple PDB files with only one chain.


Author:     D. Fastus
"""

import argparse
from Bio import PDB

arg_parser = argparse.ArgumentParser(description='Split a PDB file into multiple PDB files with only one chain.')

arg_parser.add_argument('pdb_file', type=str, help='The input PDB file.')

arg_parser.add_argument('pdb_name', type=str, help='The input PDB file.')

arg_parser.add_argument('output_dir', type=str, help='The output directory for the split PDB files.')

args = arg_parser.parse_args()


def split_pdb_chains(pdb_file, pdb_name, output_dir):
    parser = PDB.PDBParser()
    structure = parser.get_structure('PDB', pdb_file)
    io = PDB.PDBIO()

    for chain in structure.get_chains():
        io.set_structure(chain)
        chain_id = chain.get_id()
        output_file = f"{output_dir}/{pdb_name}_{chain_id}.pdb"
        io.save(output_file)

if __name__ == "__main__":
    split_pdb_chains(args.pdb_file, args.pdb_name, args.output_dir)
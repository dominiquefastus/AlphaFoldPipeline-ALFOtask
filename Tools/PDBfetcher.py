#!/usr/bin/env python

""" 
Small program to quickly fetch data from the PDB database, including:
- .pdb
- .fasta
- .mtz

All files needed to run the molecular replacement pipeline and comparison.

"""
from urllib.request import urlretrieve
from urllib.error import HTTPError

import os
import argparse
import pathlib


# construct command line interface
parser = argparse.ArgumentParser(prog="PDBfetcher.py",
                                 description="""
                                 CLI to quickly download files in fasta, pdb, mtz format
                                 of a protein using its entry id.\n\n 
                                 It's possible to download a list of proteins 
                                 and their files simultaneously. 
                                 """,
                                 epilog="Changes to PDB url should be taken into account!")

parser.add_argument("entryID", nargs='*', 
                    help="PDB entryID of the protein with the length of 4 characters")

parser.add_argument("-a", "--all", dest="all_files", action="store_true", 
                    help="download pdb, mtz and fasta file")

parser.add_argument("-p", "--pdb", dest="pdb_file", action="store_true", 
                    help="download structural data in pdb format of the protein")

parser.add_argument("-m", "--mtz", dest="mtz_file", action="store_true", 
                    help="download reflection data in mtz format of the protein")

parser.add_argument("-f", "--fasta", dest="fasta_file", action="store_true", 
                    help="download sequence data in fasta format of the protein")

parser.add_argument("-o", "--output", dest="output_path", default=os.getcwd(), type=pathlib.Path,
                    help="download sequence data in fasta format of the protein")

args = parser.parse_args()

def create_folder(name, path):
    dir_path = f'{path}/{name}_downloads'
    
    if not os.path.exists(dir_path):
        os.mkdir(f'{path}/{name}_downloads')
    
    return dir_path


for name in args.entryID:
    name = name.lower()
    
    try:
        folder_path = create_folder(name=name, path=args.output_path)
        print(f"Files for {name} will be downloaded to: {folder_path}")
    except Exception as e:
        print(f"Error processing {name}: {str(e)}")
    
    print('----------------------------')
    try:
        if args.pdb_file or args.all_files:
            if not os.path.exists(f"{folder_path}/{name}.pdb"):
                urlretrieve(f'https://files.rcsb.org/download/{name}.pdb',
                            f"{folder_path}/{name}.pdb")
                print(f'{name} pdb downloaded!')
            else:
                print(f'{name} pdb already exists!')
    except HTTPError as e:
        if e.code == 404:
            print(f'{name} pdb failed to download!')
    else:
        try:
            if args.mtz_file or args.all_files:
                if not os.path.exists(f"{folder_path}/{name}.mtz"):
                    urlretrieve(f'https://edmaps.rcsb.org/coefficients/{name}.mtz',
                                f"{folder_path}/{name}.mtz")
                    print(f'{name} mtz downloaded!')
                else:
                    print(f'{name} pdb already exists!')
        except HTTPError as e:
            if e.code == 404:
                print(f'{name} mtz failed to download!')
        else:
            try:
                if args.fasta_file or args.all_files:
                    if not os.path.exists(f"{folder_path}/{name}.fasta"):
                        urlretrieve(f'https://www.rcsb.org/fasta/entry/{name}/download', 
                                    f"{folder_path}/{name}.fasta")
                        print(f'{name} fasta downloaded!')
                    else:
                        print(f'{name} pdb already exists!')
            except HTTPError as e:
                if e.code == 404:
                    print(f'{name} pdb failed to download!')
                    
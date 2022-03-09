"""
Copyright (C) Iktos - All Rights Reserved
Unauthorized copying of this file, via any medium is strictly prohibited
Proprietary and confidential
"""

import os
from pathlib import Path
from typing import Optional

import iktos.logger as logging
from iktos.structure_utils.toolkits.rdkit import read_mol, write_mol
from rdkit.Chem import MolToSmiles

LOGGER = logging.getLogger(__name__)


def clean_directory(output_directory: str) -> None:
    """Removes all files from output directory. If the directory
    does not exist, creates it. If the directory contains directories,
    they will be ignored. If the function fails to delete one of the files,
    it will crash."""
    LOGGER.info(f'Cleaning (or creating) directory `{output_directory}`')
    directory = Path(output_directory)
    directory.mkdir(parents=True, exist_ok=True)
    file_list = directory.glob('**/*')

    for file_path in file_list:
        if file_path.is_file():
            os.remove(file_path)
        elif file_path.is_dir():
            LOGGER.debug(f'Ignoring `{file_path}` (directory)')
        else:
            LOGGER.warning(f'Path `{file_path}` is not a file nor a directory')


def get_achiral_smiles(smiles: str) -> Optional[str]:
    """Returns the corresponding achiral smiles."""
    mol = read_mol(smiles, fmt='smi')
    if mol is not None:
        return MolToSmiles(mol, isomericSmiles=False)


def get_canonical_smiles(smiles: str) -> Optional[str]:
    """Returns the corresponding canonical (RDKit) smiles."""
    mol = read_mol(smiles, fmt='smi')
    if mol is not None:
        canonical_smiles = write_mol(mol, fmt='smi')
        return canonical_smiles

import os
from typing import Dict, TypedDict

import networkx as nx
import pandas as pd
from rdkit import Chem
from rdkit.Chem.rdchem import Mol


class CompoundDict(TypedDict):
    """Dictionary containing the desired data."""

    ac_50_log: float
    graph: nx.Graph


def load_raw_data() -> pd.DataFrame:
    """Load the raw data from the csv file.

    Returns:
        df (pd.DataFrame): The raw data.

    """
    script_dir = os.path.dirname(os.path.abspath(__file__))
    file_path = os.path.join(script_dir, os.pardir, "data", "AID_450_datatable_all.csv")
    df = pd.read_csv(
        file_path,
        usecols=[
            "PUBCHEM_RESULT_TAG",
            "PUBCHEM_EXT_DATASOURCE_SMILES",
            "Qualified AC50",
        ],
        dtype={
            "PUBCHEM_RESULT_TAG": int,
            "PUBCHEM_EXT_DATASOURCE_SMILES": str,
            "Qualified AC50": float,
        },
    )
    return df


def get_molecule_by_smile(iso_smile: str) -> Mol:
    """Get the molecule object from the SMILES string.

    Args:
        iso_smile (str): The SMILES string.

    Returns:
        mol (Mol): The molecule object.

    """
    # pylint: disable=no-member
    mol = Chem.MolFromSmiles(iso_smile)
    return mol


def molecule_to_graph(mol: Mol) -> nx.Graph:
    """Convert the molecule object to a graph.

    Args:
        mol (Mol): The molecule object.

    Returns:
        graph (nx.Graph): The graph representation of the molecule.

    """
    graph = nx.Graph()
    for atom in mol.GetAtoms():  # type: ignore
        graph.add_node(
            atom.GetIdx(),
            atomic_num=atom.GetAtomicNum(),
            is_aromatic=atom.GetIsAromatic(),
        )
    for bond in mol.GetBonds():  # type: ignore
        graph.add_edge(
            bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), bond_type=bond.GetBondType()
        )
    return graph


def run_processing() -> Dict[int, CompoundDict]:
    """Run the data processing pipeline.

    Returns:
        compound_graphs (Dict[int, CompoundDict]): The processed data.

    """
    compound_graphs: Dict[int, CompoundDict] = {}
    raw_data = load_raw_data()
    raw_data.dropna(inplace=True)
    for _, row in raw_data.iterrows():
        smile = row["PUBCHEM_EXT_DATASOURCE_SMILES"]
        ac_50_log = row["Qualified AC50"]
        result_tag = row["PUBCHEM_RESULT_TAG"]
        mol = get_molecule_by_smile(smile)
        graph = molecule_to_graph(mol)
        compound_graphs[result_tag] = {"ac_50_log": ac_50_log, "graph": graph}
    return compound_graphs

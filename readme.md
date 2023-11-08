<a href="./readme.CN.md">简体中文</a> | <a href="./readme.md">English</a>

![License](https://img.shields.io/badge/license-MIT-yellowgreen)  ![Language](https://img.shields.io/badge/language-python-blue)

# ProtIon

Molecular Dynamics Simulation of Protein and Ion Interactions

## Installation

```bash
# Clone this repository using git, and then, after entering the repository directory, run the following commands

conda install -yc conda-forge mamba
mamba install -yc conda-forge rdkit openmm openmmforcefields pdbfixer mdtraj openff-toolkit numpy nglview biopython
pip install .
```

## Usage
### Adding Ions and Running Molecular Dynamics Simulation

You can refer to the ``demo/run_test.py`` script for usage examples.

```python
from protion import ProtIon

input_pdb = f"demo/inputs/folding_result.pdb"
outdir = f"demo/outputs"
# Provide target ions and their respective concentrations (unit: mol / L)
# Target ions can be specified by element name, SMILES, or with a custom name in the format "{element/SMILES}::{up to 3 characters}" separated by "::".
# For example:
ions_concentration = {
    "Ca": 8e-3, # Target ion specified by element name
    "Na": 1.6e-2,
    "Mg": 5e-2,
    "C(=O)(O)[O-]::BCR": 1.6e-2, # Target ion with SMILES "C(=O)(O)[O-]" and named "BCR"
    "Cl": 8e-3 * 2 + 5e-2 * 2,
}
pi = ProtIon(outdir)
modeller, forcefield = pi.preprocess(
    input_pdb = input_pdb, 
    ions_concentratsion = ions_concentration,
    ph = 7.
)
pi.run_md(
    modeller,
    forcefield,
    step_time: float = 4,
    total_time: float = 5e6,
    skip_time: int = 1e6,
    save_interval_time: int = 1e3,
)
```
### Ion Enrichment Capability Scoring
```python
# Assuming you have added ions as described above and completed the molecular dynamics simulation,
# we will now calculate the residue enrichment score for the Ca ion concentration.
from protion import ion_enrichment_score
import mdtraj as md

traj = md.load(
    "demo/outputs/trajectory_nowat.xtc", 
    top = "demo/outputs/_result_nowat.pdb"
)

Ca_enrichment_score = ion_enrichment_score(
    traj: md.core.trajectory.Trajectory,
    ion_resname: str,
    distance_cutoff: float = 0.4, # nm
    scheme: str = "closest",
)

```
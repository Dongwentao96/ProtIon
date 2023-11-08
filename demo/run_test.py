from pathlib import Path
from protion import ProtIon
import openmm.app as app

curdir = Path(__file__).parent
OXYGEN = app.element.oxygen

def main():
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("--run-md", action="store_true", help="Run MD simulation")
    if_run_md = parser.parse_args().run_md
    input_pdb = f"{curdir}/inputs/folding_result.pdb"
    outdir = f"{curdir}/outputs"
    ions_concentration = {
        "Ca": 8e-3,
        "Na": 1.6e-2,
        "Mg": 5e-2,
        "C(=O)(O)[O-]::BCR": 1.6e-2,
        "Cl": 8e-3 * 2 + 5e-2 * 2,
    }
    pi = ProtIon(outdir)
    modeller, forcefield = pi.preprocess(
        input_pdb = input_pdb, 
        ions_concentration = ions_concentration
    )
    
    if if_run_md:
        pi.run_md(
            modeller,
            forcefield
        )
    
if __name__ == '__main__':
    main()
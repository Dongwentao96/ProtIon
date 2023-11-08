from pathlib import Path
from protion import ProtIon
import openmm.app as app
import json

OXYGEN = app.element.oxygen

def main():
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("--input_pdb", type=str)
    parser.add_argument("--outdir", type=str)
    parser.add_argument("--ions_concentration", type=str)
    parser.add_argument("--step_time", type=float, default=4, help="Duration of each simulation step (in femtoseconds)")
    parser.add_argument("--total_time", type=float, default=5e6, help="Total simulation time (in femtoseconds)")
    parser.add_argument("--skip_time", type=float, default=1e6, help=" Time (in femtoseconds) to skip at the beginning of analysis")
    parser.add_argument("--save_interval_time", type=float, default=1e3, help="Save results every x femtoseconds")
    args = parser.parse_args()
    input_pdb = args.input_pdb
    outdir = args.outdir
    ions_concentration = json.load(open(args.ions_concentration, "r"))
    pi = ProtIon(outdir)
    modeller, forcefield = pi.preprocess(
        input_pdb = input_pdb, 
        ions_concentration = ions_concentration
    )
    pi.run_md(
        modeller,
        forcefield,
        step_time = args.step_time,
        total_time = args.total_time,
        skip_time = args.skip_time,
        save_interval_time = args.save_interval_time,
    )
    
if __name__ == '__main__':
    main()
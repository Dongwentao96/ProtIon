import numpy as np
import openmm as mm
import openmm.app as app
import openmm.unit as unit
import pandas as pd
from openff.toolkit.topology import Molecule
import os
import mdtraj as md
from pathlib import Path
from tqdm import tqdm
from openmmforcefields.generators import GAFFTemplateGenerator
from .utils import get_ion_num, get_rmsd, ion_enrichment_score

OXYGEN = app.element.oxygen


class ProtIon:
    def __init__(self, outdir: str):
        self.outdir = outdir

    def preprocess(
        self,
        input_pdb: str,
        ions_concentration: dict,
        ph: float = 7.0,
        neutralize: bool = True,
    ):
        fixed_pdb = str(Path(self.outdir, Path(input_pdb).name))
        Path(self.outdir).mkdir(parents=True, exist_ok=True)
        result_pdb_savepath = f"{self.outdir}/{Path(input_pdb).stem}_solvent_ions.pdb"

        print(
            "======Begin processing protein and ligand structures to generate topology:====="
        )
        os.system(
            f"pdbfixer {input_pdb} \
                --output={fixed_pdb} \
                --add-atoms=all \
                --keep-heterogens=none \
                --replace-nonstandard \
                --add-residues \
                --verbose \
                --ph={ph}"
        )
        forcefield = app.ForceField(
            "amber/protein.ff14SB.xml",
            "amber/tip3p_standard.xml",
            "amber/tip3p_HFE_multivalent.xml",
        )
        pdb = app.PDBFile(fixed_pdb)
        pos = pdb.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
        # center position
        pos_min = np.min(pos, axis=0)
        pos_max = np.max(pos, axis=0)
        pos_box = (pos_max + 0.2) - (pos_min - 0.2)
        box_volumn = np.prod(pos_box)
        print(f"======Box volumn: {box_volumn} nm^3.=====")
        pos_mean = pos_box / 2.0
        pos = pos - pos.mean(axis=0) + pos_mean
        modeller = app.Modeller(pdb.topology, pos * unit.nanometer)

        print("======Add solvent and ions.=====")
        # add solvent and ions
        modeller.addSolvent(
            forcefield,
            padding=0.6 * unit.nanometers,
            model="tip3p",
            neutralize=neutralize,
        )
        pdb_with_solvent = f"{self.outdir}/{Path(input_pdb).stem}_solvent.pdb"
        with open(pdb_with_solvent, "w") as f:
            app.PDBFile.writeFile(modeller.topology, modeller.positions, f)

        all_atoms = np.array(
            list(
                atom
                for chain in modeller.topology.chains()
                for residue in chain.residues()
                for atom in residue.atoms()
            )
        )
        num_total_waters = sum(
            1 for i in modeller.topology.residues() if i.name == "HOH"
        )
        print(f"======Total water number: {num_total_waters}=====")
        ions_num = {
            k: get_ion_num(v, num_total_waters) for k, v in ions_concentration.items()
        }
        print(f"======Ions number:=====\n{ions_num}")
        all_pos = np.vstack(modeller.getPositions().value_in_unit(unit.nanometer))
        all_pos_min = all_pos.min(axis=0)
        all_pos_max = all_pos.max(axis=0)
        all_pos_box = all_pos_max - all_pos_min

        def add_supported_ion(ion_name: str, ion_num: int):
            nonlocal all_pos, modeller
            num_added = 0
            pos_added = []
            new_top = app.Topology()
            new_chain = new_top.addChain()
            while num_added < ion_num:
                ion_pos = np.random.rand(3) * all_pos_box
                # check the shortest distance between ion_pos and pos
                dist = np.linalg.norm(all_pos - ion_pos, axis=1).min()
                if dist < 0.4 or any(
                    np.linalg.norm(i - ion_pos) < 0.4 for i in pos_added
                ):
                    continue
                new_res = new_top.addResidue(ion_name, new_chain)
                new_top.addAtom(ion_name, app.Element.getBySymbol(ion_name), new_res)
                pos_added.append(ion_pos)
                num_added += 1
            modeller.add(new_top, np.array(pos_added) * unit.nanometer)

        def add_molecule(ion_name: str, ion_num: int, three_letter_code: str = "UNK"):
            nonlocal all_pos, all_atoms, modeller, forcefield
            molecule = Molecule.from_smiles(ion_name, name=ion_name)
            gaff = GAFFTemplateGenerator(molecules=molecule)
            forcefield.registerTemplateGenerator(gaff.generator)
            molecule_top = molecule.to_topology().to_openmm()
            molecule_residue = list(molecule_top.residues())[0]
            molecule.generate_conformers(n_conformers=1)
            mol_conf = molecule.conformers[0]
            mol_conf.ito("nm")
            mol_pos = mol_conf.m
            lower_bound = all_pos_min - mol_pos.min(axis=0)
            upper_bound = all_pos_max - mol_pos.max(axis=0)

            new_top = app.Topology()
            new_chain = new_top.addChain()
            num_added = 0
            pos_added = []
            while num_added < ion_num:
                random_offset = np.random.uniform(lower_bound, upper_bound)
                new_mol_pos = mol_pos + random_offset
                assert ((new_mol_pos > all_pos_min) & (new_mol_pos < all_pos_max)).all()
                distance_matrix = np.sqrt(
                    np.sum((all_pos[:, np.newaxis] - new_mol_pos) ** 2, axis=2)
                )
                neighbor_atoms = all_atoms[np.where(distance_matrix < 0.4)[0]]
                neighbor_residue = list(set(i.residue for i in neighbor_atoms))
                if all(i.name == "HOH" for i in neighbor_residue):
                    modeller.delete(neighbor_residue)
                    new_res = new_top.addResidue(three_letter_code, new_chain)
                    for atom in molecule_residue.atoms():
                        new_top.addAtom(atom.name, atom.element, new_res)
                    new_atoms = list(new_res.atoms())
                    for bond in molecule_residue.bonds():
                        new_top.addBond(
                            new_atoms[bond.atom1.index],
                            new_atoms[bond.atom2.index],
                            bond.type,
                            bond.order,
                        )
                    pos_added.append(new_mol_pos)
                    num_added += 1
                    all_pos = np.vstack(
                        modeller.getPositions().value_in_unit(unit.nanometer)
                    )
                    all_atoms = np.array(
                        list(
                            atom
                            for chain in modeller.topology.chains()
                            for residue in chain.residues()
                            for atom in residue.atoms()
                        )
                    )

            modeller.add(new_top, np.vstack(pos_added) * unit.nanometer)

        def add_ion(ion_name: str, ion_num: int, three_letter_code: str = "UNK"):
            if ion_name.strip().upper() in app.Element._elements_by_symbol:
                add_supported_ion(ion_name, ion_num)
            else:
                add_molecule(ion_name, ion_num, three_letter_code)

        with tqdm(total=len(ions_num), desc="Adding ions") as pbar:
            for k, v in ions_num.items():
                pbar.set_postfix({"adding": k, "num": v})
                if "::" in k:
                    ion_name, ion_code = k.split("::")
                    add_ion(ion_name, v, ion_code)
                else:
                    add_ion(k, v)
                pbar.update(1)
            with open(result_pdb_savepath, "w") as f:
                app.PDBFile.writeFile(modeller.topology, modeller.positions, f)
        modeller.deleteWater()
        modeller.addSolvent(
            forcefield, numAdded=num_total_waters, model="tip3p", neutralize=neutralize
        )
        return modeller, forcefield

    def run_md(
        self,
        modeller: app.modeller.Modeller,
        forcefield: app.forcefield.ForceField,
        step_time: float = 4,
        total_time: float = 5e6,
        skip_time: int = 1e6,
        save_interval_time: int = 1e3,
    ):
        """_summary_

        Args:
            modeller (app.modeller.Modeller): _description_
            step_time (float, optional): Duration of each simulation step (in femtoseconds). Defaults to 4.
            total_time (float, optional): Total simulation time (in femtoseconds). Defaults to 5e6.
            skip_time (int, optional): Time (in femtoseconds) to skip at the beginning of analysis. Defaults to 1e6.
            save_interval_time (int, optional): Save results every x femtoseconds. Defaults to 1e3.
        """
        assert total_time > skip_time, "skip_time should be shorter than total_time"
        save_interval_step = int(save_interval_time / step_time)
        skip_step = int(skip_time / step_time)
        skip_frame = int(skip_step / save_interval_step)
        trajectory_savepath = f"{self.outdir}/trajectory.dcd"
        total_steps = int(total_time / step_time)
        print(
            f"total_steps: {total_steps}\nskip_step: {skip_step}\nsave_interval_step: {save_interval_step}"
        )
        water_mask = np.array(
            [i.residue.name != "HOH" for i in modeller.topology.atoms()]
        )
        traj_save_idx = np.arange(modeller.topology.getNumAtoms())[water_mask]
        system = forcefield.createSystem(
            modeller.topology,
            nonbondedMethod=app.PME,
            nonbondedCutoff=0.9 * unit.nanometers,
            constraints=app.HBonds,
            rigidWater=True,
            hydrogenMass=4 * unit.amu,
            ewaldErrorTolerance=0.0005,
        )
        print("======System build done; Start Simulation:=====")
        system.addForce(
            mm.MonteCarloBarostat(1 * unit.atmospheres, 300 * unit.kelvin, 20)
        )
        integrator = mm.LangevinMiddleIntegrator(
            300 * unit.kelvin, 1.0 / unit.picoseconds, step_time * unit.femtoseconds
        )
        # set the simulation temperature (300 K), the friction coefficient (1 ps-1), and the step size (0.004 ps)
        simulation = app.Simulation(modeller.topology, system, integrator)
        simulation.context.setPositions(modeller.positions)
        simulation.minimizeEnergy(maxIterations=1000)
        simulation.context.setVelocitiesToTemperature(300 * unit.kelvin)
        simulation.reporters.append(
            md.reporters.DCDReporter(
                trajectory_savepath, save_interval_step, atomSubset=traj_save_idx
            )
        )
        simulation.reporters.append(
            app.StateDataReporter(
                f"{self.outdir}/log.txt",
                save_interval_step,
                step=True,
                potentialEnergy=True,
                temperature=True,
                density=True,
                speed=True,
                remainingTime=True,
                totalSteps=total_steps,
            )
        )
        # every 250 time step write structure
        simulation.step(total_steps)
        print("======Simulation done; Start analysis:=====")
        modeller.deleteWater()
        no_water_savepath = f"{self.outdir}/_result_nowat.pdb"
        with open(no_water_savepath, "w") as f:
            app.PDBFile.writeFile(
                modeller.topology, modeller.positions, f, keepIds=True
            )
        traj = md.load(trajectory_savepath, top=no_water_savepath)[skip_frame:]
        traj[0].save(f"{self.outdir}/test.pdb")
        # protein index
        protein_index = np.array(
            [atom.index for atom in traj.topology.atoms if atom.residue.is_protein]
        )
        # chain A atoms
        chain_A = [
            atom for atom in traj.topology.atoms if atom.residue.chain.index == 0
        ]
        # not water index
        not_water_index = np.array(
            [
                atom.index
                for atom in traj.topology.atoms
                if atom.residue.name not in ["HOH", "WAT", "H2O"]
            ]
        )
        # sub traj
        traj = traj.atom_slice(not_water_index)
        # unwrap molecule
        traj = traj.image_molecules(anchor_molecules=[chain_A])
        # align traj
        traj = traj.superpose(traj, atom_indices=protein_index)
        # save traj
        traj.save(f"{self.outdir}/trajectory_nowat.xtc")
        traj[0].save(f"{self.outdir}/trajectory_nowat.pdb")
        rmsd_data, rmsd_time_points = get_rmsd(
            trajpath=f"{self.outdir}/trajectory_nowat.xtc",
            no_water_savepath=f"{self.outdir}/trajectory_nowat.pdb"
        )
        pd.DataFrame(
            zip(rmsd_data, rmsd_time_points), 
            columns=["rmsd", "time"]
        ).to_csv(
            f"{self.outdir}/rmsd.csv", index=False
        )
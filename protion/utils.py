import mdtraj as md
from pathlib import Path
import numpy as np


def get_ion_num(target_ion_concentration: float, h20_num: int) -> int:
    """
    Obtain the number of ions corresponding to a specific ion concentration within a certain number of water molecules.
    Args:
        target_ion_concentration (float): Concentration of the target ion (mol / L)
        h20_num (int): Number of water molecules

    Returns:
        ion_num (int): Approximate number of target ions corresponding to the target ion concentration and the number of water molecules
    """
    return int(target_ion_concentration * h20_num * 1.8e-2)


def get_rmsd(traj: md.core.trajectory.Trajectory) -> (np.ndarray, np.ndarray):
    """
    Get the RMSD for the first 100 ns of the trajectory file.

    Args:
        traj (md.core.trajectory.Trajectory): Trajectory file.
    Returns:
        rmsd_data (np.ndarray): RMSD values corresponding to each time point.
        time_points (np.ndarray): Time points corresponding to the RMSD values.
    """
    traj = traj[:1000]
    sim_time = 100  # ns
    topo = traj.topology
    backbone_atoms = topo.select("backbone")
    first_frame = traj[0]
    rmsd_data = md.rmsd(traj, first_frame, 0, atom_indices=backbone_atoms)
    time_points = np.arange(0, sim_time, sim_time / len(rmsd_data))
    return rmsd_data, time_points

def ion_enrichment_score(
    traj: md.core.trajectory.Trajectory,
    ion_resname: str,
    distance_cutoff: float = 0.4, # nm
    scheme: str = "closest",
    ) -> np.ndarray:
    """
    Calculate the enrichment score for each residue with respect to a specific ion from a trajectory file.

    Args:
        traj (md.core.trajectory.Trajectory): Trajectory file.
        ion_resname (str): Ion residue name.
        distance_cutoff (float): Distance cutoff (nm). Default: 0.4 nm.
        scheme (str): Contact scheme. Default: "closest". choices: {"ca", "closest", "closest-heavy", "sidechain", "sidechain-heavy"}. see: https://mdtraj.org/1.9.4/api/generated/mdtraj.compute_contacts.html for more details.
    Returns:
        ion_enrichment_score (np.ndarray): Ion enrichment score for each residue.
    """
    topo = traj.topology
    prot_idxs = [res.index for res in topo.chain(0).residues]
    ion_resids = []
    for res in topo.residues:
        resname = res.name
        if resname == ion_resname:
            ion_resids.append(res.index)
    if not ion_resids:
        return (np.zeros(len(prot_idxs)), 0)
    contacts = list(product(prot_idxs, ion_resids))
    distances, pairs = md.compute_contacts(
        traj=traj,
        contacts=contacts,
        scheme=scheme,
    ) if contacts else (None, contacts)
    interact_count = np.zeros(len(prot_idxs))
    for i in distances < distance_cutoff:
        for resid in pairs[i][:, 0]:
            interact_count[resid] += 1
    ion_enrichment_score = interact_count / traj.n_frames / len(ion_resids)
    return ion_enrichment_score
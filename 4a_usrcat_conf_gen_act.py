import logging as log
import os
import numpy as np
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdMolDescriptors import GetUSRScore, GetUSRCAT

# controls error logging level
log.basicConfig(level=log.DEBUG)

class NoConformersError(Exception):
    pass


def subd(dir_name):
    cwd = os.getcwd()
    path = os.path.join(cwd, dir_name)
    try:
        os.makedirs(path, exist_ok = True)
        print("Directory '%s' created successfully" % dir_name)
    except OSError as error:
        print("Directory '%s' can not be created" % dir_name)


def _get_conf_gen_params(num_threads):
    # this includes additional small ring torsion potentials and
    # macrocycle ring torsion potentials and macrocycle-specific handles
    params = AllChem.ETKDGv3()
    params.useSmallRingTorsions = True
    params.numThreads = num_threads
    params.randomSeed = 0xf00d
    return params

# Function to determine the number of conformers
def _determine_conf_num(mol):
    # citation: dx.doi.org/10.1021/ci2004658
    rot_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if rot_bonds <= 7:
        return 50
    elif rot_bonds <= 12:
        return 200
    else:
        return 300


# Function to filter and sort the conformers based on the energies
def _filter_and_sort_conformers(mol, energies, remove_unconverged):
    return_mol = Chem.Mol(mol) # due to bug #3817 we have to generate a copy of the mol
    # you cannot remove the conformers of a molecules and add them back again
    confs = list(mol.GetConformers())  # generate conformers
    if not confs:
        raise NoConformersError(f"no conformers generated for molecule {mol.GetName('_Name')}")
    energies_idx = np.argsort([x[1] for x in energies])  # sort the indices by energy values
    filtered_energies = []
    return_mol.RemoveAllConformers()
    for eidx in energies_idx:
        unconverged, e = energies[eidx]
        if remove_unconverged:
            if not unconverged:
                return_mol.AddConformer(confs[eidx], assignId=True)  # add the conformer to the molecule, ordered by energy
                filtered_energies.append(e)
        else:
            return_mol.AddConformer(confs[eidx], assignId=True)  # add the conformer to the molecule, ordered by energy
            filtered_energies.append(e)

    if remove_unconverged:
        log.debug(f"Started with {len(confs)} conformers, ended with {len(return_mol.GetConformers())}")

    if len(filtered_energies) != len(return_mol.GetConformers()):
        # a check, whereby number of filtered energies and number of conformers must match
        raise RuntimeError(f"conformers ({len(mol.GetConformers())}) and energies ({len(filtered_energies)}) lengths mismatch")

    return (return_mol, filtered_energies)


# Function to perform energy minimisation
def _energy_minimize(mol, max_tries=10, remove_unconverged=True, lowest_e=False):
    mol_name = mol.GetProp("_Name")
    conf_num = len(mol.GetConformers())
    if not conf_num:
        raise NoConformersError(f"no conformers generated for molecule {mol_name}")
    not_converged = conf_num  # at the start all confs. are not converged
    retries = 0
    energies = []
    log.debug(f"{mol_name=} {conf_num=} {not_converged=} {retries=}")
    while not_converged > 0 and retries <= max_tries:  # retry until all converged or
        energies = AllChem.MMFFOptimizeMoleculeConfs(mol)  # this returns a tuple (not_converged_flag, energy)
        not_converged = sum([e[0] for e in energies])
        log.debug(f"{not_converged} conformers not converged (try #{retries})")
        retries += 1

    mol_ordered_e, e = _filter_and_sort_conformers(mol, energies, remove_unconverged)

    if lowest_e:  # keep only one conformer, the first and lowest energy
        confs_ids = [ conf.GetId() for conf in mol_ordered_e.GetConformers()]
        for c_idx in sorted(confs_ids[1:], reverse = True): # start deleting from the biggest index
            mol_ordered_e.RemoveConformer(c_idx)
            del e[c_idx]

    return mol_ordered_e, e


# Function to perform standardisation of molecules, if enabled
def _standardise(mol):
    # removeHs, disconnect metal atoms, normalize the molecule, reionize the molecule
    clean_mol = rdMolStandardize.Cleanup(mol)

    # if many fragments, get the "parent" (the actual mol we are interested in)
    parent_clean_mol = rdMolStandardize.FragmentParent(clean_mol)

    # try to neutralize molecule
    uncharger = rdMolStandardize.Uncharger()  # annoying, but necessary as no convenience method exists
    uncharged_parent_clean_mol = uncharger.uncharge(parent_clean_mol)

    # note that no attempt is made at reionization at this step
    # nor at ionization at some pH (rdkit has no pKa caculator)
    # the main aim is to represent all molecules from different sources
    # in a (single) standard way, for use in ML, catalogue, etc.

    te = rdMolStandardize.TautomerEnumerator()
    taut_uncharged_parent_clean_mol = te.Canonicalize(uncharged_parent_clean_mol)

    return taut_uncharged_parent_clean_mol


# Function to generate conformers
def _generate_conformers_from_mol(std_mol, name=None, lowest_e_only=True, num_threads=12):

    if name:
        std_mol.SetProp("_Name", name)
    else:
        name = std_mol.GetProp("_Name")

    conf_num = _determine_conf_num(std_mol)

    # add hydrogens, critical for good conformer generation
    std_mol_h = Chem.AddHs(std_mol)

    # run ETKDG (by default)
    cids = AllChem.EmbedMultipleConfs(std_mol_h, conf_num, params=_get_conf_gen_params(num_threads))
    log.debug(f"confs generated: {len(cids)}")
    if not cids:
        raise NoConformersError(f"no conformers generated for molecule {name}")

    # return the lowest energy conformer
    mol_e_min, e = _energy_minimize(std_mol_h, lowest_e=lowest_e_only)

    return Chem.RemoveHs(mol_e_min)


def conf_gen(file, std=False):
    if std is True:
        print("Standardise is True")
        # Generate conformers
        mols = [mol for mol in Chem.SmilesMolSupplier(file, delimiter="\t", sanitize=True, titleLine=False)]
        conf_mols = []
        for i, mol in enumerate(mols):
            if mol is not None:
                try:
                   std_mol = _standardise(mol)
                   m = _generate_conformers_from_mol(std_mol)
                   conf_mols.append(m)
                except NoConformersError as e:
                    log.error(str(e))
            else:
                log.warning(f"molecule is null {i}")
        return conf_mols
    else:
        # Generate conformers
        mols = [mol for mol in Chem.SmilesMolSupplier(file, delimiter="\t", sanitize=False, titleLine=False)]
        conf_mols = []
        for i, mol in enumerate(mols):
            if mol is not None:
                try:
                   m = _generate_conformers_from_mol(mol)
                   conf_mols.append(m)
                except NoConformersError as e:
                    log.error(str(e))
            else:
                log.warning(f"molecule is null {i}")
        return conf_mols


# Function to generate USRCAT for each molecule
def usrcat_gen(mol_list):
    usrcat_list = []
    for mol in mol_list:
        try:
            nam_usr = []
            m = GetUSRCAT(mol)
            name = mol.GetProp("_Name")
            nam_usr.append(name)
            nam_usr.append(m)
            usrcat_list.append(nam_usr)
        except Exception as e:
            log.error(str(e))
    return usrcat_list


# Function to perform the similarity search of USRCAT values
def similarity_usr(act_mol, usr_list):
    sim_list = []

    for i in range(len(usr_list)):
        sim = []
        score = GetUSRScore(act_mol[1], usr_list[i][1])
        sim.append(score)
        sim.append(usr_list[i][0])
        sim.append(i)
        sim_list.append(sim)
        sim_list.sort(key=lambda x: x[0], reverse=True)
    return sim_list


valid = False
while not valid:
    try:
        choice = int(input("Press 1 for Drug-like filter, 2 for Ghose, 3 for Lipinski, 4 for QED at 0.6, 5 for REOS or 6 for Veber\nDefault is Drug-like\n"))
        valid = True
    except ValueError:
        print('Please only input digits')

file_act_list = ["act_drug_like_usrcat.txt", "act_ghose.txt", "act_lipinski_relaxed.txt", "act_qed_mid.txt", "act_reos.txt", "act_veber.txt"]
file_list = ["picks_drug_like.txt_act_23-06-22.txt", "picks_ghose.txt_act_23-06-22.txt", "picks_lipinski_relaxed.txt_act_23-06-22.txt", "picks_qed_mid.txt_act_23-06-22.txt", "picks_reos.txt_act_23-06-22.txt", "picks_veber.txt_act_23-06-22.txt"]
out = ["drug_like", "ghose", "lipinski", "qed_mid", "reos", "veber"]
out_f = ["usr_drg", "usr_ghs", "usr_lpn", "usr_qed", "usr_reo", "usr_vbr"]

fld = "data"
sbd = out_f[choice-1]  # sub-folder to output the results
os.chdir(fld)
subd(sbd)

file = (f"clst/{file_list[choice-1]}")
file_act = file_act_list[choice-1]
sim_usr = (f"{sbd}/sim_usr_act_{out[choice-1]}.txt")  # results from the similarity search
sim_usr_top = (f"{sbd}/sim_usr_top_act_{out[choice-1]}.txt")  # top 15 results per active of the similarity search
sim_usrcat = (f"{sbd}/seeded_usrcat_act_{out[choice-1]}.txt")  # lowest energy conformer for each molecule
sim_usrcat_act = (f"{sbd}/usrcat_act_{out[choice-1]}.txt")  # lowest energy conformer for each active molecule

conf_mols = conf_gen(file)
conf_act_mols = conf_gen(file_act)

usrcats_chm = usrcat_gen(conf_mols)
usrcats_act = usrcat_gen(conf_act_mols)

with open(sim_usrcat, "w") as simusr, open(sim_usrcat_act, "w") as simusract:
    for i in usrcats_chm:
        simusr.write(f"{i[0]} {' '.join('{:.5f}'.format(f) for f in i[1])}\n")
    for j in usrcats_act:
        simusract.write(f"{j[0]} {' '.join('{:.5f}'.format(f) for f in j[1])}\n")

with open(sim_usr, "w") as sim_d, open(sim_usr_top, "w") as top_d, open(file, "r") as rd_file, open(file_act, "r") as rd_act:
    src_list = rd_file.readlines()
    act_list = rd_act.readlines()
    act_ind = 0
    for i in range(len(usrcats_act)):
        sim_matrix_usr = []
        sim_matrix_usr.append(similarity_usr(usrcats_act[i], usrcats_chm))
        sim_d.write(f"< {act_list[act_ind]}")
        top_d.write(f"< {act_list[act_ind]}")
        sim_matrix_top = sim_matrix_usr[0][:15]

        for j in sim_matrix_usr:
            for k in j:
                src = src_list[k[2]].rstrip()
                sim_d.write(f"{src}\t{k[0]}\t{k[1]}\t{k[2]}\n")

        for l in sim_matrix_top:
            src = src_list[l[2]].rstrip()
            top_d.write(f"{src}\t{l[0]}\t{l[1]}\t{l[2]}\n")

        act_ind += 1

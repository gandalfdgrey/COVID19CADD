import psycopg2
import os
import rdkit
from rdkit import Chem
from utils import standardize
from rdkit.Chem import Descriptors
from rdkit.Chem import QED
from rdkit.Chem.rdMolDescriptors import GetMorganFingerprint


# Connect to an existing database
connection = psycopg2.connect(user="gandalfdgrey",
                              password="",
                              host="127.0.0.1",
                              port="5432",
                              database="disco_db")

# Function to create a sub-folder if it does not exist
def subd(dir_name):
    cwd = os.getcwd()
    path = os.path.join(cwd, dir_name)
    try:
        os.makedirs(path, exist_ok = True)
        print("Directory '%s' created successfully" % dir_name)
    except OSError as error:
        print("Directory '%s' can not be created" % dir_name)
    os.chdir(dir_name)

# Create a cursor to mols and molecule ids specific compounds
def get_mols(table):
    cursor = connection.cursor()
    cursor.execute(f"select rdkitmol, molecule_id from {table} where molecule_id in (select molecule_id from study_findings where protein_target = 'mpro')")
    print("The nunmber of molecules: ", cursor.rowcount)
    mol_list = cursor.fetchall()
    return mol_list

# Create a cursor to zinc table to retrieve molecules
def get_zinc_limit(table, limit):
    cursor = connection.cursor()
    cursor.execute(f"select rdkitmol, zinc_id from {table} limit {limit}")
    print("The nunmber of molecules: ", cursor.rowcount)
    mol_list = cursor.fetchall()
    return mol_list


def get_zinc(table):
    cursor = connection.cursor()
    cursor.execute(f"select rdkitmol, zinc_id from {table}")
    print("The nunmber of molecules: ", cursor.rowcount)
    mol_list = cursor.fetchall()
    return mol_list

# create sub-dir "data" to store output
subd("data")

actives_out = "actives.txt" # output file containing the active molecules identified from literature
lipinski = "lipinski.txt" # output file containing results from Lipinkski rule of 5 filter
lipinski_relaxed = "lipinski_relaxed.txt" # output file containing results from Lipinkski rule of 5 filter where adherance to 4 rules is mandatory
veber = "veber.txt" # output file containing results from the Veber filter
ghose = "ghose.txt" # output file containing results from the Ghose filter
ro_three = "ro3.txt" # output file containing results from the rule of 3 filter
reos = "reos.txt" # output file containing results from the REOS filter
drug_like = "drug_like.txt" # output file containing results from the drug like compounds filter
qed_top = "qed_high.txt" # QED score higher than 0.9
qed_mid = "qed_mid.txt" # QED score higher than 0.6

# Obtain the molecules in the 'molecule' table and apply filters
tab = "molecule"
actives = get_mols(tab)
actives_mol = []
for m in actives:
    n = list(m)
    mol = Chem.MolFromSmiles(n[0])  # create mol objects from SMILES
    stnd_mol = standardize.standardize(mol)  # standardise the molecules
    sanit_smiles = Chem.MolToSmiles(stnd_mol)
    mol_weight = Descriptors.ExactMolWt(stnd_mol)
    logp = Descriptors.MolLogP(stnd_mol)
    h_bond_donor = Descriptors.NumHDonors(stnd_mol)
    h_bond_acceptors = Chem.Descriptors.NumHAcceptors(stnd_mol)
    rotatable_bonds = Chem.Descriptors.NumRotatableBonds(stnd_mol)
    number_of_atoms = Chem.rdchem.Mol.GetNumAtoms(stnd_mol)
    molar_refractivity = Chem.Crippen.MolMR(stnd_mol)
    topological_surface_area_mapping = Chem.QED.properties(stnd_mol).PSA
    formal_charge = Chem.rdmolops.GetFormalCharge(stnd_mol)
    heavy_atoms = Chem.rdchem.Mol.GetNumHeavyAtoms(stnd_mol)
    num_of_rings = Chem.rdMolDescriptors.CalcNumRings(stnd_mol)
    qed_score = QED.default(stnd_mol)
    n.append(sanit_smiles)
    n.append(mol_weight)
    n.append(logp)
    n.append(h_bond_donor)
    n.append(h_bond_acceptors)
    n.append(rotatable_bonds)
    n.append(number_of_atoms)
    n.append(molar_refractivity)
    n.append(topological_surface_area_mapping)
    n.append(formal_charge)
    n.append(heavy_atoms)
    n.append(num_of_rings)
    n.append(qed_score)
    actives_mol.append(n)


# Write the active molecules to actives.txt
with open(actives_out, "w") as act:
    for l in actives_mol:
        act.write(f"{l[2]}\t{l[1]}\n")

lip = []
lip_r = []
veb = []
gho = []
rot = []
reo = []
drg = []
qt = []
qm = []

tab2 = 'zinc'
compare = get_zinc(tab2)
zinc_all = actives_mol.copy()
zinc = []

err_mols = 0

for m in compare:
    n = list(m)
    mol = Chem.MolFromSmiles(n[0])
    err = Chem.SanitizeMol(mol, catchErrors=True)
    if err != Chem.SanitizeFlags.SANITIZE_NONE:  # reports identified errors
        err_mols += 1
        continue
    try:
        stnd_mol_zinc = standardize.standardize(mol)
        sanit_smiles = Chem.MolToSmiles(stnd_mol_zinc)
        mol_weight = Descriptors.ExactMolWt(stnd_mol_zinc)
        logp = Descriptors.MolLogP(stnd_mol_zinc)
        h_bond_donor = Descriptors.NumHDonors(stnd_mol_zinc)
        h_bond_acceptors = Chem.Descriptors.NumHAcceptors(stnd_mol_zinc)
        rotatable_bonds = Chem.Descriptors.NumRotatableBonds(stnd_mol_zinc)
        number_of_atoms = Chem.rdchem.Mol.GetNumAtoms(stnd_mol_zinc)
        molar_refractivity = Chem.Crippen.MolMR(stnd_mol_zinc)
        topological_surface_area_mapping = Chem.QED.properties(stnd_mol_zinc).PSA
        formal_charge = Chem.rdmolops.GetFormalCharge(stnd_mol_zinc)
        heavy_atoms = Chem.rdchem.Mol.GetNumHeavyAtoms(stnd_mol_zinc)
        num_of_rings = Chem.rdMolDescriptors.CalcNumRings(stnd_mol_zinc)
        qed_score = QED.default(stnd_mol_zinc)
        n.append(sanit_smiles)
        n.append(mol_weight)
        n.append(logp)
        n.append(h_bond_donor)
        n.append(h_bond_acceptors)
        n.append(rotatable_bonds)
        n.append(number_of_atoms)
        n.append(molar_refractivity)
        n.append(topological_surface_area_mapping)
        n.append(formal_charge)
        n.append(heavy_atoms)
        n.append(num_of_rings)
        n.append(qed_score)
        zinc.append(n)
        zinc_all.append(n)

        # Lipinski rule of five - The Rule of 5 filter was specified by Chris Lipinski, and
        # is often referred to as the “Lipinski filter”. It is a drug-like filter , and is
        # described in Lipinski et al., Adv. Drug Deliv. Rev., 23, 3 (2001)
        if mol_weight <= 500 and logp <= 5 and h_bond_donor <= 5 and h_bond_acceptors <= 10 and rotatable_bonds <= 5:
            lp=[]
            lp.append(n[2])
            lp.append(n[1])
            lip.append(lp)

        # Lipinski rule of five - The Rule of 5 filter was specified by Chris Lipinski, and
        # is often referred to as the “Lipinski filter”. It is a drug-like filter , and is
        # described in Lipinski et al., Adv. Drug Deliv. Rev., 23, 3 (2001)
        score = 0
        if mol_weight > 500:
            score += 1
        if logp > 5:
            score += 1
        if h_bond_donor > 5:
            score += 1
        if h_bond_acceptors > 10:
            score += 1
        if rotatable_bonds > 5:
            score += 1
        if score < 2:
            lpr = []
            lpr.append(n[2])
            lpr.append(n[1])
            lip_r.append(lpr)

        # The Veber filter is a rule of thumb filter for orally active drugs described in
        # Veber et. al., J Med Chem. 2002; 45(12): 2615-23.
        if rotatable_bonds <= 10 and topological_surface_area_mapping <= 140:
            vb = []
            vb.append(n[2])
            vb.append(n[1])
            veb.append(vb)

        # The Ghose filter is a drug-like filter described in Ghose et al., J. Comb.
        # Chem., 1, 55, (1999).
        if mol_weight >= 160 and mol_weight <= 480 and logp >= float(0 - 0.4) and logp <= 5.6 and number_of_atoms >= 20 and number_of_atoms <= 70 and molar_refractivity >= 40 and molar_refractivity <= 130:
            gh = []
            gh.append(n[2])
            gh.append(n[1])
            gho.append(gh)

        # The rule of 3 filter was devised by Astex for fragment libraries as an equivalent
        # to the Lipinski rule of 5 filter for drug like libraries, and is described in
        # Congreve et al., Drug Discov. Today. 8 (19): 876–7, (2003).
        if mol_weight <= 300 and logp <= 3 and h_bond_donor <= 3 and h_bond_acceptors <= 3 and rotatable_bonds <= 3:
            rt = []
            rt.append(n[2])
            rt.append(n[1])
            rot.append(rt)

        # Rapid Elimination Of Swill filter - The REOS filter is a filter designed to
        # filter out unuseful compounds from HTS screening results, and is described in
        # Waters & Namchuk, Nature Reviews Drug Discovery 2, 259-266 (2003).
        if mol_weight >= 200 and mol_weight <= 500 and logp >= int(0-5) and logp <= 5 and h_bond_donor >= 0 and h_bond_donor <= 5 and h_bond_acceptors >= 0 and h_bond_acceptors <= 10 and formal_charge >= int(0-2) and formal_charge <= 2 and rotatable_bonds >= 0 and rotatable_bonds <= 8 and heavy_atoms >= 15 and heavy_atoms <= 50:
            res = []
            res.append(n[2])
            res.append(n[1])
            reo.append(res)

        # Drug like filter
        if mol_weight < 400 and num_of_rings > 0 and rotatable_bonds < 5 and h_bond_donor <= 5 and h_bond_acceptors <= 10 and logp < 5:
            dg = []
            dg.append(n[2])
            dg.append(n[1])
            drg.append(dg)

        # QED - Quantitative Estimate of Drug-likeness - Bickerton, G.R., Paolini, G.V.,
        # Besnard, J., Muresan, S. and Hopkins, A.L., 2012. Quantifying the chemical beauty
        # of drugs. Nature chemistry, 4(2), pp.90-98.
        if qed_score >= 0.9:
            qd = []
            qd.append(n[2])
            qd.append(n[1])
            qt.append(qd)

        if qed_score >= 0.6:
            qd = []
            qd.append(n[2])
            qd.append(n[1])
            qm.append(qd)
    except:
        pass

# Output the results of the ZINC table into individual .txt files for the different filters
with open(qed_top, "w") as q_t, open(qed_mid, "w") as q_m, open(lipinski, "w") as lpsk, open(lipinski_relaxed, "w") as lpsk_r, open(veber, "w") as vebr, open(ghose, "w") as ghos, open(ro_three, "w") as ro_t, open(reos, "w") as re_os, open(drug_like, "w") as drg_l:
    for l in lip:
        lpsk.write(f"{l[0]}\t{l[1]}\n")

    for l in lip_r:
        lpsk_r.write(f"{l[0]}\t{l[1]}\n")

    for l in veb:
        vebr.write(f"{l[0]}\t{l[1]}\n")

    for l in gho:
        ghos.write(f"{l[0]}\t{l[1]}\n")

    for l in rot:
        ro_t.write(f"{l[0]}\t{l[1]}\n")

    for l in reo:
        re_os.write(f"{l[0]}\t{l[1]}\n")

    for l in drg:
        drg_l.write(f"{l[0]}\t{l[1]}\n")

    for l in qt:
        q_t.write(f"{l[0]}\t{l[1]}\n")

    for l in qm:
        q_m.write(f"{l[0]}\t{l[1]}\n")

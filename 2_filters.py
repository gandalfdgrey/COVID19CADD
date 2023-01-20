import os
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import QED

master_key = "actives.txt"
lipinski = "act_lipinski.txt" # output file containing results from Lipinkski rule of 5 filter
lipinski_relaxed = "act_lipinski_relaxed.txt" # output file containing results from Lipinkski rule of 5 filter where adherance to 4 rules is mandatory
veber = "act_veber.txt" # output file containing results from the Veber filter
ghose = "act_ghose.txt" # output file containing results from the Ghose filter
ro_three = "act_ro3.txt" # output file containing results from the rule of 3 filter
reos = "act_reos.txt" # output file containing results from the REOS filter
drug_like = "act_drug_like.txt" # output file containing results from the drug like compounds filter
qed_top = "act_qed_high.txt" # QED score higher than 0.9
qed_mid = "act_qed_mid.txt" # QED score higher than 0.6

os.chdir("data")

# Read the actives.txt file to apply drug filters, and output the results for each filter into separate files
with open(master_key, "r") as master, open(qed_top, "w") as qt, open(qed_mid, "w") as qm, open(lipinski, "w") as lip, open(lipinski_relaxed, "w") as lip_r, open(veber, "w") as veb, open(ghose, "w") as ghos, open(ro_three, "w") as rot, open(reos, "w") as reo, open(drug_like, "w") as drg:
    for index, line in enumerate(master):
        line = master.readline().strip().split("\t")
        mol = Chem.MolFromSmiles(line[0])
        mol_weight = Descriptors.ExactMolWt(mol)
        logp = Descriptors.MolLogP(mol)
        h_bond_donor = Descriptors.NumHDonors(mol)
        h_bond_acceptors = Chem.Descriptors.NumHAcceptors(mol)
        rotatable_bonds = Chem.Descriptors.NumRotatableBonds(mol)
        number_of_atoms = Chem.rdchem.Mol.GetNumAtoms(mol)
        molar_refractivity = Chem.Crippen.MolMR(mol)
        topological_surface_area_mapping = Chem.QED.properties(mol).PSA
        formal_charge = Chem.rdmolops.GetFormalCharge(mol)
        heavy_atoms = Chem.rdchem.Mol.GetNumHeavyAtoms(mol)
        num_of_rings = Chem.rdMolDescriptors.CalcNumRings(mol)
        qed_score = QED.default(mol)

        # Lipinski rule of five - The Rule of 5 filter was specified by Chris Lipinski, and
        # is often referred to as the “Lipinski filter”. It is a drug-like filter , and is
        # described in Lipinski et al., Adv. Drug Deliv. Rev., 23, 3 (2001)
        if mol_weight <= 500 and logp <= 5 and h_bond_donor <= 5 and h_bond_acceptors <= 10 and rotatable_bonds <= 5:
            lip.write(f"{line[0]}\t{line[1]}\n")

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
            lip_r.write(f"{line[0]}\t{line[1]}\n")

        # The Veber filter is a rule of thumb filter for orally active drugs described in
        # Veber et. al., J Med Chem. 2002; 45(12): 2615-23.
        if rotatable_bonds <= 10 and topological_surface_area_mapping <= 140:
            veb.write(f"{line[0]}\t{line[1]}\n")

        # The Ghose filter is a drug-like filter described in Ghose et al., J. Comb.
        # Chem., 1, 55, (1999).
        if mol_weight >= 160 and mol_weight <= 480 and logp >= float(0 - 0.4) and logp <= 5.6 and number_of_atoms >= 20 and number_of_atoms <= 70 and molar_refractivity >= 40 and molar_refractivity <= 130:
            ghos.write(f"{line[0]}\t{line[1]}\n")

        # The rule of 3 filter was devised by Astex for fragment libraries as an equivalent
        # to the Lipinski rule of 5 filter for drug like libraries, and is described in
        # Congreve et al., Drug Discov. Today. 8 (19): 876–7, (2003).
        if mol_weight <= 300 and logp <= 3 and h_bond_donor <= 3 and h_bond_acceptors <= 3 and rotatable_bonds <= 3:
            rot.write(f"{line[0]}\t{line[1]}\n")

        # Rapid Elimination Of Swill filter - The REOS filter is a filter designed to
        # filter out unuseful compounds from HTS screening results, and is described in
        # Waters & Namchuk, Nature Reviews Drug Discovery 2, 259-266 (2003).
        if mol_weight >= 200 and mol_weight <= 500 and logp >= int(0-5) and logp <= 5 and h_bond_donor >= 0 and h_bond_donor <= 5 and h_bond_acceptors >= 0 and h_bond_acceptors <= 10 and formal_charge >= int(0-2) and formal_charge <= 2 and rotatable_bonds >= 0 and rotatable_bonds <= 8 and heavy_atoms >= 15 and heavy_atoms <= 50:
            reo.write(f"{line[0]}\t{line[1]}\n")

        # Drug like filter
        if mol_weight < 400 and num_of_rings > 0 and rotatable_bonds < 5 and h_bond_donor <= 5 and h_bond_acceptors <= 10 and logp < 5:
            drg.write(f"{line[0]}\t{line[1]}\n")

        # QED - Quantitative Estimate of Drug-likeness - Bickerton, G.R., Paolini, G.V.,
        # Besnard, J., Muresan, S. and Hopkins, A.L., 2012. Quantifying the chemical beauty
        # of drugs. Nature chemistry, 4(2), pp.90-98.
        if qed_score >= 0.9:
            qt.write(f"{line[0]}\t{line[1]}\n")
        if qed_score >= 0.6:
            qm.write(f"{line[0]}\t{line[1]}\n")

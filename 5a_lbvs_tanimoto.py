import os
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit import DataStructs


# Function to create directory
def subd(dir_name):
    cwd = os.getcwd()
    path = os.path.join(cwd, dir_name)
    try:
        os.makedirs(path, exist_ok=True)
        print("Directory '%s' created successfully" % dir_name)
    except OSError as error:
        print("Directory '%s' can not be created" % dir_name)


# Function to obtain the Dice similarity coefficient between an active and a molecule
def fp_similarity_dice(act_fp, fp_list):
    sim_list = []

    for i, v in enumerate(fp_list):
        sim = []
        dice = DataStructs.DiceSimilarity(act_fp, v)
        sim.append(dice)
        sim.append(i)
        sim_list.append(sim)
        sim_list.sort(key=lambda x: x[0], reverse=True)
    return sim_list


# Function to obtain the Tanimoto similarity coefficient between an active and a molecule
def fp_similarity_tani(act_fp, fp_list):
    sim_list = []

    for i, v in enumerate(fp_list):
        sim = []
        tani = DataStructs.TanimotoSimilarity(act_fp, v)
        sim.append(tani)
        sim.append(i)
        sim_list.append(sim)
        sim_list.sort(key=lambda x: x[0], reverse=True)
    return sim_list


valid = False
while not valid:
    try:
        choice = int(input("Press 1 for Drug-like filter, 2 for Ghose, 3 for Lipinski, 4 for QED at 0.6, 5 for REOS or 6 for Veber\n"))
        valid = True
    except ValueError:
        print('Please only input digits')

file_act_list = ["act_drug_like.txt", "act_ghose.txt", "act_lipinski_relaxed.txt", "act_qed_mid.txt", "act_reos.txt", "act_veber.txt"]
file_list = ["picks_drug_like.txt_act_23-06-22.txt", "picks_ghose.txt_act_23-06-22.txt", "picks_lipinski_relaxed.txt_act_23-06-22.txt", "picks_qed_mid.txt_act_23-06-22.txt", "picks_reos.txt_act_23-06-22.txt", "picks_veber.txt_act_23-06-22.txt"]
out = ["drg_act_sim", "ghs_act_sim", "lpn_act_sim", "qed_act_sim", "reos_act_sim", "vbr_act_sim"]
out_f = ["sim_drg", "sim_ghs", "sim_lpn", "sim_qed", "sim_reos", "sim_vbr"]

fld = "data"
sbd = out_f[choice-1]  # sub-folder to output the results
os.chdir(fld)
subd(sbd)

file = (f"clst/{file_list[choice-1]}")
file_act = file_act_list[choice-1]
sim_act_tani = (f"{sbd}/{out[choice-1]}_act_tani.txt")  # Tanimoto similarity coefficient between respective active and filtered molecules
sim_act_dice = (f"{sbd}/{out[choice-1]}_act_dice.txt")  # Dice similarity coefficient between respective active and filtered molecules
sim_dice_top = (f"{sbd}/{out[choice-1]}_dice_top.txt")  # Top 15 Dice similarity coefficients between respective active and filtered molecules
sim_tani_top = (f"{sbd}/{out[choice-1]}_tani_top.txt")  # Top 15 Tanimoto similarity coefficient between respective active and filtered molecules

# Generate the morgan fingerprints from SMILES
mol = Chem.SmilesMolSupplier(file, delimiter='\t', titleLine=False)
fps = [rdMolDescriptors.GetMorganFingerprintAsBitVect(m, 2, 2048) for m in mol if m is not None]

act_mol = Chem.SmilesMolSupplier(file_act, delimiter='\t', titleLine=False)
act_fps = [rdMolDescriptors.GetMorganFingerprintAsBitVect(m, 2, 2048) for m in act_mol if m is not None]

with open(sim_act_dice, "w") as sim_d, open(sim_dice_top, "w") as top_d, open(file, "r") as rd_file, open(file_act, "r") as rd_act:
    src_list = rd_file.readlines()
    act_list = rd_act.readlines()

    for i, v in enumerate(act_list):
        sim_matrix_dice = []
        sim_matrix_dice.append(fp_similarity_dice(act_fps[i], fps))
        sim_d.write(f"< {act_list[i]}")
        top_d.write(f"< {act_list[i]}")
        sim_matrix_dice_top = sim_matrix_dice[0][:15]
        print(f"Dice {((i+1)/len(act_list))*100}% completion")

        for j in sim_matrix_dice:
            for k in j:
                src = src_list[k[1]].rstrip()
                sim_d.write(f"{src}\t{k}\n")

        for l in sim_matrix_dice_top:
            src = src_list[l[1]].rstrip()
            top_d.write(f"{src}\t{l}\n")

sim_matrices_tani = []
with open(sim_act_tani, "w") as sim_t, open(sim_tani_top, "w") as top_t, open(file, "r") as rd_file, open(file_act, "r") as rd_act:
    src_list = rd_file.readlines()
    act_list = rd_act.readlines()
    for i, v in enumerate(act_list):
        sim_matrix_tani = []
        sim_matrix_tani.append(fp_similarity_tani(act_fps[i], fps))
        sim_matrices_tani.append(sim_matrix_tani)
        sim_t.write(f"< {act_list[i]}")
        top_t.write(f"< {act_list[i]}")
        sim_matrix_tani_top = sim_matrix_tani[0][:15]
        print(f"Tanimoto {((i+1)/len(act_list))*100}% completion")

        for j in sim_matrix_tani:
            for k in j:
                src = src_list[k[1]].rstrip()
                sim_t.write(f"{src}\t{k}\n")

        for l in sim_matrix_tani_top:
            src = src_list[l[1]].rstrip()
            top_t.write(f"{src}\t{l}\n")

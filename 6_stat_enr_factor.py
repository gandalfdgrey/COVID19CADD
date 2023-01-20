import pandas as pd
from statistics import mean
import os


def subd(dir_name):
    cwd = os.getcwd()
    path = os.path.join(cwd, dir_name)
    try:
        os.makedirs(path, exist_ok=True)
        print("Directory '%s' created successfully" % dir_name)
    except OSError as error:
        print("Directory '%s' can not be created" % dir_name)


# Function to split a string and remove header symbol
def _splitter(lst):
    fin_list = []
    ind_list = []
    for ind, i in enumerate(lst):
        if "< " in i:
            ind_list.append(ind)
        j = i.rstrip()
        k = j.split("\t")
        fin_list.append(k)
    return fin_list, ind_list


def _string_list(str):
    chr_rem = "[]"
    for chr in chr_rem:
        str = str.replace(chr, "")
    x = str.split(", ")

    return x


def df_gen(fil):
    with open(fil, "r") as fle:
        fle_list = fle.readlines()
        fle_split = _splitter(fle_list)

    # create an Empty DataFrame object
    df_mol = pd.DataFrame()

    for ind, i in enumerate(fle_split[1]):
        if ind+1 < len(fle_split[1]):
            m = []
            while i < fle_split[1][ind+1]:
                if len(fle_split[0][i]) <= 2:
                    pass
                elif len(fle_split[0][i]) == 3:
                    x = _string_list(fle_split[0][i][1])
                    m.append(x[0])
                elif len(fle_split[0][i]) > 3:
                    x = _string_list(fle_split[0][i][1])
                    m.append(x[0])
                i += 1
            df_mol[ind+1] = m

        else:
            n = []
            j = fle_split[1][-1]

            while j < len(fle_split[0]):
                if len(fle_split[0][j]) <= 2:
                    pass
                elif len(fle_split[0][j]) == 3:
                    y = _string_list(fle_split[0][j][1])
                    n.append(y[0])
                elif len(fle_split[0][j]) > 3:
                    y = _string_list((fle_split[0][j][1]))
                    n.append(y[0])
                j += 1
            df_mol[len(fle_split[1])] = n

    return df_mol.astype(str)


# Function to select the given percentage cutoff
def sample_first_prows(data, perc):

    return data.head(int(len(data)*(perc)))


# Function to calculate the enrichment factor using the provided percentage cutoff
def enrichment_factor(df, per):

    df_cnt = pd.DataFrame()

    for act in act_list:
        df_cnt[f"{act}_found"] = (df[act].isin(act_list).astype(int))
    c_all = (len(df_cnt.index))

    ef = []
    for act in act_list:
        a_all = df_cnt[f"{act}_found"].sum(axis=0, skipna=True)
        sr = sample_first_prows(df_cnt[f"{act}_found"], per)
        cx = len(sr)
        ax = sr.sum()
        en_fa = (ax/cx)/(a_all/c_all)
        ef.append(en_fa)

    return ef


valid = False
while not valid:
    try:
        choice_one = int(input("Press 1 for Fingerprint Similarity Search, Press 2 for USRCAT\n"))
        if choice_one == 1:
            choice_two = int(input("Press 1 for Drug-like filter, 2 for Ghose, 3 for Lipinski, 4 for QED at 0.6, 5 for REOS or 6 for Veber\n"))
        elif choice_one == 2:
            choice_two = int(input("Press 1 for Drug-like filter. Only the drug-like filter was used for USRCAT\n"))
        per_cutoff = float(input("Provide a percentage cutoff to calculate the enrichment factor:\n"))
        valid = True
    except ValueError:
        print('Please only input digits')


file_act = [["drg_act_sim_act_tani.txt", "ghs_act_sim_act_tani.txt", "lpn_act_sim_act_tani.txt", "qed_act_sim_act_tani.txt", "reos_act_sim_act_tani.txt", "vbr_act_sim_act_tani.txt"],["sim_usr_act_drug_like.txt"]]
act = [["act_drug_like.txt", "act_ghose.txt", "act_lipinski_relaxed.txt", "act_qed_mid.txt", "act_reos.txt", "act_veber.txt"], ["act_drug_like_usrcat.txt"]]
out = ["drg_sim", "ghs_sim", "lpn_sim", "qed_sim", "reos_sim", "vbr_sim"]
out_f = [["sim_drg", "sim_ghs", "sim_lpn", "sim_qed", "sim_reos", "sim_vbr"], ["usr_drg"]]

fld = "data"

sbd = f"{out_f[choice_one-1][choice_two-1]}/stats/"  # sub-folder to output the results
os.chdir(fld)
subd(sbd)

cutoff = per_cutoff/100

file_act = (f"{out_f[choice_one-1][choice_two-1]}/{file_act[choice_one-1][choice_two-1]}")
act = (f"{act[choice_one - 1][choice_two - 1]}")
nw = (f"{out_f[choice_one-1][choice_two-1]}/stats/enrichment_factor_{cutoff}.txt")

df_mol = df_gen(file_act)

with open(act, "r") as fle:
    fle_list = fle.readlines()
    act_list = []
    for line in fle_list:
        l = line.split("\t")
        k = l[1].rstrip()
        act_list.append(k)

df = df_mol.iloc[1:, :]
df.columns = act_list

z = enrichment_factor(df, cutoff)  # calculate the enrichment factor for the provided cutoff

print(f"The individual enrichment factors are:\n{z}\n")
print(f"The mean enrichment factor is {round(mean(z), 2)}")

with open(nw, "w") as ef_out:
    ef_out.write(f"The individual enrichment factors are:\n{z}\n")
    ef_out.write(f"The mean enrichment factor is {round(mean(z), 2)}")

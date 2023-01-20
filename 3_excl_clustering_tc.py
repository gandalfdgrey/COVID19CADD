import os
import collections
from rdkit.Chem import rdDepictor, rdMolDescriptors
from rdkit.SimDivFilters import rdSimDivPickers
from rdkit import Chem
from rdkit import DataStructs
import numpy as np
import datetime

# Function to assign compounds to the nearest cluster centroid
def assignPointsToClusters(picks, fp, fps):

    sims = np.zeros(len(picks))
    for i in range(len(picks)):
        pick = picks[i]
        if fps[pick] == fp:
            pass
        else:
            sims[i] = DataStructs.TanimotoSimilarity(fps[pick], fp)
    best = np.argmax(sims)

    pybest = best.item()
    result = picks[pybest]
    return result


def subd(dir_name):
    cwd = os.getcwd()
    path = os.path.join(cwd, dir_name)
    try:
        os.makedirs(path, exist_ok=True)
        print("Directory '%s' created successfully" % dir_name)
    except OSError as error:
        print("Directory '%s' can not be created" % dir_name)


fld = "data"
sbd = "clst"  # sub-folder containing the clustering data
os.chdir(fld)  # change directory to the folder containing the separate filtered ZINC molecule files
subd(sbd)

today = datetime.datetime.now()
date_time = today.strftime("%d-%m-%y")

file_list = ["drug_like.txt", "ghose.txt", "lipinski_relaxed.txt", "qed_high.txt", "qed_mid.txt", "reos.txt", "veber.txt"]
valid = False
while not valid:
    try:
        choice = int(input("Press 1 for Drug-like filter, 2 for Ghose, 3 for Lipinski, 4 for QED at 0.9, 5 for QED at 0.6, 6 for REOS or 7 for Veber"))
        valid = True
    except ValueError:
        print('Please only input digits')

file = file_list[choice-1]
clst = (f"{sbd}/clusters_{file}_{date_time}.txt")  # Output file for the molecules in each cluster, with the filename containing the origin file and date and time stamp
pks = (f"{sbd}/picks_{file}_{date_time}.txt")  # Output file for the cluster centroids, with the filename containing the origin file and date and time stamp

# Generate the morgan fingerprints from SMILES
mol = Chem.SmilesMolSupplier(file, delimiter='\t', titleLine=False)
fps = [rdMolDescriptors.GetMorganFingerprintAsBitVect(m, 2, 2048) for m in mol if m is not None]

# Use the sphere exclusion function provided by RD-Kit to identify diverse compounds and use the cutoff as a
# minimum distance allowed between the compounds picked and to return a set of compounds satisfying that constraint
lp = rdSimDivPickers.LeaderPicker()
thresh = 0.65  # minimum distance between cluster centroids (similarity threshold)
picks = lp.LazyBitVectorPick(fps, len(fps), thresh)
print(f"The number of picks is {len(picks)}")

lines = []
clusters = collections.defaultdict(list)
for ind, m in enumerate(picks):
    clusters[ind].append(m)

with open(file, "r") as r:
    for line in r:
        ln = line.rstrip()
        lines.append(ln)

with open(pks, "w") as pk:
    for k, v in clusters.items():
        pk.write(f"{lines[v[0]]}\tC{k} {v[0]}\n")

    for i in range(len(fps)):
        if i in picks:
            continue
        else:
            res = assignPointsToClusters(picks, fps[i], fps)
            for key, value in clusters.items():
                if res in value:
                    clusters[key].append(i)

with open(clst, "w") as clt:
    for k, v in clusters.items():
        for m in v:
            clt.write(f"{lines[m]}\tC{k} {m}\n")

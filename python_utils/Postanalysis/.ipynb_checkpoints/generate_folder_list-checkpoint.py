
PATH = "/scratch/flore0a/Preproposal/Scalability/Test1"
# importa tu función real
from utils import Export_folders_from_results

dirs = Export_folders_from_results(PATH)
dirs = [d + "/Simulation2D" for d in dirs][1:]   # mismo slicing que tú

with open("dirs.txt", "w") as f:
    for d in dirs:
        f.write(d + "\n")

print(f"Wrote {len(dirs)} paths to dirs.txt")
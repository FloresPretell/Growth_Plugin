#!/usr/bin/env python3
import argparse
import os

from utils import Export_folders_from_results

def main():
    # inputs
    ap = argparse.ArgumentParser()
    ap.add_argument("--path", required=True, help="Folder where all sims are located")
    args = ap.parse_args()

    PATH = os.path.abspath(args.path)

    # function 
    dirs = Export_folders_from_results(PATH)
    dirs = [os.path.join(d, "Simulation2D/merged") for d in dirs]  # keep your [1:] rule

    with open("dirs.txt", "w") as f:
        for d in dirs:
            f.write(d + "\n")

    print(f"Wrote {len(dirs)} paths to dirs.txt")

if __name__ == "__main__":
    main()
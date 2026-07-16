#!/usr/bin/env python3
import argparse
import os

from utils import Export_folders_from_results


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--path", required=True, help="Folder where all sims are located")
    args = ap.parse_args()

    PATH = os.path.abspath(args.path)

    dirs = Export_folders_from_results(PATH)

    # Keep your [1:] rule, but do NOT append Simulation2D
    dirs = dirs[1:]

    with open("dirs.txt", "w") as f:
        for d in dirs:
            f.write(d + "\n")

    print(f"Wrote {len(dirs)} case paths to dirs.txt")


if __name__ == "__main__":
    main()
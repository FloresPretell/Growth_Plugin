#!/usr/bin/env python3
import argparse

from util_postanalisis_export_png import export_one

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--simdir", required=True)
    ap.add_argument("--variable_input", default="t")
    args = ap.parse_args()

    export_one(args.simdir, variable_input=args.variable_input)

if __name__ == "__main__":
    main()
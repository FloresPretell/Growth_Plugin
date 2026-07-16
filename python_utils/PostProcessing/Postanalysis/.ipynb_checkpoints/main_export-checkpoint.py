#!/usr/bin/env python3
import argparse

from util_postanalisis_export_png import export_one

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--simdir", required=True)
    ap.add_argument("--fields_to_export", nargs="+", default=["u", "p", "t", "b", "ca_cyt"])

    args = ap.parse_args()

    export_one(args.simdir, fields_to_export=args.fields_to_export)

if __name__ == "__main__":
    main()
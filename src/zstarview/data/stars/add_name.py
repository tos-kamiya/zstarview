#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Add IAU proper names to stars_hip.csv based on HIP or Bayer match.
Input: 
  - stars_hip.csv (from hip_to_stars_csv.py, with HIP and optionally BayerASCII col)
  - IAU-Catalog of Star Names CSV from exopla.net
Output:
  - stars_hip_named.csv (Name filled if IAU match found)
"""

import argparse
import csv

def load_iau(file):
    hip_map = {}
    bayer_map = {}
    with open(file, newline='', encoding='utf-8') as f:
        r = csv.DictReader(f)
        for row in r:
            name = row["proper names"].strip()
            if not name:
                continue
            hip = row["HIP"].strip()
            bayer = row["Bayer ID"].strip()
            if hip.isdigit():
                hip_map[int(hip)] = name
            if bayer:
                bayer_map[bayer] = name
    return hip_map, bayer_map

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stars", default="stars_hip.csv", help="stars_hip.csv from Hipparcos")
    ap.add_argument("--iau", default="IAU-Catalog of Star Names (always up to date).csv", help="IAU-Catalog of Star Names CSV from exopla.net")  # https://exopla.net/star-names/modern-iau-star-names/
    ap.add_argument("-o", "--output", default="stars.csv")
    args = ap.parse_args()

    hip_map, bayer_map = load_iau(args.iau)

    with open(args.stars, newline='', encoding='utf-8') as fin, \
         open(args.output, 'w', newline='', encoding='utf-8') as fout:
        r = csv.DictReader(fin)
        fieldnames = r.fieldnames
        w = csv.DictWriter(fout, fieldnames=fieldnames)
        w.writeheader()
        for row in r:
            try:
                hip_val = int(row["HIP"])
            except ValueError:
                hip_val = None
            name = ""
            if hip_val and hip_val in hip_map:
                name = hip_map[hip_val]
            elif "BayerASCII" in row and row["BayerASCII"] in bayer_map:
                name = bayer_map[row["BayerASCII"]]
            if name:
                row["Name"] = name
            w.writerow(row)

    print(f"Done. Output written to {args.output}")

if __name__ == "__main__":
    main()

#!/usr/bin/env python3
# hip_to_stars_csv.py
import argparse, csv
from pathlib import Path

def fnum(s): s=s.strip(); return None if not s or s=='?' else float(s)

def parse(line:str):
    if not line or line[0]!="H": return None
    hip   = line[2:14].strip()            # HIP (cols 3-14, 1-based)
    vmag  = fnum(line[41:46])             # Hp or V depending on flag; acceptable近似
    radeg = fnum(line[51:63])
    decdeg= fnum(line[64:76])
    bmv   = fnum(line[245:251])
    if radeg is None or decdeg is None: return None
    return dict(HIP=int(hip), RA=radeg/15.0, Dec=decdeg, Vmag=vmag, BmV=bmv)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("input", default="hip_main.dat", nargs="?", help="Path to hip_main.dat (I/239)")   # hip_main.dat
    ap.add_argument("-o","--output",default="stars_hip.csv")
    args=ap.parse_args()
    with open(args.input,"r",encoding="ascii",errors="ignore") as fin, \
         open(args.output,"w",newline="",encoding="utf-8") as fout:
        w=csv.writer(fout)
        w.writerow(["HIP","HR","Name","RAh","Dec","Vmag","B-V"])
        n=0
        for L in fin:
            rec=parse(L)
            if not rec: continue
            w.writerow([rec["HIP"], "", "", f'{rec["RA"]:.6f}', f'{rec["Dec"]:.6f}',
                        f'{rec["Vmag"]:.2f}' if rec["Vmag"] is not None else "",
                        f'{rec["BmV"]:.2f}' if rec["BmV"] is not None else ""])
            n+=1
    print("wrote", n, "rows to", args.output)

if __name__=="__main__": main()

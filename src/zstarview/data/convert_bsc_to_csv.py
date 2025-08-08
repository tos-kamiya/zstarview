# convert_bsc_to_csv.py
# https://cdsarc.cds.unistra.fr/viz-bin/cat/V/50


def hms_to_decimal(h, m, s):
    return float(h) + float(m) / 60 + float(s) / 3600


def dms_to_decimal(sign, d, m, s):
    sign_mult = -1 if sign == "-" else 1
    return sign_mult * (float(d) + float(m) / 60 + float(s) / 3600)


input_file = "catalog"
output_file = "stars.csv"

with open(input_file, "r") as infile, open(output_file, "w") as outfile:
    outfile.write("HR,Name,RA,Dec,Vmag,B-V\n")
    for line in infile:
        try:
            hr = int(line[0:4])

            name = line[4:14].strip()
            assert "," not in name

            rah = int(line[75:77].strip())
            ram = int(line[77:79].strip())
            ras = float(line[79:83].strip())

            des = line[83].strip()
            ded = int(line[84:86].strip())
            dem = int(line[86:88].strip())
            desec = int(line[88:90].strip())

            vmag = float(line[102:107].strip())
            bv = float(line[109:114].strip())

            ra_deg = hms_to_decimal(rah, ram, ras)
            dec_deg = dms_to_decimal(des, ded, dem, desec)

            outfile.write(f"{hr},{name},{ra_deg:.6f},{dec_deg:.6f},{vmag:.2f},{bv:.2f}\n")
        except Exception:
            continue  # skip malformed lines

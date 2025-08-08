from astroquery.vizier import Vizier
Vizier.ROW_LIMIT = -1  # no limit
Vizier.columns = ["TYC1", "TYC2", "TYC3", "RAmdeg", "DEmdeg", "VTmag", "BTmag", "HIP"]

catalog_list = Vizier.get_catalogs("I/259/tyc2")
tycho2 = catalog_list[0]
tycho2.write("tycho2.csv", format="csv", overwrite=True)
print(f"Saved CSV with {len(tycho2)} rows: tycho2.csv")

# ruční úprava nových metadat do správného AmtDB input formátu (df upravene_meta)

# new_samples - vzorky v původním formátu od Doda
# upravene_meta - vzorky ve formátu, který používá Jirka na import do databáze

import pandas as pd

upravene_meta = pd.DataFrame(columns=puvodni_sloupce)
upravene_meta["id"] = new_samples["identifier"]
upravene_meta["id_alt"] = new_samples["alternative_identifiers"]
upravene_meta["country"] = new_samples["country"]
upravene_meta["continent"] = new_samples["continent"]
upravene_meta["geo_group"] = new_samples["region"]
upravene_meta["culture"] = new_samples["culture"]
upravene_meta["epoch"] = new_samples["epoch"]
upravene_meta["group"] = new_samples["group"]
upravene_meta["comment"] = new_samples["comment"]
upravene_meta["latitude"] = new_samples["latitude"]
upravene_meta["longitude"] = new_samples["longitude"]
upravene_meta["sex"] = new_samples["sex"]
upravene_meta["site"] = new_samples["site"]
upravene_meta["site_detail"] = new_samples["site_detail"]
upravene_meta["mt_hg"] = new_samples["mt_hg"]
upravene_meta["ychr_hg"] = new_samples["ychr_hg"]
upravene_meta["ychr_snps"] = new_samples["ychr_snps"]
upravene_meta["year_from"] = new_samples["year_from"]
upravene_meta["year_to"] = new_samples["year_to"]
upravene_meta["date_detail"] = new_samples["date_detail"]
upravene_meta["bp"] = new_samples["bp"]
upravene_meta["c14_lab_code"] = new_samples["c14_lab_code"]
upravene_meta["reference_name"] = new_samples["reference_name"]
upravene_meta["reference_link"] = new_samples["reference_link"]
upravene_meta["data_link"] = new_samples["data_link"]
upravene_meta["c14_sample_tag"] = new_samples["c14_sample_tag"]
upravene_meta["c14_layer_tag"] = new_samples["c14_layer_tag"]
upravene_meta["full_mt_tag"] = 1
upravene_meta["trusted_tag"] = 1
upravene_meta["avg_coverage"] = new_samples["avg_coverage"]
upravene_meta["amtdb_version"] = "v1.009"


upravene_meta
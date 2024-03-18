# (c) Edvard Ehler, PhD, 2020 and the https://amtdb.org/ team


import pandas as pd
import argparse


def check_hsd_line_details(line, patho_db):
    """
    input: line from HSD (haplogrep) file in a following format: SampleID   Range   Haplogroup  Polymorphisms
    output: polymorphisms with mito-pathological annotation details in the for of (sampleID,"allele1,disease1,status1;allele2,disease2,status2;allele3,disease3,status3;...")
    """
    
    patho_polys = []  # polymorphisms that have pathological annotation
    
    sample_id, sample_range, sample_hg, *polymorphisms = line.strip().split()
    
    
    for poly in polymorphisms:
        if not poly.endswith("N"):
            if poly in patho_db.index:
                #print(sample_id, poly)
                found = mitopatho_db.loc[poly, ["Disease", "Status"]].tolist()
                result = [poly, ] + found

                patho_polys.append(",".join(result))
            
            
    return sample_id, ";".join(patho_polys)


print("""
###########################################################################
# MitoPathoPy - Tool for annotating pathological mutations in human mtDNA #
###########################################################################
Author: Edvard Ehler, PhD (edvard.ehler@img.cas.cz)
Year: 2020
Source: https://github.com/EdaEhler/MitoPathoPy

This tool will search for the (potential) pathological \
mutations in you mtDNA samples. To get the right input format, please, use the Haplogrep 2 \
(https://haplogrep.i-med.ac.at/app/index.html) tool to turn your sequences (fastas) into HSD \
format file. Haplogrep can be also downloaded and run localy \
(https://github.com/seppinho/haplogrep-cmd). Pathological mtDNA mutations are taken from \
data published at https://www.mitomap.org/MITOMAP.

We invite you to check our ancient mtDNA database at: https://amtdb.org/

If you use this tool in you research, please, consider citing our article:
Ehler E, Novotný J, Juras A, Chyleński M, Moravčík O, Pačes J. AmtDB: a database of \
ancient human mitochondrial genomes. Nucleic acids research. 2019 Jan 8;47(D1):D29-32.

###########################################################################
    """)



# initiate parser, more arg groups created
parser = argparse.ArgumentParser(description="Just input and output file, nothing else is needed.")
                                 
                                 #formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument("-i", "--input", help="input file, HSD format, tab-delimited table", required=True)
parser.add_argument("-o", "--output", help="output file, tab-delimited table, one line per sample, \
    following format: sample_name<TAB>allele1,disease1,status1;allele2,disease2,status2;allele3,disease3,status3;...", required=True)

# parsing of the sys.argv
args = parser.parse_args()


mitopatho_db = pd.read_pickle("mitopatho/mitopatho_db_v1.pickle")  # database of mtDNA mutations
test_file = args.input  # input file
mt_patho_dict = {}  # where the info about found patho mutations will be stored


with open(test_file, "r") as hsd:
    for line in hsd:
        if not line.startswith("SampleID\t"): # first HSD line with header
            if line.startswith("#"):  # lines with comments
                pass
            
            else:
                sample, patho_list = check_hsd_line_details(line, patho_db=mitopatho_db)
                mt_patho_dict[sample] = patho_list


with open(args.output, mode="w", encoding="utf-8") as out_file: # create output
    for k, v in mt_patho_dict.items():
        out_file.write(k + "\t" + v + "\n")
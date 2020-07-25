import argparse
import pandas as pd


def get_containing_file(a):
    return "{}.xlsx".format(a[:4])


def print_allele_freqs(mhc):

    summary = pd.read_excel(get_containing_file(mhc), sheet_name='Summary')
    row = summary.loc[summary['Allele'] == mhc]
    black = float(row['Black'].values[0]) * 100.00
    asian = float(row['Asian'].values[0]) * 100.00
    hispanic = float(row['Hispanic'].values[0]) * 100.00
    caucasian = float(row['Caucasian'].values[0]) * 100.00

    print("{}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}".format(mhc, black, asian, hispanic, caucasian))


parser = argparse.ArgumentParser()
parser.add_argument("alleles", help="Comma-separated sequence of MHC II allele, formatted as e.g. DRB1*04:02,DQA1*03:04")
args = parser.parse_args()

alleles = args.alleles.split(",")

print("MHC II Allele\tBlack\tAsian\tHispanic\tCaucasian")
for allele in alleles:
    print_allele_freqs(allele)

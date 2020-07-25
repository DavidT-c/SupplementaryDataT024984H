import pandas as pd

netmhcii = "/home/t024984h/netMHCIIpan-4.0/netMHCIIpan"
families = ["DPA1", "DPB1", "DQA1", "DQB1", "DRB1", "DRB3", "DRB4", "DRB5", ]
ethnicities = ["Black", "Asian", "Hispanic", "Caucasian"]
threshold = 1
factor = 2.5

def print_formatted(family, alleles):
    print("Identified {} {} alleles with a minimum frequency of {}% in at least one ethnicity".format(len(alleles), family, threshold))
    sample_size = 20
    start = 0
    end = sample_size
    while start < len(alleles):
        subset = [a.replace("*", "_").replace(":","") for a in alleles[start:end]]
        print(",".join(subset))
        start += sample_size
        end += sample_size


def as_allele(a, family):
    if family.startswith("DR"):
        a = a.replace("*", "_").replace(":", "")
    elif family.startswith("DPA"):
        a = a.replace("*", "").replace(":", "")
        a = "HLA-{}-DPB10101".format(a)
    elif family.startswith("DPB"):
        a = a.replace("*", "").replace(":", "")
        a = "HLA-DPA10103-{}".format(a)
    elif family.startswith("DQA"):
            a = a.replace("*", "").replace(":", "")
            a = "HLA-{}-DQB10201".format(a)
    elif family.startswith("DQB"):
            a = a.replace("*", "").replace(":", "")
            a = "HLA-DQA10101-{}".format(a)
    return a


def print_netmhcii_command_line(family, alleles, fasta_file):
    sample_size = 20
    start = 0
    end = sample_size
    count = 1
    while start < len(alleles):
        subset = ",".join([as_allele(a, family) for a in alleles[start:end]])
        print("{0} -f SARS_CoV_2_{1}_protein.fasta -a {2} -s -BA -xlsfile {4}_{3}.xls > {4}_{3}_{1}.out".format(netmhcii, fasta_file, subset, str(count), family))
        start += sample_size
        end += sample_size
        count += 1

def get_ratios(a, vals):
    return [a/v if v > 0.0 else 100.0 for v in vals]

def is_large_difference(vals, gt_all=True):
    """
    Looks for a large ratio between values.
    If gt_all is True, then one value has to be greater than all the others by thespeiifed amount,
    otherwise it has to be that much greater than just one.
    """
    if all(x < 0.05 for x in vals):
        return False
    ratios = []
    for i in range(0, len(vals)):
        ratios = get_ratios(vals[i], vals[:i] + vals[i + 1:])
        if gt_all:
            if all(r > factor for r in ratios):
                return True
        else:
            if any(r > factor for r in ratios):
                return True
    return False
        #ratios += get_ratios(vals[i], vals[:i] + vals[i+1:])
    #return any(r > factor for r in ratios)


def find_variable_prevalence(data, alleles, gt_all=True):
    """
    Determines whether the AF for any ethnicity is greater than a defined multiple of any other
    :param data:
    :param alleles:
    :param difference:
    :return:
    """
    results = []
    for a in alleles:
        row = data.loc[data['Allele'] == a]
        black = row['Black'].values[0]
        asian = row['Asian'].values[0]
        hispanic = row['Hispanic'].values[0]
        caucasian = row['Caucasian'].values[0]
        if alleles[0].startswith("DP"):
            vals = [black, asian, caucasian]
        else:
            vals = [black, asian, hispanic, caucasian]
        if is_large_difference(vals, gt_all):
            results.append([a, black*100, asian*100, hispanic*100, caucasian*100])
    return results


def get_all_alleles(ethnicity):
    alleles = []
    frequencies = []
    for family in families:
        data = pd.read_excel('{}.xlsx'.format(family), sheet_name=ethnicity)
        alleles += data['Allele'].tolist()
        frequencies += data['Overall Allele Frequency'].tolist()
    af = list(zip(alleles, frequencies))
    return sorted(af, key=lambda x: x[1])[::-1][:1]


significant = {}
one_gt_all = []
one_gt_at_least_one = []
for family in families:
    summary = pd.read_excel('{}.xlsx'.format(family), sheet_name='Summary')
    alleles = summary['Significant'].tolist()
    found = [a for a in alleles if not pd.isna(a)]
    significant[family] = found
    one_gt_all += find_variable_prevalence(summary, found, True)
    one_gt_at_least_one += find_variable_prevalence(summary, found, False)

alleles_by_ethnicity = {}
for ethnicity in ethnicities:
    alleles_by_ethnicity[ethnicity] = get_all_alleles(ethnicity)

for c in significant:
    print_formatted(c, significant[c])

for protein in ["spike", "envelope", "membrane", "nucleocapsid"]:
    print("NetMHCIIpan command line sequences for {} :".format(protein))
    for c in significant:
        print_netmhcii_command_line(c, significant[c], protein)

for f in significant:
    print("{}\t{}".format(f, len(significant[f])))

print("\nAlleles where AF for one ethnicity is at least {} times greater than all others:".format(factor))
print("{}\t{}\t{}\t{}\t{}".format("Allele", "Black", "Asian", "Hispanic", "Caucasian"))
for v in one_gt_all:
    print("\t".join([str(l) for l in v]))
print(len(one_gt_all))

print("\nAlleles where AF for one ethnicity is at least {} times greater than at least one other:".format(factor))
print("{}\t{}\t{}\t{}\t{}".format("Allele", "Black", "Asian", "Hispanic", "Caucasian"))
for v in one_gt_at_least_one:
    print("\t".join([str(l) for l in v]))
print(len(one_gt_at_least_one))

print("\n Most common MHC II allele for each ethnicity")
for ethnicity in alleles_by_ethnicity:
    for a in alleles_by_ethnicity[ethnicity]:
        print("{}\t{}\t{:.2f}".format(ethnicity, a[0], a[1]*100))
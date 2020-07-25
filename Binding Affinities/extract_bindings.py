import argparse
import glob

protein_length = {"spike": 1273, "envelope": 75, "membrane": 222, "nucleocapsid": 419}
affinity_data = {}


def reform_allele(mhc):
    """
    The NetMHCIIpan-4.0 tool accepts alleles in the form <family><number>_<major><minor>
    e.g. DRB1_1501 but for the report alleles are formatted as DRB1*15:01.
    This method applies the report format
    :param mhc: MHC II allele in the form accepted by NetMHCIIpan
    :return: Reformatted allele name
    """
    if "_" in mhc:
        mhc = mhc.replace("_", "*")
        return "{}:{}".format(mhc[:7], mhc[7:])
    return "{}*{}:{}".format(mhc[:4], mhc[4:6], mhc[6:])


def get_strongest_allele_epitope(protein):

    highest_el_score = 0.0
    strongest_mhc = ""
    strongest_data = {}
    for mhc in affinity_data[protein]:
        # Get the most strongly binding epitope for each allele
        strongest_for_mhc = affinity_data[protein][mhc]["predictions"][0]
        if strongest_for_mhc["el_score"] > highest_el_score:
            strongest_mhc = mhc
            strongest_data = strongest_for_mhc
            highest_el_score = strongest_data["el_score"]
    outstr = "{}\t{}\t{}\t{}\t{}\t{}".format(
        protein,
        reform_allele(strongest_mhc),
        strongest_data["peptide"],
        strongest_data["core"],
        strongest_data["position"],
        highest_el_score,
    )
    return outstr


def get_weakest_allele_epitope(protein):
    lowest_el_score = 100.0
    weakest_mhc = ""
    weakest_data = {}
    for mhc in affinity_data[protein]:
        # Get the most weakly binding epitope for each allele
        weakest_for_mhc = affinity_data[protein][mhc]["predictions"][0]
        if weakest_for_mhc["el_score"] < lowest_el_score:
            weakest_mhc = mhc
            weakest_data = weakest_for_mhc
            lowest_el_score = weakest_data["el_score"]
    outstr = "{}\t{}\t{}\t{}\t{}\t{}".format(
        protein,
        reform_allele(weakest_mhc),
        weakest_data["peptide"],
        weakest_data["core"],
        weakest_data["position"],
        lowest_el_score,
    )
    return outstr


def get_affinity_data(data):
    table = {}
    mhc = "None"
    family = ""
    for line in data:
        line = line.strip()
        # Look for 'Family:' and extract the group name (needed to work out which DQ or DP allele we are looking at)
        # Look for '# Allele' and extract the allele name
        if line.startswith("Family:"):
            family = line.split()[1]

        elif line.startswith("# Allele:"):
            allele = line.split()[2]

            if allele.startswith(family):
                mhc = allele
            else:
                components = allele.split("-")
                for c in components:
                    if c.startswith(family):
                        mhc = c
            table[mhc] = {}
            table[mhc]["predictions"] = []

        # Look for a line with "binders" in it, and sxtract the strong & weak values
        elif "binders" in line:
            fields = line.strip().split()
            table[mhc]["strong_binders"] = int(fields[4])
            table[mhc]["weak_binders"] = int(fields[9])
        elif mhc in line:
            prediction = {}
            fields = line.split()
            prediction["peptide"] = fields[2]
            prediction["core"] = fields[4]
            prediction["el_score"] = float(fields[7])
            prediction["position"] = int(fields[0])
            table[mhc]["predictions"].append(prediction)

    return table


def order_alleles_by_strong_binders(protein, best):
    strong_binders = []
    for mhc in affinity_data[protein]:
        strong_binders.append([mhc, affinity_data[protein][mhc]["strong_binders"]])

    sorted_by_strongest = sorted(strong_binders, key=lambda x: x[1])[::-1][:best]
    first = True
    outstr = ""
    for m in sorted_by_strongest:
        if first:
            pr = "{} ({})".format(protein, protein_length[protein])
            first = False
        else:
            pr = ""
        outstr += "{}\t{}\t{}\t{:.2f}\n".format(
            pr, reform_allele(m[0]), m[1], (m[1] / protein_length[protein]) * 100
        )
    return outstr


def order_alleles_by_strong_and_weak_binders(protein, best):
    all_binders = []
    for mhc in affinity_data[protein]:
        strong = affinity_data[protein][mhc]["strong_binders"]
        weak = affinity_data[protein][mhc]["weak_binders"]
        both = strong + weak
        all_binders.append([mhc, both, strong, weak])

    sorted_by_total = sorted(all_binders, key=lambda x: x[1])[::-1][:best]
    first = True
    outstr = ""
    for m in sorted_by_total:
        if first:
            pr = "{} ({})".format(protein, protein_length[protein])
            first = False
        else:
            pr = ""
        outstr += "{}\t{}\t{} ({} strong; {} weak)\t{:.2f}\n".format(
            pr,
            reform_allele(m[0]),
            m[1],
            m[2],
            m[3],
            (m[1] / protein_length[protein]) * 100,
        )
    return outstr


def parse_result_files(protein):
    """
    Reads in the tabulated affinity prediction data output by NetMHCIIpan-4.0
    :param protein: Set of predictions to load
    :return: Table of binding affinities plus the single strongest-binding eipitope from the protein
    """
    files = glob.glob("*{}*.out".format(protein))
    lines = []
    for f in files:
        family = f.split("_")[0]
        lines.append("Family: {}".format(family))
        with open(f, "r") as infile:
            lines += infile.readlines()

    data = get_affinity_data(lines)

    return data


parser = argparse.ArgumentParser()

protein_options = ["spike", "envelope", "membrane", "nucleocapsid"]
parser.add_argument(
    "-p",
    "--protein",
    required=True,
    choices=protein_options + ["all"],
    help="Protein to check",
)
parser.add_argument(
    "-t",
    "--top",
    action="store_true",
    default=False,
    help="Output the MHC II allele with the strongest binding affinity for a single epitope",
)
parser.add_argument(
    "-s",
    "--strongly_binding",
    action="store_true",
    default=False,
    help="Lists the MHC II alleles calculated as having the highest number of strongly-binding epitopes",
)
parser.add_argument(
    "-b",
    "--both_weak_and_strong",
    action="store_true",
    default=False,
    help="Lists the MHC II alleles calculated as having the highest number of strongly- and weakly-binding epitopes",
)
parser.add_argument(
    "-w",
    "--weakest",
    action="store_true",
    default=False,
    help="Lists the alleles weakest predicted binding affinity",
)

args = parser.parse_args()

if args.protein == "all":
    proteins = protein_options
else:
    proteins = [args.protein]

for p in proteins:
    affinity_data[p] = parse_result_files(p)

if args.top:
    print("Protein\tMHC II Allele\tPeptide\tCore\tPosition\tEL Score")
    for p in proteins:
        print(get_strongest_allele_epitope(p))

if args.weakest:
    print("Protein\tMHC II Allele\tPeptide\tCore\tPosition\tEL Score")
    for p in proteins:
        print(get_weakest_allele_epitope(p))


if args.strongly_binding:
    print("\t\tStrong Binders")
    print("Protein (Length)\tMHC II Allele\tTotal\tper 100 Amino Acids")
    for p in proteins:
        print(order_alleles_by_strong_binders(p, 5), end="")

if args.both_weak_and_strong:
    print("\t\tTotal (Strong + Weak) Binders")
    print("Protein (Length)\tMHC II Allele\tTotal\tper 100 Amino Acids")
    for p in proteins:
        print(order_alleles_by_strong_and_weak_binders(p, 5), end="")

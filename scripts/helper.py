#!/usr/bin/python
import json
import re
import sys

from dendropy import DnaCharacterMatrix


def create_sub_files(
    alignment_file,
    dates_file,
    subtree_file,
    subtree_dates_file,
    subfasta_file,
    new_dates_file,
):
    dates_dic = read_dates(dates_file)

    # clean up comments and add dates to end of taxon names
    with open(subtree_file, "r") as fp:
        content = fp.read().replace("None", "")
        content = re.sub("NODE_\d+", "", content)
        for taxon, date in dates_dic.items():
            content = content.replace(taxon, taxon + "_" + date)

    with open(subtree_dates_file, "w") as fp:
        fp.write(content)

    # add dates to end of sequence names
    sub_aln_dic = {}
    dna = DnaCharacterMatrix.get(path=alignment_file, schema="fasta")

    for taxon, date in dates_dic.items():
        t = dna.taxon_namespace.get_taxon(label=taxon)
        new_taxon_name = taxon + "_" + date
        sub_aln_dic[new_taxon_name] = str(dna[t])
    sub_dna = DnaCharacterMatrix.from_dict(sub_aln_dic)
    sub_dna.write(path=subfasta_file, schema="fasta")

    with open(new_dates_file, "w") as fp:
        fp.write(str(len(dates_dic)))
        for taxon, date in dates_dic.items():
            fp.write("\n" + taxon + "_" + date + "\t" + date)


def prepare_physher(
    subfasta_file, subtree_file, dates_file, json_template_file, json_file, iterations
):
    dates = read_dates(dates_file)
    with open(json_template_file, "r") as fp:
        content = fp.read()
    dates_json = ",\n".join(
        ['"{}_{}":{}'.format(taxon, date, date) for taxon, date in dates.items()]
    )
    dates_json = dates_json.rstrip(",\n")
    content = (
        content.replace("FILE_TEMPLATE", subfasta_file)
        .replace("DATES_TEMPLATE", dates_json)
        .replace("TREE_TEMPLATE", subtree_file)
        .replace("DIM_TEMPLATE", str(len(dates) - 1))
        .replace("TEMPLATE_ITER", iterations)
    )
    with open(json_file, "w") as fp:
        fp.write(content)


def prepare_phylotorch(
    subfasta_file,
    subtree_file,
    dates_file,
    json_template_file,
    json_file,
    iterations,
    bito,
):
    dates = read_dates(dates_file)
    taxa = []
    datess = list(map(float, dates.values()))
    root_shift = max(datess) - min(datess)
    for taxon, date in dates.items():
        taxa.append(
            {
                "id": "{}_{}".format(taxon, date),
                "type": "Taxon",
                "attributes": {"date": float(date)},
            }
        )

    with open(json_template_file, "r") as fp:
        content = fp.read()

    content = (
        content.replace("TAXA_TEMPLATE", json.dumps(taxa))
        .replace("ITERATION_TEMPLATE", iterations)
        .replace("ROOT_SHIFT_TEMPLATE", str(root_shift))
        .replace("DIM_TEMPLATE", str(len(datess) - 2))
    )
    if bito.lower() == "true":
        content = (
            content.replace("SEQUENCES_TEMPLATE", '"' + subfasta_file + '"')
            .replace("TREE_TEMPLATE", '"' + subtree_file + '"')
            .replace('"newick"', '"file"')
            .replace('"sequences"', '"file"')
        )
    else:
        with open(subtree_file, "r") as fp:
            newick = fp.read().strip()
        alignment = DnaCharacterMatrix.get(path=subfasta_file, schema="fasta")
        sequences = []
        for name in alignment:
            sequences.append(
                {"taxon": str(name).strip("'"), "sequence": str(alignment[name])}
            )

        content = content.replace("SEQUENCES_TEMPLATE", json.dumps(sequences)).replace(
            "TREE_TEMPLATE", '"' + newick + '"'
        )

    with open(json_file, "w") as fp:
        fp.write(content)


def read_dates(dates_file):
    dates_dic = {}
    with open(dates_file, "r") as fp:
        for line in fp:
            l = line.rstrip().split("\t")
            if len(l) == 2:
                dates_dic[l[0]] = l[1]
    return dates_dic


def convert_lsd_nexus_to_newick(nexus_file, newick_file):
    with open(nexus_file, "r") as fp:
        for line in fp:
            if line.startswith("tree 1 = "):
                newick = line.replace("tree 1 = ", "").strip()
    newick = re.sub('\[&date="\d+\.?\d*"]', "", newick)
    with open(newick_file, "w") as fp:
        fp.write(newick)


if sys.argv[1] == "0":
    create_sub_files(*sys.argv[2:])
elif sys.argv[1] == "1":
    prepare_physher(*sys.argv[2:])
elif sys.argv[1] == "2":
    convert_lsd_nexus_to_newick(*sys.argv[2:])
elif sys.argv[1] == "3":
    prepare_phylotorch(*sys.argv[2:])

#!/usr/bin/env python

import re

import click
from dendropy import DnaCharacterMatrix


@click.command(help="Create sub files (alignment, tree, dates)")
@click.option("--input", type=click.UNPROCESSED, required=True, help="alignment file")
@click.option("--tree", type=click.UNPROCESSED, required=True, help="tree file")
@click.option("--dates", type=click.UNPROCESSED, required=True, help="date file")
@click.option(
    "--out_fasta", type=click.UNPROCESSED, required=True, help="alignment subset"
)
@click.option("--out_tree", type=click.UNPROCESSED, required=True, help="tree subset")
@click.option("--out_dates", type=click.UNPROCESSED, required=True, help="date subset")
def subfiles(
    input,
    tree,
    dates,
    out_fasta,
    out_tree,
    out_dates,
):
    dates_dic = read_dates(dates)

    # clean up comments and add dates to end of taxon names
    with open(tree, "r") as fp:
        content = fp.read().replace("None", "")
        content = re.sub("NODE_\d+", "", content)
        for taxon, date in dates_dic.items():
            content = content.replace(taxon, taxon + "_" + date)

    with open(out_tree, "w") as fp:
        fp.write(content)

    # add dates to end of sequence names
    sub_aln_dic = {}
    dna = DnaCharacterMatrix.get(path=input, schema="fasta")

    for taxon, date in dates_dic.items():
        t = dna.taxon_namespace.get_taxon(label=taxon)
        new_taxon_name = taxon + "_" + date
        sub_aln_dic[new_taxon_name] = str(dna[t])
    sub_dna = DnaCharacterMatrix.from_dict(sub_aln_dic)
    sub_dna.write(path=out_fasta, schema="fasta")

    with open(out_dates, "w") as fp:
        fp.write(str(len(dates_dic)))
        for taxon, date in dates_dic.items():
            fp.write("\n" + taxon + "_" + date + "\t" + date)


@click.command(help="Generate JSON file from template for physher")
@click.option("--input", type=click.UNPROCESSED, required=True, help="alignment file")
@click.option("--tree", type=click.UNPROCESSED, required=True, help="tree file")
@click.option("--dates", type=click.UNPROCESSED, required=True, help="date file")
@click.option("--template", type=click.UNPROCESSED, required=True, help="template file")
@click.option("--output", type=click.UNPROCESSED, required=True, help="JSON file")
@click.option("--iterations", type=int, required=True, help="number of iterations")
@click.option("--rate", type=float, required=True, help="substitution rate")
@click.option("--lr", type=float, required=True, help="learning rate")
@click.option("--tol", type=float, required=True, help="tolerance")
def physher(input, tree, dates, template, output, iterations, rate, lr, tol):
    dates = read_dates(dates)
    with open(template, "r") as fp:
        content = fp.read()
    dates_json = ",\n".join(
        ['"{}":{}'.format(taxon, date) for taxon, date in dates.items()]
    )
    dates_json = dates_json.rstrip(",\n")
    content = (
        content.replace("FILE_TEMPLATE", input)
        .replace("DATES_TEMPLATE", dates_json)
        .replace("TREE_TEMPLATE", tree)
        .replace("DIM_TEMPLATE", str(len(dates) - 1))
        .replace("TEMPLATE_ITER", str(iterations))
        .replace("RATE_TEMPLATE", str(rate))
        .replace("LR_TEMPLATE", str(lr))
        .replace("TOL_TEMPLATE", str(tol))
    )
    with open(output, "w") as fp:
        fp.write(content)


def read_dates(dates_file):
    dates_dic = {}
    with open(dates_file, "r") as fp:
        for line in fp:
            l = line.rstrip().split("\t")
            if len(l) == 2:
                dates_dic[l[0]] = l[1]
    return dates_dic


@click.command(help="Convert nexus file created by lsd to newick format")
@click.option("--input", type=click.UNPROCESSED, required=True, help="nexus file")
@click.option("--output", type=click.UNPROCESSED, required=True, help="newick file")
def convert(input, output):
    with open(input, "r") as fp:
        for line in fp:
            if line.startswith("tree 1 = "):
                newick = line.replace("tree 1 = ", "").strip()
    newick = re.sub(r'\[&date="\d+\.?\d*"]', "", newick)
    with open(output, "w") as fp:
        fp.write(newick)


@click.group(help="CLI tool to convert and generate files")
def cli():
    pass


cli.add_command(physher)
cli.add_command(subfiles)
cli.add_command(convert)

if __name__ == "__main__":
    cli()

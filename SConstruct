import os
import re
from os.path import join
from subprocess import check_call
from timeit import default_timer as timer
import json

from dendropy import DnaCharacterMatrix

config_file = "scons.conf"
opts = Variables(config_file)
# From https://github.com/GalSim-developers/GalSim/blob/releases/2.2/SConstruct
opts.Add(
    PathVariable(
        "DYLD_LIBRARY_PATH",
        "Set the DYLD_LIBRARY_PATH inside of SCons.  "
        + "Particularly useful on El Capitan (and later), since Apple strips out "
        + "DYLD_LIBRARY_PATH from the environment that SCons sees, so if you need it, "
        + "this option enables SCons to set it back in for you by doing "
        + "`scons DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH`.",
        "",
        PathVariable.PathAccept,
    )
)

env = Environment(ENV=os.environ)
opts.Update(env)
if "DYLD_LIBRARY_PATH" in env:
    os.environ["DYLD_LIBRARY_PATH"] = env["DYLD_LIBRARY_PATH"]

dataset_size = [20, 50, 100, 200 , 500, 750, 1000, 1250, 1500, 2000]
iterations = ARGUMENTS.get("replicates", "5000")

treetime_base = "treetime_validation"
flu_H3N2_dataset = join(treetime_base, "flu_H3N2", "subtree_samples", "dataset")
dates_dir = join(flu_H3N2_dataset, "LSD_out")
subtrees_dir = join(flu_H3N2_dataset, "subtrees")

alignment_file = join(treetime_base, "resources", "flu_H3N2", "H3N2_HA_2011_2013.fasta")
tree_file = join(treetime_base, "resources", "flu_H3N2", "H3N2_HA_2011_2013.nwk")

subdata_dir = join("flu_H3N2", "results")
json_file = join("flu_H3N2", "H3N2_HA_2011_2013.json")
phylotorch_json_file = join("flu_H3N2", "phylotorch.json")
phylotorch_libsbn_json_file = join("flu_H3N2", "phylotorch-libsbn.json")

dna = DnaCharacterMatrix.get(path=alignment_file, schema="fasta")

if not os.path.lexists(subdata_dir):
    os.makedirs(subdata_dir, exist_ok=True)


def prepare_physher(target, source, env):
    subfasta_file, subtree_file, dates_file = map(str, source)
    output_tree = subfasta_file.replace(".fasta", ".physher.tree")
    dates = read_dates(dates_file)
    with open(json_file, "r") as fp:
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
        .replace("OUTPUT_TEMPLATE", output_tree)
        .replace("TEMPLATE_ITER", iterations)
    )
    with open(str(target[0]), "w") as fp:
        fp.write(content)


def prepare_phylotorch(target, source, env):
    subfasta_file, subtree_file, dates_file, json_template = map(str, source)
    dates = read_dates(dates_file)
    taxa = []
    datess = list(map(float, dates.values()))
    root_shift = max(datess) - min(datess)
    for taxon, date in dates.items():
        taxa.append(
            {
                "id": "{}_{}".format(taxon, date),
                "type": "phylotorch.evolution.taxa.Taxon",
                "attributes": {"date": float(date)},
            }
        )

    with open(json_template, "r") as fp:
        content = fp.read()

    if "libsbn" in env:
        content = (
            content.replace("SEQUENCES_TEMPLATE", subfasta_file)
            .replace("TAXA_TEMPLATE", json.dumps(taxa))
            .replace("TREE_TEMPLATE", subtree_file)
            .replace("ITERATION_TEMPLATE", iterations)
            .replace("ROOT_SHIFT_TEMPLATE", str(root_shift))
            .replace("DIM_TEMPLATE", str(len(datess) - 2))
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

        content = (
            content.replace("SEQUENCES_TEMPLATE", json.dumps(sequences))
            .replace("TAXA_TEMPLATE", json.dumps(taxa))
            .replace("TREE_TEMPLATE", '"' + newick + '"')
            .replace("ITERATION_TEMPLATE", iterations)
            .replace("ROOT_SHIFT_TEMPLATE", str(root_shift))
            .replace("DIM_TEMPLATE", str(len(datess) - 2))
        )

    with open(str(target[0]), "w") as fp:
        fp.write(content)


def timeit(method):
    def timed(target, source, env):
        start = timer()
        f = open(str(target[0]), "w")
        method(target, source, env, **{"stdout": f})
        end = timer()
        total_time = end - start
        f.write("TIME: {}".format(total_time))
        f.close()

    return timed


@timeit
def run_physher(target, source, env, **kwargs):
    check_call(
        ["/Users/mathieu/Downloads/physher-master/Release/physher", str(source[0])],
        stdout=kwargs.get("stdout"),
    )


@timeit
def run_phylostan(target, source, env, **kwargs):
    seq_file, tree_file, script = map(str, source)
    cmd = [
        "phylostan",
        "run",
        "-i",
        seq_file,
        "-t",
        tree_file,
        "-s",
        script,
        "-o",
        str(target[1]),
        "--iter",
        env["iter"],
        "--eta",
        "0.00000001",
        "--elbo_samples",
        "1",
        "--samples",
        "1",
    ]
    cmd.extend(env["model"].split())
    print(" ".join(cmd))
    check_call(cmd, stdout=kwargs.get("stdout"))


@timeit
def run_phylotorch(target, source, env, **kwargs):
    check_call(["phylotorch", str(source[0])], stdout=kwargs.get("stdout"))


@timeit
def run_lsd(target, source, env, **kwargs):
    cmd = [
        join(treetime_base, "bin", "lsd2_mac"),
        "-i",
        str(source[0]),
        "-d",
        str(source[1]),
        "-o",
        str(source[1]).replace("_dates.txt", ".out"),
        "-s",
        "1701",
        "-c",
    ]
    check_call(cmd, stdout=kwargs.get("stdout"))


def convert_lsd_nexus_to_newick(target, source, env):
    with open(str(source[0]), "r") as fp:
        for line in fp:
            if line.startswith("tree 1 = "):
                newick = line.replace("tree 1 = ", "").strip()
    newick = re.sub('\[&date="\d+\.?\d*"]', "", newick)
    with open(str(target[0]), "w") as fp:
        fp.write(newick)


def read_dates(dates_file):
    dates_dic = {}
    with open(dates_file, "r") as fp:
        for line in fp:
            l = line.rstrip().split("\t")
            if len(l) == 2:
                dates_dic[l[0]] = l[1]
    return dates_dic


def parse_results(target, source, env):
    pattern_time = re.compile(r"TIME: (\d+\.\d+)")
    csvp = open(str(target[0]), "w")
    csvp.write("method,replicate,taxa,time\n")

    for infile in source:
        total_time = -1
        stem, method = str(infile).replace(".txt", "").split(".")
        temp = stem.split("_")
        with open(str(infile), "r") as fp:
            for line in fp:
                if line.startswith("#"):
                    continue
                line = line.rstrip("\n").rstrip("\r")
                mtime = pattern_time.match(line)
                if mtime:
                    total_time = mtime.group(1)
                    break

        csvp.write("{},{},{},{}\n".format(method, temp[-1], temp[-2], total_time))
    csvp.close()


def create_sub_files(target, source, env):
    dates_file, subtree_file = map(str, source)
    subtree_dates_file, subfasta_file, new_dates_file = map(str, target)
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


stan_script_file = join(subdata_dir, "H3N2_HA_2011_2013.stan")
stan_pkl = stan_script_file.replace(".stan", ".pkl")

env.Command(
    [stan_script_file, stan_pkl],
    None,
    "phylostan build -s "
    + join(subdata_dir, "H3N2_HA_2011_2013.stan")
    + " -m JC69 --heterochronous --estimate_rate --clock strict -c constant --compile",
)

outputs = []

for n_taxa in dataset_size:

    stem = join(subdata_dir, "H3N2_HA_2011_2013_{}_0".format(n_taxa))

    dates_file = join(dates_dir, "H3N2_HA_2011_2013_{}_0.lsd_dates.txt".format(n_taxa))
    subtree_file = join(subtrees_dir, "H3N2_HA_2011_2013_{}_0.nwk".format(n_taxa))

    subtree_dates_file = stem + ".nwk"  # branch=subst
    seq_file = stem + ".fasta"
    new_dates_file = stem + ".lsd_dates.txt"

    lsd_stem = stem + ".lsd.out"
    lsd_tree_dated = lsd_stem + ".date.nexus"  # branch=time
    lsd_tree = lsd_stem + ".nexus"  # branch=subst
    lsd_tree_newick = lsd_stem + ".nwk"  # branch=subst

    lsd_tree_dated_newick = lsd_stem + ".date.nwk"  # branch=time

    physher_tree_dated = stem + ".physher.tree"

    subjson_file = stem + ".json"
    phylotorch_subjson_file = stem + "-phylotorch.json"
    phylotorch_libsbn_subjson_file = stem + "-phylotorch-libsbn.json"

    env.Command(
        target=[subtree_dates_file, seq_file, new_dates_file],
        source=[dates_file, subtree_file],
        action=create_sub_files,
    )

    env.Command(
        target=[stem + ".lsd", lsd_tree_dated, lsd_tree, lsd_stem, lsd_tree_newick],
        source=[subtree_dates_file, new_dates_file],
        action=run_lsd,
    )

    env.Command(
        target=lsd_tree_dated_newick,
        source=lsd_tree_dated,
        action=convert_lsd_nexus_to_newick,
    )

    # physher
    env.Command(
        target=subjson_file,
        source=[seq_file, lsd_tree_dated, dates_file],
        action=prepare_physher,
    )

    env.Command(
        target=[stem + ".physher.txt", physher_tree_dated],
        source=subjson_file,
        action=run_physher,
    )

    # phylotorch
    env.Command(
        target=phylotorch_subjson_file,
        source=[seq_file, physher_tree_dated, dates_file, phylotorch_json_file],
        action=prepare_phylotorch,
    )

    env.Command(
        target=stem + ".phylotorch.txt",
        source=phylotorch_subjson_file,
        action=run_phylotorch,
    )

    # phylotorch libsbn
    env.Command(
        target=phylotorch_libsbn_subjson_file,
        source=[seq_file, physher_tree_dated, dates_file, phylotorch_libsbn_json_file],
        action=prepare_phylotorch,
        libsbn=True,
    )

    env.Command(
        target=stem + ".libsbn.txt",
        source=phylotorch_libsbn_subjson_file,
        action=run_phylotorch,
    )

    # phylostan
    env.Command(
        target=[stem + ".phylostan.txt", stem, stem + ".diag", stem + ".trees"],
        source=[seq_file, physher_tree_dated, stan_script_file],
        action=run_phylostan,
        iter=iterations,
        model="-m JC69 --heterochronous --estimate_rate --clock strict -c constant",
    )

    outputs.extend(
        [
            "{}.{}.txt".format(stem, m)
            for m in ("physher", "phylotorch", "libsbn", "phylostan")
        ]
    )

env.Command(join("flu_H3N2", "H3N2_HA_2011_2013.csv"), outputs, parse_results)

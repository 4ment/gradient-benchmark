import os
import re
from os.path import join
from subprocess import check_call
from timeit import default_timer as timer

from dendropy import DnaCharacterMatrix


env = Environment(ENV=os.environ)

N_leaves_array = [20, 50, 100, 200, 500, 750, 1000, 1250, 1500, 2000]
iterations = ARGUMENTS.get('replicates', '5000')

treetime_base = 'treetime_validation'
flu_H3N2_dataset = join(treetime_base, 'flu_H3N2', 'subtree_samples', 'dataset')
dates_dir = join(flu_H3N2_dataset, 'LSD_out')
subtrees_dir = join(flu_H3N2_dataset, 'subtrees')


alignment_file = join(treetime_base, 'resources', 'flu_H3N2', 'H3N2_HA_2011_2013.fasta')
tree_file = join(treetime_base, 'resources', 'flu_H3N2', 'H3N2_HA_2011_2013.nwk')

subdata_dir = join('flu_H3N2', 'results')
json_file = join('flu_H3N2', 'H3N2_HA_2011_2013.json')

dna = DnaCharacterMatrix.get(path=alignment_file, schema='fasta')

if not os.path.lexists(subdata_dir):
    os.makedirs(subdata_dir, exist_ok=True)


def prepare_physher(target, source, env):
    subfasta_file, subtree_file, dates_file = map(str, source)
    output_tree = subfasta_file.replace('.fasta', '.physher.tree')
    dates = read_dates(dates_file)
    with open(json_file, 'r') as fp:
        content = fp.read()
    dates_json = ',\n'.join(['"{}_{}":{}'.format(taxon, date, date) for taxon, date in dates.items()])
    dates_json = dates_json.rstrip(',\n')
    content = content.replace('FILE_TEMPLATE', subfasta_file).replace('DATES_TEMPLATE', dates_json).replace(
        'TREE_TEMPLATE', subtree_file).replace('DIM_TEMPLATE', str(len(dates) - 1)).replace('OUTPUT_TEMPLATE',
                                                                                            output_tree).replace(
        'TEMPLATE_ITER', iterations)
    with open(str(target[0]), 'w') as fp:
        fp.write(content)


def run_physher(target, source, env):
    start = timer()
    f = open(str(target[0]), 'w')
    check_call(['physher', str(source[0])], stdout=f)
    end = timer()
    total_time = end - start
    f.write('TIME: {}'.format(total_time))
    f.close()


def run_phylostan(target, source, env):
    seq_file, tree_file, script = map(str, source)
    start = timer()
    f = open(str(target[0]), 'w')
    check_call(['phylostan', 'run', '-i', seq_file, '-t', tree_file, '-s', script, '-o', seq_file + '.stan', '--iter',
                env['iter'], '--eta', '0.00000001', '-m', 'JC69', '--heterochronous', '--estimate_rate', '--clock',
                'strict', '-c', 'constant', '--elbo_samples', '1', '--samples', '1'], stdout=f)
    end = timer()
    total_time = end - start
    f.write('TIME: {}'.format(total_time))
    f.close()


def run_phylotorch(target, source, env):
    start = timer()
    f = open(str(target[0]), 'w')
    cmd = ['phylotorch', '-i', str(source[0]), '-t', str(source[1]), '--iter', env['iter'], '--eta', '0.00000001', '-m',
           'JC69', '--heterochronous', '-c', 'constant', '--elbo_samples', '1']
    if 'libsbn' in env:
        cmd.append('--libsbn')
    check_call(cmd, stdout=f)
    end = timer()
    total_time = end - start
    f.write('TIME: {}'.format(total_time))
    f.close()


def run_phylojax(target, source, env):
    start = timer()
    f = open(str(target[0]), 'w')
    check_call(['phylojax', '-i', str(source[0]), '-t', str(source[1]), '--iter', env['iter'], '--eta', '0.00000001',
                '--elbo_samples', '1'])
    end = timer()
    total_time = end - start
    f.write('TIME: {}'.format(total_time))
    f.close()


def run_lsd(target, source, env):
    start = timer()
    f = open(str(target[0]), 'w')
    cmd = ['lsd', '-i', str(source[0]), '-d', str(source[1]), '-o',
           str(source[1]).replace('_dates.txt', '.out'), '-s', '1701', '-c']
    check_call(cmd)
    end = timer()
    total_time = end - start
    f.write('TIME: {}'.format(total_time))
    f.close()


def read_dates(dates_file):
    dates_dic = {}
    with open(dates_file, 'r') as fp:
        for line in fp:
            l = line.rstrip().split('\t')
            if len(l) == 2:
                dates_dic[l[0]] = l[1]
    return dates_dic


def parse_results(target, source, env):
    pattern_time = re.compile(r'TIME: (\d+\.\d+)')
    csvp = open(str(target[0]), 'w')
    csvp.write('method,replicate,taxa,time\n')

    for infile in source:
        total_time = -1
        stem, method = str(infile).split('.')
        temp = stem.split('_')
        with open(str(infile), "r") as fp:
            for line in fp:
                if line.startswith("#"):
                    continue
                line = line.rstrip('\n').rstrip('\r')
                mtime = pattern_time.match(line)
                if mtime:
                    total_time = mtime.group(1)
                    break

        csvp.write('{},{},{},{}\n'.format(method, temp[-1], temp[-2], total_time))
    csvp.close()


for n_taxa in N_leaves_array:
    dates_file = join(dates_dir, 'H3N2_HA_2011_2013_{}_0.lsd_dates.txt'.format(n_taxa))
    subtree_file = join(subtrees_dir, 'H3N2_HA_2011_2013_{}_0.nwk'.format(n_taxa))

    dates_dic = read_dates(dates_file)

    # clean up comments and add dates to end of taxon names
    with open(subtree_file, 'r') as fp:
        content = fp.read().replace('None', '')
        content = re.sub('NODE_\d+', '', content)
        for taxon, date in dates_dic.items():
            content = content.replace(taxon, taxon + '_' + date)

    subtree_dates_file = join(subdata_dir, 'H3N2_HA_2011_2013_{}_0.nwk'.format(n_taxa))
    subfasta_file = join(subdata_dir, 'H3N2_HA_2011_2013_{}_0.fasta'.format(n_taxa))

    with open(subtree_dates_file, 'w') as fp:
        fp.write(content)

    # add dates to end of sequence names
    sub_aln_dic = {}
    for taxon, date in dates_dic.items():
        t = dna.taxon_namespace.get_taxon(label=taxon)
        new_taxon_name = taxon + '_' + date
        sub_aln_dic[new_taxon_name] = str(dna[t])
    sub_dna = DnaCharacterMatrix.from_dict(sub_aln_dic)
    sub_dna.write(path=subfasta_file, schema="fasta")

    new_dates_file = os.path.join(subdata_dir, 'H3N2_HA_2011_2013_{}_0.lsd_dates.txt'.format(n_taxa))
    with open(new_dates_file, 'w') as fp:
        fp.write(str(len(dates_dic)))
        for taxon, date in dates_dic.items():
            fp.write('\n' + taxon + '_' + date + '\t' + date)

stan_script_file = join(subdata_dir, 'H3N2_HA_2011_2013.stan'.format(n_taxa))
stan_pkl = stan_script_file.replace('.stan', '.pkl')

env.Command([stan_script_file, stan_pkl], None,
            'phylostan build -s $TARGET -m JC69 --heterochronous --estimate_rate --clock strict -c constant --compile')

outputs = []

for n_taxa in N_leaves_array:
    stem = join(subdata_dir, 'H3N2_HA_2011_2013_{}_0'.format(n_taxa))
    dates_file = join(dates_dir, 'H3N2_HA_2011_2013_{}_0.lsd_dates.txt'.format(n_taxa))
    subtree_dates_file = stem + '.nwk'  # branch=subst

    lsd_stem = stem + '.lsd.out'
    lsd_tree_dated = lsd_stem + '.date.nexus'  # branch=time
    lsd_dates = stem + '.lsd_dates.txt'

    physher_tree_dated = stem + '.physher.tree'

    seq_file = stem + '.fasta'
    subjson_file = stem + '.json'

    env.Command(target=[stem + '.lsd', lsd_tree_dated],
                source=[subtree_dates_file, lsd_dates],
                action=run_lsd)

    env.Command(target=subjson_file,
                source=[seq_file, lsd_tree_dated, dates_file],
                action=prepare_physher)

    env.Command(target=[stem + '.physher', physher_tree_dated],
                source=subjson_file,
                action=run_physher, iter=iterations)

    env.Command(target=stem + '.phylotorch',
                source=[seq_file, lsd_tree_dated],
                action=run_phylotorch, iter=iterations)

    env.Command(target=stem + '.libsbn',
                source=[seq_file, physher_tree_dated],
                action=run_phylotorch, iter=iterations, libsbn=True)

    env.Command(target=stem + '.phylojax',
                source = [seq_file, physher_tree_dated],
                action = run_phylojax, iter=iterations)

    env.Command(target=stem + '.phylostan',
                source=[seq_file, physher_tree_dated, stan_script_file],
                action=run_phylostan, iter=iterations)

    outputs.extend([stem + '.' + m for m in ('physher', 'phylotorch', 'libsbn', 'phylojax', 'phylostan')])

env.Command(join('flu_H3N2', 'H3N2_HA_2011_2013.csv'), outputs, parse_results)

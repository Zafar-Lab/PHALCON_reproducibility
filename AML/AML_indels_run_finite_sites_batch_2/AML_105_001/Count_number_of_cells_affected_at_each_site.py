
import re
from collections import defaultdict

############################################################
# PARSE TREE
############################################################

def parse_gv_tree(gv_file):
    node_labels = {}
    children = defaultdict(list)
    parent = {}

    with open(gv_file) as f:
        lines = f.readlines()

    for line in lines:
        line = line.strip()

        # node labels
        m = re.match(r'(\d+) \[label="(.*)"\]', line)
        if m:
            node = int(m.group(1))
            label = m.group(2).strip()
            node_labels[node] = label

        # edges
        m = re.match(r'(\d+) -> (\d+)', line)
        if m:
            p = int(m.group(1))
            c = int(m.group(2))
            children[p].append(c)
            parent[c] = p

    all_nodes = set(node_labels.keys())

    root = list(all_nodes - set(parent.keys()))[0]

    return root, node_labels, children

############################################################
# PARSE VCF
############################################################

def parse_vcf(vcf_file):

    mutation_order = []
    cell_genotypes = {}

    with open(vcf_file) as f:

        for line in f:
            line = line.strip()

            if line.startswith('#CHROM'):
                header = line.split('\t')
                cells = header[9:]

                for c in cells:
                    cell_genotypes[c] = {}

            elif not line.startswith('#'):

                fields = line.split('\t')

                chrom = fields[0]
                pos = fields[1]

                chrom = chrom.replace("chr", "")
                mutation = f"{chrom}:{pos}"
                mutation_order.append(mutation)

                samples = fields[9:]

                for cell, gt_field in zip(cells, samples):

                    gt = gt_field.split(':')[0]

                    mutated = 0

                    if gt in ['0/1', '1/0', '1/1']:
                        mutated = 1

                    cell_genotypes[cell][mutation] = mutated

    return mutation_order, cell_genotypes

############################################################
# COMPUTE LEAVES
############################################################

def compute_descendant_leaves(children):

    memo = {}

    def dfs(node):

        if node not in children or len(children[node]) == 0:
            memo[node] = 1
            return 1

        total = 0

        for ch in children[node]:
            total += dfs(ch)

        memo[node] = total
        return total

    roots = set(children.keys())

    for r in roots:
        dfs(r)

    return memo

############################################################
# TREE TRAVERSAL
############################################################

def normalize_mutation(x):
    return x.replace('-', '')


def traverse_tree(
    node,
    node_labels,
    children,
    active_mutations,
    node_configurations
):

    current = set(active_mutations)

    label = node_labels.get(node, '')

    if label != '':

        muts = [x.strip() for x in label.split(',') if x.strip()]

        for m in muts:

            # back mutation
            if m.endswith('-'):
                current.discard(normalize_mutation(m))

            else:
                current.add(m)

    node_configurations[node] = set(current)

    for ch in children.get(node, []):
        traverse_tree(
            ch,
            node_labels,
            children,
            current,
            node_configurations
        )

############################################################
# MATCH CELLS TO CONFIGURATIONS
############################################################

def assign_cells_exactly(
    node_configurations,
    cell_genotypes
):

    tree_mutations = set()

    for config in node_configurations.values():
        tree_mutations.update(config)

    node_counts = {n: 0 for n in node_configurations}

    for cell, geno in cell_genotypes.items():

        cell_config = {
            m for m,v in geno.items()
            if v == 1 and m in tree_mutations
        }

        for node, config in node_configurations.items():

            if cell_config == config:

                node_counts[node] += 1
                break

    return node_counts
############################################################
# MAIN
############################################################

gv_file = 'AML_105_001_indels_dbsnp_nz_finite_sites_nonempty_inferred_tree_exonic.gv'
vcf_file = 'AML_105_001_indels_nonzero_dbsnp_finite_sites_outputInference.vcf'

# parse inputs
root, node_labels, children = parse_gv_tree(gv_file)
mutation_order, cell_genotypes = parse_vcf(vcf_file)

# reconstruct mutation configurations
node_configurations = {}

traverse_tree(
    root,
    node_labels,
    children,
    set(),
    node_configurations
)

# count leaves
leaf_counts = compute_descendant_leaves(children)

# match cells
node_counts = assign_cells_exactly(
    node_configurations,
    cell_genotypes
)

for node,count in node_counts.items():
    print(node, count)

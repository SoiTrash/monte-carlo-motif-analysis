import random
import statistics
import matplotlib.pyplot as plt

from collections import defaultdict

########################################################
# SETTINGS
########################################################

N_SIMULATIONS = 1000

GENOMES = [

    {
        "name": "E. coli K12",
        "fasta": r"C:\Users\soi\Desktop\motif_project\k-mer_analysis\bacteria-procariotes\ecoli_str.K-12_substr.MG1655\GCF_000005845.2_ASM584v2_genomic.fna",
        "motif": "CGCATCCGGC",
        "observed": 150
    },

    {
        "name": "E. coli Sakai",
        "fasta": r"C:\Users\soi\Desktop\motif_project\k-mer_analysis\bacteria-procariotes\ecoli_str.Sakai_substr.RIMD 0509952\GCF_000008865.2_ASM886v2_genomic.fna",
        "motif": "CAGCGCCAGC",
        "observed": 158
    },

    {
        "name": "Mycoplasma",
        "fasta": r"C:\Users\soi\Desktop\motif_project\k-mer_analysis\bacteria-procariotes\Mycoplasmoides_genitalium_str.G37\GCF_040556925.1_ASM4055692v1_genomic.fna",
        "motif": "TTTTTTAAAA",
        "observed": 67
    },

    {
        "name": "Streptomyces",
        "fasta": r"C:\Users\soi\Desktop\motif_project\k-mer_analysis\bacteria-procariotes\streptomyces_coelicolor_str.CFB_NBC_0001\GCF_008931305.1_ASM893130v1_genomic.fna",
        "motif": "GCCGCCGCCG",
        "observed": 1030
    },

    {
        "name": "S. cerevisiae",
        "fasta": r"C:\Users\soi\Desktop\motif_project\k-mer_analysis\fungi-eucariot\Saccharomyces_cervisiae_Str.S288C\GCF_000146045.2_R64_genomic.fna",
        "motif": "AAAAAAAAAA",
        "observed": 4194
    }

]

########################################################
# READ FASTA
########################################################

def read_fasta(path):

    seq = []

    with open(path) as f:

        for line in f:

            if line.startswith(">"):
                continue

            seq.append(line.strip().upper())

    return "".join(seq)

########################################################
# BUILD MARKOV
########################################################

def build_markov(seq):

    transitions = defaultdict(
        lambda: defaultdict(int)
    )

    for i in range(len(seq)-1):

        a = seq[i]
        b = seq[i+1]

        transitions[a][b] += 1

    probs = {}

    for base in transitions:

        total = sum(
            transitions[base].values()
        )

        probs[base] = {}

        for nxt in transitions[base]:

            probs[base][nxt] = (
                transitions[base][nxt] / total
            )

    return probs

########################################################
# MONTE CARLO
########################################################

def monte_carlo_distribution(
    seq,
    motif,
    n_sims
):

    k = len(motif)

    probs = build_markov(seq)

    bases_lookup = {}
    weights_lookup = {}

    for base in probs:

        bases_lookup[base] = list(
            probs[base].keys()
        )

        weights_lookup[base] = list(
            probs[base].values()
        )

    genome_length = len(seq)

    counts = []

    for sim in range(n_sims):

        current = random.choice(
            ["A","T","G","C"]
        )

        window = current

        motif_count = 0

        for _ in range(genome_length - 1):

            current = random.choices(
                bases_lookup[current],
                weights=weights_lookup[current],
                k=1
            )[0]

            window += current

            if len(window) > k:

                window = window[-k:]

            if window == motif:

                motif_count += 1

        counts.append(motif_count)

    return counts

########################################################
# FIGURE
########################################################

fig, axes = plt.subplots(
    3,
    2,
    figsize=(14,12)
)

axes = axes.flatten()

########################################################
# RUN
########################################################

for idx, genome in enumerate(GENOMES):

    print("\n====================")
    print(genome["name"])
    print("====================")

    seq = read_fasta(
        genome["fasta"]
    )

    dist = monte_carlo_distribution(
        seq,
        genome["motif"],
        N_SIMULATIONS
    )

    mean_null = statistics.mean(dist)

    ax = axes[idx]

    ax.hist(
        dist,
        bins=30,
        density=True
    )

    ax.axvline(
        mean_null,
        linestyle="--",
        linewidth=2,
        label=f"Mean={mean_null:.1f}"
    )

    ax.axvline(
        genome["observed"],
        linewidth=3,
        label=f"Obs={genome['observed']}"
    )

    ax.set_title(
        genome["name"]
    )

    ax.legend(fontsize=8)

    # save distribution

    outfile = (
        genome["name"]
        .replace(" ","_")
        + "_distribution.txt"
    )

    with open(outfile,"w") as out:

        for x in dist:

            out.write(
                str(x) + "\n"
            )

# remove empty subplot

fig.delaxes(
    axes[-1]
)

plt.tight_layout()

plt.savefig(
    "five_genome_montecarlo_panel.png",
    dpi=300
)

plt.show()
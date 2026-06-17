import random
import statistics
from collections import defaultdict, Counter

##################################################
# SETTINGS
##################################################

FASTA_FILE = "GCF_000146045.2_R64_genomic.fna"

MOTIFS = [
    "TATTATTATTATTATTATTATTATTATTAT",
    "CAAGGATTGATAATGTAATAGGATCAATGA",
    "GGATTGATAATGTAATAGGATCAATGAATA"
]

N_SIMULATIONS = 1000

##################################################
# READ GENOME
##################################################

sequence = ""

with open(FASTA_FILE, "r") as f:

    for line in f:
        if line.startswith(">"):
            continue

        sequence += line.strip().upper()

print("Genome length:", len(sequence))

##################################################
# OBSERVED COUNTS
##################################################

k = len(MOTIFS[0])

observed_counts = {}

for motif in MOTIFS:

    count = 0

    for i in range(len(sequence) - k + 1):

        if sequence[i:i+k] == motif:
            count += 1

    observed_counts[motif] = count

##################################################
# BUILD MARKOV MODEL
##################################################

transitions = defaultdict(lambda: defaultdict(int))

for i in range(len(sequence)-1):

    current = sequence[i]
    nxt = sequence[i+1]

    transitions[current][nxt] += 1

transition_probs = {}

for base in transitions:

    total = sum(transitions[base].values())

    transition_probs[base] = {}

    for nxt in transitions[base]:

        transition_probs[base][nxt] = (
            transitions[base][nxt] / total
        )

##################################################
# INITIAL BASE DISTRIBUTION
##################################################

base_counts = Counter(sequence)

bases = ["A","T","G","C"]

base_weights = [
    base_counts["A"],
    base_counts["T"],
    base_counts["G"],
    base_counts["C"]
]

##################################################
# LOOKUP TABLES
##################################################

lookup_bases = {}
lookup_weights = {}

for base in transition_probs:

    lookup_bases[base] = list(
        transition_probs[base].keys()
    )

    lookup_weights[base] = list(
        transition_probs[base].values()
    )

##################################################
# STORAGE
##################################################

null_results = {}

for motif in MOTIFS:

    null_results[motif] = []

##################################################
# MONTE CARLO
##################################################

print("\nRunning Monte Carlo simulations...\n")

genome_length = len(sequence)

motif_set = set(MOTIFS)

for sim in range(N_SIMULATIONS):

    current = random.choices(
        bases,
        weights=base_weights,
        k=1
    )[0]

    window = current

    sim_counts = {
        motif: 0
        for motif in MOTIFS
    }

    for _ in range(genome_length - 1):

        current = random.choices(
            lookup_bases[current],
            weights=lookup_weights[current],
            k=1
        )[0]

        window += current

        if len(window) > k:
            window = window[-k:]

        if len(window) == k:

            if window in motif_set:
                sim_counts[window] += 1

    for motif in MOTIFS:

        null_results[motif].append(
            sim_counts[motif]
        )

    if (sim + 1) % 100 == 0:

        print(
            f"{sim+1}/{N_SIMULATIONS} completed"
        )

##################################################
# RESULTS
##################################################

print("\n")
print("="*90)

header = (
    f"{'Motif':<12}"
    f"{'Observed':>10}"
    f"{'Mean':>12}"
    f"{'SD':>12}"
    f"{'Lower95':>12}"
    f"{'Upper95':>12}"
    f"{'Z':>12}"
    f"{'p-value':>12}"
)

print(header)
print("="*90)

results_lines = []

for motif in MOTIFS:

    values = null_results[motif]

    mean_null = statistics.mean(values)

    sd_null = statistics.stdev(values)

    lower95 = mean_null - 1.96 * sd_null
    upper95 = mean_null + 1.96 * sd_null

    z_score = (
        observed_counts[motif]
        - mean_null
    ) / sd_null

    r = sum(
        x >= observed_counts[motif]
        for x in values
    )

    p_value = (
        r + 1
    ) / (
        N_SIMULATIONS + 1
    )

    line = (
        f"{motif:<12}"
        f"{observed_counts[motif]:>10}"
        f"{mean_null:>12.2f}"
        f"{sd_null:>12.2f}"
        f"{lower95:>12.2f}"
        f"{upper95:>12.2f}"
        f"{z_score:>12.2f}"
        f"{p_value:>12.6f}"
    )

    print(line)

    results_lines.append(line)

##################################################
# SAVE RESULTS
##################################################

with open(
    "motif_results_table.txt",
    "w"
) as out:

    out.write(header + "\n")

    for line in results_lines:

        out.write(line + "\n")

print(
    "\nSaved results to motif_results_table.txt"
)
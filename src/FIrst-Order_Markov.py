import random
from collections import Counter, defaultdict

##################################################
# READ GENOME
##################################################

sequence = ""

with open("GCF_000146045.2_R64_genomic.fna", "r") as f:
    for line in f:
        line = line.strip()

        if line.startswith(">"):
            continue

        sequence += line.upper()

print("Genome length:", len(sequence))

##################################################
# REAL GENOME GC CONTENT
##################################################

real_counts = Counter(sequence)

real_gc = (
    real_counts["G"] +
    real_counts["C"]
) / len(sequence)

print("Real GC content:", real_gc)

##################################################
# BUILD TRANSITION MATRIX
##################################################

def build_transition_matrix(seq):

    transitions = defaultdict(lambda: defaultdict(int))

    for i in range(len(seq) - 1):

        current = seq[i]
        nxt = seq[i + 1]

        transitions[current][nxt] += 1

    probs = {}

    for current in transitions:

        total = sum(transitions[current].values())

        probs[current] = {}

        for nxt in transitions[current]:

            probs[current][nxt] = (
                transitions[current][nxt] / total
            )

    return probs

transition_probs = build_transition_matrix(sequence)

##################################################
# PRINT TRANSITION MATRIX
##################################################

print("\nTransition probabilities:\n")

for base in transition_probs:

    print(base)

    for nxt in transition_probs[base]:

        print(
            f"   {base}->{nxt}: "
            f"{transition_probs[base][nxt]:.4f}"
        )

##################################################
# GENERATE MARKOV GENOME
##################################################

def generate_markov_sequence(length, transition_probs):

    current = random.choice(["A", "T", "G", "C"])

    result = [current]

    for _ in range(length - 1):

        next_bases = list(
            transition_probs[current].keys()
        )

        weights = list(
            transition_probs[current].values()
        )

        current = random.choices(
            next_bases,
            weights=weights,
            k=1
        )[0]

        result.append(current)

    return "".join(result)

print("\nGenerating Markov genome...")

markov_seq = generate_markov_sequence(
    len(sequence),
    transition_probs
)

##################################################
# CHECK LENGTH
##################################################

print(
    "Generated genome length:",
    len(markov_seq)
)

##################################################
# CHECK GC CONTENT
##################################################

markov_counts = Counter(markov_seq)

markov_gc = (
    markov_counts["G"] +
    markov_counts["C"]
) / len(markov_seq)

print(
    "Generated GC content:",
    markov_gc
)

##################################################
# COUNT 8-MERS
##################################################

k = 8

kmer_counts = {}

for i in range(len(sequence) - k + 1):

    kmer = sequence[i:i+k]

    kmer_counts[kmer] = (
        kmer_counts.get(kmer, 0) + 1
    )

top = sorted(
    kmer_counts.items(),
    key=lambda x: x[1],
    reverse=True
)

top_motifs = [
    motif
    for motif, count in top[:20]
]

##################################################
# COUNT SAME MOTIFS IN MARKOV GENOME
##################################################

def count_targets(seq, targets, k):

    counts = {m: 0 for m in targets}

    for i in range(len(seq) - k + 1):

        kmer = seq[i:i+k]

        if kmer in counts:

            counts[kmer] += 1

    return counts

markov_motif_counts = count_targets(
    markov_seq,
    top_motifs,
    k
)

##################################################
# COMPARE
##################################################

print("\nTop motif comparison:\n")

for motif in top_motifs:

    print(
        motif,
        "Real:",
        kmer_counts[motif],
        "Markov:",
        markov_motif_counts[motif]
    )
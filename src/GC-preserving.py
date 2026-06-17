import random
from collections import Counter

############################################
# STEP 1 — Read Genome
############################################

sequence = ""

with open("GCF_000005845.2_ASM584v2_genomic.fna", "r") as f:
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            continue
        sequence += line.upper()

print("Genome length:", len(sequence))

############################################
# STEP 2 — Basic Composition
############################################

counts = Counter(sequence)
total = len(sequence)
gc = counts['G'] + counts['C']

print("GC content:", gc / total)

############################################
# STEP 3 — Count All 8-mers
############################################

k = 8
kmer_counts = {}

for i in range(len(sequence) - k + 1):
    kmer = sequence[i:i+k]
    if "N" in kmer:
        continue
    kmer_counts[kmer] = kmer_counts.get(kmer, 0) + 1

print("Total unique 7-mers:", len(kmer_counts))

# Sort by frequency
top = sorted(kmer_counts.items(), key=lambda x: x[1], reverse=True)

print("\nTop 10 observed motifs:")
for motif, count in top[:10]:
    print(motif, count)

############################################
# STEP 4 — Select Top 200 Motifs
############################################

top_motifs = [m for m, c in top[:200]]

############################################
# STEP 5 — Shuffle Function
############################################

def shuffle_sequence(seq):
    seq_list = list(seq)
    random.shuffle(seq_list)
    return "".join(seq_list)

############################################
# STEP 6 — Count Only Target Motifs
############################################

def count_kmers(seq, k, targets):
    counts = {m: 0 for m in targets}
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i+k]
        if kmer in counts:
            counts[kmer] += 1
    return counts

############################################
# STEP 7 — Monte Carlo Simulation
############################################

N = 1000  # increase if desired
sim_results = {m: [] for m in top_motifs}

print("\nRunning Monte Carlo simulations...")

for sim in range(N):
    shuf = shuffle_sequence(sequence)
    counts = count_kmers(shuf, k, top_motifs)

    for m in top_motifs:
        sim_results[m].append(counts[m])

    if (sim + 1) % 100 == 0:
        print(f"Completed {sim+1} simulations")

print("Simulation complete.\n")

############################################
# STEP 8 — Compute Empirical P-values
############################################

results = []

for m in top_motifs:
    observed = kmer_counts[m]
    simulated = sim_results[m]

    p_value = (sum(c >= observed for c in simulated) + 1) / (N + 1)
    sim_mean = sum(simulated) / len(simulated)

    enrichment = observed / sim_mean if sim_mean > 0 else float('inf')

    results.append((m, observed, sim_mean, enrichment, p_value))

# Sort by p-value
results.sort(key=lambda x: x[4])

############################################
# STEP 9 — Print Top Significant Motifs
############################################

print("Top 10 most significant motifs:\n")
for r in results[:10]:
    print(f"Motif: {r[0]}")
    print(f"Observed: {r[1]}")
    print(f"Simulated mean: {r[2]:.2f}")
    print(f"Enrichment: {r[3]:.2f}")
    print(f"Empirical p-value: {r[4]:.6f}")
    print("------")

############################################
# STEP 10 — Save Results to File
############################################

with open("motif_results.txt", "w") as out:
    out.write("Motif\tObserved\tSimMean\tEnrichment\tPvalue\n")
    for r in results:
        out.write(f"{r[0]}\t{r[1]}\t{r[2]:.2f}\t{r[3]:.2f}\t{r[4]:.6f}\n")

print("\nResults saved to motif_results.txt")

sequence = ""

with open("GCF_001695515.1_ASM169551v1_genomic.fna") as f:
    for line in f:
        if not line.startswith(">"):
            sequence += line.strip().upper()

motif = "CGCTGGCG"

positions = []

for i in range(len(sequence)-len(motif)+1):
    if sequence[i:i+len(motif)] == motif:
        positions.append(i)

with open("motif_context_ASM169551v1.txt","w") as out:

    for p in positions:

        start = max(0, p-20)
        end = min(len(sequence), p+len(motif)+20)

        context = sequence[start:end]

        out.write(context + "\n")

print("Saved", len(positions), "contexts")
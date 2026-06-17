TARGET_MOTIF = "TAGTAGTAGTAGTAG"

FASTA_FILE = "GCF_040556925.1_ASM4055692v1_genomic.fna"

##################################################
# READ GENOME
##################################################

sequence = ""

with open(FASTA_FILE, "r") as f:

    for line in f:

        if line.startswith(">"):
            continue

        sequence += line.strip().upper()

##################################################
# FIND POSITIONS
##################################################

positions = []

k = len(TARGET_MOTIF)

for i in range(len(sequence) - k + 1):

    if sequence[i:i+k] == TARGET_MOTIF:

        positions.append(i)

##################################################
# OUTPUT
##################################################

print("Genome length:", len(sequence))

print("Motif:", TARGET_MOTIF)

print("Occurrences:", len(positions))

print("\nPositions:\n")

for p in positions:

    print(p)

##################################################
# SAVE
##################################################

with open(
    "motif_positions.txt",
    "w"
) as out:

    for p in positions:

        out.write(
            str(p) + "\n"
        )

print(
    "\nSaved positions to motif_positions.txt"
)
import os
import re
import math
import matplotlib.pyplot as plt

########################################################
# ROOT DIRECTORY
########################################################

ROOT = r"C:\Users\soi\Desktop\motif_project\k-mer_analysis"

########################################################
# ORGANISMS
########################################################

ORGANISMS = {
    "ecoli_str.K-12_substr.MG1655":
        "E. coli K12",

    "ecoli_str.Sakai_substr.RIMD 0509952":
        "E. coli Sakai",

    "Mycoplasmoides_genitalium_str.G37":
        "Mycoplasma",

    "streptomyces_coelicolor_str.CFB_NBC_0001":
        "Streptomyces",

    "Saccharomyces_cervisiae_Str.S288C":
        "Saccharomyces"
}

########################################################
# K VALUES
########################################################

K_VALUES = [8, 10, 12, 15]

########################################################
# PARSE TOP MOTIF
########################################################

def parse_top_result(filepath):

    with open(filepath, "r", encoding="utf-8") as f:
        text = f.read()

    pattern = r"([ATGC]+)\s+Real:\s+(\d+)\s+Markov:\s+(\d+)"

    match = re.search(pattern, text)

    if not match:
        return None

    motif = match.group(1)
    real = int(match.group(2))
    markov = int(match.group(3))

    return motif, real, markov

########################################################
# COLLECT RESULTS
########################################################

results = {}

for root, dirs, files in os.walk(ROOT):

    folder_name = os.path.basename(root)

    if folder_name not in ORGANISMS:
        continue

    organism = ORGANISMS[folder_name]

    results[organism] = {}

    print("\n===================================")
    print(organism)
    print("===================================")

    for k in K_VALUES:

        filename = f"{k}-mer.txt"
        filepath = os.path.join(root, filename)

        if not os.path.exists(filepath):

            print(f"{filename} nerastas")
            continue

        parsed = parse_top_result(filepath)

        if parsed is None:

            print(f"Nepavyko perskaityti {filename}")
            continue

        motif, real, markov = parsed

        ################################################
        # LOG ENRICHMENT
        ################################################

        enrichment = math.log10(
            (real + 1) /
            (markov + 1)
        )

        results[organism][k] = enrichment

        print(
            f"k={k} "
            f"{motif} "
            f"Real={real} "
            f"Markov={markov} "
            f"log10(E/O)={enrichment:.3f}"
        )

########################################################
# PLOT
########################################################

plt.figure(figsize=(10, 6))

for organism in results:

    xs = []
    ys = []

    for k in sorted(results[organism]):

        xs.append(k)
        ys.append(results[organism][k])

    plt.plot(
        xs,
        ys,
        marker="o",
        linewidth=2,
        markersize=7,
        label=organism
    )

########################################################
# LABELS
########################################################

plt.xlabel(
    "Motyvo ilgis (k)",
    fontsize=12
)

plt.ylabel(
    r"$\log_{10}\left(\frac{Observed+1}{Markov+1}\right)$",
    fontsize=12
)

plt.title(
    "Motyvo ilgio įtaka statistiniam praturtėjimui",
    fontsize=14
)

plt.xticks(K_VALUES)

plt.grid(
    alpha=0.3
)

plt.legend()

plt.tight_layout()

########################################################
# SAVE
########################################################

plt.savefig(
    "motyvo_ilgio_itaka.png",
    dpi=300
)

plt.show()
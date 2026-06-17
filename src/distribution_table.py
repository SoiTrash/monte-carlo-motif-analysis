import numpy as np

##################################################
# INPUT
##################################################

DATA = [

    {
        "organism": r"\textit{E. coli} K12",
        "motif": "CGCATCCGGC",
        "observed": 150,
        "file": "E._coli_K12_distribution.txt"
    },

    {
        "organism": r"\textit{E. coli} Sakai",
        "motif": "CAGCGCCAGC",
        "observed": 158,
        "file": "E._coli_Sakai_distribution.txt"
    },

    {
        "organism": r"\textit{M. genitalium}",
        "motif": "TTTTTTAAAA",
        "observed": 67,
        "file": "Mycoplasma_distribution.txt"
    },

    {
        "organism": r"\textit{S. coelicolor}",
        "motif": "GCCGCCGCCG",
        "observed": 1030,
        "file": "Streptomyces_distribution.txt"
    },

    {
        "organism": r"\textit{S. cerevisiae}",
        "motif": "AAAAAAAAAA",
        "observed": 4194,
        "file": "S._cerevisiae_distribution.txt"
    }

]

##################################################
# HEADER
##################################################

print(r"\begin{table}[H]")
print(r"\centering")
print(r"\caption{Dažniausių 10-mer motyvų statistinis reikšmingumas pagal Monte Karlo simuliacijas}")
print(r"\label{tab:motif_significance}")
print(r"\begin{tabular}{l l r r r r r}")
print(r"\hline")
print(r"Organizmas & Motyvas & Stebėta & MC vidurkis & 95\% PI & Z & p \\")
print(r"\hline")

##################################################
# CALCULATIONS
##################################################

for item in DATA:

    dist = np.loadtxt(item["file"])

    observed = item["observed"]

    mean = np.mean(dist)

    sd = np.std(
        dist,
        ddof=1
    )

    lower95 = np.percentile(
        dist,
        2.5
    )

    upper95 = np.percentile(
        dist,
        97.5
    )

    z = (
        observed - mean
    ) / sd

    extreme = np.sum(
        dist >= observed
    )

    p = (
        extreme + 1
    ) / (
        len(dist) + 1
    )

    if p <= 0.001:

        p_string = r"$<0.001$"

    else:

        p_string = f"{p:.4f}"

    print(
        f"{item['organism']} & "
        f"{item['motif']} & "
        f"{observed} & "
        f"{mean:.1f} & "
        f"[{lower95:.0f}, {upper95:.0f}] & "
        f"{z:.1f} & "
        f"{p_string} \\\\"
    )

##################################################
# FOOTER
##################################################

print(r"\hline")
print(r"\end{tabular}")
print(r"\end{table}")
'''

Each data file is treated as one sub-sample. For every wild-type/mutant pair we
compute the difference in mean count per event in each file, combine the files
into one result, and take the statistical uncertainty from the spread between files.

Usage:
    python3 analysis.py
Edit FILES below if your data files have other names or are in another folder.
'''

FILES = [f"data/output-Set{i}.txt" for i in range(11)]   # output-Set0 ... output-Set10

STRAINS = {
    211: "E. coli WT",
    -211: "E. coli mutant",
    321: "B. subtilis WT",
    -321: "B. subtilis mutant",
    2212: "P. aeruginosa WT",
    -2212: "P. aeruginosa antibiotic-resistant",
    3122: "S. pneumoniae",
    -3122: "S. pneumoniae capsule-deficient",
    3312: "M. tuberculosis",
    -3312: "M. tuberculosis drug-resistant",
    3334: "S. enterica",
    -3334: "S. enterica mutant",
}
PAIRS = [(211, -211), (321, -321), (2212, -2212),
         (3122, -3122), (3312, -3312), (3334, -3334)]


def count_file(path):
    """
    Returns the number of events and the total count of each strain.
    A header line has 2 numbers (event ID, number of bacteria);
    a bacterium line has 4 (px, py, pz, strain ID).
    """
    totals = {code: 0 for code in STRAINS}
    n_events = 0
    with open(path) as f:
        for line in f:
            parts = line.split()
            if len(parts) == 2:
                n_events += 1
            elif len(parts) == 4:
                code = int(parts[3])
                if code in totals:
                    totals[code] += 1
    return n_events, totals


def combine(values, weights):
    """Weighted mean of the sub-sample results, plus their spread.

    spread = standard deviation of the sub-sample results (uncertainty of ONE sub-sample)
    error  = spread / sqrt(N)  (uncertainty of the combined, full-sample result)
    """
    n = len(values)
    mean = sum(v * w for v, w in zip(values, weights)) / sum(weights)
    if n < 2:
        return mean, float("nan"), float("nan")
    avg = sum(values) / n
    spread = (sum((v - avg) ** 2 for v in values) / (n - 1)) ** 0.5
    return mean, spread, spread / n ** 0.5


def main(paths):
    events, per_file = [], []            # per_file[i][code] = mean count per event in file i
    for p in paths:
        n, totals = count_file(p)
        events.append(n)
        per_file.append({code: totals[code] / n for code in STRAINS})
        print(f"  read {p}: {n:,} events")

    n_sub, n_total = len(paths), sum(events)
    lines = [f"Sub-samples: {n_sub}  |  total events: {n_total:,}", ""]

    lines += ["### Wild type vs. mutant", "",
              "| Strain pair | Average difference | Statistical uncertainty |",
              "|---|---:|---:|"]
    for wt, mut in PAIRS:
        diff = [f[wt] - f[mut] for f in per_file]
        dm, _, de = combine(diff, events)
        label = f"{wt} vs {mut}"
        lines.append(f"| {label:<13} | {dm:.8f} | {de:.8f} |")

    text = "\n".join(lines)
    print(text)
    with open("results.md", "w") as out:
        out.write(text + "\n")
    print("\nTable also written to results.md")


main(FILES)

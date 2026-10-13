'''
Each data file is treated as one sub-sample. For every wild-type/mutant pair we
compute, in each file, (a) the difference in mean count per event and (b) the relative
asymmetry A = (N_WT - N_mut) / (N_WT + N_mut), overall and per momentum bin.
The files are combined into one result, and the statistical uncertainty is the
spread between files divided by sqrt(N).

Momentum bins are chosen per pair so that each bin holds roughly the same number
of cells (quantile bins, edges taken from one file). The last bin is open-ended,
so no cell is discarded.

Usage:
    python3 analysis.py
Edit FILES below if your data files have other names or are in another folder.
'''

FILES = [f"data/output-Set{i}.txt" for i in range(1, 11)]   # output-Set1 ... output-Set10 # output-Set0 ... output-Set10
N_BINS = 8                      # momentum bins per pair
EDGE_FILE = FILES[len(FILES) // 2]   # file used to set the bin edges (any full file)
EDGE_STEP = 20                  # use every 20th cell of a pair to find the edges

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


def momentum(parts):
    """Total momentum p = sqrt(px^2 + py^2 + pz^2) from a split bacterium line."""
    px, py, pz = float(parts[0]), float(parts[1]), float(parts[2])
    return (px * px + py * py + pz * pz) ** 0.5


def find_edges(path):
    """
    Bin edges per pair (keyed by the WT code) so that each of the N_BINS bins
    holds about the same number of cells. WT and mutant share the same edges.
    """
    seen = {wt: 0 for wt, _ in PAIRS}
    sample = {wt: [] for wt, _ in PAIRS}
    with open(path) as f:
        for line in f:
            parts = line.split()
            if len(parts) == 4:
                wt = abs(int(parts[3]))
                if wt in sample:
                    seen[wt] += 1
                    if seen[wt] % EDGE_STEP == 0:
                        sample[wt].append(momentum(parts))
    edges = {}
    for wt, values in sample.items():
        values.sort()
        inner = [values[k * len(values) // N_BINS] for k in range(1, N_BINS)]
        edges[wt] = [0.0] + inner + [float("inf")]
    return edges


def count_file(path, edges):
    """
    One pass over a file. Returns the number of events, the total count of each
    strain, and for each pair the WT and mutant counts per momentum bin.
    A header line has 2 numbers (event ID, number of bacteria);
    a bacterium line has 4 (px, py, pz, strain ID).
    """
    totals = {code: 0 for code in STRAINS}
    bins = {wt: ([0] * N_BINS, [0] * N_BINS) for wt, _ in PAIRS}
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
                    p = momentum(parts)
                    e = edges[abs(code)]
                    i = 0
                    while p >= e[i + 1]:      # last edge is infinity, so this stops
                        i += 1
                    bins[abs(code)][0 if code > 0 else 1][i] += 1
    return n_events, totals, bins


def asymmetry(n_wt, n_mut):
    """A = (N_WT - N_mut) / (N_WT + N_mut); None if the bin is empty."""
    total = n_wt + n_mut
    return (n_wt - n_mut) / total if total else None


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


def sig(mean, error):
    return abs(mean) / error if error and error == error and error > 0 else float("nan")


def main(paths):
    print(f"Finding momentum bin edges from {EDGE_FILE} ...")
    edges = find_edges(EDGE_FILE)

    events, per_file, per_file_bins = [], [], []
    for p in paths:
        n, totals, bins = count_file(p, edges)
        events.append(n)
        per_file.append({code: totals[code] / n for code in STRAINS})
        per_file_bins.append(bins)
        print(f"  read {p}: {n:,} events")

    n_sub, n_total = len(paths), sum(events)
    lines = [f"Sub-samples: {n_sub}  |  total events: {n_total:,}", ""]

    # --- overall difference and relative asymmetry ---
    lines += ["### Wild type vs. mutant (all momenta)", "",
              "| Strain pair | Average difference | Uncertainty | Significance | Asymmetry A | Uncertainty | Significance |",
              "|---|---:|---:|---:|---:|---:|---:|"]
    for wt, mut in PAIRS:
        diff = [f[wt] - f[mut] for f in per_file]
        dm, _, de = combine(diff, events)
        asym = [asymmetry(f[wt], f[mut]) for f in per_file]
        am, _, ae = combine(asym, events)
        label = f"{wt} vs {mut}"
        lines.append(f"| {label:<13} | {dm:.8f} | {de:.8f} | {sig(dm, de):.1f}σ "
                     f"| {am:.6f} | {ae:.6f} | {sig(am, ae):.1f}σ |")

    # --- asymmetry per momentum bin ---
    for wt, mut in PAIRS:
        e = edges[wt]
        lines += ["", f"### Asymmetry vs. momentum: {STRAINS[wt]} / {STRAINS[mut]}", "",
                  "| Momentum bin | Cells per file (WT+mut, mean) | Asymmetry A | Uncertainty | Significance |",
                  "|---|---:|---:|---:|---:|"]
        for b in range(N_BINS):
            vals, wts, cells = [], [], []
            for i in range(n_sub):
                w, m = per_file_bins[i][wt][0][b], per_file_bins[i][wt][1][b]
                cells.append(w + m)
                a = asymmetry(w, m)
                if a is not None:
                    vals.append(a)
                    wts.append(events[i])
            high = "inf" if e[b + 1] == float("inf") else f"{e[b + 1]:.3f}"
            label = f"{e[b]:.3f} - {high}"
            avg_cells = sum(cells) / len(cells)
            if len(vals) < 2:
                lines.append(f"| {label} | {avg_cells:.0f} | n/a | n/a | n/a |")
                continue
            am, _, ae = combine(vals, wts)
            lines.append(f"| {label} | {avg_cells:.0f} | {am:.5f} | {ae:.5f} | {sig(am, ae):.1f}σ |")

    text = "\n".join(lines)
    print(text)
    with open("results.md", "w") as out:
        out.write(text + "\n")
    print("\nTable also written to results.md")


if __name__ == "__main__":
    main(FILES)

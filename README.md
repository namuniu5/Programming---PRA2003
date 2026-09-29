## PRA2003-Programming Małgorzata Fic(i6362471)

Statistical analysis of 5 million simulated bacterial-tracking experiments: do wild-type and mutant bacterial strains occur in equal numbers per experiment?

---

## 1. Context and scientific question

Mutations such as antibiotic resistance or loss of a capsule often carry a **fitness cost**: the mutant grows or survives less well than its wild type when no selective pressure is present. If that is true under the simulated conditions here, mutant strains should appear *less often* per experiment than their wild-type counterparts.

A single experiment contains only a few dozen cells, so counts fluctuate strongly from run to run. Any real WT/mutant difference is small compared to that noise, which is why the analysis uses the full sample of 5M experiments and quotes every result with a statistical uncertainty.

**Research questions**

1. What is the average difference in the number of wild-type and mutant cells per experiment, and what is its statistical uncertainty?
2. Is that asymmetry statistically significant for each strain?
3. *(Planned)* Does any asymmetry depend on the momentum of the bacteria?

## 2. Data

Each input file holds 500,000 experiments (~770 MB). Every experiment starts with a header line, followed by one line per tracked bacterium:

```
<event ID> <number of bacteria>
<px> <py> <pz> <strain ID>      # momentum in 10⁻²⁰ kg·m/s
...
```

| ID | Strain | ID | Strain |
|---:|---|---:|---|
| 211 | *E. coli* WT | −211 | *E. coli* mutant |
| 321 | *B. subtilis* WT | −321 | *B. subtilis* mutant |
| 2212 | *P. aeruginosa* WT | −2212 | *P. aeruginosa* antibiotic-resistant |
| 3122 | *S. pneumoniae* | −3122 | *S. pneumoniae* capsule-deficient |
| 3312 | *M. tuberculosis* | −3312 | *M. tuberculosis* drug-resistant |
| 3334 | *S. enterica* | −3334 | *S. enterica* mutant |

## 3. Method

- **Sub-sampling.** Each data file is one sub-sample (N = 11 files, `output-Set0.txt` to `output-Set10.txt`).
- **Per sub-sample.** For every strain: mean count per experiment = total cells of that strain ÷ number of experiments in the file.
- **Central value.** Weighted average of the sub-sample means (weight = number of experiments in the file).
- **Statistical uncertainty.** Standard deviation σ of the sub-sample results divided by √N. σ describes how much the result of a *single* sub-sample scatters; dividing by √N gives the uncertainty of the combined full-sample result (central limit theorem).
- **WT vs. mutant.** The difference ⟨N<sub>WT</sub>⟩ − ⟨N<sub>mutant</sub>⟩ per experiment is computed *inside each sub-sample* and then combined the same way, so correlations between the two strains are handled automatically. A positive difference means the wild type is more abundant.
- **Significance.** Difference ÷ uncertainty, expressed in σ. Below 2σ is treated as consistent with no asymmetry; 3σ or more as significant.
- Systematic uncertainties are not evaluated.

## 4. Installation

Requirements: **Python 3.8+**. The script uses no imports at all, so nothing needs to be installed.

1. Download this repository: click the green **Code** button at the top of this page, then **Download ZIP**, and unzip it. If you use git, clone it with the link from the same button.
2. Check your Python version:

```bash
python3 --version        # should be 3.8 or newer
```
3. Download the 11 data files (`output-Set0.txt` to `output-Set10.txt`) from SURFdrive (link on Canvas) and put them in a folder called `data/` inside the repository:

```
<repo>/
├── analysis.py
├── README.md
└── data/
    ├── output-Set0.txt
    ├── ...
    └── output-Set10.txt
```

The data files are not included in this repository because of their size (~770 MB each).

## 5. Usage

```bash
python3 analysis.py
```

The list of input files is set in the `FILES` variable at the top of `analysis.py`. Edit it if your files have different names or are in a different folder.

The script reads the files line by line, so memory use stays low. It prints progress for each file, then the result table.

## 6. Results

| Strain pair | Average difference | Statistical uncertainty | Significance |
|---|---:|---:|---:|
| 211 vs -211   | 0.03229671 | 0.00451870 | 7.1σ |
| 321 vs -321   | 0.00569100 | 0.00329367 | 1.7σ |
| 2212 vs -2212 | 0.02387314 | 0.00237277 | 10.1σ |
| 3122 vs -3122 | 0.00490256 | 0.00058316 | 8.4σ |
| 3312 vs -3312 | 0.00044067 | 0.00048793 | 0.9σ |
| 3334 vs -3334 | 0.00003551 | 0.00006854 | 0.5σ |

### Interpretation

All six differences are positive: in every pair the wild type is at least as abundant as its mutant.

- **Significant asymmetry (> 5σ):** *E. coli* (211 / −211), *P. aeruginosa* (2212 / −2212) and *S. pneumoniae* (3122 / −3122). Their mutants are clearly depleted, consistent with a fitness cost of the mutation.
- **No significant asymmetry:** *B. subtilis* (1.7σ) shows at most a weak hint. *M. tuberculosis* (0.9σ) and *S. enterica* (0.5σ) are consistent with zero, but these strains are so rare per experiment that the sample cannot resolve a difference of the size seen in the other pairs.

The hypothesis of a fitness cost is supported for three of the six strains and cannot be tested conclusively for the other three at this sample size.


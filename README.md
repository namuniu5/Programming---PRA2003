## PRA2003-Programming Małgorzata Fic(i6362471)

Statistical analysis of 5 million simulated bacterial-tracking experiments: do wild-type and mutant bacterial strains occur in equal numbers per experiment, and does any difference depend on the momentum of the bacteria?

---

## 1. Context and scientific question

Mutations such as antibiotic resistance or loss of a capsule often carry a **fitness cost**: the mutant grows or survives less well than its wild type when no selective pressure is present. If that is true under the simulated conditions here, mutant strains should appear *less often* per experiment than their wild-type counterparts.

A single experiment contains only a few dozen cells, so counts fluctuate strongly from run to run. Any real WT/mutant difference is small compared to that noise, which is why the analysis uses the full sample of 5 million experiments and quotes every result with a statistical uncertainty.

**Research questions**

1. What is the average difference in the number of wild-type and mutant cells per experiment, and what is its statistical uncertainty?
2. Is that asymmetry statistically significant for each strain?
3. Does the asymmetry depend on the momentum of the bacteria, $p = \sqrt{p_x^2 + p_y^2 + p_z^2}$?

## 2. Data

Each input file holds 500,000 experiments (about 0.8 GB). Every experiment starts with a header line, followed by one line per tracked bacterium:

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

The analysis uses `output-Set1.txt` to `output-Set10.txt`: **ten files, 5,000,000 experiments in total**. `output-Set0.txt` is a test file with a single experiment (29 bacteria) and is not used.

## 3. Method

### 3.1 Quantities

For each of the six WT/mutant pairs two quantities are computed:

- **Difference** ⟨N<sub>WT</sub>⟩ − ⟨N<sub>mut</sub>⟩: the average number of WT cells per experiment minus the average number of mutant cells per experiment (Research Questions 1 and 2).
- **Relative asymmetry** $A = \dfrac{N_\text{WT} - N_\text{mut}}{N_\text{WT} + N_\text{mut}}$: the same difference divided by the total number of cells of the pair. $A$ is dimensionless, runs from −1 (only mutants) to +1 (only WT), and is 0 for equal numbers. The raw difference scales with how common a strain is, so it cannot be compared between abundant strains (*E. coli*, about 18 cells per experiment) and rare ones (*S. enterica*, about 0.001). $A$ can, and it is also the only sensible quantity inside a momentum bin, where the number of cells varies strongly from bin to bin. For $A \ll 1$, the WT is more abundant than the mutant by a factor of about $1 + 2A$.

A positive value means the wild type is more abundant.

### 3.2 Sub-sampling and statistical uncertainty

The uncertainties come from the **sub-sampling method**, which does not assume a particular distribution of the counts (such as Poisson statistics). It works as follows.

1. **Sub-samples.** Each of the ten files is analysed on its own. Each file is one sub-sample of $n = 500{,}000$ independent experiments.
2. **One result per file.** In file $k$, the total count of each strain divided by $n$ gives its mean count per experiment. From these the difference $d_k$ and the asymmetry $A_k$ of each pair are computed *inside the file*. This gives ten values $x_1, \dots, x_{10}$ of each quantity $x$.
3. **Central value.** The mean of the ten values, $\bar{x} = \frac{1}{N}\sum_k x_k$ with $N = 10$. All files have the same number of experiments, so no weighting is needed (the code weights by the number of experiments, and the weights are all equal).
4. **Spread.** The sample standard deviation of the ten values (with $N-1$ in the denominator):
   $$\sigma = \sqrt{\frac{1}{N-1}\sum_{k=1}^{N}(x_k - \bar{x})^2}$$
   $\sigma$ describes how much the result of **one** sub-sample (500,000 experiments) fluctuates.
5. **Statistical uncertainty of the result.**
   $$\delta = \frac{\sigma}{\sqrt{N}}$$
   The central value is the average of $N$ independent sub-samples, so it fluctuates less than a single sub-sample by a factor $\sqrt{N}$ (the standard error of the mean). $\delta$ is the uncertainty of the result for the full sample of 5,000,000 experiments, and it is the uncertainty quoted everywhere in this README. For example, the ten values of $A$ for *E. coli* have $\sigma \approx 1.1 \times 10^{-4}$, so $\delta = \sigma/\sqrt{10} \approx 3.6 \times 10^{-5}$.
6. **Significance.** $|\bar{x}|/\delta$, the distance of the result from 0 in units of its uncertainty. The convention used here: **≥ 3σ significant, 2σ to 3σ inconclusive, below 2σ not significant** (consistent with no asymmetry).

**Why the difference is computed inside each file.** WT and mutant cells of a pair come from the same experiments, so their counts fluctuate together (an experiment with many cells tends to contain more of both). Taking the difference or $A$ within each file before combining lets these correlated fluctuations cancel. Combining separate uncertainties of the WT and mutant counts would ignore this and overstate the uncertainty.

**Assumptions and limits of the method**

- It assumes the ten files are independent and statistically identical. This holds well for the whole-sample results, but not for the highest momentum bins (Section 3.3).
- With only ten sub-samples, $\sigma$ itself is only known to about $1/\sqrt{2(N-1)} \approx 24\%$. A significance of 2.9σ and one of 3.1σ are therefore not meaningfully different.
- Systematic uncertainties (for example from the simulation) are not evaluated.

**Which uncertainty to compare with.** If the plain standard deviation $\sigma$ of the ten files is quoted as the uncertainty, all uncertainties in this README become $\sqrt{10} \approx 3.16$ times larger and all significances 3.16 times smaller. The asymmetries themselves do not change. Under that convention *B. subtilis* would be 1.7σ and *M. tuberculosis* 0.9σ instead of 5.5σ and 2.9σ. Those two verdicts therefore depend on this choice. $\sigma/\sqrt{N}$ is used here because the result reported is the average over all ten files.

### 3.3 Momentum dependence

For every cell of a pair the total momentum $p = \sqrt{p_x^2 + p_y^2 + p_z^2}$ is calculated, and the cell is placed in one of eight momentum bins. In every bin, $A$ is calculated separately **for each file**, and the ten values are combined exactly as in Section 3.2 (mean, $\sigma$, $\delta = \sigma/\sqrt{10}$). A file with no cells of the pair in a bin is skipped for that bin.

**Bin edges.** The bins have different widths, because the edges are the 1/8, 2/8, … 7/8 quantiles of the momentum of the pair in `output-Set6.txt` (found by sampling every 20th cell of the pair). This gives bins that hold about the same number of cells, so every bin has comparable statistics. The last bin is open-ended, so no cell is discarded. WT and mutant of a pair share the same edges. Momentum is in units of 10⁻²⁰ kg·m/s.

| Pair | Inner bin edges (the first bin starts at 0, the last is open-ended) |
|---|---|
| *E. coli* | 0.450, 0.724, 1.110, 1.710, 2.696, 4.379, 7.859 |
| *B. subtilis* | 0.674, 1.093, 1.661, 2.498, 3.814, 6.016, 10.620 |
| *P. aeruginosa* | 0.909, 1.481, 2.225, 3.290, 4.936, 7.746, 13.511 |
| *S. pneumoniae* | 1.002, 1.654, 2.452, 3.663, 5.490, 8.595, 15.022 |
| *M. tuberculosis* | 1.199, 1.882, 2.781, 4.141, 6.213, 9.683, 17.067 |
| *S. enterica* | 1.300, 1.737, 2.296, 3.826, 5.838, 12.008, 25.991 |

**Check for unreliable bins.** For each bin the uncertainty $\delta$ is compared with the uncertainty expected from counting statistics alone, $\sqrt{(1-A^2)/N_\text{cells}}$, where $N_\text{cells}$ is the number of cells in the bin over all ten files (estimated as ten times the mean number of cells per file). If $\delta$ is more than twice this value, the ten files differ from each other in that bin by more than random counting allows, so the assumption of identical sub-samples does not hold and $\delta$ is not a reliable measure of the uncertainty. Such bins are marked with **†** in the tables and drawn as hollow markers in the figures, and they are not interpreted. This applies to the highest one or two bins of five of the six pairs (and to the lowest bin of *E. coli*). The cause was not investigated, because the per-file counts per bin were not saved.

## 4. Installation

Requirements: **Python 3.8+**. `analysis.py` uses no imports at all, so nothing needs to be installed for it. `plot_from_results.py` needs **matplotlib**:

```bash
python3 -m pip install matplotlib
```

1. Download this repository: click the green **Code** button at the top of this page, then **Download ZIP**, and unzip it. If you use git, clone it with the link from the same button.
2. Check your Python version:

```bash
python3 --version        # should be 3.8 or newer
```
3. Download the data files `output-Set1.txt` to `output-Set10.txt` from SURFdrive (link on Canvas) and put them in a folder called `data/` inside the repository:

```
<repo>/
├── analysis.py
├── plot_from_results.py
├── README.md
├── asymmetry_*.png          (figures used below)
└── data/
    ├── output-Set1.txt
    ├── ...
    └── output-Set10.txt
```

The data files are not included in this repository because of their size (about 0.8 GB each).

## 5. Usage

Run both scripts from a terminal inside the repository folder.

**Step 1: analysis.** Reads the data and writes the result tables.

```bash
python3 analysis.py
```

The list of input files is set in the `FILES` variable at the top of `analysis.py` (`output-Set1.txt` to `output-Set10.txt`). Edit it if your files have different names or are in another folder. The script reads each file line by line, so memory use stays low, and counts the strains and fills the momentum bins in the same pass. It first finds the bin edges from `output-Set6.txt`, then reads all ten files (this takes a while: tens of minutes with plain Python), prints progress for each file, and writes all tables to `results.md`. `N_BINS` (default 8) and `EDGE_STEP` (default 20) at the top of the script set the number of momentum bins and the sampling used to find the edges.

**Step 2: figures.** Reads `results.md` (not the data) and makes the figures.
```bash
python3 plot_from_results.py results.md
```

It writes `asymmetry_<WT ID>.png` for each pair and `asymmetry_all_pairs.png` (all six pairs in one figure) next to `results.md`. Without an argument it reads `~/Desktop/results.md`.

## 6. Results

Average over the full sample of 5,000,000 experiments, from the ten sub-samples. Uncertainties are $\sigma/\sqrt{10}$ (Section 3.2). The "WT excess" is the factor by which the WT outnumbers the mutant, about $2A$.

| Strain pair | Average difference per experiment | Significance | Asymmetry $A$ | Significance | WT excess |
|---|---:|---:|---:|---:|---:|
| 211 vs −211 (*E. coli*) | 0.02983 ± 0.00132 | 22.6σ | 0.000810 ± 0.000036 | 22.6σ | 0.16 % |
| 321 vs −321 (*B. subtilis*) | 0.00526 ± 0.00096 | 5.5σ | 0.001136 ± 0.000208 | 5.5σ | 0.23 % |
| 2212 vs −2212 (*P. aeruginosa*) | 0.02205 ± 0.00069 | 31.8σ | 0.009980 ± 0.000315 | 31.7σ | 2.0 % |
| 3122 vs −3122 (*S. pneumoniae*) | 0.00453 ± 0.00017 | 26.6σ | 0.008941 ± 0.000335 | 26.7σ | 1.8 % |
| 3312 vs −3312 (*M. tuberculosis*) | 0.00041 ± 0.00014 | 2.9σ | 0.005630 ± 0.001972 | 2.9σ | 1.1 % |
| 3334 vs −3334 (*S. enterica*) | 0.000033 ± 0.000020 | 1.6σ | 0.015341 ± 0.009266 | 1.7σ | 3.1 % |

The difference and $A$ give the same significance to within 0.1σ. The full momentum tables are in the appendix (Section 9).

### Interpretation

**Research Questions 1 and 2.** All six differences are positive: in every pair the wild type is on average more abundant than its mutant.

- **Significant (≥ 3σ):** *E. coli* (22.6σ), *P. aeruginosa* (31.7σ), *S. pneumoniae* (26.7σ) and *B. subtilis* (5.5σ). Their mutants are depleted, consistent with a fitness cost of the mutation. The size of the effect differs strongly: about 2 % more WT for *P. aeruginosa* and *S. pneumoniae*, but only 0.16 % for *E. coli* and 0.23 % for *B. subtilis*. The *E. coli* effect is tiny, and it can only be measured because *E. coli* is by far the most common strain, so its uncertainty is very small.
- **Inconclusive (2σ to 3σ):** *M. tuberculosis* (2.9σ). The asymmetry is positive and of a similar size to the one in *P. aeruginosa*, but the data cannot confirm it.
- **Not significant (< 2σ):** *S. enterica* (1.7σ). The strain is about 1,000 times rarer than *P. aeruginosa*, so the uncertainty (±1.9 %) is as large as the effects seen in the other pairs, and a difference of that size could not have been detected.

The fitness-cost hypothesis is supported for four of the six strains and cannot be tested conclusively for the other two at this sample size. As noted in Section 3.2, the verdicts for *B. subtilis* and *M. tuberculosis* depend on quoting $\sigma/\sqrt{N}$ and not $\sigma$.

**Research Question 3.** The figures below show $A$ per momentum bin. No pair has a clearly negative bin: the few negative values are less than 1.5σ from zero or lie in unreliable bins.

- ***P. aeruginosa*** is the only pair with a hint of momentum dependence. $A$ is flat at 0.0075–0.0094 in the first five bins (each 6σ to 11σ from zero) and rises to 0.0129 ± 0.0016 in the 4.9–7.7 bin, about 2.9σ above the lowest bin. This is moderate evidence for an increase with momentum, not proof: eight bins were examined, and the two highest bins are unreliable (†).
- ***E. coli*** rises slightly, from 0.00062 ± 0.00020 (1.7–2.7) to 0.00133 ± 0.00011 (4.4–7.9), a difference of about 3σ, but the effect is very small in absolute terms.
- ***S. pneumoniae*** is flat at 0.0055–0.0097 in seven bins, with no trend.
- ***B. subtilis*, *M. tuberculosis* and *S. enterica*:** no bin is significantly different from the overall result, and no trend can be seen. For *S. enterica* only 67 to 209 cells per file fall in each bin, so only a very large change of $A$ with momentum could have been detected.

**Limitations.** In five of the six pairs the highest one or two bins are unreliable (†), and so is the lowest *E. coli* bin: the ten files scatter much more there than counting statistics predict (for example 6 times more for *P. aeruginosa* 7.7–13.5 and about 26 times more for its top bin), so the quoted uncertainty is dominated by differences between files and not by counting. Nothing can be concluded from those bins. Apart from that, the main limit is statistical: for the rare strains there are too few cells to resolve an asymmetry of about 1 %.

## 7. Figures

All figures show the relative asymmetry $A = (N_\text{WT} - N_\text{mut})/(N_\text{WT} + N_\text{mut})$ (dimensionless) against the momentum bin, with the bin edges in units of 10⁻²⁰ kg·m/s on the x-axis. The bins have roughly equal numbers of cells, so they are drawn with equal spacing. Each point is the mean of the ten files, and the error bars are $\sigma/\sqrt{10}$ of the ten per-file values of $A$ (Section 3.2). The blue band is the overall $A$ of the pair $\pm$ its uncertainty, and the grey line marks $A = 0$. Hollow markers are unreliable bins (†, Section 3.3). A triangle at the edge of the plot marks an unreliable bin whose value lies outside the plotted range, with its value written next to it. The vertical range of each plot is set by the reliable bins.

<img width="1600" height="1000" alt="image" src="https://github.com/user-attachments/assets/823bbe81-0604-4ace-aebf-e242d94cb552" />
**Figure 1.** *P. aeruginosa* WT vs antibiotic-resistant. $A$ is positive in all bins and flat at 0.0075–0.0094 up to $p \approx 4.9$ (each 6σ to 11σ from zero), then rises to 0.0129 ± 0.0016 in the 4.9–7.7 bin, about 2.9σ above the lowest bin. The two highest bins (0.0178 ± 0.0054 and 0.038 ± 0.024, triangles) are unreliable, so the rise above $p \approx 8$ cannot be confirmed.

<img width="1600" height="1000" alt="asymmetry_211" src="https://github.com/user-attachments/assets/03045e6d-c2db-4718-aa07-4413fe2c8824" />
**Figure 2.** *E. coli* WT vs mutant. $A$ is positive in every bin and small (0.0001 to 0.0013 in the reliable bins). It increases slightly with momentum, with the highest reliable bin (0.00133 ± 0.00011) about 3σ above the 1.7–2.7 bin. The lowest bin (hollow) and the highest bin (triangle, 0.0049 ± 0.0031) are unreliable.


<img width="1600" height="1000" alt="image" src="https://github.com/user-attachments/assets/3a7e83a2-17bc-4537-bb4d-f17d4d5673e4" />
**Figure 3.** *S. pneumoniae* WT vs capsule-deficient. The seven reliable bins scatter around the overall value (0.0089) between 0.0055 and 0.0097, with no trend with momentum. The highest bin (−0.0045 ± 0.018, triangle) is unreliable.


<img width="1600" height="1000" alt="image" src="https://github.com/user-attachments/assets/74d5a1e3-6f8e-40f2-9256-dbdaa7a628c1" />
**Figure 4.** *B. subtilis* WT vs mutant. The overall asymmetry is small (0.0011 ± 0.0002) and the individual bins are mostly less than 1.5σ from zero. Only the 3.8–6.0 bin (0.0021 ± 0.0006, 3.4σ) stands out, and the 1.1–1.7 bin is slightly negative. No trend with momentum can be identified. The two highest bins are unreliable.


<img width="1600" height="1000" alt="image" src="https://github.com/user-attachments/assets/9556bf83-c0fe-4fd8-abc0-6e0b57352294" />
**Figure 5.** *M. tuberculosis* WT vs drug-resistant. Every reliable bin has a large error bar (about 0.004 to 0.006) and is at most 2σ from zero, consistent with the small positive overall value (0.0056 ± 0.0020). The two highest bins (about 4,000 cells per file) are unreliable; the highest (0.104 ± 0.100) is off the scale.


<img width="1600" height="1000" alt="asymmetry_3334" src="https://github.com/user-attachments/assets/829b42e1-f611-4d1b-8afb-d9b18a9c34ee" />
**Figure 6.** *S. enterica* WT vs mutant. With only 67 to 209 cells per bin and file, the error bars (0.02 to 0.05) are as large as or larger than the overall asymmetry (0.015 ± 0.009). All bins are consistent with it and with zero (at most 1.8σ), and the values scatter without a pattern. Only a very large dependence on momentum could have been detected.

All six pairs together are in `asymmetry_all_pairs.png`.

## 8. Use of AI

AI assisted with extending the analysis and structuring this README. I ran everything on the course data myself, and all numbers and figures here come from those runs.

## 9. Appendix: momentum tables

† marks an unreliable bin (Section 3.3). Momentum bins are in units of 10⁻²⁰ kg·m/s, "cells" is the mean number of WT + mutant cells per file in the bin, and the uncertainty is $\sigma/\sqrt{10}$.

<details>
<summary><i>E. coli</i> (211 / −211)</summary>

| Momentum bin | Cells per file | Asymmetry $A$ | Significance |
|---|---:|---:|---:|
| 0.000 – 0.450 † | 2,075,968 | 0.00084 ± 0.00069 | 1.2σ |
| 0.450 – 0.724 | 2,093,007 | 0.00011 ± 0.00023 | 0.5σ |
| 0.724 – 1.110 | 2,141,995 | 0.00084 ± 0.00022 | 3.9σ |
| 1.110 – 1.710 | 2,261,418 | 0.00067 ± 0.00012 | 5.5σ |
| 1.710 – 2.696 | 2,579,492 | 0.00062 ± 0.00020 | 3.1σ |
| 2.696 – 4.379 | 3,016,434 | 0.00099 ± 0.00017 | 5.8σ |
| 4.379 – 7.859 | 2,166,849 | 0.00133 ± 0.00011 | 11.7σ |
| 7.859 – ∞ † | 2,075,260 | 0.00490 ± 0.00314 | 1.6σ |

</details>

<details>
<summary><i>B. subtilis</i> (321 / −321)</summary>

| Momentum bin | Cells per file | Asymmetry $A$ | Significance |
|---|---:|---:|---:|
| 0.000 – 0.674 | 261,995 | 0.00071 ± 0.00059 | 1.2σ |
| 0.674 – 1.093 | 268,361 | 0.00025 ± 0.00054 | 0.5σ |
| 1.093 – 1.661 | 279,621 | −0.00101 ± 0.00066 | 1.5σ |
| 1.661 – 2.498 | 307,066 | 0.00072 ± 0.00086 | 0.8σ |
| 2.498 – 3.814 | 367,198 | 0.00082 ± 0.00066 | 1.2σ |
| 3.814 – 6.016 | 305,579 | 0.00207 ± 0.00061 | 3.4σ |
| 6.016 – 10.620 † | 260,854 | −0.00007 ± 0.00133 | 0.0σ |
| 10.620 – ∞ † | 264,144 | 0.00386 ± 0.00133 | 2.9σ |

</details>

<details>
<summary><i>P. aeruginosa</i> (2212 / −2212)</summary>

| Momentum bin | Cells per file | Asymmetry $A$ | Significance |
|---|---:|---:|---:|
| 0.000 – 0.909 | 126,328 | 0.00791 ± 0.00074 | 10.7σ |
| 0.909 – 1.481 | 128,826 | 0.00813 ± 0.00134 | 6.1σ |
| 1.481 – 2.225 | 140,160 | 0.00882 ± 0.00101 | 8.8σ |
| 2.225 – 3.290 | 163,620 | 0.00747 ± 0.00076 | 9.8σ |
| 3.290 – 4.936 | 171,401 | 0.00943 ± 0.00115 | 8.2σ |
| 4.936 – 7.746 | 127,275 | 0.01294 ± 0.00155 | 8.3σ |
| 7.746 – 13.511 † | 124,735 | 0.01781 ± 0.00538 | 3.3σ |
| 13.511 – ∞ † | 122,368 | 0.03828 ± 0.02353 | 1.6σ |

</details>

<details>
<summary><i>S. pneumoniae</i> (3122 / −3122)</summary>

| Momentum bin | Cells per file | Asymmetry $A$ | Significance |
|---|---:|---:|---:|
| 0.000 – 1.002 | 28,342 | 0.00862 ± 0.00189 | 4.6σ |
| 1.002 – 1.654 | 29,925 | 0.00546 ± 0.00172 | 3.2σ |
| 1.654 – 2.452 | 32,058 | 0.00674 ± 0.00120 | 5.6σ |
| 2.452 – 3.663 | 40,684 | 0.00725 ± 0.00165 | 4.4σ |
| 3.663 – 5.490 | 36,809 | 0.00974 ± 0.00140 | 7.0σ |
| 5.490 – 8.595 | 28,798 | 0.00725 ± 0.00186 | 3.9σ |
| 8.595 – 15.022 | 28,405 | 0.00970 ± 0.00221 | 4.4σ |
| 15.022 – ∞ † | 28,181 | −0.00452 ± 0.01809 | 0.2σ |

</details>

<details>
<summary><i>M. tuberculosis</i> (3312 / −3312)</summary>

| Momentum bin | Cells per file | Asymmetry $A$ | Significance |
|---|---:|---:|---:|
| 0.000 – 1.199 | 4,400 | 0.00545 ± 0.00517 | 1.1σ |
| 1.199 – 1.882 | 3,962 | 0.00730 ± 0.00595 | 1.2σ |
| 1.882 – 2.781 | 4,778 | −0.00199 ± 0.00438 | 0.5σ |
| 2.781 – 4.141 | 6,110 | 0.00569 ± 0.00458 | 1.2σ |
| 4.141 – 6.213 | 4,789 | 0.00649 ± 0.00360 | 1.8σ |
| 6.213 – 9.683 | 4,096 | 0.01289 ± 0.00646 | 2.0σ |
| 9.683 – 17.067 † | 4,081 | −0.00759 ± 0.01796 | 0.4σ |
| 17.067 – ∞ † | 4,008 | 0.10420 ± 0.09973 | 1.0σ |

</details>

<details>
<summary><i>S. enterica</i> (3334 / −3334)</summary>

| Momentum bin | Cells per file | Asymmetry $A$ | Significance |
|---|---:|---:|---:|
| 0.000 – 1.300 | 117 | 0.04481 ± 0.02821 | 1.6σ |
| 1.300 – 1.737 | 67 | 0.02130 ± 0.04638 | 0.5σ |
| 1.737 – 2.296 | 84 | 0.01083 ± 0.03844 | 0.3σ |
| 2.296 – 3.826 | 209 | 0.03556 ± 0.02022 | 1.8σ |
| 3.826 – 5.838 | 172 | 0.00685 ± 0.02223 | 0.3σ |
| 5.838 – 12.008 | 204 | 0.03005 ± 0.02773 | 1.1σ |
| 12.008 – 25.991 | 148 | −0.01714 ± 0.03147 | 0.5σ |
| 25.991 – ∞ | 79 | 0.01924 ± 0.04058 | 0.5σ |

</details>

---

**Course:** Programming (PRA2003) | Maastricht University
**Author:** Małgorzata Fic (i6362471)
**Date:** October 2026

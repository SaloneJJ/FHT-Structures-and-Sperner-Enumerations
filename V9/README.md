# Fractal Hyper-Trees — programs and logs of the trilogy, version 9

Archived version (Zenodo): https://doi.org/10.5281/zenodo.24296051 — the program `fht.py` at the root of the repository implements an earlier definition (v4.0.0); the code of the trilogy is in this folder.

This folder contains the programs, logs and pre-registered criteria cited in the three articles
(version 9, November 2026), signed Jean-Jacques Salone:

- **I.** Combinatorial foundations: definition, enumeration, invariants, morphisms and the algebra C ⊕ H ⊕ M3(C);
- **II.** The finite spectral triple of the Standard Model from flat FHTs: K26, the descent from the top algebra, closures, flavour and masses;
- **III.** The continuum from above: the three-sphere and the holonomy theorem, space as a Gromov boundary, chirality, native spin and Dirac operator, masses on the boundary and the spectral action.

## Contents
- `code/`: one entry per name cited in the articles (scripts, folders of scripts with their logs, data files), plus their local dependencies. The table `INVENTAIRE.md` gives, for each name cited in an article, the corresponding file here and its place in the author's working tree.
- `code/FHT_Enum_results.json`: the counts f(v) of Article I for v ≤ 9 (with the counts by depth), as written by `FHT_Enumeration.py`.
- `code/V52_paysage_v5_f30.json.xz`: catalogue of the flat FHTs of the landscape (xz-compressed; `xz -d` to restore), which contains the class members K_27_1021, …, K_36_44 of Article II.
- `code/verifications_FHT29/` (`verif_EDSU`, `tests`, `calculs`): the confrontations with NuFIT 6.1 and the numerical tests of Article III, with their criteria written before the computations. The NuFIT 6.1 tables are not redistributed; download them from http://www.nu-fit.org.

## Running
Python 3.9 with numpy and scipy (`requirements.txt`); `pynauty` for some automorphism computations.
**Caveat:** many scripts were written inside the author's working tree and refer to it by relative or absolute paths
(e.g. `../FHT_28/...`). They are deposited as run, for traceability; to rerun one, place its dependencies as indicated
in `INVENTAIRE.md` or adapt the path at the top of the script. Logs are the outputs of the runs cited in the articles.

## Method
Each computation cited in the articles has its success criterion written before it was run; any change made after a
first result is declared in the corresponding criteria or results file. Statements in the articles carry an explicit
status (proved, computed, conjecture; derived, posed, borrowed…).

## Use of generative AI
The programs were written jointly by the author and large language models of the Claude family (Anthropic: Claude
Sonnet 5.5, Claude Opus 5.5, Claude Fable 5.2); the author checked their algorithms and ran the long computations.
See the declaration in each article.

## Licences
Code: MIT (`LICENSE`). Logs, criteria, data and documentation: CC BY 4.0 (`LICENSE-docs.md`).

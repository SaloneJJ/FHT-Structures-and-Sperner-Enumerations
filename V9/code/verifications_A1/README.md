# Verifications for Article I (Fractal Hyper-Trees I. Combinatorial foundations), revised version

Independent checks added during the revision of Article I (November 2026). Python 3 with numpy and scipy.
Each script prints its results; the `.log` files are the outputs of the runs cited in the article.

| script | what it checks (Article I) | run |
|---|---|---|
| `test_axioms_literal.py` | literal test of the axioms (A1)–(A5) over all edge families and all level functions with values in {1,…,v−1}, v ≤ 4: 1, 4, 23 isomorphism classes (Section "Enumeration", validation) | `python3 test_axioms_literal.py` (a few seconds) |
| `twisted_decomp.py` | numerical Wedderburn decomposition of real twisted group algebras R^η[G] (centre, blocks of a generic central element, type R/C/H by the signature of the trace form); checks on R^sgn[S4] = A_F and R^sgn[S5]; decompositions for S4×S3 and S4×C2 and all their sign characters (census of S4, sectors) | `python3 twisted_decomp.py` |
| `valpha54_check.py` | recomputes, for the representative groups of the 54 classes (v,α), v ≤ 6 (generators read from `valpha54_table_rows.tex`, extracted from the table of the article), the type (r,c,h) of R[G] and the decomposition of R^η[G] for every sign character; lists the sectors A_F | `python3 valpha54_check.py` (needs `twisted_decomp.py`) |
| `gl23_check.py` | R^det[GL(2,3)] = A_F ⊕ H ⊕ H ⊕ M2(H); det = sign of the permutation of the 9 points of F_3^2; counts of g² = ±I used in the proof of the spinorial sector | `python3 gl23_check.py` (needs `twisted_decomp.py`) |
| `norms_check.py` | character norms: Λ² of the standard representation of S5 and spin 3/2 restricted to the binary octahedral group 2O are irreducible; degrees 2, 2, 4 of the faithful representations | `python3 norms_check.py` |
| `sep_signature.py` | separation of the 260 FHTs with at most five vertices by numerical signatures (i)–(vii), recomputed from the tables of the appendices (`A1_appendix_tables.tex`); variant with the flat types of the hosts | `python3 sep_signature.py A1_appendix_tables.tex` |

| `fht_moves.py` | independent implementation of the moves of Section 5 (fusion, fission, explosion, aggregation, leaf addition and removal), literal axioms, reading of a realization, canonical form | library |
| `check_moves.py` | checks of Section 8 on the moves, exhaustive for v ≤ 6: classes by closure from the star (1, 1, 4, 23, 231, 16947), validity and inverses, Ω1, edge toggles, components of the quotient graphs, restricted descents (seed 2026) | `python3 check_moves.py --vmax 6 --seed 2026` (about 17 min) |
| `census_s4.py` | census of the FHTs whose automorphism group contains S4 (faithful actions, invariant generation, canonical form, orders of Aut), v ≤ 7 | `python3 census_s4.py 4 5 6 7` (about 2 min) |
| `census_s4_big.py` | same census for v = 8 (actions 0–8) and v = 9 (action 8+1), one JSON checkpoint per action; \|Aut\| by the exact sequence of thm:autrec (validated on v = 6, 7). Result: 10 618 (v=8) and 3 104 (v=9, 8+1). Action 4+1+1+1+1 at v=8: out of memory (8 GB), not computed. Log: `census_v8_v9.log` | `python3 census_s4_big.py 8 0 1 2 3 4 5 6 7 8` ; `python3 census_s4_big.py 9 0` (several hours) |

The programs fht_moves.py, check_moves.py, census_s4.py and census_s4_big.py are exact (integers and sets); they replace the programs of the first version, which were not kept.
The algebra checks are in floating point; they are checks, not exact certifications.
Written with the assistance of Claude (Anthropic), checked and run by the author. Code: MIT; logs: CC BY 4.0.

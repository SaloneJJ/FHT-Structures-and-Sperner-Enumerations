# AF26 — covariance de AF23 et signes des transports (14/11/2026, soir)
* `D3.py` : cordes ou murs (critères 02_theoreme_holonomies/CRITERES_D3_cordes_ou_murs.md) — contrôle échoué, voir EXAMEN_RACINES §6.
* `cov.py`, `cov2.py` : test de covariance de H_or sous des symétries du tesséracte (N = 1) ; localisation : seuls les blocs entre facettes cassent.
* `cov3.py`, `cov4.py` : existe-t-il des signes de transports rendant H_or covariant ? (système linéaire sur ℤ/2, par générateur).
* `D4.py` : sensibilité du spectre aux signes (critères CRITERES_D4_signes_transports.md). N = 2, ≈ 6 min.
Lancement : `cd ~/Desktop/2026FHT/FHT_27/04_code/AF26_covariance && python3 D4.py | tee D4.log` (et de même pour les autres).
* `D5b.py` (Kramers Θ, N = 1 or 2), `D5c.py` (J by facet), `D5d.py` (J and Θ by cell), `D5e.py` (unitary chirality, control) : criteria CRITERES_D5_chiralite_et_Theta.md ;
  usage `python3 D5d.py 1`. Results in EXAMEN_RACINES §6.8.

#!/bin/bash
# Propositions P1 (branche du jumeau a N = 4) et P4 (flot C3' et pas en theta) -- criteres : 05_audit/AUDIT_FABLE_PROPOSITIONS.md (13/11/2026, 20 h 30).
# A lancer SUR LE MAC, venv active (numpy + scipy), depuis n'importe ou :
#     bash ~/Desktop/2026FHT/FHT_26/03_code/CalculCodes/AUDIT_FABLE/lancer_P1_P4.sh            (tout, en sequence)
#     bash .../lancer_P1_P4.sh P1        ou        bash .../lancer_P1_P4.sh P4                    (une seule proposition)
# Duree (ordre de grandeur, d'apres V338_logs.txt) : P1 = 2 x ~2 h ; P4 = ~25 min (N = 3, 256 pas) + ~4 h (N = 4, 128 pas). Memoire : ~6 Go par calcul N = 4.
# Les calculs tournent en sequence (jamais deux N = 4 a la fois : 2 x 6 Go approche la limite de 14 Go).
# Journaux : AUDIT_FABLE/P1_P4_logs/ ; a la fin, evaluer_P1_P4.py applique les criteres et ecrit P1_P4_logs/BILAN.txt.
# Pour laisser tourner la nuit sans que le Mac s'endorme :  caffeinate -i bash .../lancer_P1_P4.sh
set -u
AF="$(cd "$(dirname "$0")" && pwd)"; CC="$(cd "$AF/.." && pwd)"; LOG="$AF/P1_P4_logs"; mkdir -p "$LOG"
cd "$CC"                                   # V336/V338 lisent V162_FHT_R3_jauge.py dans le dossier courant
export PYTHONPATH="$CC/V331_relecture${PYTHONPATH:+:$PYTHONPATH}"   # module rvA_common
python3 -c "import numpy, scipy" || { echo "numpy et scipy sont requis (activer le venv)"; exit 1; }
QUOI="${1:-tout}"
run () {   # nom, commande...
  local nom="$1"; shift
  if [ -s "$LOG/$nom.txt" ] && grep -q 'flot net par couche' "$LOG/$nom.txt"; then echo "deja fait : $nom (journal present)"; return; fi
  echo "$(date '+%d/%m %H:%M') debut $nom : $*" | tee -a "$LOG/deroulement.txt"
  "$@" > "$LOG/$nom.txt" 2>&1
  echo "$(date '+%d/%m %H:%M') fin   $nom (code $?)" | tee -a "$LOG/deroulement.txt"
}
if [ "$QUOI" = tout ] || [ "$QUOI" = P4 ]; then
  run P4_N3_256 python3 V336_FHT_flot.py 3 0 1 1 256 1 unif 1.0
fi
if [ "$QUOI" = tout ] || [ "$QUOI" = P1 ]; then
  run P1_N4_eps0.0 python3 V338_FHT_flot.py 4 0 1 1 64 1 unif 0.0
  run P1_N4_eps0.2 python3 V338_FHT_flot.py 4 0 1 1 64 1 unif 0.2
fi
if [ "$QUOI" = tout ] || [ "$QUOI" = P4 ]; then
  run P4_N4_128 python3 V338_FHT_flot.py 4 0 1 1 128 1 unif 1.0
fi
python3 "$AF/evaluer_P1_P4.py" "$LOG" | tee "$LOG/BILAN.txt"

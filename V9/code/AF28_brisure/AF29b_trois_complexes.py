#!/usr/bin/env python3
"""AF29b -- un seul code : calcule l'un apres l'autre X_LR', X_PS', X_grand (solveur version 2), relit les 5 complexes deja faits
(K26, K26_V25, X_comp, X_chi, X_SM : option --reprise) et imprime le bilan final sur les 8 complexes.
Usage : python3 AF29b_trois_complexes.py | tee AF29b_trois_complexes.log"""
import os, sys
H_ = os.path.dirname(os.path.abspath(__file__))
main = os.path.join(H_, 'AF29b_action_spectrale_descente.py')
sys.argv = [main, '--reprise', 'K26', 'K26_V25', 'X_comp', 'X_chi', 'X_SM', "X_LR'", "X_PS'", 'X_grand']
exec(compile(open(main).read(), main, 'exec'), {'__file__': main, '__name__': '__main__'})

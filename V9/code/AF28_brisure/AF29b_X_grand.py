#!/usr/bin/env python3
"""AF29b -- lanceur pour le seul complexe X_grand (meme code que AF29b_action_spectrale_descente.py ; sorties AF29_X_grand.npz et AF29_X_grand_bilan.json).
Usage : python3 AF29b_X_grand.py | tee AF29b_X_grand.log"""
import os, sys
H_ = os.path.dirname(os.path.abspath(__file__))
sys.argv = [os.path.join(H_, 'AF29b_action_spectrale_descente.py'), "X_grand"]
exec(compile(open(sys.argv[0]).read(), sys.argv[0], 'exec'), {'__file__': sys.argv[0], '__name__': '__main__'})

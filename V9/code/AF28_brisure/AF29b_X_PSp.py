#!/usr/bin/env python3
"""AF29b -- lanceur pour le seul complexe X_PS' (meme code que AF29b_action_spectrale_descente.py ; sorties AF29_X_PSp.npz et AF29_X_PSp_bilan.json).
Usage : python3 AF29b_X_PSp.py | tee AF29b_X_PSp.log"""
import os, sys
H_ = os.path.dirname(os.path.abspath(__file__))
sys.argv = [os.path.join(H_, 'AF29b_action_spectrale_descente.py'), "X_PS'"]
exec(compile(open(sys.argv[0]).read(), sys.argv[0], 'exec'), {'__file__': sys.argv[0], '__name__': '__main__'})

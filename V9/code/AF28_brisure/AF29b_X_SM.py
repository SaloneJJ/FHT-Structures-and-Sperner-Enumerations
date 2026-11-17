#!/usr/bin/env python3
"""AF29b -- lanceur pour le seul complexe X_SM (meme code que AF29b_action_spectrale_descente.py ; sorties AF29_X_SM.npz et AF29_X_SM_bilan.json).
Usage : python3 AF29b_X_SM.py | tee AF29b_X_SM.log"""
import os, sys
H_ = os.path.dirname(os.path.abspath(__file__))
sys.argv = [os.path.join(H_, 'AF29b_action_spectrale_descente.py'), "X_SM"]
exec(compile(open(sys.argv[0]).read(), sys.argv[0], 'exec'), {'__file__': sys.argv[0], '__name__': '__main__'})

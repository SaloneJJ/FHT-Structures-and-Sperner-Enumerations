"""Outils d analyse des balayages : integrande M1 par couche."""
import numpy as np
def integrand(z, eps=None):
    eps = z['eps'] if eps is None else eps
    E, S, W = z['E'], z['S'], z['W']
    m = (np.abs(E) < eps) & ~np.isnan(E)
    return np.einsum('tnk,tk->tn', W, np.where(m, S, 0.0)) / (2 * eps)   # f_n(theta)

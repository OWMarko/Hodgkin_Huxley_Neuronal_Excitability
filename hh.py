import math

import numpy as np
from scipy.integrate import solve_ivp

# Paramètres classiques utilisés par le mécanisme hh de NEURON
CM = 1.0 # µF/cm**2
G_NA, G_K, G_L = 120.0, 36.0, 0.3 # mS/cm**2
E_NA, E_K, E_L = 50.0, -77.0, -54.3  # mV


def phi(x):
    if abs(x) < 1e-7:
        return 1.0 - x / 2.0 + x * x / 12.0
    
    return x / np.expm1(x) #prolongée par 1 en x=0


def taux(v):

    am = phi(-(v + 40.0) /10.0)
    bm = 4.0 * np.exp(-(v+65.0) /18.0)
    ah = 0.07 * np.exp(-(v+65.0) /20.0)
    bh = 1.0 / (1.0 + np.exp(-(v+35.0) /10.0))
    an = 0.1 * phi(-(v+55.0)/10.0)
    bn = 0.125 * np.exp(-(v+65.0) / 80.0) #1/ms

    return am, bm, ah, bh, an, bn


def etat_initial():

    am, bm, ah, bh, an, bn = taux(-65.0) #Potentiel nominal -65 mV

    return np.array([-65.0, am / (am+bm), ah / (ah+bh), an / (an+bn)])


def derivee(etat, courant, g_na=G_NA, g_k=G_K):

    v, m, h, n = etat
    am, bm, ah, bh, an, bn = taux(v)
    i_na = g_na * m**3 * h * (v-E_NA)
    i_k = g_k * n**4 * (v-E_K)
    i_l = G_L * (v-E_L)
    
    return np.array([
        (courant - i_na - i_k - i_l) / CM,
        am * (1.0-m) - bm*m,
        ah * (1.0-h) - bh*h,
        an * (1.0-n) - bn*n,])


def simuler(amplitude, debut, fin, duree, g_na=G_NA, g_k=G_K, pas=0.05):

    if (not np.isfinite((debut, fin, duree, pas)).all()
            or not (0 < debut < fin < duree and pas > 0)):
        raise ValueError("Temps finis et ordonnés, pas fini strictement positif requis.")

    etat = etat_initial()
    morceaux_t, morceaux_y = [], []

    for a, b, courant in ((0.0, debut, 0.0), (debut, fin, amplitude), (fin, duree, 0.0)):

        n = max(1, math.ceil((b-a) / pas))
        temps = np.linspace(a, b, n + 1) # uniformes sur chaque intervalle

        sol = solve_ivp(lambda t, y: derivee(y, courant, g_na, g_k),
                        (a, b), etat, method="DOP853", t_eval=temps,
                        rtol=1e-9, atol=(1e-11, 1e-12, 1e-12, 1e-12),
                        max_step=0.1) # pas internes sont adaptatifs mais limités à 0.1ms
        
        if not sol.success or not np.isfinite(sol.y).all():
            raise RuntimeError("Échec de l'intégration numérique")
        
        debut_index = 0 if not morceaux_t else 1 # La première valeur des segments suivants est déjà présente
        morceaux_t.append(sol.t[debut_index:])
        morceaux_y.append(sol.y[:, debut_index:])
        etat = sol.y[:, -1]

    return np.concatenate(morceaux_t), np.concatenate(morceaux_y, axis=1)


def compter_pics(t, v, debut=100.0, fin=250.0): #Compte les passages ascendants à 0 mV dans [debut, fin)

    indices = np.flatnonzero((v[:-1] < 0) & (v[1:] >= 0))

    temps_pics = t[indices] - v[indices] * (t[indices+1] - t[indices]) / (v[indices+1] - v[indices]) # t doit être croissant et v fini

    return int(np.count_nonzero((temps_pics >= debut) & (temps_pics < fin)))

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from hh import compter_pics, simuler


def experience_impulsion(dossier):

    resultats = {i: simuler(i, 20, 21, 60) for i in (0, 1, 10)}
    fig, axes = plt.subplots(3, 1, figsize=(8, 8), layout="constrained")

    for i, (t, y) in resultats.items():
        axes[0].plot(t, y[0], label=f"I = {i} µA/cm²")

    axes[0].set(xlim=(18, 50), ylabel="V (mV)", title="Réponse à une impulsion de courant")
    axes[0].legend()

    for i in (0, 1):
        t, y = resultats[i]
        axes[1].plot(t, y[0], label=f"I = {i} µA/cm²")

    axes[1].set(xlim=(18, 40), ylabel="V (mV)", title="Agrandissement : réponse faible")
    axes[1].legend()

    t, y = resultats[10]
    for k, nom in enumerate(("m · activation Na", "h · inactivation Na", "n · activation K"), 1):
        axes[2].plot(t, y[k], label=nom)

    axes[2].set(xlim=(18, 40), ylim=(0, 1), xlabel="Temps (ms)", ylabel="Fraction (sans dimension)", title="Variables des canaux")
    axes[2].legend()

    for ax in axes:
        ax.axvspan(20, 21, color="grey", alpha=0.15)

    fig.savefig(dossier / "impulsion.png", dpi=150)
    plt.close(fig)
    print(f"Impulsion forte : Vmax ≈ {resultats[10][1][0].max():.3f} mV")


def experience_conductances(dossier):

    cas = {
        "Référence": (120.0, 36.0),
        "Na −20 %": (96.0, 36.0),
        "K +20 %": (120.0, 43.2),}
    
    fig, axes = plt.subplots(2, 1, figsize=(8, 7), layout="constrained")
    noms, comptes = [], []

    for nom, (g_na, g_k) in cas.items():
        t, y = simuler(14, 50, 250, 300, g_na, g_k)
        n_pics = compter_pics(t, y[0])
        noms.append(nom)
        comptes.append(n_pics)
        axes[0].plot(t, y[0], label=nom)
        print(f"{nom} : {n_pics} événements sur [100, 250) ms")

    axes[0].set(xlim=(100, 145), ylim=(-85, 50), ylabel="V (mV)", xlabel="Temps (ms)", title="Courant soutenu : I = 14 µA/cm²")
    axes[0].legend()

    axes[1].bar(noms, comptes)
    axes[1].set(ylabel="Nombre de passages à 0 mV",title="Comptage sur [100, 250) ms")
    axes[1].set_ylim(0, max(comptes) + 3)
    
    fig.savefig(dossier / "conductances.png", dpi=150)
    plt.close(fig)


def main():
    dossier = Path(__file__).resolve().parent / "figures"
    dossier.mkdir(exist_ok=True)
    experience_impulsion(dossier)
    experience_conductances(dossier)


if __name__ == "__main__":
    main()

import numpy as np
import pytest

from hh import compter_pics, derivee, etat_initial, phi, simuler, taux


def test_limites_des_taux():
    assert taux(-40)[0] == pytest.approx(1.0)
    assert taux(-55)[4] == pytest.approx(0.1)
    assert phi(0) == 1.0


def test_portes_stationnaires_a_moins_65():
    x = etat_initial()
    dx = derivee(x, 0)
    assert np.allclose(dx[1:], 0, atol=1e-14)
    assert dx[0] == pytest.approx(0.03032370918249644, abs=1e-11)


def test_taux_positifs():
    for v in np.linspace(-100, 50, 31):
        assert np.all(np.array(taux(v)) > 0)


def test_reponse_impulsion():
    t0, y0 = simuler(0, 20, 21, 60)
    t1, y1 = simuler(1, 20, 21, 60)
    t10, y10 = simuler(10, 20, 21, 60)
    assert np.array_equal(t0, t1) and np.array_equal(t1, t10)
    assert y0[0].max() < -60
    assert y1[0].max() < 0
    assert y10[0].max() > 30
    assert np.all((y10[1:] >= -1e-9) & (y10[1:] <= 1 + 1e-9))


def test_temps_croissants_et_commutations():
    t, y = simuler(10, 20, 21, 60)
    assert np.all(np.diff(t) > 0)
    assert np.sum(t == 20) == 1 and np.sum(t == 21) == 1
    assert y.shape == (4, t.size)


def test_conductances_et_comptages():
    cas = [(120, 36, 12), (96, 36, 9), (120, 43.2, 10)]
    for g_na, g_k, attendu in cas:
        t, y = simuler(14, 50, 250, 300, g_na, g_k)
        assert compter_pics(t, y[0]) == attendu


def test_detection_sur_signal_simple():
    t = np.array([0, 1, 2, 3, 4], dtype=float)
    v = np.array([-30, 10, -30, 10, -30], dtype=float)
    assert compter_pics(t, v, debut=0, fin=4) == 2
    assert compter_pics(t, v, debut=0, fin=1) == 1


def test_comptage_stable_quand_on_affine_lechantillonnage():
    t1, y1 = simuler(14, 50, 250, 300, pas=0.1)
    t2, y2 = simuler(14, 50, 250, 300, pas=0.05)
    assert compter_pics(t1, y1[0]) == compter_pics(t2, y2[0]) == 12


def test_pas_maximal_sur_segments_non_divisibles():
    pas = 0.3
    t, y = simuler(0, 0.2, 0.7, 1.05, pas=pas)
    marge = 8 * np.finfo(float).eps * max(1.0, t[-1])
    assert np.all(np.diff(t) > 0)
    assert np.max(np.diff(t)) <= pas + marge
    for borne in (0.0, 0.2, 0.7, 1.05):
        assert np.count_nonzero(t == borne) == 1
    assert np.array_equal(y[:, 0], etat_initial())
    assert np.isfinite(y).all()


def test_temps_et_pas_non_finis_rejetes():
    for nom in ("debut", "fin", "duree", "pas"):
        for valeur in (np.nan, np.inf, -np.inf):
            arguments = dict(debut=0.2, fin=0.7, duree=1.05, pas=0.3)
            arguments[nom] = valeur
            with pytest.raises(ValueError):
                simuler(0, **arguments)

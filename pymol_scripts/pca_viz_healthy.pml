load /home/user/Documents/Stijn_Jasper/CasusD/fold_healthy_egfr/fold_healthy_egfr_model_0-MD.gro, Healthy
load_traj /home/user/Documents/Stijn_Jasper/CasusD/fold_healthy_egfr/whol.xtc, Healthy

remove solvent
remove name na+cl

intra_fit Healthy

remove solvent
remove name na+cl

import numpy as np

X = np.array([cmd.get_model("Healthy", state=idx + 1).get_coord_list() for idx in range(0, cmd.count_states("Healthy"))])
mean = X.mean(axis=0)
X = (X - mean).reshape(len(X), -1)

vals, vecs = np.linalg.eigh(X @ (X.T / len(X)))
loadings = X.T @ vecs[:, ::-1]
loadings /= (loadings ** 2).sum(axis=0, keepdims=True) ** 0.5
scores = X @ loadings

L1 = loadings[:, 0].reshape((-1, 3))
S1 = scores[:, 0]

xmin, xmax = S1.min() * L1 + mean, S1.max() * L1 + mean


# create pca movements
cmd.load_cgo([u for (xs, ys, zs), (xe, ye, ze) in zip(xmin, xmax) for u in [9.0, xs, ys, zs, xe, ye, ze, 0.1, 0, 0, 1, 1, 0, 0]], 'pc1')

show stick, Healthy

color yellow, Healthy
color white, resi 164

set opaque_background, on

set_view (\
    -0.680888057,   -0.639156938,    0.357364565,\
     0.714106381,   -0.471343905,    0.517577291,\
    -0.162363753,    0.607610345,    0.777351856,\
     0.000000000,    0.000000000,  -47.846450806,\
    81.214691162,   73.242485046,   36.959102631,\
    40.350788116,   55.342113495,  -20.000000000 )

png healthy_movement_site.png


# create mean model
M = cmd.get_model("Healthy")
for idx, mx in enumerate(mean): M.atom[idx].coord = list(mx)
cmd.load_model(M, 'mean')

hide everything, Healthy
hide everything, pc1

set_view (\
    -0.680888057,   -0.639156938,    0.357364565,\
     0.714106381,   -0.471343905,    0.517577291,\
    -0.162363753,    0.607610345,    0.777351856,\
     0.000000000,    0.000000000,  -26.282503128,\
    81.214683533,   73.242477417,   36.959106445,\
    -0.313288718,   52.878246307,  -20.000000000 )

png mean_sticks_healthy.png


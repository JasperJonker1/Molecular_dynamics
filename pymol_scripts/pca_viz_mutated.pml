load /home/user/Documents/Stijn_Jasper/CasusD/fold_tumor_leu858arg/fold_tumor_leu858arg_model_0-MD.gro, Tumor
load_traj /home/user/Documents/Stijn_Jasper/CasusD/fold_tumor_leu858arg/disease.xtc, Tumor

remove solvent
remove name na+cl

intra_fit Tumor
import numpy as np

X = np.array([cmd.get_model("Tumor", state=idx + 1).get_coord_list() for idx in range(0, cmd.count_states("Tumor"))])
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

show stick, Tumor

color yellow, Tumor
color white, resi 164

set opaque_background, on


set_view (\
     0.362789661,    0.103574730,    0.926011443,\
     0.218060195,   -0.975647151,    0.023696518,\
     0.905907273,    0.193328321,   -0.376544356,\
     0.000000000,    0.000000000,  -42.855102539,\
    63.013862610,   69.613998413,   21.331027985,\
    31.426670074,   54.283546448,  -20.000000000 )

png mutation_movement_site.png


hide everything, Tumor
hide everything, pc1

# create mean model
M = cmd.get_model("Tumor")
for idx, mx in enumerate(mean): M.atom[idx].coord = list(mx)
cmd.load_model(M, 'mean')

set_view (\
    -0.817288637,    0.347239733,    0.459667951,\
    -0.225106403,    0.542005181,   -0.809663296,\
    -0.530284345,   -0.765206158,   -0.364814341,\
     0.000000000,    0.000000000,  -25.917985916,\
    63.013847351,   69.613990784,   21.331027985,\
    -1.366039991,   53.202045441,  -20.000000000 )

png mean_sticks_mutated.png




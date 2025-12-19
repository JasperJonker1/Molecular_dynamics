# NOTE replace "/home/user/Documents/" with "/students/2025-2026/MvK" in final script to make it work in BIN

# load files with trajectory
load /home/user/Documents/Stijn_Jasper/CasusD/fold_healthy_egfr/fold_healthy_egfr_model_0-MD.gro, Healthy
load /home/user/Documents/Stijn_Jasper/CasusD/fold_tumor_leu858arg/fold_tumor_leu858arg_model_0-MD.gro, Tumor

load_traj /home/user/Documents/Stijn_Jasper/CasusD/fold_healthy_egfr/whol.xtc, Healthy
load_traj /home/user/Documents/Stijn_Jasper/CasusD/fold_tumor_leu858arg/disease.xtc, Tumor

# remove water and salt
remove solvent
select salt, name na+cl
remove salt

color blue, Healthy
color red, Tumor

# aligning
intra_fit Healthy
intra_fit Tumor

align Tumor, Healthy

# select part around the mutation (mutation = L -> R)
select normal, pepseq VKITDFGLAK
select mutation, pepseq VKITDFGRAK

color white, mutation
color orange, normal

show sticks, mutation
show sticks, normal

hide everything, hydro
util.cbaw not elem C

deselect

hide cartoon

set_view (\
     0.416936725,   -0.845008016,    0.334848881,\
    -0.081855856,    0.331992716,    0.939721465,\
    -0.905242085,   -0.419215202,    0.069252782,\
    -0.000177117,   -0.000056691,  -91.481369019,\
    74.553916931,   68.491058350,   31.426887512,\
    30.065603256,  153.054321289,  -20.000000000 )

png mutation_site.png

show cartoon, Healthy
show cartoon, Tumor

show sticks, mutation
show sticks, normal

deselect

mset 1

set_view (\
     0.998427391,   -0.031380139,    0.046456799,\
     0.025084889,    0.991147280,    0.130379260,\
    -0.050136916,   -0.129008234,    0.990375280,\
     0.000026464,    0.000074082, -290.625762939,\
    65.712379456,   65.119903564,   29.442941666,\
   229.132431030,  352.121276855,  -20.000000000 )

png tumor_healthy.png

mset 1 - 2001

movie.produce simulation_EGFR.mpg


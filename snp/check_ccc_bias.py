import numpy as np

def ccc(i, j, q):
    Fi = 1 - (fi / q)
    Fj = 1 - (fj / q)
    CCC_ij = 4.5 * Rij * Fi * Fj
    return CCC_ij

def duo(Dij, fi, fj):
    DUOij = 4 * Dij * (1 - (fi / 1.5)) * (1 - (fj / 1.5))
    return DUOij

num_occ_a = 160000
num_occ_b = 100
num_occ = 160000
num_vec = 2000000
Dij_norm = (num_vec - num_occ) / num_vec
fi = num_occ_a / num_vec
fj = num_occ_b / num_vec

result_norm = duo(Dij_norm, fi, fj)
print(result_norm)

result_one = duo(Dij_norm, fi, fj)
print(result_one)

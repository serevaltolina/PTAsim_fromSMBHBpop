from __future__ import print_function, division         

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy import stats
import scipy
import h5py

'''
    This scripts computes the average of a set of posterior distributions using KDE objects.
    N.B.: in two universes of our dataset (n. 1029 and 1055), the obtained posterior is disjoint from the others.
    This is due to peculiar feature of the GWB realisation that cause the recovery to be strongly biased 
    far away from the nominal values of the GWB parameters. Including those two realisations would reduce 
    the average to zero; thus, we decided to exclude them.


'''


home_dir = './SMBHB_datasets/univ_'

# KDE method
#----------------------------------------------------------------------------

points = 80000 # reweighted samples

u_list = np.arange(1001, 1101, 1)
A_test = np.arange(-15.5,-13.5,0.005)
gamma_test = np.arange(0,7,0.0175)
N_univ = 98 # int(len(u_list)) - disjoint posteriors
cont = 0

# ked priors
points_gamma_p = []
points_A_p = []
for i in range(points):
    points_gamma_p.append(np.random.random()*7)
    points_A_p.append(-13.5-np.random.random()*2)
points_gamma_p = np.array(points_gamma_p)
points_A_p = np.array(points_A_p)

kernel_gamma_p = stats.gaussian_kde(points_gamma_p)
kernel_A_p = stats.gaussian_kde(points_A_p)

prior_gamma = kernel_gamma_p(gamma_test) 
prior_A = kernel_A_p(A_test)

# combining posteriors for A_GWB and gamma_GWB
for univ in u_list:

    # !!!
    # some posteriors are disjoint ans need to be excluded
    # !!!

    data_dir = home_dir + str(univ) + '/'
    chain =  np.loadtxt(data_dir+'chain_1.txt')
    A_chain = chain[-points:,-5]
    gamma_chain = chain[-points:,-6]

    w_dir = data_dir + 'weights.txt'
    w = []
    with open(w_dir,'r') as f:
        lines = f.readlines()
        lines = lines[1:]
        for j in lines:
            w.append(float(j.split()[2]))
        f.close()
    w = np.array(w)

    kernel_A = stats.gaussian_kde(A_chain, weights=w)
    kernel_gamma = stats.gaussian_kde(gamma_chain, weights=w)

    if cont == 0:
        data_A = kernel_A(A_test)
        data_gamma = kernel_gamma(gamma_test)
        cont += 1
    else:
        data_temp_A = kernel_A(A_test) 
        data_temp_gamma = kernel_gamma(gamma_test) 
        data_A = data_A + data_temp_A
        data_gamma = data_gamma * data_temp_gamma

    print(univ,' done')

data_A = data_A - (N_univ-1)*prior_A 
data_gamma = data_gamma / prior_gamma**(N_univ-1) 


# plot
#----------------------------------------------------------------------------


fig = plt.figure()
plt.plot(A_test,data_A, label='realistic data set', color='darkorange')
plt.axvline(x=np.log10(2.4*10**(-15)), color='black', label='$log_{10}(A_{GWB})$ = 2.4e-15')
plt.xlabel('$log_{10}(A_{GWB})$', fontsize=13)
plt.title('Combined posterior - $log_{10}(A_{GWB})$', fontsize=15)
plt.legend(loc='best')
plt.ylim(0,45)
plt.xlim(-15.3,-14)
fig.savefig('./combined_postA.png')

fig = plt.figure()
plt.plot(gamma_test,data_gamma, label='realistic data set', color='darkorange')
plt.axvline(x=4.33, color='black', label='$\gamma_{GWB}$ = 13/3')
plt.xlabel('$\gamma_{GWB}$', fontsize=13)
plt.title('Combined posterior - $\gamma_{GWB}$', fontsize=15)
plt.legend(loc='best')
plt.xlim(3.75,5.5)
fig.savefig('./combined_postgamma.png')


h5f = h5py.File(f'./combined_post.h5', 'w')
h5f.create_dataset('/A_realinj', data=np.asarray(data_A))
h5f.create_dataset('/gamma_realinj', data=np.asarray(data_gamma))
h5f.close()


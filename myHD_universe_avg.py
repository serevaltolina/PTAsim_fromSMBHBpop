from __future__ import print_function, division

import numpy as np
import matplotlib.pyplot as plt
import h5py

from optimal_statistic_covariances import OptimalStatistic, get_HD_curve

rhos = np.zeros([10])
sigmas = np.zeros([10])
n_univ = 100

for i in range (n_univ):

    data_dir = './SMBHB_datasets/univ_' + str(1001+i) + '/'

    # before this script run OS on the single datasets and save results as indicated here
    f = h5py.File(data_dir+'OS_10bin.h5', 'r')
    avg_angles = np.array(f['avg_angles'])
    binned_os = np.array(f['binned_os'])
    bin_cov_matrix = np.array(f['bin_cov_matrix'])
    interval = np.array(f['interval']) # theoretical uncertainty

    errors = np.sqrt(np.diag(bin_cov_matrix))

    rhos += binned_os
    sigmas += errors**2

rhos = rhos/n_univ
sigmas = np.sqrt(sigmas)/n_univ

print('FINAL rhos and sigma = ', rhos, sigmas)

A_gwb = 2.4e-15
A2 = A_gwb**2

fig = plt.figure()
plt.errorbar(avg_angles, rhos/A2, sigmas/A2, capsize=10, color='deeppink', fmt='o')
angles = np.linspace(0, np.pi, num=100)
plt.plot(angles, get_HD_curve(angles), c='k', ls='--', label='HD function') 
plt.scatter(avg_angles, get_HD_curve(avg_angles)+interval/A2, color='grey', marker='x', label='theoretical variance')  
plt.scatter(avg_angles, get_HD_curve(avg_angles)-interval/A2, color='grey', marker='x') 
plot_dir = './HDfromOS_alluniv_avg.png'
plt.xlabel('$\gamma_{ab}$')
plt.ylabel('$\Gamma(\gamma_{ab})$')
plt.legend(loc='best')
fig.savefig(plot_dir)


# saving results
h5f = h5py.File(f'./OSfinal_10bin_univavg.h5', 'w')
h5f.create_dataset('/avg_angles', data=np.asarray(avg_angles) )
h5f.create_dataset('/rhos', data=np.asarray(rhos) )
h5f.create_dataset('/sigmas', data=np.asarray(sigmas) )
h5f.close()

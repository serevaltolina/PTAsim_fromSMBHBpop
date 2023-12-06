import sys, math

import numpy as np
import matplotlib.pyplot as plt

import os
import shutil
import json

import glob
import scipy.linalg as sl
from scipy import interpolate
from scipy.integrate import quad
import scipy


#--------------------------------------------------------------------------------------------------
#           PP PLOT FUNCTION
#--------------------------------------------------------------------------------------------------

def PP_plot(u_list, N_L_evaluations, PP_output):

    '''
        This function builds the PP plot from a set of mcmc posteriors saved in chain_dir and
        re-weighted accordingly to weights.txt
        
    '''

    total_count_matrix = []
    n_univ = len(u_list)

    # PP plot P_x values
    P_x = np.arange(0.025,1.0,0.025)

    # reading A_gwb samples and weights for each universe
    #--------------------------------------------------------------
    for u in u_list:

        points = 80000 # number of samples of the chain for which we computed the weights
        var = 1
        data_dir = home_dir + str(u) + '/'
        chain_dir = data_dir + 'chain_1.txt'
        chain = np.loadtxt(chain_dir)
        A_chain = chain[-points:,-5] 

        w = []
        w = np.array(w)
                
        w_dir = data_dir + 'weights.txt'
        with open(w_dir,'r') as f:
            lines = f.readlines()
            for i in range(len(lines)):
                if i == 0:
                    continue;
                else:
                    w.append(float(lines[i].split()[2]))
            f.close()
        w = np.array(w)

        # posterior fit
        #--------------------------------------------------------------

        n, bins, _ = plt.hist(A_chain,bins=80, density = True, weights=w, alpha = 0.8)
        step = float(((np.max(A_chain)-np.min(A_chain))/len(bins))/2)

        x_n = np.concatenate((bins-step,[bins[-1]+step]))
        n = np.concatenate(([0.],n,[0.]))

        f_post = interpolate.interp1d(x_n,n)
        x_test = np.arange(x_n[0],x_n[-1]+0.001,0.1)

        # CDF calculation
        #--------------------------------------------------------------

        N_L_evaluations = int(N_L_evaluations)
        Amin = np.min(A_chain)
        Amax = np.max(A_chain)
        logA_inj = np.log10(2.4*10**(-15))
        logA_test = np.linspace(Amin,Amax,N_L_evaluations)
        y_cdf = np.array([tup[0] for tup in [quad(f_post, a, b) for a, b in [(a, b) for a, b in zip(logA_test, logA_test[1:len(logA_test)])]]] + [0]).cumsum()
        y_cdf[-1] = 1.

        # inverse CDF interpolation
        f_cdf = interpolate.interp1d(y_cdf,logA_test)
    
        # PP plot y data counting
        #--------------------------------------------------------------
        cont = []
        for x in P_x:
            y = f_cdf(x)
            if y >= logA_inj:
                cont.append(1)
            else:
                cont.append(0)

        cont = np.array(cont)
        total_count_matrix.append(cont)

        print('universe '+str(u)+' done')

    total_count_matrix = np.array(total_count_matrix)
    print(total_count_matrix.shape)
    
    # PP plot y-axis data
    #--------------------------------------------------------------
    PP_y = []
    for i in range(len(P_x)):
        PP_y.append((np.sum(total_count_matrix[:,i]))/n_univ)
    PP_y = np.array(PP_y)

    P_x = np.concatenate(([0],P_x,[1.]))
    PP_y = np.concatenate(([0],PP_y,[1.]))


    # PP plot
    #--------------------------------------------------------------
    fig = plt.figure()
    plt.plot(P_x,PP_y,color='darkorange')
    x_theory = np.arange(0,1.1,0.1)
    y_theory = x_theory
    y_theory_std = np.sqrt(x_theory*(1-x_theory)/100)
    plt.plot(x_theory,y_theory,ls='--',color='royalblue',label='theoretical expect.')
    plt.fill_between(x_theory,y_theory+y_theory_std,y_theory-y_theory_std, color='lightblue', alpha=1,  label ='1 $\sigma$')
    plt.fill_between(x_theory,y_theory+3*y_theory_std,y_theory-3*y_theory_std, color='lightblue', alpha=0.4, label ='3 $\sigma$')
    plt.xlabel('$P_{x}$', fontsize=14)
    plt.ylabel('P($Agwb_{inj}$ is found within $P_{x}$)', fontsize=14)
    plt.grid()
    plt.legend(loc='best')
    
    fig.savefig(PP_output+'.png')

    # saving PP plot points
    #--------------------------------------------------------------
    
    with open(PP_output+'.txt', 'w') as ff:
        for i in range (len(P_x)):
            ff.writelines(str(P_x[i])+'\t\t\t\t'+str(PP_y[i])+'\n')
        ff.close()
    

        
#--------------------------------------------------------------------------------------------------
#          Main script
#--------------------------------------------------------------------------------------------------

home_dir = './SMBHB_datasets/univ_'
universe_list = np.arange(1001,1101,1)
PP_output = './PP_plot'
PP_plot(universe_list, 1000, PP_output)

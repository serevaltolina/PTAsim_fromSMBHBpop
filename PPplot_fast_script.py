import sys, os

import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as sl

from scipy import interpolate
from scipy.integrate import quad


#-----------------------------------------------------------------------------
#       PP PLOT FUNCTION
#-----------------------------------------------------------------------------

def PP_plot_L(u_list, logA_test):

    '''
        This function builds the PP plot from a set of likelihood evaluations saved in L_evaluations_9fbin.txt,
        computed from a pta model with A_GWB as the only free parameter.
        
    '''

    total_count_matrix = []
    n_univ = len(u_list)

    # PP plot P_x values
    P_x = np.arange(0.025,1.0,0.025)

    for u in u_list:

        # reading logL output data
        data_dir = home_dir + str(u) + 'L_evaluations_9fbin.txt'
        L_list = []
        with open(data_dir, 'r') as f:
            lines = f.readlines()
            for i in lines[1:]:
                L_list.append(float(i.split()[1]))
 
            L_list = np.array(L_list)
            max_L = np.max(L_list)
            plot_L_list = np.exp(L_list-max_L)

            f.close()

        # nominal A_GWB of the realisations
        logA_inj = np.log10(2.4*10**(-15))
        
        # likelihood plot interpolation function
        f_plot_L = interpolate.interp1d(logA_test,plot_L_list)
        norm_factor, err = quad(f_plot_L, logA_min, logA_max)
 
        # likelihood pplot normalization
        f_plot_L_n = interpolate.interp1d(logA_test,plot_L_list/norm_factor)

        # CDF calculation
        y_cdf = np.array([tup[0] for tup in [quad(f_plot_L_n, a, b) for a, b in [(a, b) for a, b in zip(logA_test, logA_test[1:len(logA_test)])]]] + [0]).cumsum()

        # inverse CDF interpolation
        f_cdf = interpolate.interp1d(y_cdf,logA_test)

        # PP plot y data counting
        cont = []
        for x in P_x:
            y = f_cdf(x)
            if y >= logA_inj:
                cont.append(1)
            else:
                cont.append(0)
        cont = np.array(cont)
        total_count_matrix.append(cont)
        print(str(u)+'_done')

    total_count_matrix = np.array(total_count_matrix)

    # PP plot y-axis data
    PP_y = []
    for i in range(len(P_x)):
        PP_y.append((np.sum(total_count_matrix[:,i]))/n_univ)
    PP_y = np.array(PP_y)

    P_x = np.concatenate(([0],P_x,[1.]))
    PP_y = np.concatenate(([0],PP_y,[1.]))

    # PP plot
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
    plot_dir = './PPplot_Agwb_fast.png'
    fig.savefig(plot_dir)

    with open('./PPplot_Agwb_fast.png.txt', 'w') as f:
        for i in range(len(P_x)):
            f.writelines(str(P_x[i])+'\t\t\t\t'+str(PP_y[i])+'\n')
        f.close()


#-----------------------------------------------------------------------------

# main script

logA_max = -13.5
logA_min = -15.5
N_L_evaluations = 1000
logA_test = np.linspace(logA_min,logA_max,N_L_evaluations)

home_dir = './SMBHB_datasets/univ_'
u_list = np.arange(1001,1101,1)
PP_plot_L(u_list, logA_test)
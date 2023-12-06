from __future__ import print_function, division         
import sys, json

import numpy as np
import matplotlib.pyplot as plt

from enterprise.pulsar import Pulsar
import enterprise.signals.parameter as parameter
from enterprise.signals import utils
from enterprise.signals import signal_base
from enterprise.signals import white_signals
from enterprise.signals import gp_signals

'''
    Likelihood evaluations sapling the logA_GWB prior, when logA_GWB is the only free parameter.
'''

univ_code = str(sys.argv[1]) # range: 1001 - 1100

home_dir = './SMBHB_datasets/univ_'+univ_code+'/'
noise_dict_dir = './complete_noise_dict.json'

def pulsar_model(psr, model_base):
    
    model = model_base
    
    # red noise
    if psr.name+'_red_noise_gamma' in noise_dict:
        log10_A = parameter.Uniform(-18,-11)
        gamma = parameter.Uniform(0,7)
        pl = utils.powerlaw(log10_A=log10_A, gamma=gamma)
        rn = gp_signals.FourierBasisGP(spectrum=pl, components=30)
        model += rn
        
    # DM variations
    if psr.name+'_dm_gp_gamma' in noise_dict:
        log10_A_dm = parameter.Uniform(-18,-8)
        gamma_dm = parameter.Uniform(0, 7)
        pl_dm = utils.powerlaw(log10_A=log10_A_dm, gamma=gamma_dm)
        dm_basis = utils.createfourierdesignmatrix_dm(nmodes=100)
        dm = gp_signals.BasisGP(priorFunction=pl_dm, basisFunction=dm_basis, name='dm')
        model += dm
        
    return model(psr)

# reading the names of the pulsars in the array
pulsar_list_dir = './pulsar_list.dat'
pulsar_names = np.loadtxt(pulsar_list_dir,dtype=str)



#**********************************************************************************
#       RECOVERY SCRIPT
#**********************************************************************************


# defining psr object
psrs = []
for name in pulsar_names:
    par_dir = home_dir + 'parfiles/' + name + '.par'
    tim_dir = home_dir + 'tims/' + name + '.tim'
    psr = Pulsar(par_dir, tim_dir)
    psrs.append(psr)


# noise dictionary reading
with open(noise_dict_dir, 'r') as fp:
    noise_dict = json.load(fp)
    fp.close()

# find the maximum time span to set GW frequency sampling
tmin = [p.stoas.min() for p in psrs]
tmax = [p.stoas.max() for p in psrs]
Tspan = np.max(tmax) - np.min(tmin)

# fixing EFAC and EQUAD
efac = parameter.Constant(1.0) 
equad = parameter.Constant(-6)

# white noise
ef = white_signals.MeasurementNoise(efac=efac)
eq = white_signals.TNEquadNoise(log10_tnequad = equad)

# timing model
tm = gp_signals.TimingModel() 

# gwb
orf = utils.hd_orf()
log10_A_gwb = parameter.Uniform(-15.5,-13.5)
gamma_gwb = parameter.Constant()
pl_gwb = utils.powerlaw(log10_A=log10_A_gwb, gamma=gamma_gwb)
gw = gp_signals.FourierBasisCommonGP(pl_gwb, orf, components=9, name='gw', Tspan=Tspan)

# model base
model_base = tm + ef + eq + gw

# PTA object definition
pta = signal_base.PTA([pulsar_model(psr,model_base) for psr in psrs])

# pta fixed parameters definition
pta.set_default_params(noise_dict) 


# LIKELIHOOD EVALUATIONS
#-------------------------------------------------------------

N_L_evaluations = 1000
A_min = -15.5
A_max = -13.5
xs_values = np.linspace(A_min,A_max,N_L_evaluations)
L_output = home_dir + 'L_evaluations_9fbin.txt'

with open(L_output, 'w') as f:
    f.writelines('log10_A_gwb_prior \t\t lnL\n')
    for i in range (N_L_evaluations):
        xs = {par.name: xs_values[i] for par in pta.params}
        L = pta.get_lnlikelihood(xs)
        f.writelines(str(xs['gw_log10_A'])+'\t\t'+str(L)+'\n')


# lnL plot
#------------------------------------------------------------

with open(L_output,'r') as f:
    lines = f.readlines()
    lnL = []
    logA = []
    cont = 0
    for i in lines:
        if cont == 0:
            cont = 1
            continue;
        else:
            logA.append(float(i.split()[0]))
            lnL.append(float(i.split()[1]))
    lnL = np.array(lnL)
    logA = np.array(logA)
    max_L = np.max(lnL)

    L = np.exp(lnL-max_L)
    f.close()

nominal_A = np.log10(2.4*10**(-15))
fig = plt.figure()
plt.plot(logA,L)
plt.axvline(nominal_A, color='orange')
plt.title('Injected log10(A_GWB) = 2.4e-15')
plt.xlabel('log10_A_GWB')
plt.ylabel('exp(lnL-max(lnL))')
plt.grid()
fig.savefig(home_dir+'L_plot9.png')

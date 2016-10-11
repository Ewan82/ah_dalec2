import data_class as dc
import mod_class as mc
import numpy as np
import random as rand
import pickle

# Background cov mat
b_cor = pickle.load(open('b_edc_cor.p', 'r'))
b_std = np.sqrt(np.diag(pickle.load(open('b_edc.p', 'r'))))
b_std[10] = 0.25*b_std[10]
b_std[1] = 0.25*b_std[1]
b_std[0:17] = 0.5*b_std[0:17]
D = np.zeros_like(b_cor)
np.fill_diagonal(D, b_std)
b = 0.6 * np.dot(np.dot(D, b_cor), D)


def var_ens(east_west, size_ens=10):
    edc_ens = pickle.load(open('xa_param_ens.p', 'r'))
    param_ens = rand.sample(edc_ens, size_ens)
    output = [run_4dvar_desroziers(pvals, east_west) for pvals in param_ens]
    f = open('misc/var_ens_out2', 'w')
    pickle.dump(output, f)
    f.close()
    return output


def perturb_obs(ob_arr, ob_err_arr):
    for ob in enumerate(ob_arr):
        ob_arr[ob[0]] = ob[1] + np.random.normal(0, ob_err_arr[ob[0]])
    return ob_arr


def run_4dvar_desroziers(pvals, east_west):
    d = dc.DalecData(2015, 2016, 'nee_day_'+east_west+', nee_night_'+east_west+', c_roo_'+east_west+','
                     ' c_woo_'+east_west+', clma, lai_'+east_west,
                     nc_file='../../alice_holt_data/ah_data_daily_test_nee.nc')
    d.B = b
    d.ob_err_dict['clma'] = 0.33 * d.ob_err_dict['clma']
    d.ob_err_dict['lai'] = 0.33 * d.ob_err_dict['lai']
    m = mc.DalecModel(d)
    m.yoblist = perturb_obs(m.yoblist, m.yerroblist)
    out = m.find_min_tnc_cvt(pvals, dispp=0)
    return m.yoblist, pvals, out


def create_ensemble_trunc_uc(d, mean, cov):
    """ Creates an ensemble of parameter values satisfying ecological constraints.
    """
    param_ensemble = []
    failed_ensemble = []
    while len(param_ensemble) < 1500:
        pvals = np.random.multivariate_normal(mean, cov)
        if pvals_test_uc(d, pvals) == True:
            param_ensemble.append(pvals)
            print '%i' %len(param_ensemble)
        else:
            failed_ensemble.append(pvals)
            continue
    return param_ensemble, failed_ensemble


def pvals_test_uc(d, pvals):
    """ Test if pvals pass constraints.
    """
    #calculate carbon fluxes
    f_auto = pvals[1]
    f_fol = (1-f_auto)*pvals[2]
    f_lab = (1-f_auto-f_fol)*pvals[12]
    f_roo = (1-f_auto-f_fol-f_lab)*pvals[3]
    f_woo = (1 - f_auto - f_fol - f_lab - f_roo)

    #universal constraint tests
    uc = [10*pvals[16] > pvals[18],
          pvals[8] < pvals[7],
          pvals[0] > pvals[8],
          pvals[5] < 1/(365.25*pvals[4]),
          pvals[6] > pvals[8]*np.exp(pvals[9]*d.t_mean[0]),
          0.2*f_roo < (f_fol+f_lab) < 5*f_roo,
          pvals[14]-pvals[11]>45]
    if all(uc) == True:
        return True
    else:
        return False
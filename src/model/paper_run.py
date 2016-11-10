import data_class as dc
import mod_class as mc
import numpy as np
import sympy as smp
import plot as p
import re
import os
import pickle


d = dc.DalecData(2015, 2016, 'clma', nc_file='../../alice_holt_data/ah_data_daily_test_nee2.nc', scale_nee=0)


def experiment_bmat(f_name):
    # Construct B
    b_cor = pickle.load(open('b_edc_cor.p', 'r'))
    b_std = np.sqrt(np.diag(pickle.load(open('b_edc.p', 'r'))))
    b_std[0:17] = b_std[0:17]*0.5
    D = np.zeros_like(b_cor)
    np.fill_diagonal(D, b_std)
    b = 0.6*np.dot(np.dot(D, b_cor), D)  #*0.6
    experiment(f_name, b)
    return 'done!'


def experiment_bmat_ceff(f_name):
    # Construct B
    b_cor = pickle.load(open('b_edc_cor.p', 'r'))
    b_std = np.sqrt(np.diag(pickle.load(open('b_edc.p', 'r'))))
    b_std[10] = 0.25*b_std[10]
    b_std[1] = 0.25*b_std[1]
    #b_std[2] = 0.1*b_std[2]
    b_std[0:17] = b_std[0:17]*0.5
    D = np.zeros_like(b_cor)
    np.fill_diagonal(D, b_std)
    b = 0.6*np.dot(np.dot(D, b_cor), D)  #*0.6
    experiment(f_name, b)
    return 'done!'


def experiment(f_name, b_mat, xb=d.xb_ew_lai_hi):
    east_west_joint_run(xb, f_name+'nee/', 'nee, clma', b_mat)
    east_west_joint_run(xb, f_name+'needn/', 'nee_day, nee_night, clma', b_mat)
    east_west_joint_run(xb, f_name+'lai/', 'lai, clma', b_mat)
    east_west_joint_run(xb, f_name+'needn_lai/', 'nee_day, nee_night, lai, clma', b_mat)
    east_west_joint_run(xb, f_name+'needn_lai_cw/', 'nee_day, nee_night, lai, clma, c_woo', b_mat)
    east_west_joint_run(xb, f_name+'needn_cw/', 'nee_day, nee_night, c_woo', b_mat)
    east_west_joint_run(xb, f_name+'needn_lai_cw_cr/', 'nee_day, nee_night, lai, clma, c_woo, c_roo', b_mat)
    east_west_joint_run(xb, f_name+'nee_needn_lai_cw_cr/', 'nee, nee_day, nee_night, lai, clma, c_woo, c_roo', b_mat)
    east_west_joint_run(xb, f_name+'nee_lai_cw_cr/', 'nee, lai, clma, c_woo, c_roo', b_mat)
    east_west_joint_run(xb, f_name+'nee_lai/', 'nee, lai, clma', b_mat)
    east_west_joint_run(xb, f_name+'neeconst_needn_lai_cw_cr/', 'nee, nee_day, nee_night, lai, clma, c_woo, c_roo',
                        b_mat, rm='nee')
    return 'done!'


def save_plots(f_name, xb, xa_east, xa_west, d_e, d_w, me, mw):
    # Plot 4dvar time series
    ax, fig = p.plot_obs_east_west('nee', xa_east, xa_west, d_e, d_w)
    fig.savefig(f_name+'_nee.png', bbox_inches='tight')
    ax, fig = p.plot_obs_east_west('nee_night', xa_east, xa_west, d_e, d_w)
    fig.savefig(f_name+'_neen.png', bbox_inches='tight')
    ax, fig = p.plot_obs_east_west('nee_day', xa_east, xa_west, d_e, d_w)
    fig.savefig(f_name+'_need.png', bbox_inches='tight')
    ax, fig = p.plot_obs_east_west('lai', xa_east, xa_west, d_e, d_w)
    fig.savefig(f_name+'_lai.png', bbox_inches='tight')
    ax, fig = p.plot_obs_east_west('c_woo', xa_east, xa_west, d_e, d_w)
    fig.savefig(f_name+'_cwoo.png', bbox_inches='tight')
    ax, fig = p.plot_obs_east_west('gpp', xa_east, xa_west, d_e, d_w)
    fig.savefig(f_name+'_gpp.png', bbox_inches='tight')
    ax, fig = p.plot_obs_east_west('rt', xa_east, xa_west, d_e, d_w)
    fig.savefig(f_name+'_rt.png', bbox_inches='tight')
    ax, fig = p.plot_obs_east_west('rh', xa_east, xa_west, d_e, d_w)
    fig.savefig(f_name+'_rh.png', bbox_inches='tight')
    ax, fig = p.plot_obs_east_west('ra', xa_east, xa_west, d_e, d_w)
    fig.savefig(f_name+'_ra.png', bbox_inches='tight')

    ax, fig = p.plot_obs_east_west_cum('nee_day', xa_east, xa_west, d_e, d_w,
                                     y_label='Cumulative NEE (g C m$^{-2}$ day$^{-1}$)')
    fig.savefig(f_name+'nee_day_cum.png', bbox_inches='tight')
    ax, fig = p.plot_obs_east_west_cum('nee', xa_east, xa_west, d_e, d_w,
                                     y_label='Cumulative NEE (g C m$^{-2}$ day$^{-1}$)')
    fig.savefig(f_name+'nee_cum.png', bbox_inches='tight')
    ax, fig = p.plot_obs_east_west_cum('rt', xa_east, xa_west, d_e, d_w,
                                     y_label='Cumulative ecosystem respiration (g C m$^{-2}$ day$^{-1}$)')
    fig.savefig(f_name+'rt_cum.png', bbox_inches='tight')
    ax, fig = p.plot_obs_east_west_cum('gpp', xa_east, xa_west, d_e, d_w,
                                     y_label='Cumulative GPP (g C m$^{-2}$ day$^{-1}$)')
    fig.savefig(f_name+'gpp_cum.png', bbox_inches='tight')

    # Plot error in analysis and background
    ax, fig = p.plot_inc_east_west(xb, xa_east, xa_west)
    fig.savefig(f_name+'_xa_inc.pdf', bbox_inches='tight')
    # Plot error cov mats
    # ax, fig = p.plot_bmat(p.cov2cor(me.dC.B))
    # fig.savefig(f_name+'_bmat.png', bbox_inches='tight')
    ax, fig = p.plot_rmat(p.cov2cor(me.rmatrix))
    fig.savefig(f_name+'_rmat_east.png', bbox_inches='tight')
    ax, fig = p.plot_rmat(p.cov2cor(mw.rmatrix))
    fig.savefig(f_name+'_rmat_west.png', bbox_inches='tight')
    return 'done'


def east_west_joint_run(xb, f_name, obs_str, b_mat, rm='None'):
    if not os.path.exists(f_name):
        os.makedirs(f_name)
    # east data
    obs_east = ob_str_east_west(obs_str, 'east', rm_obs=rm)
    de = dc.DalecData(2015, 2016, obs_east,
                      nc_file='../../alice_holt_data/ah_data_daily_test_nee2.nc', scale_nee=1)
    de.B = b_mat
    # obs err scaling
    # west data
    obs_west = ob_str_east_west(obs_str, 'west', rm_obs=rm)
    dw = dc.DalecData(2015, 2016, obs_west,
                      nc_file='../../alice_holt_data/ah_data_daily_test_nee2.nc', scale_nee=1)
    dw.B = b_mat
    # obs err scaling
    # setup model
    me = mc.DalecModel(de)
    me.rmatrix = r_mat_corr(me.yerroblist, me.ytimestep, me.y_strlst, me.rmatrix, corr=0.3, tau=2.)[1]
    mw = mc.DalecModel(dw)
    mw.rmatrix = r_mat_corr(mw.yerroblist, mw.ytimestep, mw.y_strlst, mw.rmatrix, corr=0.3, tau=2.)[1]
    # run DA scheme
    xa_e = me.find_min_tnc_cvt(xb, f_name+'east_assim')
    xa_w = mw.find_min_tnc_cvt(xb, f_name+'west_assim')
    # save plots
    save_plots(f_name, xb, xa_e[1], xa_w[1], de, dw, me, mw)
    return 'done'


# ------------------------------------------------------------------------------
# R matrix
# ------------------------------------------------------------------------------

def r_mat_corr(yerroblist, ytimestep, y_strlst, r_diag, corr=0.3, tau=1., cut_off=4.):
    """ Creates a correlated R matrix.
    """
    r_corr = np.eye(len(ytimestep)) #MAKE SURE ALL VALUES ARE FLOATS FIRST!!!!
    for i in xrange(len(ytimestep)):
        if y_strlst[i] == 'nee_day' or y_strlst[i] == 'nee_night':
            for j in xrange(len(ytimestep)):
                if y_strlst[j] == 'nee_day' or y_strlst[j] == 'nee_night':
                    if abs(ytimestep[i]-ytimestep[j]) < cut_off:
                        r_corr[i, j] = corr*np.exp(-(abs(float(ytimestep[i])-float(ytimestep[j]))**2)/float(tau)**2) \
                                      + (1-corr)*smp.KroneckerDelta(ytimestep[i], ytimestep[j])
                    if y_strlst[j] == 'nee_day' and y_strlst[i] == 'nee_night':
                        r_corr[i, j] = corr*np.exp(-(abs(float(ytimestep[i])-float(ytimestep[j]))**2)/float(tau)**2)
                    elif y_strlst[i] == 'nee_day' and y_strlst[j] == 'nee_night':
                        r_corr[i, j] = corr*np.exp(-(abs(float(ytimestep[i])-float(ytimestep[j]))**2)/float(tau)**2)
    r = np.dot(np.dot((np.sqrt(r_diag)), r_corr), np.sqrt(r_diag))
    return r_corr, r


# ------------------------------------------------------------------------------
# Ob_str
# ------------------------------------------------------------------------------

def ob_str_east_west(ob_str, east_west, rm_obs='None'):
    obs_lst = re.findall(r'[^,;\s]+', ob_str)
    ob_east_west = ['nee', 'nee_day', 'nee_night', 'lai', 'c_woo', 'c_roo']
    if rm_obs != 'None':
        ob_east_west.remove(rm_obs)
    new_ob_str = ''
    for ob in obs_lst:
        if ob in ob_east_west:
            ob = ob + '_' + east_west
        else:
            ob = ob
        new_ob_str += ob + ', '
    return new_ob_str


# ------------------------------------------------------------------------------
# Prior
# ------------------------------------------------------------------------------


def east_west_joint_run_prior(xb, f_name, obs_str, b_mat, rm='None'):
    if not os.path.exists(f_name):
        os.makedirs(f_name)
    # east data
    obs_east = ob_str_east_west(obs_str, 'east', rm_obs=rm)
    de = dc.DalecData(2012, 2014, obs_east,
                      nc_file='../../alice_holt_data/ah_data_daily_test_nee2.nc', scale_nee=1)
    de.B = b_mat
    # obs err scaling
    # west data
    obs_west = ob_str_east_west(obs_str, 'west', rm_obs=rm)
    dw = dc.DalecData(2012, 2014, obs_west,
                      nc_file='../../alice_holt_data/ah_data_daily_test_nee2.nc', scale_nee=1)
    dw.B = b_mat
    # obs err scaling
    # setup model
    me = mc.DalecModel(de)
    me.rmatrix = r_mat_corr(me.yerroblist, me.ytimestep, me.y_strlst, me.rmatrix, corr=0.3, tau=2.)[1]
    mw = mc.DalecModel(dw)
    mw.rmatrix = r_mat_corr(mw.yerroblist, mw.ytimestep, mw.y_strlst, mw.rmatrix, corr=0.3, tau=2.)[1]
    # run DA scheme
    xa_e = me.find_min_tnc_cvt(xb, f_name+'east_assim')
    xa_w = mw.find_min_tnc_cvt(xb, f_name+'west_assim')
    # save plots
    save_plots(f_name, xb, xa_e[1], xa_w[1], de, dw, me, mw)
    return 'done'


def experiment_prior(f_name, b_mat, xb=d.xb_ew_lai_hi):
    east_west_joint_run(xb, f_name+'nee/', 'nee, clma', b_mat)
    east_west_joint_run(xb, f_name+'needn/', 'nee_day, nee_night, clma', b_mat)
    east_west_joint_run(xb, f_name+'lai/', 'lai, clma', b_mat)
    east_west_joint_run(xb, f_name+'needn_lai/', 'nee_day, nee_night, lai, clma', b_mat)
    east_west_joint_run(xb, f_name+'nee_needn_lai_cw_cr/', 'nee, nee_day, nee_night, lai, clma', b_mat)
    east_west_joint_run(xb, f_name+'nee_lai/', 'nee, lai, clma', b_mat)
    east_west_joint_run(xb, f_name+'neeconst_needn_lai_cw_cr/', 'nee, nee_day, nee_night, lai, clma, c_woo, c_roo',
                        b_mat, rm='nee')
    return 'done!'


def experiment_prior_run(f_name):
    # Construct B
    b = pickle.load(open('b_edc.p', 'r'))
    experiment_prior(f_name, b, xb=d.edinburgh_mean)
    return 'done!'
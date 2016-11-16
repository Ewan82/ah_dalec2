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
    b_std[0:17] = b_std[0:17]  # *0.5
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
    #b_std[1] = 0.25*b_std[1]
    #b_std[2] = 0.25*b_std[2] # Maybe get rid of this constraint
    b_std[0:17] = b_std[0:17]  # *0.5
    D = np.zeros_like(b_cor)
    np.fill_diagonal(D, b_std)
    b = 0.6*np.dot(np.dot(D, b_cor), D)  #*0.6
    experiment(f_name, b)
    return 'done!'


def experiment_bmat_ceff_fauto(f_name):
    # Construct B
    b_cor = pickle.load(open('b_edc_cor.p', 'r'))
    b_std = np.sqrt(np.diag(pickle.load(open('b_edc.p', 'r'))))
    b_std[10] = 0.25*b_std[10]
    b_std[1] = 0.25*b_std[1]
    #b_std[2] = 0.25*b_std[2] # Maybe get rid of this constraint
    b_std[0:17] = b_std[0:17]  # *0.5
    D = np.zeros_like(b_cor)
    np.fill_diagonal(D, b_std)
    b = 0.6*np.dot(np.dot(D, b_cor), D)  #*0.6
    experiment(f_name, b)
    return 'done!'


def experiment_bmat_ceff_ffol(f_name):
    # Construct B
    b_cor = pickle.load(open('b_edc_cor.p', 'r'))
    b_std = np.sqrt(np.diag(pickle.load(open('b_edc.p', 'r'))))
    b_std[10] = 0.25*b_std[10]
    #b_std[1] = 0.25*b_std[1]
    b_std[2] = 0.25*b_std[2] # Maybe get rid of this constraint
    b_std[0:17] = b_std[0:17]  # *0.5
    D = np.zeros_like(b_cor)
    np.fill_diagonal(D, b_std)
    b = 0.6*np.dot(np.dot(D, b_cor), D)  #*0.6
    experiment(f_name, b)
    return 'done!'


def experiment_bmat_ceff_fauto_ffol(f_name):
    # Construct B
    b_cor = pickle.load(open('b_edc_cor.p', 'r'))
    b_std = np.sqrt(np.diag(pickle.load(open('b_edc.p', 'r'))))
    b_std[10] = 0.25*b_std[10]
    b_std[1] = 0.25*b_std[1]
    b_std[2] = 0.25*b_std[2] # Maybe get rid of this constraint
    b_std[0:17] = b_std[0:17]  # *0.5
    D = np.zeros_like(b_cor)
    np.fill_diagonal(D, b_std)
    b = 0.6*np.dot(np.dot(D, b_cor), D)  #*0.6
    experiment(f_name, b)
    return 'done!'


def experiment_bmat_ceff_fauto_ffol_flab(f_name):
    # Construct B
    b_cor = pickle.load(open('b_edc_cor.p', 'r'))
    b_std = np.sqrt(np.diag(pickle.load(open('b_edc.p', 'r'))))
    b_std[10] = 0.25*b_std[10]
    b_std[1] = 0.25*b_std[1]
    b_std[2] = 0.25*b_std[2] # Maybe get rid of this constraint
    b_std[12] = 0.25*b_std[12]
    b_std[0:17] = b_std[0:17]  # *0.5
    D = np.zeros_like(b_cor)
    np.fill_diagonal(D, b_std)
    b = 0.6*np.dot(np.dot(D, b_cor), D)  #*0.6
    experiment(f_name, b)
    return 'done!'


def experiment(f_name, b_mat, xb=d.xb_ew_lai_hi):
    #east_west_joint_run(xb, f_name+'nee/', 'nee, clma', b_mat)
    #east_west_joint_run(xb, f_name+'needn/', 'nee_day, nee_night, clma', b_mat)
    east_west_joint_run(xb, f_name+'nee_needn/', 'nee, nee_day, nee_night', b_mat)
    #east_west_joint_run(xb, f_name+'lai/', 'lai, clma', b_mat)
    #east_west_joint_run(xb, f_name+'cw/', 'c_woo, clma', b_mat)
    #east_west_joint_run(xb, f_name+'needn_lai/', 'nee_day, nee_night, lai, clma', b_mat)
    #east_west_joint_run(xb, f_name+'needn_lai_cw/', 'nee_day, nee_night, lai, clma, c_woo', b_mat)
    #east_west_joint_run(xb, f_name+'needn_cw/', 'nee_day, nee_night, c_woo', b_mat)
    #east_west_joint_run(xb, f_name+'needn_lai_cw_cr/', 'nee_day, nee_night, lai, clma, c_woo, c_roo', b_mat)
    east_west_joint_run(xb, f_name+'nee_needn_lai/', 'nee, nee_day, nee_night, lai, clma', b_mat)
    east_west_joint_run(xb, f_name+'nee_needn_lai_cw/', 'nee, nee_day, nee_night, lai, clma, c_woo', b_mat)
    #east_west_joint_run(xb, f_name+'nee_lai_cw_cr/', 'nee, lai, clma, c_woo, c_roo', b_mat)
    #east_west_joint_run(xb, f_name+'nee_lai/', 'nee, lai, clma', b_mat)
    #east_west_joint_run(xb, f_name+'neeconst_needn_lai_cw_cr/', 'nee, nee_day, nee_night, lai, clma, c_woo, c_roo',
    #                    b_mat, rm='nee')
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

    ax, fig, cum_east, cum_west = p.plot_obs_east_west_cum('nee_day', xa_east, xa_west, d_e, d_w,
                                     y_label='Cumulative NEE (g C m$^{-2}$ day$^{-1}$)')
    fig.savefig(f_name+'nee_day_cum.png', bbox_inches='tight')
    ax, fig, cum_east, cum_west = p.plot_obs_east_west_cum('nee', xa_east, xa_west, d_e, d_w,
                                     y_label='Cumulative NEE (g C m$^{-2}$ day$^{-1}$)')
    fig.savefig(f_name+'nee_cum.png', bbox_inches='tight')
    ax, fig, cum_east, cum_west = p.plot_obs_east_west_cum('rt', xa_east, xa_west, d_e, d_w,
                                     y_label='Cumulative ecosystem respiration (g C m$^{-2}$ day$^{-1}$)')
    fig.savefig(f_name+'rt_cum.png', bbox_inches='tight')
    ax, fig, cum_east, cum_west = p.plot_obs_east_west_cum('gpp', xa_east, xa_west, d_e, d_w,
                                     y_label='Cumulative GPP (g C m$^{-2}$ day$^{-1}$)')
    fig.savefig(f_name+'gpp_cum.png', bbox_inches='tight')

    # Plot error in analysis and background
    ax, fig = p.plot_inc_east_west(xb, xa_east, xa_west)
    fig.savefig(f_name+'_xa_inc.png', bbox_inches='tight')
    ax, fig = p.plot_table(xb, xa_east, xa_west)
    fig.savefig(f_name+'_table.png', bbox_inches='tight')
    # Plot error cov mats
    b = d_e.B
    a_east = me.acovmat(xa_east)
    a_west = mw.acovmat(xa_west)
    ax, fig = p.plot_var_red_east_west(b, a_east, a_west)
    fig.savefig(f_name+'var_red.png', bbox_inches='tight')
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
                      nc_file='../../alice_holt_data/ah_data_daily_test_nee3.nc', scale_nee=1)
    de.B = b_mat
    # obs err scaling
    # west data
    obs_west = ob_str_east_west(obs_str, 'west', rm_obs=rm)
    dw = dc.DalecData(2015, 2016, obs_west,
                      nc_file='../../alice_holt_data/ah_data_daily_test_nee3.nc', scale_nee=1)
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
        if y_strlst[i] == 'nee':
            for j in xrange(len(ytimestep)):
                if y_strlst[j] == 'nee':
                    if abs(ytimestep[i]-ytimestep[j]) < cut_off:
                        r_corr[i, j] = corr*np.exp(-(abs(float(ytimestep[i])-float(ytimestep[j]))**2)/float(tau)**2) \
                                      + (1-corr)*smp.KroneckerDelta(ytimestep[i], ytimestep[j])
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


def east_west_joint_run_prior(xb, f_name, obs_str, b_mat, rm='None', end_yr=2014, start_yr=2012):
    if not os.path.exists(f_name):
        os.makedirs(f_name)
    # east data
    obs_east = ob_str_east_west(obs_str, 'east', rm_obs=rm)
    de = dc.DalecData(start_yr, end_yr, obs_east,
                      nc_file='../../alice_holt_data/ah_data_daily_test_nee3.nc', scale_nee=0)
    de.B = b_mat
    # obs err scaling
    # west data
    obs_west = ob_str_east_west(obs_str, 'west', rm_obs=rm)
    dw = dc.DalecData(start_yr, end_yr, obs_west,
                      nc_file='../../alice_holt_data/ah_data_daily_test_nee3.nc', scale_nee=0)
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
    east_west_joint_run_prior(xb, f_name+'nee_needn/', 'nee, nee_day, nee_night, clma', b_mat)
    east_west_joint_run_prior(xb, f_name+'nee_needn_1314/', 'nee, nee_day, nee_night, clma', b_mat, end_yr=2015,
                              start_yr=2013)
    east_west_joint_run_prior(xb, f_name+'nee/', 'nee, clma', b_mat)
    east_west_joint_run_prior(xb, f_name+'needn/', 'nee_day, nee_night, clma', b_mat)
    east_west_joint_run_prior(xb, f_name+'lai/', 'lai, clma', b_mat)
    east_west_joint_run_prior(xb, f_name+'needn_lai/', 'nee_day, nee_night, lai, clma', b_mat)
    east_west_joint_run_prior(xb, f_name+'nee_needn_lai/', 'nee, nee_day, nee_night, lai, clma', b_mat)
    east_west_joint_run_prior(xb, f_name+'nee_lai/', 'nee, lai, clma', b_mat)
    east_west_joint_run_prior(xb, f_name+'neeconst_needn_lai/', 'nee, nee_day, nee_night, lai, clma,',
                        b_mat, rm='nee')
    return 'done!'


def experiment_prior_run(f_name):
    # Construct B
    b_cor = pickle.load(open('b_edc_cor.p', 'r'))
    b_std = np.sqrt(np.diag(pickle.load(open('b_edc.p', 'r'))))
    b_std[0:17] = b_std[0:17]*0.5
    D = np.zeros_like(b_cor)
    np.fill_diagonal(D, b_std)
    b = 0.6*np.dot(np.dot(D, b_cor), D)  #*0.6
    experiment_prior(f_name, b, xb=d.xb_ew_lai_hi)
    return 'done!'


def experiment_prior_run2(f_name):
    # Construct B
    #b_cor = pickle.load(open('b_edc_cor.p', 'r'))
    #b_std = np.sqrt(np.diag(pickle.load(open('b_edc.p', 'r'))))
    #b_std[0:17] = b_std[0:17]*0.5
    #D = np.zeros_like(b_cor)
    #np.fill_diagonal(D, b_std)
    #b = 0.6*np.dot(np.dot(D, b_cor), D)  #*0.6
    b = pickle.load(open('b_edc.p', 'r'))
    experiment_prior(f_name, b, xb=d.edinburgh_median)
    return 'done!'


def save_paper_plots(f_name, exp_name):
    if not os.path.exists(f_name):
        os.makedirs(f_name)
    east = pickle.load(open(exp_name+'east_assim', 'r'))
    west = pickle.load(open(exp_name+'west_assim', 'r'))
    b = east['b_mat']
    # east data
    de = dc.DalecData(2015, 2016, 'clma',
                      nc_file='../../alice_holt_data/ah_data_daily_test_nee2.nc', scale_nee=1)
    de.B = b
    de.ob_dict = east['obs']
    de.ob_err_dict = east['obs_err']
    # obs err scaling
    # west data
    dw = dc.DalecData(2015, 2016, 'clma',
                      nc_file='../../alice_holt_data/ah_data_daily_test_nee2.nc', scale_nee=1)
    dw.B = b
    dw.ob_dict = west['obs']
    dw.ob_err_dict = west['obs_err']
    # obs err scaling
    # setup model
    me = mc.DalecModel(de)
    me.rmatrix = r_mat_corr(me.yerroblist, me.ytimestep, me.y_strlst, me.rmatrix, corr=0.3, tau=2.)[1]
    mw = mc.DalecModel(dw)
    mw.rmatrix = r_mat_corr(mw.yerroblist, mw.ytimestep, mw.y_strlst, mw.rmatrix, corr=0.3, tau=2.)[1]
    # a_east = pickle.load(open('a_east.p', 'r'))
    # a_west = pickle.load(open('a_west.p', 'r'))
    a_east = me.acovmat(east['xa'])
    a_west = mw.acovmat(west['xa'])

    annual_flux_lst = []
    ax, fig = p.plot_var_red_east_west(b, a_east, a_west)
    fig.savefig(f_name+'var_red.png', bbox_inches='tight')

    ax, fig = p.plot_east_west_paper('rh', east['xa'], west['xa'], de, dw, a_east, a_west,
                                     y_label=r'Heterotrophic respiration (g C m$^{-2}$ day$^{-1}$)')
    fig.savefig(f_name+'rh.png', bbox_inches='tight')
    ax, fig = p.plot_east_west_paper('nee_day', east['xa'], west['xa'], de, dw, a_east, a_west,
                                     y_label=r'NEE$_{day}$ (g C m$^{-2}$ day$^{-1}$)', y_lim=[-15, 5])
    fig.savefig(f_name+'nee_day.png', bbox_inches='tight')
    ax, fig, cum_east, cum_west = p.plot_east_west_paper_cum('nee_day', east['xa'], west['xa'], de, dw, a_east, a_west,
                                     y_label='Cumulative NEE (g C m$^{-2}$ day$^{-1}$)')
    fig.savefig(f_name+'nee_day_cum.png', bbox_inches='tight')
    ax, fig, cum_east, cum_west = p.plot_east_west_paper_cum('nee', east['xa'], west['xa'], de, dw, a_east, a_west,
                                     y_label='Cumulative NEE (g C m$^{-2}$ day$^{-1}$)')
    fig.savefig(f_name+'nee_cum.png', bbox_inches='tight')
    annual_flux_lst.append(cum_east)
    annual_flux_lst.append(cum_west)
    ax, fig = p.plot_east_west_paper('nee_night', east['xa'], west['xa'], de, dw, a_east, a_west,
                                     y_label=r'NEE$_{night}$ (g C m$^{-2}$ day$^{-1}$)')
    fig.savefig(f_name+'nee_night.png', bbox_inches='tight')
    ax, fig = p.plot_east_west_paper('gpp', east['xa'], west['xa'], de, dw, a_east, a_west,
                                     y_label=r'Gross primary production (g C m$^{-2}$ day$^{-1}$)')
    fig.savefig(f_name+'gpp.png', bbox_inches='tight')
    ax, fig, cum_east, cum_west = p.plot_east_west_paper_cum('gpp', east['xa'], west['xa'], de, dw, a_east, a_west,
                                     y_label='Cumulative GPP (g C m$^{-2}$ day$^{-1}$)')
    fig.savefig(f_name+'gpp_cum.png', bbox_inches='tight')
    annual_flux_lst.append(cum_east)
    annual_flux_lst.append(cum_west)
    ax, fig = p.plot_east_west_paper('lai', east['xa'], west['xa'], de, dw, a_east, a_west,
                                     y_label=r'Leaf area index')
    fig.savefig(f_name+'lai.png', bbox_inches='tight')
    ax, fig = p.plot_east_west_paper('c_woo', east['xa'], west['xa'], de, dw, a_east, a_west,
                                     y_label=r'Woody biomass and coarse root carbon (g C m$^{-2}$)',
                                     y_lim=[9000, 14500])
    fig.savefig(f_name+'c_woo.png', bbox_inches='tight')
    ax, fig = p.plot_east_west_paper('ra', east['xa'], west['xa'], de, dw, a_east, a_west,
                                     y_label=r'Autotrophic respiration (g C m$^{-2}$ day$^{-1}$)')
    fig.savefig(f_name+'ra.png', bbox_inches='tight')
    ax, fig = p.plot_east_west_paper('rt', east['xa'], west['xa'], de, dw, a_east, a_west,
                                     y_label=r'Total ecosystem respiration (g C m$^{-2}$ day$^{-1}$)')
    fig.savefig(f_name+'rt.png', bbox_inches='tight')
    ax, fig, cum_east, cum_west = p.plot_east_west_paper_cum('rt', east['xa'], west['xa'], de, dw, a_east, a_west,
                                     y_label='Cumulative ecosystem respiration (g C m$^{-2}$ day$^{-1}$)')
    fig.savefig(f_name+'rt_cum.png', bbox_inches='tight')
    annual_flux_lst.append(cum_east)
    annual_flux_lst.append(cum_west)
    ax, fig = p.plot_inc_east_west(east['xb'], east['xa'], west['xa'])
    fig.savefig(f_name+'xa_inc.png', bbox_inches='tight')
    f = open(f_name+'annual_fluxes.txt', 'w')
    for item in annual_flux_lst:
        f.write("%s\n" % item)
    f.close()
    return 'done!'


def do_plots(f_name):
    for item in ['nee_needn', 'nee_needn_lai', 'nee_needn_lai_cw']:
        save_paper_plots(f_name+item+'_pp/', f_name+item+'/')
    return 'done'
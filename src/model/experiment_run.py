import data_class as dc
import mod_class as mc
import mod_class_param_choice as mc_p
import numpy as np
import plot as p
import pickle


def save_plots(f_name, xb, xa_east, xa_west, d_e, d_w):
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
    ax, fig = p.plot_obs_east_west('c_roo', xa_east, xa_west, d_e, d_w)
    fig.savefig(f_name+'_croo.png', bbox_inches='tight')
    ax, fig = p.plot_obs_east_west('gpp', xa_east, xa_west, d_e, d_w)
    fig.savefig(f_name+'_gpp.png', bbox_inches='tight')
    ax, fig = p.plot_obs_east_west('rt', xa_east, xa_west, d_e, d_w)
    fig.savefig(f_name+'_rt.png', bbox_inches='tight')
    ax, fig = p.plot_obs_east_west('rh', xa_east, xa_west, d_e, d_w)
    fig.savefig(f_name+'_rh.png', bbox_inches='tight')
    ax, fig = p.plot_obs_east_west('ra', xa_east, xa_west, d_e, d_w)
    fig.savefig(f_name+'_ra.png', bbox_inches='tight')
    ax, fig = p.plot_obs_east_west('cl', xa_east, xa_west, d_e, d_w)
    fig.savefig(f_name+'_clit.png', bbox_inches='tight')
    ax, fig = p.plot_obs_east_west('cf', xa_east, xa_west, d_e, d_w)
    fig.savefig(f_name+'_cfol.png', bbox_inches='tight')

    # Plot error in analysis and background
    ax, fig = p.plot_inc_east_west(xb, xa_east, xa_west)
    fig.savefig(f_name+'_xa_inc.png', bbox_inches='tight')
    return 'done'


def east_west_joint_run(xb, f_name):
    de = dc.DalecData(2015, 2016, 'nee_day_east, c_roo_east, c_woo_east, clma, lai_east',
                      nc_file='../../alice_holt_data/ah_data_daily_test_nee.nc')
    de.ob_err_dict['clma'] = (1./3.) * de.ob_err_dict['clma']
    de.ob_err_dict['lai'] = (1./3.) * de.ob_err_dict['lai']
    B = pickle.load(open('b_edc.p', 'r'))
    B = np.delete(B, 10, axis=0)
    B = 0.8*np.delete(B, 10, axis=1)
    de.B = B
    #de.B = np.concatenate((de.edinburgh_std[:10], de.edinburgh_std[11:]))**2*np.eye(22)
    #A = pickle.load(open('A_D.p', 'r'))
    #A = np.delete(A, (10), axis=0)
    #A = np.delete(A, (10), axis=1)
    #de.B = 1.2*A
    dw = dc.DalecData(2015, 2016, 'nee_day_west, c_roo_west, c_woo_west, clma, lai_west',
                      nc_file='../../alice_holt_data/ah_data_daily_test_nee.nc')
    dw.ob_err_dict['clma'] = (1./3.) * dw.ob_err_dict['clma']
    dw.ob_err_dict['lai'] = (1./3.) * dw.ob_err_dict['lai']
    dw.B = B
    #dw.B = np.concatenate((de.edinburgh_std[:10], de.edinburgh_std[11:]))**2*np.eye(22)
    #dw.B = A
    me = mc_p.DalecModel(de)
    mw = mc_p.DalecModel(dw)
    xa_e = me.find_min_tnc_cvt(xb, f_name+'east_assim')
    xa_w = mw.find_min_tnc_cvt(xb, f_name+'west_assim')
    xbb = np.array(me.create_ordered_lst(np.array(xb.tolist()[0], dtype=np.float)))
    save_plots(f_name, xbb, xa_e[2], xa_w[2], de, dw)
    return 'done'


def east_west_joint_run_ceff_fauto(xb, f_name):
    de = dc.DalecData(2015, 2016, 'nee_day_east, nee_night_east, c_roo_east, c_woo_east, clma, lai_east',
                      nc_file='../../alice_holt_data/ah_data_daily_test_nee.nc')
    #de.ob_err_dict['clma'] = (1./3.) * de.ob_err_dict['clma']
    #de.ob_err_dict['lai'] = (1./3.) * de.ob_err_dict['lai']
    B = pickle.load(open('b_edc.p', 'r'))
    B = np.delete(B, 10, axis=0)
    B = np.delete(B, 10, axis=1)
    B = np.delete(B, 1, axis=0)
    B = 0.8*np.delete(B, 1, axis=1)
    de.B = B
    #de.B = np.concatenate((de.edinburgh_std[:10], de.edinburgh_std[11:]))**2*np.eye(22)
    #A = pickle.load(open('A_D.p', 'r'))
    #A = np.delete(A, (10), axis=0)
    #A = np.delete(A, (10), axis=1)
    #de.B = 1.2*A
    dw = dc.DalecData(2015, 2016, 'nee_day_west, nee_night_west, c_roo_west, c_woo_west, clma, lai_west',
                      nc_file='../../alice_holt_data/ah_data_daily_test_nee.nc')
    #dw.ob_err_dict['clma'] = (1./3.) * dw.ob_err_dict['clma']
    #dw.ob_err_dict['lai'] = (1./3.) * dw.ob_err_dict['lai']
    dw.B = B
    #dw.B = np.concatenate((de.edinburgh_std[:10], de.edinburgh_std[11:]))**2*np.eye(22)
    #dw.B = A
    me = mc_p.DalecModel(de)
    mw = mc_p.DalecModel(dw)
    xa_e = me.find_min_tnc_cvt(xb, f_name+'east_assim')
    xa_w = mw.find_min_tnc_cvt(xb, f_name+'west_assim')
    xbb = np.array(me.create_ordered_lst(np.array(xb.tolist()[0], dtype=np.float)))
    save_plots(f_name, xbb, xa_e[2], xa_w[2], de, dw)
    return 'done'


def east_west_joint_run_full(xb, f_name, clma_er=1., lai_er=1., need_er=1., neen_er=1., cr_er=1., cw_er=1.):
    f_name += 'clmaer%r_laier%r_needer%r_neener%r_crer%r_cwer%r' %(clma_er, lai_er, need_er, neen_er, cr_er, cw_er)
    # Construct B
    b_cor = pickle.load(open('b_edc_cor.p', 'r'))
    b_std = np.sqrt(np.diag(pickle.load(open('b_edc.p', 'r'))))
    b_std[10] = 0.25*b_std[10]
    b_std[1] = 0.25*b_std[1]
    b_std[0:17] = 0.5*b_std[0:17]
    D = np.zeros_like(b_cor)
    np.fill_diagonal(D, b_std)
    b = 0.6 * np.dot(np.dot(D, b_cor), D)
    # east data
    de = dc.DalecData(2015, 2016, 'nee_day_east, nee_night_east, c_roo_east, c_woo_east, clma, lai_east',
                      nc_file='../../alice_holt_data/ah_data_daily_test_nee.nc')
    de.B = b
    # obs err scaling
    de.ob_err_dict['clma'] = clma_er * de.ob_err_dict['clma']
    de.ob_err_dict['lai'] = lai_er * de.ob_err_dict['lai']
    de.ob_err_dict['nee_day'] = need_er * de.ob_err_dict['nee_day']
    de.ob_err_dict['nee_night'] = neen_er * de.ob_err_dict['nee_night']
    de.ob_err_dict['c_roo'] = cr_er * de.ob_err_dict['c_roo']
    de.ob_err_dict['c_woo'] = cw_er * de.ob_err_dict['c_woo']

    # west data
    dw = dc.DalecData(2015, 2016, 'nee_day_west, nee_night_west, c_roo_west, c_woo_west, clma, lai_west',
                      nc_file='../../alice_holt_data/ah_data_daily_test_nee.nc')
    dw.B = b
    # obs err scaling
    dw.ob_err_dict['clma'] = clma_er * dw.ob_err_dict['clma']
    dw.ob_err_dict['lai'] = lai_er * dw.ob_err_dict['lai']
    dw.ob_err_dict['nee_day'] = need_er * dw.ob_err_dict['nee_day']
    dw.ob_err_dict['nee_night'] = neen_er * dw.ob_err_dict['nee_night']
    dw.ob_err_dict['c_roo'] = cr_er * dw.ob_err_dict['c_roo']
    dw.ob_err_dict['c_woo'] = cw_er * dw.ob_err_dict['c_woo']
    # setup model
    me = mc.DalecModel(de)
    mw = mc.DalecModel(dw)
    # run DA scheme
    xa_e = me.find_min_tnc_cvt(xb, f_name+'east_assim')
    xa_w = mw.find_min_tnc_cvt(xb, f_name+'west_assim')
    #save plots
    save_plots(f_name, xb, xa_e[1], xa_w[1], de, dw)
    return 'done'


def east_west_joint_run_full2(xb, f_name, clma_er=1., lai_er=1., need_er=1., neen_er=1., cr_er=1., cw_er=1.):
    f_name += 'clmaer%r_laier%r_needer%r_neener%r_crer%r_cwer%r' %(clma_er, lai_er, need_er, neen_er, cr_er, cw_er)
    # Construct B
    b_cor = pickle.load(open('b_edc_cor.p', 'r'))
    b_std = np.sqrt(np.diag(pickle.load(open('b_edc.p', 'r'))))
    b_std[10] = 0.25*b_std[10]
    b_std[1] = 0.25*b_std[1]
    b_std[0:17] = 0.5*b_std[0:17]
    D = np.zeros_like(b_cor)
    np.fill_diagonal(D, b_std)
    b = np.dot(np.dot(D, b_cor), D)
    # east data
    de = dc.DalecData(2015, 2016, 'nee_day_east, nee_night_east, c_roo_east, c_woo_east, clma, lai_east',
                      nc_file='../../alice_holt_data/ah_data_daily_test_nee.nc')
    de.B = b
    # obs err scaling
    de.ob_err_dict['clma'] = clma_er * de.ob_err_dict['clma']
    de.ob_err_dict['lai'] = lai_er * de.ob_err_dict['lai']
    de.ob_err_dict['nee_day'] = need_er * de.ob_err_dict['nee_day']
    de.ob_err_dict['nee_night'] = neen_er * de.ob_err_dict['nee_night']
    de.ob_err_dict['c_roo'] = cr_er * de.ob_err_dict['c_roo']
    de.ob_err_dict['c_woo'] = cw_er * de.ob_err_dict['c_woo']

    # west data
    dw = dc.DalecData(2015, 2016, 'nee_day_west, nee_night_west, c_roo_west, c_woo_west, clma, lai_west',
                      nc_file='../../alice_holt_data/ah_data_daily_test_nee.nc')
    dw.B = b
    # obs err scaling
    dw.ob_err_dict['clma'] = clma_er * dw.ob_err_dict['clma']
    dw.ob_err_dict['lai'] = lai_er * dw.ob_err_dict['lai']
    dw.ob_err_dict['nee_day'] = need_er * dw.ob_err_dict['nee_day']
    dw.ob_err_dict['nee_night'] = neen_er * dw.ob_err_dict['nee_night']
    dw.ob_err_dict['c_roo'] = cr_er * dw.ob_err_dict['c_roo']
    dw.ob_err_dict['c_woo'] = cw_er * dw.ob_err_dict['c_woo']
    # setup model
    me = mc.DalecModel(de)
    mw = mc.DalecModel(dw)
    # run DA scheme
    xa_e = me.find_min_tnc_cvt(xb, f_name+'east_assim')
    xa_w = mw.find_min_tnc_cvt(xb, f_name+'west_assim')
    #save plots
    save_plots(f_name, xb, xa_e[1], xa_w[1], de, dw)
    return 'done'


def east_west_joint_run_nee_err(xb, f_name, clma_er=1., lai_er=1., need_er=1., neen_er=1., cr_er=1., cw_er=1.):
    f_name += 'clmaer%r_laier%r_needer%r_neener%r_crer%r_cwer%r' %(clma_er, lai_er, need_er, neen_er, cr_er, cw_er)
    # Construct B
    b_cor = pickle.load(open('b_edc_cor.p', 'r'))
    b_std = np.sqrt(np.diag(pickle.load(open('b_edc.p', 'r'))))
    b_std[10] = 0.1*b_std[10]
    b_std[1] = 0.1*b_std[1]
    #b_std[0:17] = 0.5*b_std[0:17]
    D = np.zeros_like(b_cor)
    np.fill_diagonal(D, b_std)
    b = np.dot(np.dot(D, b_cor), D)
    # east data
    de = dc.DalecData(2015, 2016, 'nee_day_east, nee_night_east, c_roo_east, c_woo_east, clma, lai_east',
                      nc_file='../../alice_holt_data/ah_data_daily_test_nee.nc', scale_nee=1)
    de.B = b
    # obs err scaling
    de.ob_err_dict['clma'] = clma_er * de.ob_err_dict['clma']
    de.ob_err_dict['lai'] = lai_er * de.ob_err_dict['lai']
    de.ob_err_dict['nee_day'] = need_er * de.ob_err_dict['nee_day']
    de.ob_err_dict['nee_night'] = neen_er * de.ob_err_dict['nee_night']
    de.ob_err_dict['c_roo'] = cr_er * de.ob_err_dict['c_roo']
    de.ob_err_dict['c_woo'] = cw_er * de.ob_err_dict['c_woo']

    # west data
    dw = dc.DalecData(2015, 2016, 'nee_day_west, nee_night_west, c_roo_west, c_woo_west, clma, lai_west',
                      nc_file='../../alice_holt_data/ah_data_daily_test_nee.nc', scale_nee=1)
    dw.B = b
    # obs err scaling
    dw.ob_err_dict['clma'] = clma_er * dw.ob_err_dict['clma']
    dw.ob_err_dict['lai'] = lai_er * dw.ob_err_dict['lai']
    dw.ob_err_dict['nee_day'] = need_er * dw.ob_err_dict['nee_day']
    dw.ob_err_dict['nee_night'] = neen_er * dw.ob_err_dict['nee_night']
    dw.ob_err_dict['c_roo'] = cr_er * dw.ob_err_dict['c_roo']
    dw.ob_err_dict['c_woo'] = cw_er * dw.ob_err_dict['c_woo']
    # setup model
    me = mc.DalecModel(de)
    mw = mc.DalecModel(dw)
    # run DA scheme
    xa_e = me.find_min_tnc_cvt(xb, f_name+'east_assim')
    xa_w = mw.find_min_tnc_cvt(xb, f_name+'west_assim')
    #save plots
    save_plots(f_name, xb, xa_e[1], xa_w[1], de, dw)
    return 'done'

# ------------------------------------------------------------------------------
# Exp runs
# ------------------------------------------------------------------------------

def new_b_run(xb, f_name):
    east_west_joint_run_full(xb, f_name, clma_er=1, lai_er=1, need_er=1, neen_er=1, cr_er=1, cw_er=1)
    east_west_joint_run_full(xb, f_name, clma_er=1, lai_er=1, need_er=2.5, neen_er=1, cr_er=2, cw_er=2)
    east_west_joint_run_full(xb, f_name, clma_er=0.33, lai_er=0.33, need_er=1.5, neen_er=0.75, cr_er=2, cw_er=2)
    east_west_joint_run_full(xb, f_name, clma_er=0.33, lai_er=0.33, need_er=1, neen_er=0.5, cr_er=1, cw_er=1)
    east_west_joint_run_full(xb, f_name, clma_er=1.5, lai_er=1.5, need_er=3, neen_er=2, cr_er=0.5, cw_er=0.5)
    east_west_joint_run_full(xb, f_name, clma_er=0.5, lai_er=0.5, need_er=1.5, neen_er=0.5, cr_er=3, cw_er=3)
    return 'done'


def new_b_run2(xb, f_name):
    east_west_joint_run_full2(xb, f_name, clma_er=1, lai_er=1, need_er=1, neen_er=1, cr_er=1, cw_er=1)
    east_west_joint_run_full2(xb, f_name, clma_er=1, lai_er=1, need_er=2.5, neen_er=1, cr_er=2, cw_er=2)
    east_west_joint_run_full2(xb, f_name, clma_er=0.33, lai_er=0.33, need_er=1.5, neen_er=0.75, cr_er=2, cw_er=2)
    east_west_joint_run_full2(xb, f_name, clma_er=0.33, lai_er=0.33, need_er=1, neen_er=0.5, cr_er=1, cw_er=1)
    east_west_joint_run_full2(xb, f_name, clma_er=1.5, lai_er=1.5, need_er=3, neen_er=2, cr_er=0.5, cw_er=0.5)
    east_west_joint_run_full2(xb, f_name, clma_er=0.5, lai_er=0.5, need_er=1.5, neen_er=0.5, cr_er=3, cw_er=3)
    return 'done'


def new_b_run_scale_nee_err(xb, f_name):
    east_west_joint_run_full2(xb, f_name, clma_er=1, lai_er=1, need_er=1, neen_er=1, cr_er=1, cw_er=1)
    east_west_joint_run_full2(xb, f_name, clma_er=1, lai_er=1, need_er=2.5, neen_er=1, cr_er=2, cw_er=2)
    east_west_joint_run_full2(xb, f_name, clma_er=0.33, lai_er=0.33, need_er=1.5, neen_er=0.75, cr_er=2, cw_er=2)
    east_west_joint_run_full2(xb, f_name, clma_er=0.33, lai_er=0.33, need_er=1, neen_er=0.5, cr_er=1, cw_er=1)
    east_west_joint_run_full2(xb, f_name, clma_er=1.5, lai_er=1.5, need_er=3, neen_er=2, cr_er=0.5, cw_er=0.5)
    east_west_joint_run_full2(xb, f_name, clma_er=0.5, lai_er=0.5, need_er=1.5, neen_er=0.5, cr_er=3, cw_er=3)
    return 'done'

# ------------------------------------------------------------------------------
# East West run
# ------------------------------------------------------------------------------


def east_west_run(f_name, ob_list, east_west):
    ob_str = ''
    for ob in ob_list:
        if ob == 'clma':
            ob_str += ob+','
        else:
            ob_str += ob+'_'+east_west+','
    d = dc.DalecData(2015, 2016, ob_str)
    d.B = d.make_b(d.edinburgh_std)
    m = mc.DalecModel(d)
    assim_results, xa = m.find_min_tnc_cvt(d.edinburgh_mean, f_name+'_assim_res')
    # Plot 4dvar time series
    ax, fig = p.plot_4dvar('nee', d, xb=d.edinburgh_mean, xa=xa)
    fig.savefig(f_name+'_nee.png', bbox_inches='tight')
    ax, fig = p.plot_4dvar('nee_day', d, xb=d.edinburgh_mean, xa=xa)
    fig.savefig(f_name+'_need.png', bbox_inches='tight')
    ax, fig = p.plot_4dvar('nee_night', d, xb=d.edinburgh_mean, xa=xa)
    fig.savefig(f_name+'_neen.png', bbox_inches='tight')
    ax, fig = p.plot_4dvar('lai', d, xb=d.edinburgh_mean, xa=xa)
    fig.savefig(f_name+'_lai.png', bbox_inches='tight')
    ax, fig = p.plot_4dvar('c_woo', d, xb=d.edinburgh_mean, xa=xa)
    fig.savefig(f_name+'_cwoo.png', bbox_inches='tight')
    ax, fig = p.plot_4dvar('c_roo', d, xb=d.edinburgh_mean, xa=xa)
    fig.savefig(f_name+'_croo.png', bbox_inches='tight')

    # Plot scatter plots of obs
    ax, fig = p.plot_scatter('nee_day', xa, d, len(d.I), 'a')
    fig.savefig(f_name+'_need_scat.png', bbox_inches='tight')
    ax, fig = p.plot_scatter('nee_night', xa, d, len(d.I), 'a')
    fig.savefig(f_name+'_neen_scat.png', bbox_inches='tight')

    # Plot error in analysis and background
    ax, fig = p.plot_a_inc(d.edinburgh_mean, xa, east_west)
    fig.savefig(f_name+'_xa_inc.png', bbox_inches='tight')
    return 'all experimented'


# ------------------------------------------------------------------------------
# East West run new B
# ------------------------------------------------------------------------------


def east_west_run_b(f_name, east_west, net_file="None"):
    if east_west == 'east':
        obs = 'nee_day_east, nee_night_east, clma, lai_east, c_woo_east, c_roo_east'
    elif east_west == 'west':
        obs = 'nee_day_west, nee_night_west, clma, lai_west, c_woo_west, c_roo_west'
    if net_file != "None":
        d = dc.DalecData(2015, 2016, obs, nc_file=net_file)
    else:
        d = dc.DalecData(2015, 2016, obs)
    m = mc.DalecModel(d)
    assim_results, xa = m.find_min_tnc_cvt(d.edinburgh_mean, f_name+'_assim_res')
    # Plot 4dvar time series
    ax, fig = p.plot_4dvar('nee', d, xb=d.edinburgh_mean, xa=xa)
    fig.savefig(f_name+'_nee.png', bbox_inches='tight')
    ax, fig = p.plot_4dvar('nee_day', d, xb=d.edinburgh_mean, xa=xa)
    fig.savefig(f_name+'_need.png', bbox_inches='tight')
    ax, fig = p.plot_4dvar('nee_night', d, xb=d.edinburgh_mean, xa=xa)
    fig.savefig(f_name+'_neen.png', bbox_inches='tight')
    ax, fig = p.plot_4dvar('lai', d, xb=d.edinburgh_mean, xa=xa)
    fig.savefig(f_name+'_lai.png', bbox_inches='tight')
    ax, fig = p.plot_4dvar('c_woo', d, xb=d.edinburgh_mean, xa=xa)
    fig.savefig(f_name+'_cwoo.png', bbox_inches='tight')
    ax, fig = p.plot_4dvar('c_roo', d, xb=d.edinburgh_mean, xa=xa)
    fig.savefig(f_name+'_croo.png', bbox_inches='tight')

    # Plot scatter plots of obs
    ax, fig = p.plot_scatter('nee_day', xa, d, len(d.I), 'a')
    fig.savefig(f_name+'_need_scat.png', bbox_inches='tight')
    ax, fig = p.plot_scatter('nee_night', xa, d, len(d.I), 'a')
    fig.savefig(f_name+'_neen_scat.png', bbox_inches='tight')

    # Plot error in analysis and background
    ax, fig = p.plot_a_inc(d.edinburgh_mean, xa, east_west)
    fig.savefig(f_name+'_xa_inc.png', bbox_inches='tight')
    return 'all experimented'


# ------------------------------------------------------------------------------
# East West run new B
# ------------------------------------------------------------------------------
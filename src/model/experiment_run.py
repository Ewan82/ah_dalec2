import data_class as dc
import mod_class as mc
import plot as p
import pickle


def save_plots(f_name, xa_east, xa_west, d_e, d_w):
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
    ax, fig = p.plot_obs_east_west('cl', xa_east, xa_west, d_e, d_w)
    fig.savefig(f_name+'_clit.png', bbox_inches='tight')
    ax, fig = p.plot_obs_east_west('cf', xa_east, xa_west, d_e, d_w)
    fig.savefig(f_name+'_cfol.png', bbox_inches='tight')

    # Plot error in analysis and background
    ax, fig = p.plot_inc_east_west(d_e.edinburgh_mean, xa_east, xa_west)
    fig.savefig(f_name+'_xa_inc.png', bbox_inches='tight')
    return 'done'


def east_west_joint_run(xb, f_name):
    de = dc.DalecData(2015, 2016, 'nee_day_east, nee_night_east, c_roo_east, c_woo_east, clma, lai_east')
    de.B = pickle.load(open('b_edc.p', 'r'))
    dw = dc.DalecData(2015, 2016, 'nee_day_west, nee_night_west, c_roo_west, c_woo_west, clma, lai_west')
    dw.B = pickle.load(open('b_edc.p', 'r'))
    me = mc.DalecModel(de)
    mw = mc.DalecModel(dw)
    xa_e = me.find_min_tnc_cvt(xb)
    xa_w = mw.find_min_tnc_cvt(xb)
    save_plots(f_name, xa_e[1], xa_w[1], de, dw)
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

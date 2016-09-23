import data_class as dc
import mod_class as mc
import plot as p
import pickle


# ------------------------------------------------------------------------------
# East West run
# ------------------------------------------------------------------------------


def east_west_run(f_name, easy_west):
    if easy_west == 'east':
        obs = 'nee_day_east, nee_night_east, clma, lai_east, c_woo_east, c_roo_east'
    elif easy_west == 'west':
        obs = 'nee_day_west, nee_night_west, clma, lai_west, c_woo_west, c_roo_west'
    d = dc.DalecDataTwin(2015, 2016, obs)
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
    ax, fig = p.plot_a_inc(d.edinburgh_mean, xa)
    fig.savefig(f_name+'_xa_inc.png', bbox_inches='tight')
    return 'all experimented'


# ------------------------------------------------------------------------------
# East West run new B
# ------------------------------------------------------------------------------


def east_west_run_b(f_name, easy_west, net_file="None"):
    if easy_west == 'east':
        obs = 'nee_day_east, nee_night_east, clma, lai_east, c_woo_east, c_roo_east'
    elif easy_west == 'west':
        obs = 'nee_day_west, nee_night_west, clma, lai_west, c_woo_west, c_roo_west'
    if net_file != "None":
        d = dc.DalecDataTwin(2015, 2016, obs, nc_file=net_file)
    else:
        d = dc.DalecDataTwin(2015, 2016, obs)
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
    ax, fig = p.plot_a_inc(d.edinburgh_mean, xa)
    fig.savefig(f_name+'_xa_inc.png', bbox_inches='tight')
    return 'all experimented'

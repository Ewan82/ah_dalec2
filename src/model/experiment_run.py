import data_class as dc
import mod_class as mc
import plot as p
import pickle


def day_night_twin_run(start, end, obs, f_name, obs_loc):
    d = dc.DalecDataTwin(start, end, obs)
    pik_obs = pickle.load(open(obs_loc, 'r'))
    d.ob_dict = pik_obs['obs']
    d.ob_err_dict = pik_obs['obs_err']
    m = mc.DalecModel(d)
    assim_results, xa = m.find_min_tnc_cvt(d.xb, f_name+'assim_res')
    d2 = dc.DalecDataTwin(start, 2013, obs)
    # Plot 4dvar time series
    ax, fig = p.plot_4dvar_twin('nee', d2, xa=xa, obdict_a=d.ob_dict)
    fig.savefig(f_name+'_nee.pdf', bbox_inches='tight')
    ax, fig = p.plot_4dvar_twin('nee_day', d2, xa=xa, obdict_a=d.ob_dict)
    fig.savefig(f_name+'_need.pdf', bbox_inches='tight')
    ax, fig = p.plot_4dvar_twin('nee_night', d2, xa=xa, obdict_a=d.ob_dict)
    fig.savefig(f_name+'_neen.pdf', bbox_inches='tight')
    ax, fig = p.plot_4dvar_twin('lai', d2, xa=xa, obdict_a=d.ob_dict)
    fig.savefig(f_name+'_lai.pdf', bbox_inches='tight')

    # Plot scatter plots of obs
    ax, fig = p.plot_scatter_twin('nee', xa, d2, len(d.I), 'f')
    fig.savefig(f_name+'_nee_scat.pdf', bbox_inches='tight')
    ax, fig = p.plot_scatter_twin('nee_day', xa, d2, len(d.I), 'f')
    fig.savefig(f_name+'_need_scat.pdf', bbox_inches='tight')
    ax, fig = p.plot_scatter_twin('nee_night', xa, d2, len(d.I), 'f')
    fig.savefig(f_name+'_neen_scat.pdf', bbox_inches='tight')
    ax, fig = p.plot_scatter_twin('lai', xa, d2, len(d.I), 'f')
    fig.savefig(f_name+'_lai_scat.pdf', bbox_inches='tight')

    # Plot error in analysis and background
    ax, fig = p.plottwinerr(d.x_truth, d.xb, xa)
    fig.savefig(f_name+'_twin_err.pdf', bbox_inches='tight')
    return 'all experimented'

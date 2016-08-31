"""Plotting functions related to dalecv2.
"""
import numpy as np
import matplotlib.pyplot as plt
import model2 as mc

import seaborn as sns

# ------------------------------------------------------------------------------
# Plot observation time series
# ------------------------------------------------------------------------------

def plotgpp(cf, dC, start, fin):
    """Plots gpp using acm equations given a cf value, a dataClass and a start
    and finish point. NOTE cf is treated as constant in this plot
    (unrealistic).
    """
    xlist = np.arange(start, fin, 1)
    gpp = np.ones(fin - start)*-9999.
    for x in xrange(start, fin):
        gpp[x-start] = mc.acm(cf, dC.p17, dC.p11, dC, x)
    plt.plot(xlist, gpp)
    plt.show()


def plotphi(onoff, pvals, dC, start, fin):
    """Plots phi using phi equations given a string "fall" or "onset", a
    dataClass and a start and finish point. Nice check to see dynamics.
    """
    xlist = np.arange(start, fin, 1)
    phi = np.ones(fin - start)*-9999.
    for x in xrange(start, fin):
        if onoff == 'onset':
            phi[x-start] = m.phi_onset(pvals[11], pvals[13], dC, x)
        elif onoff == 'fall':
            phi[x-start] = m.phi_fall(pvals[14], pvals[15], pvals[4], dC, x)
    plt.plot(xlist, phi)
    plt.show()


def plotphi2(pvals, dC, start, fin):
    """Plots phi using phi equations given a string "fall" or "onset", a
    dataClass and a start and finish point. Nice check to see dynamics.
    """
    sns.set_context(rc={'lines.linewidth':.8, 'lines.markersize':6})
    xlist = np.arange(start, fin, 1)
    phion = np.ones(fin - start)*-9999.
    phioff = np.ones(fin - start)*-9999.
    sns.set_context(rc={'lines.linewidth':.8, 'lines.markersize':6})
    fig, ax = plt.subplots( nrows=1, ncols=1,)
    for x in xrange(start, fin):
        phion[x-start] = m.phi_onset(pvals[11], pvals[13], dC, x)
        phioff[x-start] = m.phi_fall(pvals[14], pvals[15], pvals[4], dC, x)
    ax.plot(xlist, phion)
    ax.plot(xlist, phioff)
    ax.set_xlabel('Day of year')
    ax.set_ylabel('Rate of C allocation')
    ax.set_xlim(0, 365)
    return ax, fig


def plot_drive_dat(dat, dC):
    """Plots a specified observation using obs eqn in obs module. Takes an
    observation string, a dataClass (dC) and a start and finish point.
    """
    sns.set_context(rc={'lines.linewidth': 0.8, 'lines.markersize': 6})
    fig, ax = plt.subplots(nrows=1, ncols=1)
    drive_dict = {'t_mean': dC.t_mean, 't_max': dC.t_max, 't_min': dC.t_min,
                  'I': dC.I, 't_day': dC.t_day, 't_night': dC.t_night}
    palette = sns.color_palette("colorblind", 11)

    ax.plot(dC.dates, drive_dict[dat], color=palette[0])
    ax.set_xlabel('Year')
    ax.set_ylabel(dat)
    plt.gcf().autofmt_xdate()
    return ax, fig


def plot_ob_dict(ob, dC):
    """Plots a specified observation using obs eqn in obs module. Takes an
    observation string, a dataClass (dC) and a start and finish point.
    """
    sns.set_context(rc={'lines.linewidth': 0.8, 'lines.markersize': 6})
    fig, ax = plt.subplots(nrows=1, ncols=1)

    palette = sns.color_palette("colorblind", 11)

    ax.plot(dC.dates, dC.ob_dict[ob], 'o', color=palette[0])
    ax.set_xlabel('Year')
    ax.set_ylabel(ob)
    plt.gcf().autofmt_xdate()
    return ax, fig


def plot_obs(ob, pvals, dC):
    """Plots a specified observation using obs eqn in obs module. Takes an
    observation string, a dataClass (dC) and a start and finish point.
    """
    sns.set_context(rc={'lines.linewidth':.8, 'lines.markersize':6})
    fig, ax = plt.subplots(nrows=1, ncols=1)
    m = mc.DalecModel(dC)
    mod_lst = m.mod_list(pvals)
    obs_lst = m.oblist(ob, mod_lst)

    palette = sns.color_palette("colorblind", 11)

    ax.plot(dC.dates, obs_lst, color=palette[0])
    ax.set_xlabel('Year')
    ax.set_ylabel(ob)
    plt.gcf().autofmt_xdate()
    return ax, fig


def plot_4dvar(ob, dC, xb=None, xa=None, erbars=1, awindl=None, obdict_a=None):
    """Plots a model predicted observation value for two initial states (xb,xa)
    and also the actual observations taken of the physical quantity. Takes a ob
    string, two initial states (xb,xa), a dataClass and a start and finish
    time step.
    """
    sns.set_context(rc={'lines.linewidth':.8, 'lines.markersize':6})
    fig, ax = plt.subplots(nrows=1, ncols=1)
    m = mc.DalecModel(dC)
    palette = sns.color_palette("colorblind", 11)
    if xb != None:
        mod_lst = m.mod_list(xb)
        obs_lst = m.oblist(ob, mod_lst)
        ax.plot(dC.dates, obs_lst, color=palette[0])
    if xa != None:
        mod_lst = m.mod_list(xa)
        obs_lst = m.oblist(ob, mod_lst)
        ax.plot(dC.dates, obs_lst, color=palette[1])

    ob_dict = dC.ob_dict
    ob_err_dict = dC.ob_err_dict
    if ob in ob_dict.keys():
        if erbars == True:
            ax.errorbar(dC.dates, ob_dict[ob], yerr=ob_err_dict[ob],
                         fmt='o', label=ob+'_o', color=palette[2], alpha=0.7)
        else:
            ax.plot(dC.dates, ob_dict[ob], 'o', label=ob+'_o', color=palette[2])
    if obdict_a != None:
        ax.plt.plot(dC.dates, obdict_a[ob], 'o')

    if awindl != None:
        ax.axvline(x=dC.dates[awindl], color='k', ls='dashed')

    ax.set_xlabel('Year')
    ax.set_ylabel(ob)
    plt.gcf().autofmt_xdate()

    return ax, fig


def plotscatterobs(ob, pvals, dC, awindl, bfa='a'):
    """Plots scatter plot of obs vs model predicted values. Takes an initial
    parameter set, a dataClass (must have only desired ob for comparison
    specified in dC), assimilation window length and whether a comparison of
    background 'b', forecast 'f' or analysis 'a' is desired.
    """
    sns.set_context('poster', font_scale=1.5, rc={'lines.linewidth': 1., 'lines.markersize': 6.})
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 10))
    sns.set_style('ticks')
    palette = sns.color_palette("colorblind", 11)
    m = mc.DalecModel(dC)
    mod_lst = m.mod_list(pvals)
    obs_lst = m.oblist(ob, mod_lst)
    y, yerr, ytimestep = m.obscost()
    hx = m.hxcost(mod_lst)
    y_obs = dC.ob_dict[ob]
    plt_ob_lst = (y_obs/y_obs)*obs_lst
    one_one = np.arange(int(min(min(y_obs), min(plt_ob_lst)))-1, int(max(max(y_obs), max(plt_ob_lst)))+1)
    plt.plot(one_one, one_one, color=palette[0])
    if bfa == 'b' or bfa == 'a':
        ax.plot(y_obs[0:awindl], plt_ob_lst[0:awindl], 'o', color=palette[1])
        error = np.sqrt(np.nansum((y_obs[0:awindl] - plt_ob_lst[0:awindl])**2)/len(y_obs[0:awindl]))
        yhx = np.nanmean(y_obs[0:awindl] - plt_ob_lst[0:awindl])
    elif bfa == 'f':
        ax.plot(y_obs[awindl:], plt_ob_lst[awindl:], 'o', color=palette[1])
        error = np.sqrt(np.nansum((y_obs[awindl:] - plt_ob_lst[awindl:])**2)/len(y_obs[awindl:]))
        yhx = np.nanmean(y_obs[awindl:] - plt_ob_lst[awindl:])
    else:
        raise Exception('Please check function input for bfa variable')
    plt.xlabel(ob.upper()+r' observations (g C m$^{-2}$ day$^{-1}$)')
    plt.ylabel(ob.upper()+' model (g C m$^{-2}$ day$^{-1}$)')
    #plt.title(bfa+'_error=%f, mean(y-hx)=%f' %(error,yhx))
    print bfa+'_error=%f, mean(y-hx)=%f' %(error, yhx)
    #plt.xlim((-20, 15))
    #plt.ylim((-20, 15))
    return ax, fig





import numpy as np
import matplotlib.mlab as mlab
import csv

import mod_class as mc
import data_class as ahd


def write_dict2csv(dic, f_name):
    keys = sorted(dic.keys())
    with open(f_name, "wb") as outfile:
        writer = csv.writer(outfile, delimiter="\t")
        writer.writerow(keys)
        writer.writerows(zip(*[dic[key] for key in keys]))


def run_assimilation(yr):
    d = ahd.DalecData(yr, yr+1, 'nee, nee_day, nee_night', nc_file='../../alice_holt_data/ah_data_daily_test_nee3.nc')
    m = mc.DalecModel(d)
    find_min, xa = m.find_min_tnc_cvt(d.xb_ew_lai_hi)
    result_dic = {}
    result_dic['time'] = [t.strftime('%m/%d/%Y') for t in d.dates]
    result_dic['incident_radiation'] = d.I.tolist()
    result_dic['t_mean'] = d.t_mean.tolist()
    result_dic['doy'] = d.D.tolist()
    mod_lst_xb = m.mod_list(d.xb_ew_lai_hi)
    result_dic['nee_day_xb'] = m.oblist('nee_day', mod_lst_xb).tolist()
    result_dic['nee_xb'] = m.oblist('nee', mod_lst_xb).tolist()
    mod_lst_xa = m.mod_list(xa)
    result_dic['nee_day_xa'] = m.oblist('nee_day', mod_lst_xa).tolist()
    result_dic['nee_xa'] = m.oblist('nee', mod_lst_xa).tolist()
    return result_dic


def multiple_assimilations():
    for yr in np.arange(1999, 2016, 1):
        res_dic = run_assimilation(yr)
        write_dict2csv(res_dic, 'assimilation_results_'+str(yr)+'.csv')
    return 'done'
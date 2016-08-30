import numpy as np
import netCDF4 as nC
import matplotlib.mlab as mlab
import datetime as dt


def open_csv(filename, missing_val='N/A'):
    """Opens a csv file into a recorded array.
    """
    return mlab.csv2rec(filename, missing=missing_val)


def open_netcdf(filename):
    """Opens a netCDF file
    """
    return nC.Dataset(filename, 'a')


def ah_str2date(date_str):
    """Converts string into datetime object for alice holt spreadsheet.
    """
    return dt.datetime.strptime(date_str, "%d/%m/%Y %H:%M")


def add_data2nc(nc_file, csv_file, data_title, nc_title, date_col='date_combined'):
    """ Adds data to a netCDF file
    :param nc_file: netCDF file location str
    :param csv_file: csv file location str
    :param data_title: title column for data to add as str
    :param nc_title: title of nc variable to add it to as str
    :param date_col: title of date column as str
    :return: nothing
    """
    nc_ob = open_netcdf(nc_file)
    var = nc_ob.variables[nc_title]
    times = nc_ob.variables['time']
    dat_arr = open_csv(csv_file)
    for x in xrange(len(dat_arr[date_col])):
        try:
            idx = nC.date2index(dat_arr[date_col][x], times)
        except ValueError:
            print x
        var[idx, 0, 0] = dat_arr[data_title][x]
    nc_ob.close()
    return 'data updated!'


def nc_create_var(nc_file, var_name, dims, dtype='f8'):
    """ Adds new variable to netCDF data file
    :param nc_file: netcdf data file location as str
    :param var_name: name for new variable as str
    :param dims: dimensions of new variable, tuple containing dimension strings
    :param dtype: data type of new variable
    :return:
    """
    nc_ob = open_netcdf(nc_file)
    nc_ob.createVariable(var_name, dtype, dims)
    nc_ob.close()
    return 'variable added!'


def find_date(nc_time_int, nc_time_ob):
    """ Returns the date as a datetime given any datetime object
    :param nc_time_int: date as an integer corresponding to units in nc file
    :param nc_time_ob: netcdf time object
    :return: date as a datetime object
    """
    day = nC.num2date(nc_time_int, nc_time_int.units).date()
    return dt.datetime.combine(day, dt.datetime.min.time())


def create_date_list(start_date, end_date, del_t='day'):
    """ Creates a list of daily or yearly datetime objects
    :param start_date: start date for list as datetime
    :param end_date: end date for list as datetime
    :return: datetime list
    """
    times = []
    if del_t == 'day':
        delta = dt.timedelta(hours=24)
    elif del_t == 'year':
        delta = dt.timedelta(years=1)
    date = start_date
    while date <= end_date:
        times.append(date)
        date = date + delta
    return times


def nc_doy(doy, times):
    """ Finds day of year for netcdf time dimension
    :param doy: day of year netcdf variable
    :param times: half hourly times as netcdf variable
    :return:
    """
    for t in times[:]:
        day_time = nC.num2date(t, times.units)
        idx = nC.date2index(day_time, times)
        tt = day_time.timetuple()
        doy[idx] = tt.tm_yday
    return 'yay'


def nc_mean_daily_temp(half_hourly_temps, daily_mean_temps, times, time_lst):
    """ Finds mean daily temperatures from half hourly temperature data
    :param half_hourly_temps: half hourly temperatures as netcdf variable
    :param daily_mean_temps: empty daily mean temperature netcdf variable
    :param times: half hourly times as netcdf variable
    :param time_lst: list of daily datetime objects
    :return:
    """
    for t in time_lst:
        idx = nC.date2index(t, times)
        mean_daily_temp = np.mean(half_hourly_temps[idx:idx+48, 0, 0])
        daily_mean_temps[idx, 0, 0] = mean_daily_temp
    return 'yay'


def nc_max_daily_temp(half_hourly_temps, daily_max_temps, times, time_lst):
    for t in time_lst:
        idx = nC.date2index(t, times)
        max_daily_temp = np.max(half_hourly_temps[idx:idx+48, 0, 0])
        daily_max_temps[idx, 0, 0] = max_daily_temp
    return 'yay'


def nc_min_daily_temp(half_hourly_temps, daily_min_temps, times, time_lst):
    for t in time_lst:
        idx = nC.date2index(t, times)
        min_daily_temp = np.min(half_hourly_temps[idx:idx+48, 0, 0])
        daily_min_temps[idx, 0, 0] = min_daily_temp
    return 'yay'


def nc_total_daily_rg(half_hourly_rg, daily_rg, times, time_lst):
    """ Finds total daily global radiation from half hourly global radiation data
    :param half_hourly_rg: half hourly global radiation as netcdf variable (W m-2)
    :param daily_rg: empty total daily global radiation netcdf variable (M J m-2 day-1)
    :param times: half hourly times as netcdf variable
    :param time_lst: list of daily datetime objects
    :return:
    """
    for t in time_lst:
        idx = nC.date2index(t, times)
        total_daily_rg = 30*60*1e-6*np.sum(half_hourly_rg[idx:idx+48, 0, 0])  # Convert W m-2 to M J m-2 day-1
        daily_rg[idx, 0, 0] = total_daily_rg
    return 'yay'


def nc_day_len(is_day, day_len, times, time_lst):
    """ Finds total daily global radiation from half hourly global radiation data
    :param is_day: half hourly 'is day' netcdf variable with values of 1 day or 0 night
    :param day_len: empty day length netcdf variable (hours)
    :param times: half hourly times as netcdf variable
    :param time_lst: list of daily datetime objects
    :return:
    """
    for t in time_lst:
        idx = nC.date2index(t, times)
        where_day = np.where(is_day[idx:idx+48, 0, 0] == 1)[0]
        day_len[idx, 0, 0] = len(where_day)*0.5
    return 'yay'


def nc_night_len(is_day, night_len, times, time_lst):
    """ Finds total daily global radiation from half hourly global radiation data
    :param is_day: half hourly 'is day' netcdf variable with values of 1 day or 0 night
    :param day_len: empty day length netcdf variable (hours)
    :param times: half hourly times as netcdf variable
    :param time_lst: list of daily datetime objects
    :return:
    """
    for t in time_lst:
        idx = nC.date2index(t, times)
        night_idx = idx + np.where(is_day[idx:idx+48, 0, 0] == 1)[0][-1] + 1
        night_length = 0
        while is_day[night_idx, 0, 0] == 0:
            night_length += 1  # if final night incomplete will return false night length for last day of data
            night_idx += 1
            if night_idx >= len(times):
                night_length = float('NaN')
                break
        night_len[idx, 0, 0] = night_length*0.5
    return 'yay'


def nc_day_mean_temp(is_day, hh_temp, mean_t_day, times, time_lst):
    """ Finds total daily global radiation from half hourly global radiation data
    :param is_day: half hourly 'is day' netcdf variable with values of 1 day or 0 night
    :param hh_temp: half hourly temperatures netcdf variable
    :param mean_t_day: netcdf variable to fill with mean daytime temperatures
    :param times: half hourly times as netcdf variable
    :param time_lst: list of daily datetime objects
    :return:
    """
    for t in time_lst:
        idx = nC.date2index(t, times)
        where_day = np.where(is_day[idx:idx+48, 0, 0] == 1)[0]
        mean_t_day[idx, 0, 0] = np.mean(hh_temp[idx+where_day[0]:idx+where_day[-1]])
    return 'yay'


def nc_night_mean_temp(is_day, hh_temp, mean_t_night, times, time_lst):
    """ Finds total daily global radiation from half hourly global radiation data
    :param is_day: half hourly 'is day' netcdf variable with values of 1 day or 0 night
    :param hh_temp: half hourly temperatures netcdf variable
    :param mean_t_day: netcdf variable to fill with mean nighttime temperatures
    :param times: half hourly times as netcdf variable
    :param time_lst: list of daily datetime objects
    :return:
    """
    for t in time_lst:
        idx = nC.date2index(t, times)
        night_idx1 = idx + np.where(is_day[idx:idx+48, 0, 0] == 1)[0][-1] + 1
        night_idx2 = night_idx1
        while is_day[night_idx2, 0, 0] == 0:
            night_idx2 += 1
            if night_idx2 >= len(times):
                night_idx2 = float('NaN')
                break
        if np.isnan(night_idx2) == True:
            mean_t_night[idx, 0, 0] = float('NaN')
            break
        mean_t_night[idx, 0, 0] = np.mean(hh_temp[night_idx1:night_idx2, 0, 0])
    return 'yay'


def clip_co2_flux(co2_flux, clipped_co2_flux, above_clip=50., below_clip=-60.):
    """ Clips co2_flux data and saves in clipped_co2_flux netcdf variable
    :param co2_flux: co2 flux measurements from tower as netcdf variable
    :param clipped_co2_flux: netcdf variable to save clipped measurements in
    :return:
    """
    clipped_co2_flux[:, 0, 0] = co2_flux[:, 0, 0]
    clip_nee = clipped_co2_flux[:, 0, 0]
    # Value of +/- 70 u mol m-2 s-1 chosen for upper and lower limit to clip.
    clip_nee[clip_nee > above_clip] = float('NaN')
    clip_nee[clip_nee < below_clip] = float('NaN')
    clipped_co2_flux[:, 0, 0] = clip_nee
    return 'yay'


def clip_co2_flux_mean(co2_flux, clipped_co2_flux, idx1, idx2):
    """ Clips co2_flux data and saves in clipped_co2_flux netcdf variable
    :param co2_flux: co2 flux measurements from tower as netcdf variable
    :param clipped_co2_flux: netcdf variable to save clipped measurements in
    :return:
    """
    clipped_co2_flux[idx1:idx2, 0, 0] = co2_flux[idx1:idx2, 0, 0]
    clip_nee = clipped_co2_flux[idx1:idx2, 0, 0]
    # Value of +/- 70 u mol m-2 s-1 chosen for upper and lower limit to clip.
    clip_mean_pos = np.nanmean(clip_nee[clip_nee > 0])
    clip_std_pos = np.nanstd(clip_nee[clip_nee > 0])
    clip_mean_neg = np.nanmean(clip_nee[clip_nee < 0])
    clip_std_neg = np.nanstd(clip_nee[clip_nee < 0])
    clip_nee[clip_nee > clip_mean_pos+3*clip_std_pos] = float('NaN')
    clip_nee[clip_nee < clip_mean_neg-3*clip_std_neg] = float('NaN')
    clipped_co2_flux[idx1:idx2, 0, 0] = clip_nee
    return 'yay'


def find_indices_year(times, year):
    """ Returns the first and last time index for a given year
    :param times: time netcdf variable
    :param year: year to find indices for as an integer
    :return: first index, final index
    """
    year_entries = [x for x in times[:] if nC.num2date(x, times.units).year == year]
    idx1 = np.where(times[:] == year_entries[0])[0][0]
    idx2 = np.where(times[:] == year_entries[-1])[0][0]
    return idx1, idx2


def find_indices_month(time_arr, month, time_units):
    """ Returns the first and last time index for a given month
    :param times: slice from a netcdf time variable
    :param month: month to find indices for as an integer
    :return: first index, final index
    """
    month_entries = [x for x in time_arr if nC.num2date(x, time_units).month == month]
    idx1 = np.where(time_arr == month_entries[0])[0][0]
    idx2 = np.where(time_arr == month_entries[-1])[0][0]
    return idx1, idx2


def make_year_lst(times):
    start_yr = nC.num2date(times[0], times.units).year
    end_yr = nC.num2date(times[-1], times.units).year
    year_lst = np.arange(start_yr, end_yr+1)
    return year_lst


def clip_co2_flux_wrapper(co2_flux, clipped_co2_flux, times):
    """ Clips co2 flux data by year and month using clip_co2_flux_mean fn
    :param co2_flux: co2 flux netcdf variable
    :param clipped_co2_flux: netcdf variable to put clipped flux observations
    :param times: time netcdf variable
    :return:
    """
    yr_lst = make_year_lst(times)
    for yr in yr_lst:
        yr_idx1, yr_idx2 = find_indices_year(times, yr)
        clip_co2_flux_mean(co2_flux, clipped_co2_flux, yr_idx1, yr_idx2)
        times_yr = times[yr_idx1:yr_idx2]
        for month in np.arange(1, 13):
            month_idx1, month_idx2 = find_indices_month(times_yr, month, times.units)
            clip_co2_flux_mean(co2_flux, clipped_co2_flux, yr_idx1+month_idx1, yr_idx1+month_idx2)
    return 'yay'


def quality_control_co2_flux(clipped_co2_flux, qc_co2_flux, nee, idx1, idx2, idx):
    """ Quality controls flux data and processes it to daily or half daily
    :param clipped_co2_flux: half hourly clipped co2 flux as netcdf variable
    :param qc_co2_flux: qc flags as nc variable corresponding to half hourly co2 flux
    :param nee: nc variable to fill with processed data
    :param idx1: start index to being quality control and data averaging
    :param idx: start of day index
    :return:
    """
    fill = 0
    qc_flag = 0

    for x in xrange(idx1, idx2):
        if np.isnan(clipped_co2_flux[x, 0, 0]) == True:
            fill += 1
        elif qc_co2_flux[x, 0, 0] == 2:
            qc_flag += 1
        elif qc_co2_flux[x, 0, 0] == 1:
             qc_flag += 1
        else:
            continue

    if fill > 0.:
        nee[idx, 0, 0] = float('NaN')
    elif qc_flag > 6.:
        nee[idx, 0, 0] = float('NaN')
    else:
        # u mol m-2 s-1 to g C m-2 day-1 (CHECK what units do we want day night in?)
        nee[idx, 0, 0] = 12.011*1e-6 * (idx2-idx1)*30*60 * np.mean(clipped_co2_flux[idx1:idx2, 0, 0])


def process_co2_flux_daily(clipped_co2_flux, qc_co2_flux, daily_nee, times, time_lst):
    """ Produces a daily NEE product
    :param clipped_co2_flux: half hourly clipped co2 flux as netcdf variable
    :param qc_co2_flux: qc flags as nc variable corresponding to half hourly co2 flux
    :param daily_nee: nc variable to fill with processed data
    :param times: half hourly times as netcdf variable
    :param time_lst: list of daily datetime objects
    :return:
    """
    for t in time_lst:
        idx = nC.date2index(t, times)
        quality_control_co2_flux(clipped_co2_flux, qc_co2_flux, daily_nee, idx, idx+48, idx)
    return 'yay'


def process_co2_flux_daytime(clipped_co2_flux, qc_co2_flux, nee_day, is_day, times, time_lst):
    """ Produces a daytime NEE product
    :param clipped_co2_flux: half hourly clipped co2 flux as netcdf variable
    :param qc_co2_flux: qc flags as nc variable corresponding to half hourly co2 flux
    :param nee_day: nc variable to fill with processed data
    :param times: half hourly times as netcdf variable
    :param time_lst: list of daily datetime objects
    :return:
    """
    for t in time_lst:
        idx = nC.date2index(t, times)
        where_day = np.where(is_day[idx:idx+48, 0, 0] == 1)[0]
        quality_control_co2_flux(clipped_co2_flux, qc_co2_flux, nee_day, idx+where_day[0], idx+where_day[-1], idx)
    return 'yay'


def process_co2_flux_nighttime(clipped_co2_flux, qc_co2_flux, nee_night, is_day, times, time_lst):
    """ Produces a nighttime NEE product
    :param clipped_co2_flux: half hourly clipped co2 flux as netcdf variable
    :param qc_co2_flux: qc flags as nc variable corresponding to half hourly co2 flux
    :param nee_day: nc variable to fill with processed data
    :param times: half hourly times as netcdf variable
    :param time_lst: list of daily datetime objects
    :return:
    """
    for t in time_lst:
        idx = nC.date2index(t, times)
        night_idx1 = idx + np.where(is_day[idx:idx+48, 0, 0] == 1)[0][-1] + 1
        night_idx2 = night_idx1
        while is_day[night_idx2, 0, 0] == 0:
            night_idx2 += 1
            if night_idx2 >= len(times):
                night_idx2 = float('NaN')
                break
        if np.isnan(night_idx2) == True:
            nee_night[idx, 0, 0] = float('NaN')
            break
        else:
            quality_control_co2_flux(clipped_co2_flux, qc_co2_flux, nee_night, night_idx1, night_idx2-1, idx)
            if nee_night[idx, 0, 0] < 0:
                nee_night[idx, 0, 0] = float('NaN')
            else:
                continue
    return 'yay'

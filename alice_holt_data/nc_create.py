import numpy as np
import netCDF4 as nC
import matplotlib.mlab as ml
import time as tt
import pandas as pd
import datetime as dt


def create_netcdf_dataset():
    """Creates netcdf dataset for half hourly flux measurements at Alice Holt
    """
    dataset = nC.Dataset('ah_data_half_hourly.nc', 'w', format='NETCDF4_CLASSIC')
    time = dataset.createDimension('time', None)
    lat = dataset.createDimension('lat', 1)
    lon = dataset.createDimension('lon', 1)
    times = dataset.createVariable('time', np.float64, ('time',))
    latitudes = dataset.createVariable('latitude', np.float32, ('lat',))
    latitudes[0] = 51.153526
    longitudes = dataset.createVariable('longitude', np.float32, ('lon',))
    longitudes[0] = -0.85835201

    dataset.description = 'Alice Holt Straits Inclosure half hourly data'
    dataset.history = 'Created ' + tt.ctime(tt.time())
    dataset.source = 'Ewan Pinnington, University of Reading. email: ewan.pinnington@gmail.com'

    latitudes.units = 'degrees north'
    longitudes.units = 'degrees east'
    times.units = 'minutes since 1970-01-01 00:00:00.0'
    times.calendar = 'gregorian'

    air_temp = dataset.createVariable('air_temp', 'f4', ('time', 'lat', 'lon'))
    air_temp.units = 'degC'
    air_temp.description = 'air temperature at 27m'
    air_temp.standard_name = 'air_temperature'
    soil_temp = dataset.createVariable('soil_temp', 'f4', ('time', 'lat', 'lon'))
    soil_temp.units = 'degC'
    soil_temp.description = 'soil temperature at 3cm'
    rg = dataset.createVariable('rg', 'f4', ('time', 'lat', 'lon'))
    rg.units = 'W m-2'
    rg.standard_name = 'surface_downwelling_shortwave_flux_in_air'
    co2_flux = dataset.createVariable('co2_flux', 'f4', ('time', 'lat', 'lon'))
    co2_flux.units = 'umol m-2 s-1'
    co2_flux.standard_name = 'surface_upward_mole_flux_of_carbon_dioxide'
    co2_flux.description = 'unprocessed Alice Holt flux tower record'
    qc_co2_flux = dataset.createVariable('qc_co2_flux', 'i1', ('time', 'lat', 'lon'))
    qc_co2_flux.units = 'none'
    qc_co2_flux.description ='quality control flag for half hourly co2 flux observations, 0 - good, 2 - bad'
    u_star = dataset.createVariable('u_star', 'f4', ('time', 'lat', 'lon'))
    u_star.units = 'm s-1'
    wind_dir = dataset.createVariable('wind_dir', 'f4', ('time', 'lat', 'lon'))
    wind_dir.units = 'degree'
    wind_dir.standard_name = 'wind_from_direction'
    wind_dir.description = 'wind direction degrees from north'
    foot_print = dataset.createVariable('foot_print', 'f4', ('time', 'lat', 'lon'))
    foot_print.units = 'm'
    foot_print.description = 'distance from tower where 90% of co2 flux is measured'
    return dataset


def open_netcdf(filename):
    """Opens a netCDF file
    """
    return nC.Dataset(filename, 'a')


def open_xls_sheet(filename, sheet_name):
    return pd.read_excel(filename, sheet_name)


def add_data2nc(nc_data, pd_df, data_title, nc_title, date_col='date_combined'):
    """ Adds data to a netCDF file
    :param nc_data: netCDF data set object
    :param pd_df: pandas data frame object
    :param data_title: title column for data to add as str
    :param nc_title: title of nc variable to add it to as str
    :param date_col: title of date column as str
    :return: nothing
    """
    var = nc_data.variables[nc_title]
    times = nc_data.variables['time']
    for x in xrange(len(pd_df[date_col])):
        try:
            tm = pd_df[date_col][x]  # datetime for var
            # Round datetime to nearest 10min mark
            discard = dt.timedelta(minutes=tm.minute % 10,
                             seconds=tm.second,
                             microseconds=tm.microsecond)
            tm -= discard
            if discard >= dt.timedelta(minutes=5):
                tm += dt.timedelta(minutes=10)
            # Find datetime index
            idx = nC.date2index(tm, times, select='exact')
        except TypeError:
            print x
            break
        except ValueError:
            print x
            break
        var[idx, 0, 0] = pd_df[data_title][x]


def add_excel_ah_obs(xls_file, nc_file, start_yr=1999, end_yr=2016):
    years = np.arange(start_yr, end_yr)
    nc_data = open_netcdf(nc_file)
    nc_vars = ['air_temp', 'soil_temp', 'rg', 'co2_flux', 'qc_co2_flux', 'u_star', 'wind_dir', 'foot_print']
    for yr in years:
        print yr
        pd_df = open_xls_sheet(xls_file, str(yr))
        for var_title in nc_vars:
            print var_title
            add_data2nc(nc_data, pd_df, var_title, var_title)
    nc_data.close()
    return 'net_cdf file updated!'


import numpy as np
import collections as col
import netCDF4 as Nc


class DalecData():

    def __init__(self, start_date=None, end_date=None, ob_str=None, delta_t=None):

        # 'I.C. for carbon pools gCm-2'   range
        self.clab = 75.0               # (10,1e3)
        self.cf = 0.0                  # (10,1e3)
        self.cr = 135.0                # (10,1e3)
        self.cw = 14313.0              # (3e3,3e4)
        self.cl = 70.                  # (10,1e3)
        self.cs = 18624.0              # (1e3, 1e5)

        # 'Parameters for optimization'                     range
        self.p1 = 1.1e-5  # theta_min, cl to cs decomp      (1e-5 - 1e-2)day^-1
        self.p2 = 0.45  # f_auto, fraction of GPP respired  (0.3 - 0.7)
        self.p3 = 0.01  # f_fol, frac GPP to foliage        (0.01 - 0.5)
        self.p4 = 0.457  # f_roo, frac GPP to fine roots    (0.01 - 0.5)
        self.p5 = 3.  # clspan, leaf lifespan               (1.0001 - 5)
        self.p6 = 4.8e-5  # theta_woo, wood C turnover      (2.5e-5 - 1e-3)day^-1
        self.p7 = 6.72e-3  # theta_roo, root C turnover rate(1e-4 - 1e-2)day^-1
        self.p8 = 0.024  # theta_lit, litter C turnover     (1e-4 - 1e-2)day^-1
        self.p9 = 2.4e-5  # theta_som, SOM C turnover       (1e-7 - 1e-3)day^-1
        self.p10 = 0.0193  # Theta, temp dependence exp fact(0.018 - 0.08)
        self.p11 = 90.  # ceff, canopy efficiency param     (10 - 100)
        self.p12 = 140.  # d_onset, clab release date       (1 - 365) (60,150)
        self.p13 = 0.7  # f_lab, frac GPP to clab           (0.01 - 0.5)
        self.p14 = 27.  # cronset, clab release period      (10 - 100)
        self.p15 = 308.  # d_fall, date of leaf fall        (1 - 365) (242,332)
        self.p16 = 35.  # crfall, leaf fall period          (10 - 100)
        self.p17 = 24.2  # clma, leaf mass per area         (10 - 400)gCm^-2

        self.paramdict = col.OrderedDict([('theta_min', self.p1),
                                          ('f_auto', self.p2), ('f_fol', self.p3),
                                          ('f_roo', self.p4), ('clspan', self.p5),
                                          ('theta_woo', self.p6), ('theta_roo', self.p7),
                                          ('theta_lit', self.p8), ('theta_som', self.p9),
                                          ('Theta', self.p10), ('ceff', self.p11),
                                          ('d_onset', self.p12), ('f_lab', self.p13),
                                          ('cronset', self.p14), ('d_fall', self.p15),
                                          ('crfall', self.p16), ('clma', self.p17),
                                          ('clab', self.clab), ('cf', self.cf),
                                          ('cr', self.cr), ('cw', self.cw), ('cl', self.cl),
                                          ('cs', self.cs)])
        self.paramdict2 = {'theta_min': self.p1, 'f_auto': self.p2, 'f_fol': self.p3,
                            'f_roo': self.p4, 'clspan': self.p5, 'theta_woo': self.p6, 'theta_roo': self.p7,
                            'theta_lit': self.p8, 'theta_som': self.p9, 'Theta': self.p10, 'ceff': self.p11,
                            'd_onset': self.p12, 'f_lab': self.p13, 'cronset': self.p14, 'd_fall': self.p15,
                            'crfall': self.p16, 'clma': self.p17, 'clab': self.clab, 'cf': self.cf,
                            'cr': self.cr, 'cw': self.cw, 'cl': self.cl, 'cs': self.cs}
        self.pvals = np.array(self.paramdict.values())

        self.ahpvals = np.array([9.41e-04, 4.7e-01, 2.8e-01, 2.60e-01, 1.01e+00, 2.6e-04,
                                 2.48e-03, 3.38e-03, 2.6e-06, 1.93e-02, 9.0e+01, 1.4e+02,
                                 4.629e-01, 2.7e+01, 3.08e+02, 3.5e+01, 5.2e+01, 78.,
                                 2., 134., 14257.32, 68.95, 18625.77])

        self.edinburghpvals = np.array([0.000189966332469257, 0.565343476756027,
                                        0.015313852599075, 0.229473358726997, 1.3820788381002,
                                        2.56606744808776e-05, 0.000653099081656547, 0.00635847131570823,
                                        4.32163613374937e-05, 0.0627274280370167, 66.4118798958804,
                                        122.361932206327, 0.372217324163812, 114.092521668926,
                                        308.106881011017, 63.6023224321684, 201.056970845445,
                                        201.27512854457, 98.9874539256948, 443.230119619488,
                                        20293.9092250464, 141.405866537237, 2487.84616355469])

        self.xb = self.pvals

        self.bnds = ((1e-5, 1e-2), (0.3, 0.7), (0.01, 0.5), (0.01, 0.5), (1.0001, 10.),
                     (2.5e-5, 1e-3), (1e-4, 1e-2), (1e-4, 1e-2), (1e-7, 1e-3), (0.018, 0.08),
                     (10, 100), (1, 365), (0.01, 0.5), (10, 100), (1, 365), (10, 100), (10, 400),
                     (10, 1000), (10, 1000), (10, 1000), (100, 1e5), (10, 1000), (100, 2e5))

        self.xa = None

        # Constants for ACM model
        self.acmwilliamsxls = np.array([0.0155, 1.526, 324.1, 0.2017,
                                        1.315, 2.595, 0.037, 0.2268,
                                        0.9576])
        self.acmreflex = np.array([0.0156935, 4.22273, 208.868, 0.0453194,
                                   0.37836, 7.19298, 0.011136, 2.1001,
                                   0.789798])
        self.acm = self.acmreflex  # (currently using params from REFLEX)
        self.phi_d = -2.5  # max. soil leaf water potential difference
        self.R_tot = 1.  # total plant-soil hydrolic resistance
        self.lat = 0.89133965  # latitutde of forest site in radians
                               # lat = 51.153525 deg, lon = -0.858352 deg

        # misc
        self.ca = 390.0  # atmospheric carbon
        self.radconv = 365.25 / np.pi
        self.delta_t = delta_t

        # 'Background standard deviations for carbon pools & B matrix'
        self.sigb_clab = 7.5  # 20%
        self.sigb_cf = 10.0  # 20%
        self.sigb_cw = 1000.  # 20%
        self.sigb_cr = 13.5 # 20%
        self.sigb_cl = 7.0  # 20%
        self.sigb_cs = 1500.  # 20%

        # 'Observation standard deviations for carbon pools and NEE'
        self.sigo_clab = 7.5  # 10%
        self.sigo_cf = 10.0  # 10%
        self.sigo_cw = 1000.  # 10%
        self.sigo_cr = 13.5  # 30%
        self.sigo_cl = 7.0  # 30%
        self.sigo_cs = 1500. # 30%
        self.sigo_nee = 0.71  # gCm-2day-1
        self.sigo_lf = 0.5
        self.sigo_lw = 0.5
        self.sigo_litresp = 0.5
        self.sigo_soilresp = 0.6
        self.sigo_rtot = 0.71
        self.sigo_rh = 0.6

        self.errdict = {'clab': self.sigo_clab, 'cf': self.sigo_cf,
                        'cw': self.sigo_cw, 'cl': self.sigo_cl, 'cr': self.sigo_cr,
                        'cs': self.sigo_cs, 'nee': self.sigo_nee,
                        'lf': self.sigo_lf, 'lw': self.sigo_lw,
                        'litresp': self.sigo_litresp,
                        'soilresp': self.sigo_soilresp,
                        'rtot': self.sigo_rtot,
                        'rh': self.sigo_rh}

    def extract_data(self):
        data_set = Nc.Dataset(self.filename, 'a')

        # 'Daily temperatures degC'
        self.t_mean = self.fluxdata['t_mean']
        self.t_max = self.fluxdata['t_max']
        self.t_min = self.fluxdata['t_min']
        self.t_range = np.array(self.t_max) - np.array(self.t_min)

        # 'Driving Data'
        self.I = self.fluxdata['i']  # incident radiation
        self.ca = 390.0  # atmospheric carbon
        self.D = self.fluxdata['day']  # day of year
        self.year = self.fluxdata['year']  # Year
        self.month = self.fluxdata['month']  # Month
        self.date = self.fluxdata['date']  # Date in month


#class AliceHolt(DalecData):
"""Dalecv2 model class takes a data class and then uses functions to run the
dalecv2 model.
"""
import numpy as np
import algopy


class DalecModel():

    def __init__(self, dataclass, timestep=0, strtrun=0):
        """dataClass and timestep at which to run the dalecv2 model.
        """
        self.dC = dataclass
        self.x = timestep
        if self.dC.k is None:
            self.lenrun = self.dC.lenrun
        else:
            self.lenrun = self.dC.lenrun*self.dC.k
        self.xb = self.dC.xb
        self.modcoston = 1
        #self.modobdict = {'gpp': self.gpp, 'nee': self.nee, 'rt': self.rec,
        #                  'cf': self.cf, 'clab': self.clab, 'cr': self.cr,
        #                  'cw': self.cw, 'cl': self.cl, 'cs': self.cs,
        #                  'lf': self.lf, 'lw': self.lw, 'lai': self.lai,
        #                  'litresp': self.litresp, 'soilresp': self.soilresp,
        #                  'rtot': self.rtot, 'rh': self.rh}
        self.startrun = strtrun
        self.endrun = self.lenrun
        #self.yoblist, self.yerroblist = self.obscost()
        #self.rmatrix = self.rmat(self.yerroblist)
        #self.nume = 100


# ------------------------------------------------------------------------------
# Model functions
# ------------------------------------------------------------------------------

    @staticmethod
    def fitpolynomial(ep, multfac):
        """Polynomial used to find phi_f and phi (offset terms used in
        phi_onset and phi_fall), given an evaluation point for the polynomial
        and a multiplication term.
        """
        cf = [2.359978471e-05, 0.000332730053021, 0.000901865258885,
              -0.005437736864888, -0.020836027517787, 0.126972018064287,
              -0.188459767342504]
        polyval = cf[0]*ep**6 + cf[1]*ep**5 + cf[2]*ep**4 + cf[3]*ep**3 + cf[4]*ep**2 + \
            cf[5]*ep**1 + cf[6]*ep**0
        phi = polyval*multfac
        return phi

    def temp_term(self, Theta):
        """Calculates the temperature exponent factor for carbon pool
        respiration's given a value for Theta parameter.
        """
        temp_term = np.exp(Theta*self.dC.t_mean[self.x])
        return temp_term

    def acm(self, cf, clma, ceff, acm):
        """Aggregated canopy model (ACM) function
        ------------------------------------------
        Takes a foliar carbon (cf) value, leaf mass per area (clma) and canopy
        efficiency (ceff) and returns the estimated value for Gross Primary
        Productivity (gpp) of the forest at that time.
        """
        t_range = 0.5*(self.dC.t_max[self.x] - self.dC.t_min[self.x])
        L = cf / clma
        q = acm[1] - acm[2]
        gc = (abs(self.dC.phi_d))**acm[8] / \
             (t_range + acm[4]*self.dC.R_tot)
        p = ((ceff*L) / gc)*np.exp(acm[6]*self.dC.t_max[self.x])
        ci = 0.5*(self.dC.ca + q - p + np.sqrt((self.dC.ca + q - p)**2
                  - 4*(self.dC.ca*q - p*acm[1])))
        E0 = (acm[5]*L**2) / (L**2 + acm[7])
        delta = -23.4*np.cos((360.*(self.dC.D[self.x] + 10) / 365.) *
                             (np.pi/180.))*(np.pi/180.)
        s = 24*np.arccos((- np.tan(self.dC.lat)*np.tan(delta))) / np.pi
        if s >= 24.:
            s = 24.
        elif s <= 0.:
            s = 0.
        else:
            s = s
        gpp = (E0*self.dC.I[self.x]*gc*(self.dC.ca - ci))*(acm[0]*s +
                                                           acm[3]) / \
              (E0*self.dC.I[self.x] + gc*(self.dC.ca - ci))
        return gpp

    def phi_onset(self, d_onset, cronset):
        """Leaf onset function (controls labile to foliar carbon transfer)
        takes d_onset value, cronset value and returns a value for phi_onset.
        """
        releasecoeff = np.sqrt(2.)*cronset / 2.
        magcoeff = (np.log(1.+1e-3) - np.log(1e-3)) / 2.
        offset = self.fitpolynomial(1+1e-3, releasecoeff)
        phi_onset = (2. / np.sqrt(np.pi))*(magcoeff / releasecoeff) * \
            np.exp(-(np.sin((self.dC.D[self.x] - d_onset + offset) /
                     self.dC.radconv)*(self.dC.radconv / releasecoeff))**2)
        return phi_onset

    def phi_fall(self, d_fall, crfall, clspan):
        """Leaf fall function (controls foliar to litter carbon transfer) takes
        d_fall value, crfall value, clspan value and returns a value for
        phi_fall.
        """
        releasecoeff = np.sqrt(2.)*crfall / 2.
        magcoeff = (np.log(clspan) - np.log(clspan - 1.)) / 2.
        offset = self.fitpolynomial(clspan, releasecoeff)
        phi_fall = (2. / np.sqrt(np.pi))*(magcoeff / releasecoeff) * \
            np.exp(-(np.sin((self.dC.D[self.x] - d_fall + offset) /
                   self.dC.radconv)*self.dC.radconv / releasecoeff)**2)
        return phi_fall

    def dalecv2(self, p, delta_t):
        """DALECV2 carbon balance model
        -------------------------------
        evolves carbon pools to the next time step, taking the 6 carbon pool
        values and 17 parameters at time t and evolving them to time t+1.
        Outputs both the 6 evolved C pool values and the 17 constant parameter
        values.

        phi_on = phi_onset(d_onset, cronset)
        phi_off = phi_fall(d_fall, crfall, clspan)
        gpp = acm(cf, clma, ceff)
        temp = temp_term(Theta)

        clab2 = (1 - phi_on)*clab + (1-f_auto)*(1-f_fol)*f_lab*gpp
        cf2 = (1 - phi_off)*cf + phi_on*clab + (1-f_auto)*f_fol*gpp
        cr2 = (1 - theta_roo)*cr + (1-f_auto)*(1-f_fol)*(1-f_lab)*f_roo*gpp
        cw2 = (1 - theta_woo)*cw + (1-f_auto)*(1-f_fol)*(1-f_lab)*(1-f_roo)*gpp
        cl2 = (1-(theta_lit+theta_min)*temp)*cl + theta_roo*cr + phi_off*cf
        cs2 = (1 - theta_som*temp)*cs + theta_woo*cw + theta_min*temp*cl
        """
        out = algopy.zeros(23, dtype=p)
        # ACM model
        if self.dC.isday[self.x] == 0.0:
            gpp = 0.0
        elif self.dC.isday[self.x] == 1.0:
            gpp = self.acm(p[18], p[16], p[10], self.dC.acm)
        else:
            ValueError('isday value incorrect')
        auto_resp = self.acm(p[18], p[16], p[10], self.dC.acm)*p[1]
        # Labile release and leaf fall factors
        phi_on = self.phi_onset(p[11], p[13])
        phi_off = self.phi_fall(p[14], p[15], p[4])
        # Temperature factor
        temp = self.temp_term(p[9])

        out[17] = p[17] + (-phi_on*p[17] + (1-p[2])*p[12]*gpp - (1-p[2])*p[12]*auto_resp)*delta_t
        out[18] = p[18] + (-phi_off*p[18] + phi_on*p[17] + (gpp - auto_resp)*p[2])*delta_t
        out[19] = p[19] + (-p[6]*p[19] + (1-p[2])*(1-p[12])*p[3]*gpp\
                           - (1-p[2])*(1-p[12])*p[3]*auto_resp)*delta_t
        out[20] = p[20] + (-p[5]*p[20] + (1-p[2])*(1-p[12])*(1-p[3])*gpp\
                           - (1-p[2])*(1-p[12])*(1-p[3])*auto_resp)*delta_t
        out[21] = p[21] + (-(p[7]+p[0])*temp*p[21] + p[6]*p[19] + phi_off*p[18])*delta_t
        out[22] = p[22] + (-p[8]*temp*p[22] + p[5]*p[20] + p[0]*temp*p[21])*delta_t
        out[0:17] = p[0:17]
        return out

    def dalecv2new(self, pdiff, p, delta_t, output):
        """DALECV2 carbon balance model
        -------------------------------
        evolves carbon pools to the next time step, taking the 6 carbon pool
        values and 17 parameters at time t and evolving them to time t+1.
        Outputs both the 6 evolved C pool values and the 17 constant parameter
        values.

        phi_on = phi_onset(d_onset, cronset)
        phi_off = phi_fall(d_fall, crfall, clspan)
        gpp = acm(cf, clma, ceff)
        temp = temp_term(Theta)

        clab2 = (1 - phi_on)*clab + (1-f_auto)*(1-f_fol)*f_lab*gpp
        cf2 = (1 - phi_off)*cf + phi_on*clab + (1-f_auto)*f_fol*gpp
        cr2 = (1 - theta_roo)*cr + (1-f_auto)*(1-f_fol)*(1-f_lab)*f_roo*gpp
        cw2 = (1 - theta_woo)*cw + (1-f_auto)*(1-f_fol)*(1-f_lab)*(1-f_roo)*gpp
        cl2 = (1-(theta_lit+theta_min)*temp)*cl + theta_roo*cr + phi_off*cf
        cs2 = (1 - theta_som*temp)*cs + theta_woo*cw + theta_min*temp*cl
        """
        for x in xrange(len(output)):
            p[output[x]] = pdiff[x]
        out = algopy.zeros(len(output), dtype=pdiff)

        phi_on = phi_onset(d_onset, cronset)
        phi_off = phi_fall(d_fall, crfall, clspan)
        gpp = acm(cf, clma, ceff)
        temp = temp_term(Theta)

        clab = (1 - phi_on)*clab + (1-f_auto)*(1-f_fol)*f_lab*gpp
        cf = (1 - phi_off)*cf + phi_on*clab + (1-f_auto)*f_fol*gpp
        cr = (1 - theta_roo)*cr + (1-f_auto)*(1-f_fol)*(1-f_lab)*f_roo*gpp
        cw = (1 - theta_woo)*cw + (1-f_auto)*(1-f_fol)*(1-f_lab)*(1-f_roo)*gpp
        cl = (1-(theta_lit+theta_min)*temp)*cl + theta_roo*cr + phi_off*cf
        cs = (1 - theta_som*temp)*cs + theta_woo*cw + theta_min*temp*cl

        for x in xrange(len(output)):
            if output[x] in ['clab', 'cf', 'cr', 'cw', 'cl', 'cs']:
                output[x] = output[x]+'2'
            out[x] = p[output[x]]
        return out

    def dalecv2diff(self, p, delta_t):
        """DALECV2 carbon balance model
        -------------------------------
        evolves carbon pools to the next time step, taking the 6 carbon pool
        values and 17 parameters at time t and evolving them to time t+1.
        Ouputs an array of just the 6 evolved C pool values.

        phi_on = phi_onset(d_onset, cronset)
        phi_off = phi_fall(d_fall, crfall, clspan)
        gpp = acm(cf, clma, ceff)
        temp = temp_term(Theta)

        clab2 = (1 - phi_on)*clab + (1-f_auto)*(1-f_fol)*f_lab*gpp
        cf2 = (1 - phi_off)*cf + phi_on*clab + (1-f_auto)*f_fol*gpp
        cr2 = (1 - theta_roo)*cr + (1-f_auto)*(1-f_fol)*(1-f_lab)*f_roo*gpp
        cw2 = (1 - theta_woo)*cw + (1-f_auto)*(1-f_fol)*(1-f_lab)*(1-f_roo)*gpp
        cl2 = (1-(theta_lit+theta_min)*temp)*cl + theta_roo*cr + phi_off*cf
        cs2 = (1 - theta_som*temp)*cs + theta_woo*cw + theta_min*temp*cl
        """
        out = algopy.zeros(6, dtype=p)
        # ACM model
        if self.dC.isday[self.x] == 0.0:
            gpp = 0.0
        elif self.dC.isday[self.x] == 1.0:
            gpp = self.acm(p[18], p[16], p[10], self.dC.acm)
        else:
            ValueError('isday value incorrect')
        auto_resp = self.acm(p[18], p[16], p[10], self.dC.acm)*p[1]
        # Labile release and leaf fall factors
        phi_on = self.phi_onset(p[11], p[13])
        phi_off = self.phi_fall(p[14], p[15], p[4])
        # Temperature factor
        temp = self.temp_term(p[9])

        out[0] = p[17] + (-phi_on*p[17] + (1-p[2])*p[12]*gpp - (1-p[2])*p[12]*auto_resp)*delta_t
        out[1] = p[18] + (-phi_off*p[18] + phi_on*p[17] + (gpp - auto_resp)*p[2])*delta_t
        out[2] = p[19] + (-p[6]*p[19] + (1-p[2])*(1-p[12])*p[3]*gpp\
                           - (1-p[2])*(1-p[12])*p[3]*auto_resp)*delta_t
        out[3] = p[20] + (-p[5]*p[20] + (1-p[2])*(1-p[12])*(1-p[3])*gpp\
                           - (1-p[2])*(1-p[12])*(1-p[3])*auto_resp)*delta_t
        out[4] = p[21] + (-(p[7]+p[0])*temp*p[21] + p[6]*p[19] + phi_off*p[18])*delta_t
        out[5] = p[22] + (-p[8]*temp*p[22] + p[5]*p[20] + p[0]*temp*p[21])*delta_t
        return out

    def jac_dalecv2(self, p, delta_t, vals_to_update, pvals):
        """Using algopy package calculates the jacobian for dalecv2 given a
        input vector p.
        """
        p = algopy.UTPM.init_jacobian(p)
        if len(p)==23:
            return algopy.UTPM.extract_jacobian(self.dalecv2(p, delta_t))
        else:
            for x in xrange(len(p)):
                pvals[vals_to_update[x]]=p[x]
            return algopy.UTPM.extract_jacobian(self.dalecv2(pvals, delta_t))

    def jac3_dalecv2(self, p, paramdict, delta_t, output):
        """Using algopy package calculates the jacobian for dalecv2 given a
        input vector p.
        """
        p = algopy.UTPM.init_jacobian(p)
        return algopy.UTPM.extract_jacobian(self.dalecv2new(p, paramdict, delta_t, output))

    def jac2_dalecv2(self, p, delta_t):
        """Use algopy reverse mode ad calc jac of dv2.
        """
        mat = np.ones((23, 23))*-9999.
        mat[0:17] = np.eye(17, 23)
        p = algopy.UTPM.init_jacobian(p)
        mat[17:] = algopy.UTPM.extract_jacobian(self.dalecv2diff(p, delta_t))
        return mat

    def mod_list(self, pvals):
        """Creates an array of evolving model values using dalecv2 function.
        Takes a list of initial param values.
        """
        mod_list = np.concatenate((np.array([pvals]),
                                  np.ones((self.endrun-self.startrun, len(pvals)))*-9999.))

        self.x = self.startrun
        for t in xrange(self.endrun-self.startrun):
            mod_list[(t+1)] = self.dalecv2(mod_list[t])
            self.x += 1

        self.x -= self.endrun
        return mod_list

    def linmod_list(self, pvals):
        """Creates an array of linearized models (Mi's) taking a list of
        initial param values and a run length (lenrun).
        """
        mod_list = np.concatenate((np.array([pvals]),
                                  np.ones((self.endrun-self.startrun, len(pvals)))*-9999.))
        matlist = np.ones((self.endrun-self.startrun, 23, 23))*-9999.

        self.x = self.startrun
        for t in xrange(self.endrun-self.startrun):
            mod_list[(t+1)] = self.dalecv2(mod_list[t])
            matlist[t] = self.jac2_dalecv2(mod_list[t])
            self.x += 1

        self.x -= self.endrun
        return mod_list, matlist
"""Dalecv2 model class takes a data class and then uses functions to run the
dalecv2 model.
"""
import numpy as np
import copy as cp
import algopy


class DalecModel():

    def __init__(self, dataclass, timestep=0, delta_t=1., output='all'):
        """dataClass and timestep at which to run the dalecv2 model.
        """
        self.dC = dataclass
        self.x = timestep
        self.delta_t = delta_t
        if output == 'all':
            self.output = self.dC.paramdict.keys()
        else:
            self.output = output

        if self.dC.k is None:
            self.lenrun = self.dC.lenrun
        else:
            self.lenrun = self.dC.lenrun*self.dC.k
        self.xb = self.dC.xb
        self.modcoston = 1
        self.modobdict = {'gpp': self.gpp, 'nee': self.nee, 'rt': self.rec,
                          'cf': self.cf, 'clab': self.clab, 'cr': self.cr,
                          'cw': self.cw, 'cl': self.cl, 'cs': self.cs,
                          'lf': self.lf, 'lw': self.lw, 'lai': self.lai,
                          'litresp': self.litresp, 'soilresp': self.soilresp,
                          'rtot': self.rtot, 'rh': self.rh}
        self.endrun = self.lenrun
        #self.yoblist, self.yerroblist = self.obscost()
        #self.rmatrix = self.rmat(self.yerroblist)
        #self.nume = 100


# ------------------------------------------------------------------------------
# Model functions
# ------------------------------------------------------------------------------

    @staticmethod
    def fit_polynomial(ep, multfac):
        """Polynomial used to find phi_f and phi (offset terms used in
        phi_onset and phi_fall), given an evaluation point for the polynomial
        and a multiplication term.
        """
        cf = [2.359978471e-05, 0.000332730053021, 0.000901865258885,
              -0.005437736864888, -0.020836027517787, 0.126972018064287,
              -0.188459767342504]
        poly_val = cf[0]*ep**6 + cf[1]*ep**5 + cf[2]*ep**4 + cf[3]*ep**3 + cf[4]*ep**2 + \
            cf[5]*ep**1 + cf[6]*ep**0
        phi = poly_val*multfac
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
        offset = self.fit_polynomial(1+1e-3, releasecoeff)
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
        offset = self.fit_polynomial(clspan, releasecoeff)
        phi_fall = (2. / np.sqrt(np.pi))*(magcoeff / releasecoeff) * \
            np.exp(-(np.sin((self.dC.D[self.x] - d_fall + offset) /
                   self.dC.radconv)*self.dC.radconv / releasecoeff)**2)
        return phi_fall

    def dalecv2(self, pdiff, paramdict):
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
        p = cp.copy(paramdict)
        for x in xrange(len(self.output)):
            p[self.output[x]] = pdiff[x]
        out = algopy.zeros(len(self.output), dtype=pdiff)
        # ACM model
        if self.dC.isday[self.x] == 0.0:
            gpp = 0.0
        elif self.dC.isday[self.x] == 1.0:
            gpp = self.acm(p['cf'], p['clma'], p['ceff'], self.dC.acm)
        else:
            ValueError('isday value incorrect')
        auto_resp = self.acm(p['cf'], p['clma'], p['ceff'], self.dC.acm)*p['f_auto']
        # Labile release and leaf fall factors
        phi_on = self.phi_onset(p['d_onset'], p['cronset'])
        phi_off = self.phi_fall(p['d_fall'], p['crfall'], p['clspan'])
        # Temperature factor
        temp = self.temp_term(p['Theta'])

        p['clab2'] = p['clab'] + (-phi_on*p['clab'] + (1-p['f_fol'])*p['f_lab']*gpp
                                  - (1-p['f_fol'])*p['f_lab']*auto_resp)*self.delta_t
        p['cf2'] = p['cf'] + (-phi_off*p['cf'] + phi_on*p['clab'] + (gpp - auto_resp)*p['f_fol'])*self.delta_t
        p['cr2'] = p['cr'] + (-p['theta_roo']*p['cr'] + (1-p['f_fol'])*(1-p['f_lab'])*p['f_roo']*gpp
                    - (1-p['f_fol'])*(1-p['f_lab'])*p['f_roo']*auto_resp)*self.delta_t
        p['cw2'] = p['cw'] + (-p['theta_woo']*p['cw'] + (1-p['f_fol'])*(1-p['f_lab'])*(1-p['f_roo'])*gpp
                    - (1-p['f_fol'])*(1-p['f_lab'])*(1 -p['f_roo'])*auto_resp)*self.delta_t
        p['cl2'] = p['cl'] + (-(p['theta_lit']+p['theta_min'])*temp*p['cl']
                              + p['theta_roo']*p['cr'] + phi_off*p['cf'])*self.delta_t
        p['cs2'] = p['cs'] + (-p['theta_som']*temp*p['cs'] + p['theta_woo']*p['cw']
                              + p['theta_min']*temp*p['cl'])*self.delta_t
        x_pvals = p.values() #np.array([p['clab2'], p['cf2'], p['cr2'], p['cw2'], p['cl2'], p['cs2']])
        for x in xrange(len(self.output)):
            if self.output[x] in ['clab', 'cf', 'cr', 'cw', 'cl', 'cs']:
                out[x] = p[self.output[x]+'2']
            else:
                out[x] = p[self.output[x]]
        del p
        return out, x_pvals

    def jac_dalecv2(self, p, paramdict):
        """Using algopy package calculates the jacobian for dalecv2 given a
        input vector p.
        """
        p = algopy.UTPM.init_jacobian(p)
        return algopy.UTPM.extract_jacobian(self.dalecv2(p, paramdict)[0])

    def dalecv2_test(self, pdiff, paramdict):
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
        p = cp.copy(paramdict)
        for x in xrange(len(self.output)):
            p[self.output[x]] = pdiff[x]
        out = algopy.zeros(len(self.output), dtype=pdiff)
        # ACM model
        if self.dC.isday[self.x] == 0.0:
            gpp = 0.0
        elif self.dC.isday[self.x] == 1.0:
            gpp = self.acm(p['cf'], p['clma'], p['ceff'], self.dC.acm)
        else:
            ValueError('isday value incorrect')
        auto_resp = self.acm(p['cf'], p['clma'], p['ceff'], self.dC.acm)*p['f_auto']
        # Labile release and leaf fall factors
        phi_on = self.phi_onset(p['d_onset'], p['cronset'])
        phi_off = self.phi_fall(p['d_fall'], p['crfall'], p['clspan'])
        # Temperature factor
        temp = self.temp_term(p['Theta'])

        p['clab'] = p['clab'] + (-phi_on*p['clab'] + (1-p['f_fol'])*p['f_lab']*gpp
                                  - (1-p['f_fol'])*p['f_lab']*auto_resp)*self.delta_t
        p['cf'] = p['cf'] + (-phi_off*p['cf'] + phi_on*p['clab'] + (gpp - auto_resp)*p['f_fol'])*self.delta_t
        p['cr'] = p['cr'] + (-p['theta_roo']*p['cr'] + (1-p['f_fol'])*(1-p['f_lab'])*p['f_roo']*gpp
                    - (1-p['f_fol'])*(1-p['f_lab'])*p['f_roo']*auto_resp)*self.delta_t
        p['cw'] = p['cw'] + (-p['theta_woo']*p['cw'] + (1-p['f_fol'])*(1-p['f_lab'])*(1-p['f_roo'])*gpp
                    - (1-p['f_fol'])*(1-p['f_lab'])*(1 -p['f_roo'])*auto_resp)*self.delta_t
        p['cl'] = p['cl'] + (-(p['theta_lit']+p['theta_min'])*temp*p['cl']
                              + p['theta_roo']*p['cr'] + phi_off*p['cf'])*self.delta_t
        p['cs'] = p['cs'] + (-p['theta_som']*temp*p['cs'] + p['theta_woo']*p['cw']
                              + p['theta_min']*temp*p['cl'])*self.delta_t
        x_pvals = p.values() #np.array([p['clab2'], p['cf2'], p['cr2'], p['cw2'], p['cl2'], p['cs2']])
        for x in xrange(len(self.output)):
                out[x] = p[self.output[x]]
        del p
        return out, x_pvals

    def jac_dalecv2_test(self, p, paramdict):
        """Using algopy package calculates the jacobian for dalecv2 given a
        input vector p.
        """
        p = algopy.UTPM.init_jacobian(p)
        return algopy.UTPM.extract_jacobian(self.dalecv2_test(p, paramdict)[0])

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

    @staticmethod
    def mfac(matlist, timestep):
        """matrix factorial function, takes a list of matrices and a time step,
        returns the matrix factoral.
        """
        if timestep == -1.:
            return np.eye(23)
        mat = matlist[0]
        for t in xrange(0, timestep):
            mat = np.dot(matlist[t+1], mat)
        return mat


# ------------------------------------------------------------------------------
# Observation functions
# ------------------------------------------------------------------------------

    def gpp(self, p):
        """Function calculates gross primary production (gpp).
        """
        gpp = self.acm(p[18], p[16], p[10], self.dC.acm)
        return gpp

    def rec(self, p):
        """Function calculates total ecosystem respiration (rec).
        """
        rec = p[1]*self.acm(p[18], p[16], p[10], self.dC.acm) + \
            (p[7]*p[21] + p[8]*p[22])*self.temp_term(p[9])
        return rec

    def nee(self, p):
        """Function calculates Net Ecosystem Exchange (nee).
        """
        nee = -(1. - p[1])*self.acm(p[18], p[16], p[10]) + \
            (p[7]*p[21] + p[8]*p[22])*self.temp_term(p[9])
        return nee

    def litresp(self, p):
        """Function calculates litter respiration (litresp).
        """
        litresp = p[7]*p[21]*self.temp_term(p[9])
        return litresp

    def soilresp(self, p):
        """Function calculates soil respiration (soilresp). (heterotrophic)
        """
        soilresp = p[8]*p[22]*self.temp_term(p[9])
        return soilresp

    def rh(self, p):
        """Fn calculates rh (soilresp+litrep).
        """
        rh = self.litresp(p) + self.soilresp(p)
        return rh

    def rtot(self, p):
        """Function calculates soil + root respiration (soilrootresp).
        """
        rtot = p[8]*p[22]*self.temp_term(p[9]) + 5. #Figure this out boi
        return rtot

    def lai(self, p):
        """Fn calculates leaf area index (cf/clma).
        """
        lai = p[18] / p[16]
        return lai

    def lf(self, p):
        """Fn calulates litter fall.
        """
        lf = self.phi_fall(p[14], p[15], p[4])*p[18]
        return lf

    def lw(self, p):
        """Fn calulates litter fall.
        """
        lw = p[5]*p[20]
        return lw

    def clab(self, p):
        """Fn calulates labile carbon.
        """
        clab = p[17]
        return clab

    def cf(self, p):
        """Fn calulates foliar carbon.
        """
        cf = p[18]
        return cf

    def cr(self, p):
        """Fn calulates root carbon.
        """
        cr = p[19]
        return cr

    def cw(self, p):
        """Fn calulates woody biomass carbon.
        """
        cw = p[20]
        return cw

    def cl(self, p):
        """Fn calulates litter carbon.
        """
        cl = p[21]
        return cl

    def cs(self, p):
        """Fn calulates soil organic matter carbon.
        """
        cs = p[22]
        return cs

    def linob(self, ob, pvals):
        """Function returning jacobian of observation with respect to the
        parameter list. Takes an obs string, a parameters list, a dataClass
        and a time step x.
        """
        dpvals = algopy.UTPM.init_jacobian(pvals)
        return algopy.UTPM.extract_jacobian(self.modobdict[ob](dpvals))


# ------------------------------------------------------------------------------
# Assimilation functions
# ------------------------------------------------------------------------------

    def bmat(self, corr=False, varyp=False):
        """Attempt at creating a b matrix.
        """
        pmat = np.ones((23, 1000))*9999.
        modevmat = np.ones((23, 1000))*9999.

        if varyp is False:
            for x in xrange(23):
                if x < 17.:
                    for i in xrange(1000):
                        pmat[x, i] = self.dC.pvals[x]
                elif x >= 17.:
                    pmat[x] = np.random.normal(self.dC.pvals[x],
                                               self.dC.pvals[x]*0.3, 1000)
        else:
            for x in xrange(23):
                pmat[x] = np.random.normal(self.dC.pvals[x],
                                           self.dC.pvals[x]*0.3, 1000)

        for x in xrange(1000):
            if pmat[4, x] < 1.0001:
                pmat[4, x] = 1.0001

        for x in xrange(1000):
            modevmat[:, x] = self.mod_list(pmat[:, x])[-1]

        if corr is False:
            return np.cov(modevmat)
        elif corr is True:
            return np.corrcoef(modevmat)

    def obscost(self):
        """Function returning list of observations and a list of their
        corresponding error values. Takes observation dictionary and an
        observation error dictionary.
        """
        yoblist = np.array([])
        yerrlist = np.array([])
        for t in xrange(self.startrun, self.endrun):
            for ob in self.dC.ob_dict.iterkeys():
                if np.isnan(self.dC.ob_dict[ob][t]) != True:
                    yoblist = np.append(yoblist, self.dC.ob_dict[ob][t])
                    yerrlist = np.append(yerrlist,
                                         self.dC.oberrdict[ob+'_err'][t])
        return yoblist, yerrlist

    def hxcost(self, pvallist):
        """Function returning a list of observation values as predicted by the
        DALEC model. Takes a list of model values (pvallist), an observation
        dictionary and a dataClass (dC).
        """
        hx = np.array([])
        self.x = self.startrun
        for t in xrange(self.startrun, self.endrun):
            for ob in self.dC.ob_dict.iterkeys():
                if np.isnan(self.dC.ob_dict[ob][t]) != True:
                    hx = np.append(hx,
                                   self.modobdict[ob](pvallist[t-self.startrun]))
            self.x += 1

        self.x -= self.endrun
        return hx

    def rmat(self, yerr):
        """Returns observation error covariance matrix given a list of
        observation error values.
        """
        r = (yerr**2)*np.eye(len(yerr))
        return r

    def hmat(self, pvallist, matlist):
        """Returns a list of observation values as predicted by DALEC (hx) and
        a linearzied observation error covariance matrix (hmat). Takes a list
        of model values (pvallist), a observation dictionary, a list of
        linearized models (matlist) and a dataClass (dC).
        """
        hx = np.array([])
        hmat = np.array([])
        self.x = self.startrun
        for t in xrange(self.startrun, self.endrun):
            temp = []
            for ob in self.dC.ob_dict.iterkeys():
                if np.isnan(self.dC.ob_dict[ob][t]) != True:
                    hx = np.append(hx,
                                   self.modobdict[ob](pvallist[t-self.startrun]))
                    temp.append([self.linob(ob, pvallist[t-self.startrun])])
            self.x += 1
            if len(temp) != 0.:
                hmat = np.append(hmat, np.dot(np.vstack(temp),
                                 self.mfac(matlist, t-self.startrun-1)))
            else:
                continue

        self.x -= self.endrun
        hmat = np.reshape(hmat, (len(hmat)/23, 23))
        return hx, hmat

    def modcost(self, pvals):
        """model part of cost fn.
        """
        return np.dot(np.dot((pvals-self.xb), np.linalg.inv(self.dC.B)), (pvals-self.xb).T)

    def obcost(self, pvals):
        """Observational part of cost fn.
        """
        pvallist = self.mod_list(pvals)
        hx = self.hxcost(pvallist)
        return np.dot(np.dot((self.yoblist-hx), np.linalg.inv(self.rmatrix)), (self.yoblist-hx).T)

    def cost(self, pvals):
        """4DVAR cost function to be minimized. Takes an initial state (pvals),
        an observation dictionary, observation error dictionary, a dataClass
        and a start and finish time step.
        """
        ob_cost = self.obcost(pvals)
        if self.modcoston is True:
            mod_cost = self.modcost(pvals)
        else:
            mod_cost = 0
        cost = 0.5*ob_cost + 0.5*mod_cost
        return cost

    def gradcost(self, pvals):
        """Gradient of 4DVAR cost fn to be passed to optimization routine.
        Takes an initial state (pvals), an obs dictionary, an obs error
        dictionary, a dataClass and a start and finish time step.
        """
        pvallist, matlist = self.linmod_list(pvals)
        hx, hmatrix = self.hmat(pvallist, matlist)
        obcost = np.dot(hmatrix.T, np.dot(np.linalg.inv(self.rmatrix),
                                          (self.yoblist-hx).T))
        if self.modcoston is True:
            modcost = np.dot(np.linalg.inv(self.dC.B), (pvals-self.xb).T)
        else:
            modcost = 0
        gradcost = - obcost + modcost
        return gradcost
#!/usr/bin/env python
import numpy as np

from cython_funcs import *

class MoI:
    '''Class for calculating the potential in a semi-infinite slice of neural tissue.
    Set-up:


              SALINE -> sigma_S 

    <----------------------------------------------------> z = + h
    
              TISSUE -> sigma_T


                   o -> charge_pos = [x',y',z']


    <-----------*----------------------------------------> z = 0
                 \-> elec_pos = [x,y,z] 

                 ELECTRODE GLASS PLATE -> sigma_G 
        

    Arguments:
         sigma_G: 0.0, # Conductivity below electrode
         sigma_T: 0.3, # Conductivity of tissue
         sigma_S: 1.5, # Conductivity of saline
         h: 300., # Slice thickness in um
         steps : 20, # How many steps to include of the infinite series

    The cython functions can be used independently of this class.
    This package assumes that the slice extends from 0 to h in the z-direction, however, the convention
    from -a to a (with a = h/2) is used in functions that end with '_shifted'.

    '''
    def __init__(self, sigma_G=0.0, sigma_T=0.3, sigma_S=1.5, h=300., steps=20, **kwargs):

        self.sigma_G = sigma_G
        self.sigma_T = sigma_T
        self.sigma_S = sigma_S
        self._check_for_anisotropy()
        self.h = h
        self.a = self.h/2.
        self.steps = steps

    def _anisotropic_saline_scaling(self):
        """ To make formula work in anisotropic case we scale the conductivity of the
        saline to be a scalar k times the tissue conductivity. (Wait 1990)
        """

        ratios = np.array(self.sigma_S) / np.array(self.sigma_T)

        if np.abs(ratios[0] - ratios[2]) <= 1e-15:
            scale_factor = ratios[0]
        elif np.abs(ratios[1] - ratios[2]) <= 1e-15:
            scale_factor = ratios[1]
        else:
            raise RuntimeError("Anisotropy structure not understood.")

        sigma_S_scaled = scale_factor * np.array(self.sigma_T)
        sigma_T_net = np.sqrt(self.sigma_T[0] * self.sigma_T[2])
        sigma_S_net = np.sqrt(sigma_S_scaled[0] * sigma_S_scaled[2])

        # print "Sigma_T: %s, Sigma_S: %s, Sigma_S_scaled: %s, scale factor: %g" % (
        #     self.sigma_T, self.sigma_S, sigma_S_scaled, scale_factor)
        self.anis_W = (sigma_T_net - sigma_S_net)/(sigma_T_net + sigma_S_net)

    def _check_for_anisotropy(self):
        """ Checks if input conductivities are tensors or scalars
        and sets self.is_anisotropic correspondingly
        """
        types = [type(self.sigma_T), type(self.sigma_S), type(self.sigma_G)]

        if (list in types) or (np.ndarray in types):
            self.is_anisotropic = True

            if type(self.sigma_G) in [list, np.ndarray]:
                if len(self.sigma_G) != 3:
                    raise ValueError("Conductivity vector but not with size 3")
                self.sigma_G = np.array(self.sigma_G)
            else:
                self.sigma_G = np.array([self.sigma_G, self.sigma_G, self.sigma_G])
            if type(self.sigma_T) in [list, np.ndarray]:
                if len(self.sigma_T) != 3:
                    raise ValueError("Conductivity vector but not with size 3")
                self.sigma_T = np.array(self.sigma_T)
            else:
                self.sigma_T = np.array([self.sigma_T, self.sigma_T, self.sigma_T])

            if type(self.sigma_S) in [list, np.ndarray]:
                if len(self.sigma_S) != 3:
                    raise ValueError("Conductivity vector but not with size 3")
                self.sigma_S = np.array(self.sigma_S)
            else:
                self.sigma_S = np.array([self.sigma_S, self.sigma_S, self.sigma_S])
       
            self._anisotropic_saline_scaling()
            if (self.sigma_G[0] == self.sigma_G[1] == self.sigma_G[2]) and \
               (self.sigma_T[0] == self.sigma_T[1] == self.sigma_T[2]) and \
               (self.sigma_S[0] == self.sigma_S[1] == self.sigma_S[2]):
                print "Isotropic conductivities can be given as scalars."         
        else:
            self.is_anisotropic = False
 
    def in_domain(self, elec_pos, charge_pos):
        """ Checks if elec_pos and charge_pos is within valid area.
        Otherwise raise exception."""

        # If inputs are single positions
        if (np.array(elec_pos).shape == (3,)) and (np.array(charge_pos).shape == (3,)):
            elec_pos = [elec_pos]
            charge_pos = [charge_pos]

        for epos in elec_pos:
            if not np.abs(epos[2]) < 1e-14:
                raise RuntimeError("Electrode plane not at z=0.")
        for cpos in charge_pos:
            if np.abs(cpos[2]) > self.h + 1e-9:
                raise RuntimeError("Charge not within valid range.")
        for cpos in charge_pos:
            for epos in elec_pos:
                dist = np.sqrt(np.sum((np.array(cpos) - np.array(epos))**2))
                if dist < 1e-6:
                    raise RuntimeError("Charge and electrode at same position!")

    def in_domain_shifted(self, elec_pos, charge_pos):
        """ Checks if elec_pos and charge_pos is within valid area.
        Otherwise raise exception."""

        # If inputs are single positions
        if (np.array(elec_pos).shape == (3,)) and (np.array(charge_pos).shape == (3,)):
            elec_pos = [elec_pos]
            charge_pos = [charge_pos]

        for epos in elec_pos:
            if not np.abs(epos[2] + self.a) <= 1e-14:
                raise RuntimeError("Electrode not within valid range.")
        for cpos in charge_pos:

            if np.abs(cpos[2]) >= self.a + 1e-9:
                raise RuntimeError("Charge not within valid range.")
        for cpos in charge_pos:
            for epos in elec_pos:
                dist = np.sqrt(np.sum((np.array(cpos) - np.array(epos))**2))
                if dist < 1e-6:
                    raise RuntimeError("Charge and electrode at same position!")

    def anisotropic_saline_scaling(self, charge_pos, elec_pos, imem=1):
        """ Calculate the moi point source potential with saline conductivity
        sigma_S is scaled to k * sigma_T. There is also a much faster cython version of this"""
        self.in_domain(elec_pos, charge_pos)
        x0, y0, z0 = charge_pos[:]
        x, y, z = elec_pos[:]

        def _omega(dz):
            return 1/np.sqrt(self.sigma_T[0]*self.sigma_T[2]*(y - y0)**2 +
                             self.sigma_T[0]*self.sigma_T[1]*dz**2 +
                             self.sigma_T[1]*self.sigma_T[2]*(x - x0)**2) 
        phi = _omega(-z0)
        n = 1
        while n < self.steps:
            phi += self.anis_W**n * (_omega(2*n*self.h - z0) + _omega(-2*n*self.h - z0))
            n += 1   
        phi *= 2*imem/(4*np.pi)
        return phi

    def isotropic_moi_shifted(self, charge_pos, elec_pos, imem=1):
        """ This function calculates the potential at the position elec_pos = [x, y, z]
        set up by the charge at position charge_pos = [x0, y0, z0]. To get get the potential
        from multiple charges, the contributions must be summed up.
        """
        def _omega(dz):
            return 1/np.sqrt((y - y0)**2 + (x - x0)**2 + dz**2)
        x0, y0, z0 = charge_pos[:]
        x, y, z = elec_pos[:]
        phi = _omega(z - z0)
        n = 0
        WTS = (self.sigma_T - self.sigma_S)/(self.sigma_T + self.sigma_S)
        WTG = (self.sigma_T - self.sigma_G)/(self.sigma_T + self.sigma_G)
        while n < self.steps:
            if n == 0:
                phi += WTS * _omega(z + z0 - (4*n + 2)*self.a) + WTG * _omega(z + z0 + (4*n + 2)*self.a)
            else:
                phi += (WTS*WTG)**n *(WTS * _omega(z + z0 - (4*n + 2)*self.a) + WTG * _omega(z + z0 + (4*n + 2)*self.a)
                                      + _omega(z - z0 + 4*n*self.a) + _omega(z - z0 - 4*n*self.a))
            n += 1
        phi *= imem/(4*np.pi*self.sigma_T)
        return phi

    def isotropic_moi(self, charge_pos, elec_pos, imem=1):
        """ This function calculates the potential at the position elec_pos = [x, y, z]
        set up by the charge at position charge_pos = [x0, y0, z0]. To get get the potential
        from multiple charges, the contributions must be summed up.
        """
        def _omega(dz):
            return 1/np.sqrt((y - y0)**2 + (x - x0)**2 + dz**2)
        x0, y0, z0 = charge_pos[:]
        x, y, z = elec_pos[:]
        phi = _omega(z - z0)
        n = 0
        WTS = (self.sigma_T - self.sigma_S)/(self.sigma_T + self.sigma_S)
        WTG = (self.sigma_T - self.sigma_G)/(self.sigma_T + self.sigma_G)
        while n < self.steps:
            if n == 0:
                phi += (WTS * _omega(z + z0 - 2*(n + 1)*self.h) +
                        WTG * _omega(z + z0 + 2*n*self.h))
            else:
                phi += (WTS*WTG)**n * (WTS * _omega(z + z0 - 2*(n + 1)*self.h) +
                                       WTG * _omega(z + z0 + 2*n*self.h) +
                                       _omega(z - z0 + 2*n*self.h) + _omega(z - z0 - 2*n*self.h))
            n += 1
        phi *= imem/(4*np.pi*self.sigma_T)
        return phi

    def line_source_moi_shifted(self, comp_start, comp_end, comp_length, elec_pos, imem=1):
        """ Calculate the moi line source potential at electrode plane"""
        self.in_domain_shifted(elec_pos, comp_start)
        self.in_domain_shifted(elec_pos, comp_end)
        x0, y0, z0 = comp_start[:]
        x1, y1, z1 = comp_end[:]
        x, y, z = elec_pos[:]
        dx = x1 - x0
        dy = y1 - y0
        dz = z1 - z0
        a_x = x - x0
        a_y = y - y0
        W = (self.sigma_T - self.sigma_S)/(self.sigma_T + self.sigma_S)
        def _omega(a_z):
            #See Rottman integration formula 46) page 137 for explanation
            factor_a = comp_length*comp_length
            factor_b = - a_x*dx - a_y*dy - a_z * dz
            factor_c = a_x*a_x + a_y*a_y + a_z*a_z
            b_2_ac = factor_b*factor_b - factor_a * factor_c
            if np.abs(b_2_ac) <= 1e-16:
                num = factor_a + factor_b
                den = factor_b
            else:
                num = factor_a + factor_b + comp_length*np.sqrt(factor_a + 2*factor_b + factor_c)
                den = factor_b + comp_length*np.sqrt(factor_c)
            return np.log(num/den)
        phi = _omega(-self.a - z0)
        n = 1
        while n < self.steps:
            phi += W**n * (_omega((4*n-1)*self.a - z0) + _omega(-(4*n+1)*self.a - z0))
            n += 1   
        phi *= 2*imem/(4*np.pi*self.sigma_T * comp_length)
        return phi

    def line_source_moi(self, comp_start, comp_end, comp_length, elec_pos, imem=1):
        """ Calculate the moi line source potential at electrode plane"""
        self.in_domain(elec_pos, comp_start)
        self.in_domain(elec_pos, comp_end)
        x0, y0, z0 = comp_start[:]
        x1, y1, z1 = comp_end[:]
        x, y, z = elec_pos[:]
        dx = x1 - x0
        dy = y1 - y0
        dz = z1 - z0
        a_x = x - x0
        a_y = y - y0
        W = (self.sigma_T - self.sigma_S)/(self.sigma_T + self.sigma_S)
        def _omega(a_z):
            #See Rottman integration formula 46) page 137 for explanation
            factor_a = comp_length*comp_length
            factor_b = - a_x*dx - a_y*dy - a_z * dz
            factor_c = a_x*a_x + a_y*a_y + a_z*a_z
            b_2_ac = factor_b*factor_b - factor_a * factor_c
            if np.abs(b_2_ac) <= 1e-16:
                num = factor_a + factor_b
                den = factor_b
            else:
                num = factor_a + factor_b + \
                      comp_length*np.sqrt(factor_a + 2*factor_b + factor_c)
                den = factor_b + comp_length*np.sqrt(factor_c)
            return np.log(num/den)
        phi = _omega(-z0)
        n = 1
        while n < self.steps:
            phi += W**n * (_omega(2*n*self.h - z0) + _omega(-2*n*self.h - z0))
            n += 1
        phi *= 2*imem/(4*np.pi*self.sigma_T * comp_length)
        return phi

    def point_source_moi_at_elec_plane_shifted(self, charge_pos, elec_pos, imem=1):
        """ Calculate the moi point source potential assuming electrode at MEA electrode plane (elec_z = -a)"""
        self.in_domain_shifted(elec_pos, charge_pos)
        x0, y0, z0 = charge_pos[:]
        x, y, z = elec_pos[:]
        W = (self.sigma_T - self.sigma_S)/(self.sigma_T + self.sigma_S)
        def _omega(dz):
            return 1/np.sqrt((y - y0)**2 + (x - x0)**2 + dz**2)
        phi = _omega(-self.a - z0)
        n = 1
        while n < self.steps:
            phi += W**n * (_omega((4*n-1)*self.a - z0) + _omega(-(4*n+1)*self.a - z0))
            n += 1   
        return 2*phi*imem/(4*np.pi*self.sigma_T)

    def point_source_moi_at_elec_plane(self, charge_pos, elec_pos, imem=1):
        """ Calculate the moi point source potential assuming electrode at MEA electrode plane (elec_z = 0)"""
        self.in_domain(elec_pos, charge_pos)
        x0, y0, z0 = charge_pos[:]
        x, y, z = elec_pos[:]
        W = (self.sigma_T - self.sigma_S)/(self.sigma_T + self.sigma_S)
        def _omega(dz):
            return 1/np.sqrt((y - y0)**2 + (x - x0)**2 + dz**2)
        phi = _omega(z0)
        n = 1
        while n < self.steps:
            phi += W**n * (_omega(2*n*self.h - z0) + _omega(-2*n*self.h - z0))
            n += 1
        return 2*phi*imem/(4*np.pi*self.sigma_T)

    def potential_at_elec_big_average(self, elec_pos, r, n_avrg_points, function, func_args):
        """ Calculate the potential at electrode 'elec_index' with n_avrg_points points"""
        phi = 0.
        for pt in xrange(n_avrg_points):
            pt_pos = np.array([(np.random.rand() - 0.5) * 2 * r,
                               (np.random.rand() - 0.5) * 2 * r])
            # If outside electrode
            while np.sum(pt_pos**2) > r**2:
                pt_pos = np.array([(np.random.rand() - 0.5) * 2 * r,
                                   (np.random.rand() - 0.5) * 2 * r])
            avrg_point_pos = [elec_pos[0] + pt_pos[0],
                              elec_pos[1] + pt_pos[1],
                              elec_pos[2]]
            phi += function(*func_args, elec_pos=avrg_point_pos)
        return phi/n_avrg_points

    def make_mapping(self, ext_sim_dict, xmid=None, ymid=None, zmid=None,
                                xstart=None, ystart=None, zstart=None,
                                xend=None, yend=None, zend=None):
        """ Make a mapping given two arrays of electrode positions"""
        elec_x = ext_sim_dict['elec_x'] # Array
        elec_y = ext_sim_dict['elec_y'] # Array
        elec_z = ext_sim_dict['elec_z'] # Scalar

        n_elecs = len(elec_x)
        if ext_sim_dict['include_elec']:
            n_avrg_points = ext_sim_dict['n_avrg_points']
        if ext_sim_dict['use_line_source']:
            function = self.line_source_moi
            n_compartments = len(xstart)
        else:
            function = self.point_source_moi_at_elec_plane
            n_compartments = len(xmid)

        mapping = np.zeros((n_elecs, n_compartments))
        for comp in xrange(n_compartments):
            for elec in xrange(n_elecs):
                elec_pos = [elec_x[elec], elec_y[elec], elec_z]
                if ext_sim_dict['use_line_source']:
                    comp_start = np.array([xstart[comp], ystart[comp], zstart[comp]])
                    comp_end = np.array([xend[comp], yend[comp], zend[comp]])
                    comp_length = np.sqrt(np.sum((comp_end - comp_start)**2))
                    func_args = [comp_start, comp_end, comp_length]
                else:
                    charge_pos = [xmid[comp], ymid[comp], zmid[comp]]
                    func_args = [charge_pos]
                if ext_sim_dict['include_elec']:
                    mapping[elec, comp] = self.potential_at_elec_big_average(elec_pos, ext_sim_dict['elec_radius'],
                                                                             n_avrg_points, function, func_args)
                else:
                    mapping[elec, comp] = function(*func_args, elec_pos=elec_pos)
        return mapping

    def make_mapping_free_z(self, ext_sim_dict, xmid=None, ymid=None, zmid=None,
                                xstart=None, ystart=None, zstart=None,
                                xend=None, yend=None, zend=None):
        """ Make a mapping given three arrays of electrode positions.
        Electrode z-position is not confined to electrode plane"""
        elec_x = ext_sim_dict['elec_x']  # Array
        elec_y = ext_sim_dict['elec_y']  # Array
        elec_z = ext_sim_dict['elec_z']  # Scalar or array

        if type(elec_z) in [int, float]:
            # If elec_z is a single number, we make an array out of it.
            elec_z = elec_z * np.ones(len(elec_x))

        n_elecs = len(elec_x)
        if ext_sim_dict['include_elec']:
            raise NotImplementedError()
        if ext_sim_dict['use_line_source']:
            raise NotImplementedError()

        mapping = np.zeros((n_elecs, len(xmid)))
        for comp in xrange(len(xmid)):
            for elec in xrange(n_elecs):
                elec_pos = [elec_x[elec], elec_y[elec], elec_z[elec]]
                charge_pos = [xmid[comp], ymid[comp], zmid[comp]]
                mapping[elec, comp] = self.isotropic_moi(charge_pos, elec_pos)
        return mapping


    def make_mapping_cython(self, ext_sim_dict, xmid=None, ymid=None, zmid=None,
                            xstart=None, ystart=None, zstart=None,
                            xend=None, yend=None, zend=None):
        """ Make a mapping given two arrays of electrode positions.
        Electrode z position in confined to the electrode plane

        Parameters:
        ---------------
        ext_sim_dict: {'include_elec': Boolean value, specifying if the result should be averaged over an electrode
                       'n_avrg_points': If 'include_elec' == True, this is the number of averaging points
                       'elec_radius': If 'include_elec' == True, this is the electrode radii
                       'use_line_source': Boolean value deciding whether the line- or point-source formula is used.
                       'elec_x': Array with electrode x-positions
                       'elec_y': Array with electrode y-positions
                       }
        Point source positions
        (x, y, z)mid: Middle point (for point source formula)
        (x, y, z)start: Starting point (for line source formula)
        (x, y, z)end: Ending point (for line source formula)

        Returns:
        -----------------
        mapping : A (number of electrodes x number of sources) matrix mapping from source-current to potentials.

        Phi(x, t) = np.dot(mapping, source_currents(t))

        """
        print "Making mapping. Cython style"
        elec_x = ext_sim_dict['elec_x']
        elec_y = ext_sim_dict['elec_y']
        elec_z = ext_sim_dict['elec_z']

        if xmid != None:
            xmid = np.array(xmid, order='C')
            ymid = np.array(ymid, order='C')
            zmid = np.array(zmid, order='C')

        if xstart != None:
            xend = np.array(xend, order='C')
            yend = np.array(yend, order='C')
            zend = np.array(zend, order='C')
            xstart = np.array(xstart, order='C')
            ystart = np.array(ystart, order='C')
            zstart = np.array(zstart, order='C')
            
        if ext_sim_dict['include_elec']:
            n_avrg_points = ext_sim_dict['n_avrg_points']
            elec_r = ext_sim_dict['elec_radius']
            if ext_sim_dict['use_line_source']:
                function = LS_with_elec_mapping
                func_args = [self.sigma_T, self.sigma_S, self.h, self.steps, n_avrg_points, elec_r,
                             elec_x, elec_y, xstart, ystart, zstart, xend, yend, zend]
            else:
                function = PS_with_elec_mapping
                func_args = [self.sigma_T, self.sigma_S, self.h, self.steps, n_avrg_points,
                             elec_r, elec_x, elec_y, xmid, ymid, zmid]
        else:
            if ext_sim_dict['use_line_source']:
                function = LS_without_elec_mapping
                func_args = [self.sigma_T, self.sigma_S, self.h, self.steps, elec_x, elec_y,
                             xstart, ystart, zstart, xend, yend, zend]
            else:
                function = PS_without_elec_mapping
                func_args = [self.sigma_T, self.sigma_S, self.h, self.steps, elec_x, elec_y, xmid, ymid, zmid]
        mapping = function(*func_args)
        return mapping


def plot_decay_with_distance_example():
    # Source positions
    xmid = np.array([0.])
    ymid = np.array([0.])
    zmid = np.array([100.])

    # Electrode positions
    elec_x = np.linspace(0, 1000, 100)  # um
    elec_y = np.zeros(len(elec_x))  # um
    elec_z = 0.  # um
    elec_clrs = ['r', 'b', 'g']

    # Source currents
    t = np.linspace(0, 1, 1)  # ms
    imem = np.array([[1.]])  # nA

    h = 200 # slice thickness [um]

    ext_sim_dict = {'include_elec': False,
                    'elec_x': elec_x,
                    'elec_y': elec_y,
                    'elec_z': elec_z,
                    'n_avrg_points': 100,
                    'elec_radius': 10,
                    'use_line_source': False,
                    }

    mapping = MoI(h=h, sigma_S=0.3, sigma_T=0.3).make_mapping_cython(ext_sim_dict, xmid=xmid, ymid=ymid, zmid=zmid)
    phi_homo = 1000 * np.dot(mapping, imem)  # uV

    mapping = MoI(h=h, sigma_S=2.0, sigma_T=0.3).make_mapping_cython(ext_sim_dict, xmid=xmid, ymid=ymid, zmid=zmid)
    phi_hetro = 1000 * np.dot(mapping, imem)  # uV

    import pylab as plt
    plt.close('all')
    plt.subplots_adjust(hspace=0.5)
    plt.subplot(211, ylim=[elec_z - 20, h + 20], xlim=[-600, 600],
                ylabel='z [$\mu m$]', xlabel='x [$\mu m$]', aspect='equal')
    plt.plot(xmid, zmid, 'y*', label='Sources', ms=15)
    plt.scatter(elec_x, np.ones(len(elec_x))*elec_z, c=elec_clrs, marker='s', s=35, label='Electrodes')
    plt.plot([-10000, 10000], [elec_z, elec_z], 'k-')
    plt.plot([-10000, 10000], [h, h], 'k-')
    plt.legend()
    plt.subplot(223, title='Measured potentials', ylabel='$\Phi(t)$ [$\mu V$]', xlabel='Distance')
    plt.plot(elec_x, phi_homo[:, 0], '-', c='k', lw=2)
    plt.plot(elec_x, phi_hetro[:, 0], '-', c='g', lw=2)

    plt.subplot(224, title='Normalized\nmeasured potentials', ylabel='$\Phi(t)$', xlabel='Distance')
    l1, = plt.plot(elec_x, phi_homo[:, 0] / phi_homo[0, 0], '-', c='k', lw=2, label='Infinite slice')
    l2, = plt.plot(elec_x, phi_hetro[:, 0] / phi_hetro[0, 0], '-', c='g', lw=2, label='Saline immersion')

    l3, = plt.plot(elec_x, np.sqrt(elec_x[0]**2 + 100**2) / np.sqrt(elec_x**2 + 100**2), ':', lw=2, c='r', label='1/r')
    plt.legend(frameon=False)

    plt.show()

def plot_sinusoidal_current_example():
    # Source positions
    xmid = np.array([-100., 100.])
    ymid = np.array([0., 0.])
    zmid = np.array([80., 110.])

    # Electrode positions
    elec_x = np.array([-100., 0., 100.])  # um
    elec_y = np.array([0., 0., 0.])  # um
    elec_z = 0.  # um
    elec_clrs = ['r', 'b', 'g']

    # Source currents
    t = np.linspace(0, 1, 100)  # ms
    imem = np.array([np.sin(2 * np.pi * 100 * t + theta) for theta in [0, np.pi]])  # nA

    h = 200 # slice thickness [um]

    ext_sim_dict = {'include_elec': True,
                    'elec_x': elec_x,
                    'elec_y': elec_y,
                    'elec_z': elec_z,
                    'n_avrg_points': 100,
                    'elec_radius': 10,
                    'use_line_source': False,
                    }
    mapping = MoI(h=h).make_mapping_cython(ext_sim_dict, xmid=xmid, ymid=ymid, zmid=zmid)
    phi_homo = 1000 * np.dot(mapping, imem)  # uV

    import pylab as plt
    plt.close('all')
    plt.subplots_adjust(hspace=0.5)
    plt.subplot(211, ylim=[elec_z - 20, h + 20], xlim=[-600, 600],
                ylabel='z [$\mu m$]', xlabel='x [$\mu m$]', aspect='equal')
    plt.plot(xmid, zmid, 'y*', label='Sources', ms=15)
    plt.scatter(elec_x, np.ones(len(elec_x))*elec_z, c=elec_clrs, marker='s', s=35, label='Electrodes')
    plt.plot([-10000, 10000], [elec_z, elec_z], 'k-')
    plt.plot([-10000, 10000], [h, h], 'k-')
    plt.legend()
    plt.subplot(212, title='Measured potentials', ylabel='$\Phi(t)$ [$\mu V$]', xlabel='Time [ms]')
    [plt.plot(t, phi_homo[idx, :], c=elec_clrs[idx]) for idx in xrange(len(elec_x))]

    plt.show()

if __name__ == '__main__':
    plot_sinusoidal_current_example()
    plot_decay_with_distance_example()

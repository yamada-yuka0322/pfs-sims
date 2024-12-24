'''

module for HOD to populate halos with galaxies 


'''
import numpy as np
from scipy.stats import norm
from scipy.special import erf

# halotools
from halotools.empirical_models import HodModelFactory, OccupationComponent
from halotools.sim_manager import UserSuppliedHaloCatalog
from halotools.empirical_models import NFWPhaseSpace


class CustomCentrals(OccupationComponent):
    r"""HOD model for the central galaxy occupation. Populates central galaxies depending on the halo mass.
    """
    def __init__(self, threshold=-23, **params):
        r"""
        Parameters
        ----------
        threshold : float, optional
            Logarithm of the primary galaxy property threshold. If the primary
            galaxy property is luminosity, it is given in h=1 solar luminosity
            units.

        params : dictionary
            Dictionary of HOD parameters used for the gaussian HOD model.
            Must include 'logMc','sigma_M' and 'Ac'
            
        Example
        --------
        >>> cen_model = CustomCentrals(**params_dict)

        """
        super(CustomCentrals, self).__init__(gal_type='centrals', threshold=threshold, upper_occupation_bound=1)
        self.param_dict = params

    def mean_occupation(self, table):
        r"""Expected number of central galaxies in a halo. Derived from the Gaussian HOD models and the given HOD parameters.

         Parameters
         ----------
         table : object
             Data table storing halo catalog.

         Returns
         -------
         mean_ncen : array
             Mean number of central galaxies in halos of the input mass.
        """
        halo_mass = table['halo_mvir'][:]
        logMc = self.param_dict['logMc']
        sigma_M = self.param_dict['sigma_M']
        Ac = self.param_dict['Ac']
        
        return Ac/(2*np.pi)**0.5/sigma_M * np.exp(-(np.log10(halo_mass) - logMc)**2/2/sigma_M**2)

class CustomSatellites(OccupationComponent):
    r"""Power-law model for the satellite galaxy occupation. 
    Occupates satellite galaxies depending on the halo mass.
    """
    def __init__(self, threshold=-23, **params):
        r"""
        Parameters
        ----------
        threshold : float, optional
            Logarithm of the primary galaxy property threshold. If the primary
            galaxy property is luminosity, it is given in h=1 solar luminosity
            units.

        params : dictionary
            Dictionary of HOD parameters used for the gaussian HOD model.
            Must include 'As','logM0' and 'alpha'
            
        Example
        --------
        >>> sat_model = CustomSatellites(**params_dict)

        """
        super(CustomSatellites, self).__init__(gal_type='satellites', threshold=threshold, upper_occupation_bound=float('inf'), **kwargs)
        self.param_dict = params

    def mean_occupation(self, table):
        r"""Expected number of central galaxies in a halo. Derived from the Gaussian HOD models and the given HOD parameters.

         Parameters
         ----------
         table : object
             Data table storing halo catalog.

         Returns
         -------
         mean_ncen : array
             Mean number of satellite galaxies in halos of the input mass.
        """
        halo_mass = table['halo_mvir'][:]
        As = self.param_dict['As']
        logM0 = self.param_dict['logM0']
        alpha = self.param_dict['alpha']

        logM1 = 13
        dist = As*((halo_mass-10**logM0)/10**logM1)**alpha

        return np.nan_to_num(np.maximum(dist, 0) * (halo_mass > 10**logM0))


class Velocity(object):
    r"""Populate satellite galaxy velocity with gaussian distribution aroung the halo velocity.
    The dispersion is calculated as halo particle velocity dispersion times galaxy velocity bias fv.
    """
    def __init__(self, gal_type, fv):
        r"""
        Parameters
        ----------
        gal_type : string
            Type of galaxy

        fv : float
            galaxy velocity bias.
            
        Example
        --------
        >>> sat_velocity = Velocity(gal_type = 'satellites', fv = 1.23)

        """
        self.gal_type = gal_type
        self._mock_generation_calling_sequence = ['assign_velocity']
        self._galprop_dtypes_to_allocate = np.dtype([('vx', 'f8'),('vy', 'f8'), ('vz', 'f8')])
        self.list_of_haloprops_needed = ['halo_vx', 'halo_vy', 'halo_vz', 'halo_sigmav']
        self.fv = fv

    def assign_velocity(self, table, seed):
        table['vx'][:] = table['halo_vx'][:]+norm.rvs(loc=0.0, scale=1.0, size=len(table["halo_vz"][:]))*self.fv/3**0.5*table['halo_sigmav'][:]
        table['vy'][:] = table['halo_vy'][:]+norm.rvs(loc=0.0, scale=1.0, size=len(table["halo_vz"][:]))*self.fv/3**0.5*table['halo_sigmav'][:]
        table['vz'][:] = table['halo_vz'][:]+norm.rvs(loc=0.0, scale=1.0, size=len(table["halo_vz"][:]))*self.fv/3**0.5*table['halo_sigmav'][:]
        

def PopulateGalaxies(halo, **cent_params, **sat_params):
    """ Populate 
    
    Parameters
    ----------
    halo : object
        Halo class instance, which holds the information of halos.
        Must include
        redshift, boxsize, particle mass, halo position, halo velocity, virial radius, scale radius and the halo particle velocity dispersion

    cent_params : Dictionary
        Dictionary of HOD parameters used for the central galaxy occupation.
        Must include 'logMc','sigma_M' and 'Ac'.

    sat_params : Dictionary
        Dictionary of HOD parameters used for the satellite galaxy occupation.
        Must include 'As','logM0' and 'alpha'.
        
    Returns
    --------
    table of galaxy properties

    """
    redshift = halo.redshift

    custom_phase_space = NFWPhaseSpace(conc_mass_model='direct_from_halo_catalog', redshift=redshift)
    custom_cens = CustomCentrals(**cent_params, redshift = redshift)
    custom_sats = CustomSatellites(**sat_params, redshift = redshift)
    sat_velocity = Velocity('satellites')

    position = halo.position
    velocity = halo.velocity
    halocat = UserSuppliedHaloCatalog(redshift=halo.redshift, Lbox=halo.Lbox, particle_mass=halo.particle_mass, 
                                      halo_x=position[:,0]+halo.Lbox/2.0, halo_y=position[:,1]+halo.Lbox/2.0, halo_z=position[:,2]+halo.Lbox/2.0, 
                                      halo_id=halo.ids, halo_mvir=halo.mass, halo_vx = velocity[:,0], halo_vy = velocity[:,1], 
                                      halo_vz = velocity[:,2], halo_rvir = halo.rvir, halo_nfw_conc= halo.rvir/halo.rs, 
                                      halo_sigmav = halo.sigmav, halo_upid = halo.upid, halo_hostid = halo.pid)
    
    custom_hod_model  = HodModelFactory(centrals_occupation=custom_cens, 
                                    satellites_occupation=custom_sats, 
                                    satellites_profile=custom_phase_space, 
                                    satellites_velocity=sat_velocity,
                                    redshift=redshift)# Populate the halo catalog with galaxiescustom_hod_model.populate_mock(halocat)
    Num_ptcl_requirement = 150
    custom_hod_model.populate_mock(halocat, Num_ptcl_requirement= Num_ptcl_requirement)

    galaxies = custom_hod_model.mock.galaxy_table
    return galaxies

'''

module for HOD to populate halos with galaxies 


'''
#HOD model
from scipy.stats import norm
from halotools.empirical_models import HodModelFactory, OccupationComponent
from halotools.sim_manager import UserSuppliedHaloCatalog
from halotools.empirical_models import NFWPhaseSpace, NFWProfile
import numpy as np
import math
from scipy.special import erf, erfc
from abacusnbody.data.compaso_halo_catalog import CompaSOHaloCatalog
import asdf

#Central galaxy occupation
class CustomCentrals(OccupationComponent):
    def __init__(self, threshold=-23, **kwargs):
        # Call the parent constructor
        super(CustomCentrals, self).__init__(gal_type='centrals', threshold=threshold, upper_occupation_bound=1, **kwargs)

        # Define custom parameters
        self.param_dict = {
            'logMc': 11.89,
            'sigma_M': 0.11,
            #'Ac' : 0.0628,
            'Ac': 0.08,
            'gamma' : 7.06
        }

    def mean_occupation(self, **kwargs):
        # Custom central occupation function (example)
        halo_mass = kwargs['table']['halo_mvir']
        logMc = self.param_dict['logMc']
        sigma_M = self.param_dict['sigma_M']
        Ac = self.param_dict['Ac']
        gamma = self.param_dict['gamma']
        
        # Central occupation based on a Gaussian function
        #return Ac/(2*math.pi)**0.5/sigma_M * np.exp(-(np.log10(halo_mass) - logMc)**2/2/sigma_M**2)*(1+erf(gamma*(np.log10(halo_mass) - logMc)/2/sigma_M**0.5))
        return Ac/(2*math.pi)**0.5/sigma_M * np.exp(-(np.log10(halo_mass) - logMc)**2/2/sigma_M**2)

#Satellite galaxy occupation
class CustomSatellites(OccupationComponent):
    def __init__(self, threshold=-23, **kwargs):
        super(CustomSatellites, self).__init__(gal_type='satellites', threshold=threshold, upper_occupation_bound=float('inf'), **kwargs)

        # Define custom parameters
        self.param_dict = {
            #'As': 0.0064,
            'As': 0.00502,
            'logM0': 11.72,
            'logM1_': 6.12,
            'alpha': -0.31
        }

    def mean_occupation(self, **kwargs):
        halo_mass = kwargs['table']['halo_mvir']
        As = self.param_dict['As']
        logM0 = self.param_dict['logM0']
        logM1_ = self.param_dict['logM1_']
        alpha = self.param_dict['alpha']

        #logM1 = logM1_ + 1/alpha*np.log10(As)
        logM1 = 13
        dist = As*((halo_mass-10**logM0)/10**logM1)**alpha
        
        # Power-law model for satellite occupation
        return np.nan_to_num(np.maximum(dist, 0) * (halo_mass > 10**logM0))

#velocity of satellite galaxies
class Velocity(object):

    def __init__(self, gal_type):

        self.gal_type = gal_type
        self._mock_generation_calling_sequence = ['assign_velocity']
        self._galprop_dtypes_to_allocate = np.dtype([('vx', 'f8'),('vy', 'f8'), ('vz', 'f8')])
        self.list_of_haloprops_needed = ['halo_vx', 'halo_vy', 'halo_vz', 'halo_sigmav']

    def assign_velocity(self, **kwargs):
        table = kwargs['table']
        table['vx'][:] = table['halo_vx']+norm.rvs(loc=0.0, scale=1.0, size=len(table["halo_vz"]))*1.23/3**0.5*table['halo_sigmav']
        table['vy'][:] = table['halo_vy']+norm.rvs(loc=0.0, scale=1.0, size=len(table["halo_vz"]))*1.23/3**0.5*table['halo_sigmav']
        table['vz'][:] = table['halo_vz']+norm.rvs(loc=0.0, scale=1.0, size=len(table["halo_vz"]))*1.23/3**0.5*table['halo_sigmav']

def LoadHalos(path, filename):
    directory = path
    filename = filename
    af = asdf.open(directory+filename)
    cat = CompaSOHaloCatalog(directory, 
                         fields = ['v_com', 'x_com', 'N', 'id', 'r98_com', 'r25_com', 'sigmav3d_com'],
                         cleaned=False,
                         halo_lc=False)
    
    velocity  = cat.halos['v_com'].data                           #center of mass velocity in km/s 
    position = cat.halos['x_com'].data                            #center of mass coordimate in Mpc/h
    mass = cat.halos['N'].data*af['header']['ParticleMassHMsun']  #total mass of DM particles in halos detected with CompaSO in M_solar/h
    ids = cat.halos['id'].data                                    #id of halos
    r98 = cat.halos['r98_com'].data                               #radius that includes 98% of the mass in Mpc/h
    r25 = cat.halos['r25_com'].data                               #radius that includes 25% of the mass in Mpc/h
    sigmav = cat.halos['sigmav3d_com'].data                       #velocity dispersion of DM particles within the halo km/s

    x = position[:,0]
    y= position[:,1]
    z = position[:,2]

    vx = velocity[:,0]
    vy = velocity[:,1]
    vz = velocity[:,2]

    rvir = r98                                                    #We will use r98 as the virial radius
    rs = r25                                                      #We will use r25 as the scale radius
    
    #By default, Halotools can use the halo hierarchy to populate galaxies. 
    #For host halos, halo_upid = -1, whereas for subhalos halo_upid is a long integer equal to the halo_id of the host halo.
    #The halo_pid of an order-N subhalo stores the halo_id of the associated order-N-1 subhalo hosting it.
    #However, CompaSO used in Abacus Summit cannot detect subhalos, so all of the halos in the halo catalog will be the parent halo
    
    upid = np.ones(len(x))*(-1)
    pid = upid
    return x, y, z, vx, vy, vz, rvir, rs, upid, pid

def PopulateGalaxies(path, filename):
    directory = path
    filename = filename
    af = asdf.open(directory+filename)
    redshift = af['header']['Redshift']

    custom_phase_space = NFWPhaseSpace(conc_mass_model='direct_from_halo_catalog', redshift=redshift)
    custom_cens = CustomCentrals(redshift=redshift)
    custom_sats = CustomSatellites(redshift=redshift)
    sat_velocity = Velocity('satellites')

    x, y, z, vx, vy, vz = LoadHalos(path, filename)

    halocat =UserSuppliedHaloCatalog(redshift=redshift, Lbox=Lbox, particle_mass=particle_mass, 
                                     halo_x=x+af['header']['BoxSize']/2.0, halo_y=y+af['header']['BoxSize']/2.0, 
                                     halo_z=z+af['header']['BoxSize']/2.0, halo_id=ids, halo_mvir=mass, halo_upid=upid,
                                     halo_vx = vx, halo_vy = vy, halo_vz =vz, 
                                     halo_hostid = pid, halo_rvir = rvir, halo_nfw_conc= rvir/rs, halo_sigmav = sigmav)
    
    custom_hod_model  = HodModelFactory(centrals_occupation=custom_cens, 
                                    satellites_occupation=custom_sats, 
                                    satellites_profile=custom_phase_space, 
                                    satellites_velocity=sat_velocity,
                                    redshift=redshift)# Populate the halo catalog with galaxiescustom_hod_model.populate_mock(halocat)
    Num_ptcl_requirement = 150
    # Populate the mock galaxy catalog
    custom_hod_model.populate_mock(halocat, Num_ptcl_requirement= Num_ptcl_requirement)

    # Access the galaxy catalog
    galaxies = custom_hod_model.mock.galaxy_table
    return galaxies
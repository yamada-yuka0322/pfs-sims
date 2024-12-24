'''

module to interface with halo catalogs from N-body simulations.

to-do: 
    - [ ] implement halo class object for general halo catalog
    - [ ] implement i/o for existing halo catalogs (abacus, uchuu)
    - [ ] documentation for where to find halo catalogs. 
'''
import os 
import numpy as np 
from astropy.table import Table


class Halos(object):
    def __init__(self):
        """ Flexible class object for general halo catalogs. 
        
        """
        # attributes of the halo catalog
        self.redshift = None
        self.Lbox = None              # in units of Mpc/h
        self.particle_mass = None     # in Msun/h
        # add additional attributes here as we go.
        
        # catalog with halo properties such as position, velocity, mass, virial radius, etc
        # in an astropy.table.Table 
        self.catalog = Table()

    def _check_catalog(self): 
        ''' check that the minimum set of halo catalog attributes and halo properties are included.
        '''
        # check attributes
        if self.redshift is None: raise ValueError("specify redshift of halo catalog")
        if self.Lbox is None: raise ValueError("specify box size (Lbox) in Mpc/h")
        if self.particle_mass is None: raise ValueError("specify particle mass in Msun/h")

        # check catalog
        min_attrs = ['position', 'velocity', 'mass']
        for attr in min_attrs: 
            if attr not in self.catalog.columns(): raise ValueError("%s not specified in catalog" % attr)
        return None


def load_Abacus(path):
    """Load the Abacus halo catalog given path of the asdf data files

    [add links to abacus documentation here]

    
    """
    import asdf 
    try:
        af = asdf.open(path + '/halo_info_000.asdf')
        cat = CompaSOHaloCatalog(path, 
                                 fields = ['v_L2com', 'x_L2com', 'N', 'id', 'r98_L2com', 'r25_L2com', 'sigmav3d_L2com'],
                                 cleaned=True,
                                 halo_lc=False)
    except Exception as e:
        raise ValueError(f"Error loading halo catalog: {e}")

    halos = Halos()
    # set catalog attributes here
    halos.redshift = 
    halos.Lbox = 
    halos.particle_mass = af['header']['ParticleMassHMsun'] # CHECK UNITS OF THIS

    # set halo properties here
    halos.catalog['position'] = cat.halos['x_L2com'].data
    halos.catalog['velocity'] = cat.halos['v_L2com'].data
    halos.catalog['mass']     = cat.halos['N'].data * af['header']['ParticleMassHMsun']
    
    # check that the minimum catalog details are set
    halos._check_catalog()
    
    return halos

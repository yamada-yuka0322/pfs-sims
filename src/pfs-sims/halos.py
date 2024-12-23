'''

module to interface with halo catalogs from N-body simulations.

to-do: 
    - [ ] implement halo class object for general halo catalog
    - [ ] implement i/o for existing halo catalogs (abacus, uchuu)
    - [ ] documentation for where to find halo catalogs. 
'''

class halo():
    def __init__(self, path):
        """
        Initialize the class and load the halo catalog.

        Parameters:
        - path: path to the directory holding the halo_info files.
        """
        self.redshift = None
        self.Lbox = None
        self.particle_mass = None

        self.position = None
        self.velocity = None
        self.rvir = None
        self.rs = None

        self.sigmav = None
        self.mass = None
        self.ids = None

        self.pid = None
        self.upid = None

        self.path = path
        self.load_Abacus()

    def load_Abacus(self):
        """Loads the halo catalog and sets halo properties."""
        try:
            af = asdf.open(self.path + '/halo_info_000.asdf')
            cat = CompaSOHaloCatalog(self.path, 
                                     fields = ['v_L2com', 'x_L2com', 'N', 'id', 'r98_L2com', 'r25_L2com', 'sigmav3d_L2com'],
                                     cleaned=True,
                                     halo_lc=False)
            
            self.velocity  = cat.halos['v_L2com'].data
            self.position = cat.halos['x_L2com'].data
            self.mass = cat.halos['N'].data*af['header']['ParticleMassHMsun']
            self.ids = cat.halos['id'].data
            self.rvir = cat.halos['r98_L2com'].data
            self.rs = cat.halos['r25_L2com'].data
            self.sigmav = cat.halos['sigmav3d_L2com'].data

            self.upid = np.ones(len(self.ids))*(-1)
            self.pid = self.upid
        except Exception as e:
            print(f"Error loading halo catalog: {e}")
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.coordinates as coord
import astropy.units as u

def plot_screens_3D(sources,screens,bounds=[-0.5,0.5,-0.5,0.5,-0.5,0.5],include_legend=False,cm="magma"):
    """
    Plot the 3D locations and velocity vectors of scintillation screens.
    
    Parameters
    ----------
    sources: a pandas data frame or dictionary containing the Galactic longitude, Galactic latitude,
        and distance of the source (typically pulsar)
    screens: a pandas data frame or dictionary containing the source name, distance, orientation,
        and velocity of scintillation screen(s)
    bounds: an array of plotting bounds in kpc, [xmin,xmax,ymin,ymax,zmin,zmax]
    include_legend: a boolean to set whether legend of sources is plotted
    cm : color map to use, either string for named colour map or matplotlib colormap
        
    Returns
    -------
    None
    """
    
    if type(cm) is str:
        cmap = plt.get_cmap("magma")
    else:
        cmap = cm
    slicedCM = cmap(np.linspace(0, 1, len(sources)))
    ax = plt.figure(figsize=(7,7)).add_subplot(projection='3d')
    psr_coords = SkyCoord(l=sources['GL']<<u.deg,b=sources['GB']<<u.deg,distance = sources['DIST']<<u.kpc,frame='galactic')
    i = 0
    for psr_name in sources['NAME']:
        ax.plot3D([0,psr_coords[i].cartesian.x.value],
              [0,psr_coords[i].cartesian.y.value],
              [0,psr_coords[i].cartesian.z.value],
              label=psr_name,color=slicedCM[i])
    
        picrs = psr_coords[i].icrs
   
        scr = screens[screens['PULSAR']==psr_name]
        if np.isnan(scr['v_meas'].values[0]):
            scr_coords = SkyCoord(ra = np.tile(picrs.ra,len(scr)),dec = np.tile(picrs.dec,len(scr)),
                      distance = scr['DIST'].values*u.kpc,
                      frame='icrs').galactic.cartesian
            ax.scatter3D(scr_coords.x,scr_coords.y,scr_coords.z,s=50,color=slicedCM[i])
        else:
            mu_ra = scr['v_meas']*np.sin(scr['theta']*u.deg)/scr['DIST']
            mu_dec = scr['v_meas']*np.cos(scr['theta']*u.deg)/scr['DIST']
            scr_coords = SkyCoord(ra = np.tile(picrs.ra,len(scr)),dec = np.tile(picrs.dec,len(scr)),
                      distance = scr['DIST'].values*u.kpc,
                      pm_ra_cosdec=mu_ra.values*u.km*u.rad/(u.s*u.kpc), 
                      pm_dec = mu_ra.values*u.km*u.rad/(u.s*u.kpc),
                      frame='icrs').galactic.cartesian
            scr_vec = scr_coords.differentials['s']
            ax.quiver(scr_coords.x,scr_coords.y,scr_coords.z,scr_vec.d_x,scr_vec.d_y,scr_vec.d_z,
                  length=0.1,normalize=True,color=slicedCM[i])
            ax.scatter3D(scr_coords.x,scr_coords.y, scr_coords.z, s=50, color=slicedCM[i])
    
        i += 1
    
    ax.set_xlim(bounds[0],bounds[1])
    ax.set_ylim(bounds[2],bounds[3])
    ax.set_zlim(bounds[4],bounds[5])

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    
    if include_legend:
        plt.legend()
        
def summary_plots(sources,screens):
    """
    Plot some histograms and comparisons of screens data
    
    Parameters
    ----------
    sources: a pandas data frame or dictionary containing the distance, number of screens, and dispersion measure
        of radio sources (typically pulsars)
    screens: a pandas data frame or dictionary containing the source name, and fractional screen distances
        of scintillation screen(s)
        
    Returns
    -------
    None
    
    """
    plt.figure(figsize=(7.5,7.5))
    plt.subplot(221)
    plt.hist(screens['s'])
    plt.xlabel('Fractional screen distance, s')
    plt.ylabel('Counts')
    plt.subplot(222)
    plt.hist(sources['NSCREENS'],bins=np.arange(0.5,np.max(sources['NSCREENS'])+1.5,1),rwidth=0.5)
    plt.xlabel('Number of screens')
    plt.ylabel('Counts')
    plt.subplot(223)
    plt.scatter(sources['DM'],sources['NSCREENS'])
    plt.xlabel(r'Dispersion Measure (pc cm$^{-3}$)')
    plt.ylabel('Number of screens')
    plt.yticks(np.arange(1,np.max(sources['NSCREENS'])+1,1))
    plt.subplot(224)
    plt.scatter(sources['DIST'],sources['NSCREENS'])
    plt.xlabel(r'Source Distance, $d_p$ (kpc)')
    plt.ylabel('Number of screens')
    plt.yticks(np.arange(1,np.max(sources['NSCREENS'])+1,1))
    

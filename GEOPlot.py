"""
Plot the positions of the geosynchronous objects over a given range of
time, relative to a field of observation
"""

import argparse as ap
import datetime
from astropy.time import Time
from astropy import units as u
from astropy.coordinates import Longitude, Latitude, EarthLocation
from skyfield.api import load, Topos, utc
from skyfield.sgp4lib import EarthSatellite
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# currently set for La Palma
SITE_LATITUDE = 28.7603135
SITE_LONGITUDE = -17.8796168
SITE_ELEVATION = 2387

SITE_LOCATION = EarthLocation(lat=SITE_LATITUDE*u.deg,
                              lon=SITE_LONGITUDE*u.deg,
                              height=SITE_ELEVATION*u.m)

def argParse():
    """
    Argument parser settings
    """
    
    parser = ap.ArgumentParser()
    
    parser.add_argument('tle_path',
                        help='path to relevant tle file',
                        type=str)
    
    parser.add_argument('out_path',
                        help='path to output directory',
                        type=str)
    
    parser.add_argument('start_utc',
                        help='start of night [utc], format: "yyyy-MM-ddThh:mm:ss"',
                        type=str)
    
    parser.add_argument('field_coords',
                        help='coordinates of tracked field, format: "(hh:mm:ss,dd:mm:ss)"',
                        type=str)
    
    parser.add_argument('field_size',
                        help='coordinates of tracked field, format: "(height_deg,width_deg)"',
                        type=str)
    
    parser.add_argument('n_step',
                        help='number of timesteps to plot',
                        type=int)
    
    parser.add_argument('time_step',
                        help='time resolution of output plots',
                        type=float)
    
    return parser.parse_args()

def parseInput(epoch, field_coords, field_size):
    """
    Read in epoch and requested field coordinates
    """
    
    start_utc = datetime.datetime.strptime(epoch, '%Y-%m-%dT%H:%M:%S')
    start_utc = start_utc.replace(tzinfo=utc)
    
    # unpack coords - avoids argparse issues
    field_ra = Longitude(field_coords.split(',')[0][1:], u.hourangle)
    field_dec = Latitude(field_coords.split(',')[1][:-1], u.deg)
    
    # unpack field dimensions
    field_height = Longitude(field_size.split(',')[0][1:], u.deg)
    field_width = Latitude(field_size.split(',')[1][:-1], u.deg)
    
    return start_utc, field_ra, field_dec, field_height, field_width

def generateEphemeris(name, tle1, tle2, time):
    """
    Generate ephemeris for the given satellite at given epoch
    
    Parameters
    ----------
    name : str
        Name of the satellite (line 0 of 3le)
    tle1 : str
        First line of the satellite's tle (line 1 of 3le)
    tle2 : str
        Second line of the satellite's tle (line 2 of 3le)
    date_str : str
        Desired epoch in format YYYY-mm-ddTHH:MM:SS
    
    Returns
    -------
    ra : 
    
    Raises
    ------
    None
    """
    
    observer = Topos(SITE_LATITUDE, 
                     SITE_LONGITUDE, 
                     elevation_m=SITE_ELEVATION)
    
    target = EarthSatellite(tle1, 
                            tle2, 
                            name)
    
    ra, dec, distance = (target - observer).at(time).radec()

    return [name, ra.hours, dec.degrees]

def getField(ha, dec, height, width):
    """
    Obtain rectangular patch corresponding to field of view
    """
    
    x = ha - width / 2
    y = dec - height / 2
    
    r = Rectangle(xy=(x.hourangle, y.deg),
                  width=width.hourangle,
                  height=height.deg)
    
    r.set_facecolor('none')
    r.set_edgecolor('red')
    
    return r

if __name__ == "__main__":
    
    args = argParse()
    print(args.field_coords)
    
    # load timescale for skyfield
    ts = load.timescale()
    
    # read inputs
    start_utc, field_ra, field_dec, field_height, field_width = parseInput(args.start_utc,
                                                                           args.field_coords,
                                                                           args.field_size)
    
    tle_file = 'geo_3le_' + start_utc.strftime('%Y%m%d') + '.txt'

    count = 0
    while count < args.n_step:
        
        date = start_utc + datetime.timedelta(minutes=count * args.time_step)
        date = date.replace(tzinfo=utc)
        time = ts.utc(date)
        
        print('Processing ', str(date))
        
        name = []
        tle1 = []
        tle2 = []
        with open(args.tle_path + tle_file) as f:
            for line in f:
                line = line.rstrip('\n')
                if line[0] == '0':
                    name.append(line[2:]) # remove line number for name
                elif line[0] == '1':
                    tle1.append(line)
                else:
                    tle2.append(line)

        ha = []
        dec = []
        for (i, j, k) in zip(name, tle1, tle2):
            
            ephemeris = generateEphemeris(i, j, k, time)
            
            ra = Longitude(ephemeris[1], u.hourangle)
            lst = Time(date, scale='utc', location=SITE_LOCATION).sidereal_time('apparent')
            
            ha.append((lst - ra).wrap_at(12*u.hourangle).hourangle)
            dec.append(ephemeris[2])
        
        field_ha = (lst - field_ra).wrap_at(12*u.hourangle)
        
        # obtain field for plot
        field_ax1 = getField(field_ha,
                             field_dec,
                             field_height,
                             field_width)
        
        field_ax2 = getField(field_ha,
                             field_dec,
                             field_height,
                             field_width)
        
        ax1 = plt.subplot2grid((3,1), (0,0), rowspan=2)
        ax2 = plt.subplot2grid((3,1), (2,0))
        
        ax1.plot(ha, dec, 'k.', ms=5)
        
        ax1.add_artist(field_ax1)
        
        ax1.set_title(str(date)) # label plots with utc time
        
        ax1.set_ylabel('Declination / $^\circ$', y=0.25)
        
        ax1.set_xlim(-5,5)
        ax1.set_ylim(-30-3,30+3)
        
        ax1.set_xticks([])

        ax2.plot(ha, dec, 'k.', ms=5)
        
        ax2.add_artist(field_ax2)
        
        ax2.set_xlabel('Hour angle / hr')
        
        ax2.set_xlim(-5,5)
        ax2.set_ylim(-6,-3)
        
        plt.savefig(args.out_path + 'geo_tle_' + str(date) + '.png')
        plt.close()
        
        count += 1

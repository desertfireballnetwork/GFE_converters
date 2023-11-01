# imports
# system packages
import os

# science packages
#import numpy as np
from astropy.table import Table
from astropy.time import Time
from astropy.coordinates import EarthLocation, AltAz, ICRS
import astropy.units as u


def RA_Dec_to_Alt_Az():
    # TODO
    pass


def Alt_Az_to_RA_Dec(az_array,
                     alt_array,
                     datetime_list,
                     obs_latitude,
                     obs_longitude,
                     obs_elevation):
    
    station = EarthLocation(lat=obs_latitude*u.deg,
                        lon=obs_longitude*u.deg,
                        height=obs_elevation*u.meter)
    
    # Datetime is ISO 8601 UTC format
    datetime_array = Time(datetime_list)

    # Convert horizontal to equatorial coordinates
    horizontal_data = AltAz(az_array*u.deg,alt_array*u.deg,
                            obstime=datetime_array,location=station)
    icrs_instance = ICRS()
    equatorial_data = horizontal_data.transform_to(icrs_instance)
    
    ra = equatorial_data.ra.value
    dec = equatorial_data.dec.value
    return ra, dec

def main(ifile):
    t = Table.read(ifile)
    
    ra, dec = Alt_Az_to_RA_Dec(t['azimuth'],
                     t['altitude'],
                     t['datetime'],
                     t.meta['obs_latitude'],
                     t.meta['obs_longitude'],
                     t.meta['obs_elevation'])
    
    t['ra'] = ra
    t['dec'] = dec
    
    for c in ['JD']:
        if c in t.colnames:
            t.remove_column(c)
    
    # Fix METADATA
    t.meta['camera_id'] = t.meta['dfn_camera_codename']
    t.meta['cx'] = t.meta['NAXIS1']
    t.meta['cy'] = t.meta['NAXIS2']

    t.meta['obs_az'] = 0.
    t.meta['obs_ev'] = 90.
    t.meta['fov_vert'] = 90.
    t.meta['fov_horiz'] = 360.

    t.meta['file_standard'] = 'GFE_1.2'
    
    for k in ['NAXIS1', 'NAXIS2', 'dfn_camera_codename']:
        if k in t.meta:
            t.meta.pop(k)
            
            
    o_fname = '_'.join([(t['datetime'][0][:19]).replace(':','_'),'DFN', t.meta['location'].replace('_','')]) + '.ecsv'
    o_fname = os.path.join(os.path.dirname(ifile), o_fname)
    
    
    # write the file to disk
    t.write(o_fname, format='ascii.ecsv', delimiter=',', overwrite=False)
    print(f'table has been written to {o_fname}')

if __name__ == "__main__":
    import argparse
    # parse arguments
    parser = argparse.ArgumentParser(description='Convert DFN astrometry ECSV to GFE ECSV standard')
    parser.add_argument("-i", "--ifile", type=str, required=True, help="input filename")
    
    args = parser.parse_args()
    
    main(args.ifile)


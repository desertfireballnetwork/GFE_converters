
# system packages
import os

#file handlers
import xmltodict  #XML  - for UFOAnalyzer

# date handling
from datetime import datetime, timedelta

# numerical packages
import numpy as np
import astropy.units as u
from astropy.table import Table
from astropy.time import Time

# definitions of constants:
ISO_FORMAT = "%Y-%m-%dT%H:%M:%S.%f"  #defines a consistent iso date format

def isoStr(iso_datetime_string):
     return datetime.strptime(iso_datetime_string, ISO_FORMAT)

def ufo_to_std(ufo_1):
    # This function takes a UFO file, passed as a nested dictionary, and returns a table in DFN format.
    # UFO format is a deeply nested XML file. The nesting is:
    # The whole XML is in the dictionary ufo_1
    # ufoanalyzer_record  "  "  ufo_2 - a dictionary of station data
    # ua2_objects         "  "  ufo_3 - intermediate, suppressed
    # ua2_object          "  "  ufo_4 - a dictionary of observation metadata
    # ua2_objpath         "  "  ufo_5 - intermediate, suppressed
    # ua2_fdata2          "  "  ufo_6 - the dictionary of trajectory data
    
    # Note on UFO capture algorithm:
    # Assuming head=30 and the video is interlaced 25 fps (so 
    # effectively 50 fps), the capture algorithm seems to be:
    # 1. Event detected at time X.  This is used as the 
    #    timestamp and is recorded in the file.
    # 2. Save the framestack from time X plus 30 full 
    #    frames (60 interlaced half-frames) beforehand
    # 3. Now treat each of your half-frames as frames.  So 
    #    time X is frame 61.
    # 4. Examine each frame from frame 1 to the end of the 
    #    frame stack to see whether the event started earlier 
    #    or later than you thought.  
    # 5. Rather than frame 61, it can sometimes be frame 54, 
    #    or 59, or 64 when the first real event is detected.  
    #    Save this as “fs”.
    # 6. List all of the frames where you think you know what 
    #    happened, starting with fno=fs, and skipping frames 
    #    that can’t be analysed.  
        
     
    ttt_list = []
    ufo_2=ufo_1['ufoanalyzer_record']

    ufo_4=ufo_2['ua2_objects']['ua2_object']
    meteor_count = len(ufo_4)
    if meteor_count > 10 :
        # if ufo4 has 59 elements then it's a single meteor but the data is less nested
        meteor_count = 1

    # Now get metadata from ufo_2
    
    #location
    obs_latitude = float(ufo_2['@lat'])
    obs_longitude = float(ufo_2['@lng'])
    obs_elevation = float(ufo_2['@alt'])
    
    #camera station site name
    origin = "UFOAnalyzer" + '_Ver_'+ ufo_2['@u2']  # or other formal network names
    location = ufo_2['@lid']
    telescope = ufo_2['@sid']     #no spaces or special characters
    camera_id = location + '_' + telescope

    #observer and instrument
    observer = ufo_2['@observer']
    instrument = ufo_2['@cam']
    comment = ufo_2['@memo']
    cx = int(ufo_2['@cx'])
    cy = int(ufo_2['@cy'])

    image_file = ufo_2['@clip_name']+'.AVI'
    astrometry_number_stars = int(float(ufo_2['@rstar']))
    lens = ufo_2['@lens']
 
    # calculate event timings - file timestamp
    timestamp_str = ufo_2['@y'] + '-' + ufo_2['@mo'] + '-' + ufo_2['@d']
    timestamp_str += 'T' + ufo_2['@h'] + ':' + ufo_2['@m'] + ':' + ufo_2['@s']
    timestamp = Time(timestamp_str)

    # frame rate and beginning and middle of clip
    multiplier = 1 + int(ufo_2['@interlaced'])
    head = int(ufo_2['@head']) * multiplier
    tail = int(ufo_2['@tail']) * multiplier
    frame_rate = float(ufo_2['@fps']) * multiplier
    

    # now loop through each meteor
    
    for k in range(meteor_count):
    
        if meteor_count == 1 :
            ufo_5 = ufo_4
        else:
            ufo_5 = ufo_4[k]
        
        sec = float(ufo_5['@sec'])
        nlines = int(ufo_5['@sN'])
        ufo_6=ufo_5['ua2_objpath']['ua2_fdata2']
    
        no_frames = int(ufo_5['@fN'])+ head + tail
        fs = int(ufo_5['@fs'])
        exposure_time = (no_frames-1.0)/frame_rate
        AVI_start_sec = -float(head)/frame_rate  
        AVI_mid_sec =  exposure_time * 0.5 + AVI_start_sec
        AVI_start_time = str(timestamp + timedelta(seconds=AVI_start_sec))
        AVI_mid_time = str(timestamp + timedelta(seconds=AVI_mid_sec))  
        timestamp_frame = head + 1      # The file timestamp is the first frame after the "head"
        exposure_time = (no_frames - 1.0) / frame_rate 
    
        fov_vert = 0.0
        fov_horiz =float(ufo_2['@vx'])
        if cx > 0:
            fov_vert = fov_horiz * cy / cx
    
        # construction of the metadata dictionary
        meta_dic = {'obs_latitude': obs_latitude,
               'obs_longitude': obs_longitude,
               'obs_elevation': obs_elevation,
               'origin': origin,
               'location': location,
               'telescope': telescope,
               'camera_id': camera_id,
               'observer': observer,
               'comment': comment,
               'instrument': instrument,
               'lens': lens,
               'cx' : cx,     
               'cy' : cy,     
               'photometric_band' : 'Unknown',     
               'image_file' : image_file,
               'isodate_start_obs': AVI_start_time,   
               'isodate_calib': AVI_mid_time,   
               'exposure_time': exposure_time,   
               'astrometry_number_stars' : astrometry_number_stars,
               # 'photometric_zero_point': float(ufo_2['@mimMag']),    
               # 'photometric_zero_point_uncertainty': 0.0,
               'mag_label': 'mag',     
               'no_frags': 1,
               'obs_az': float(ufo_2['@az']),     
               'obs_ev': float(ufo_2['@ev']),     
               'obs_rot': float(ufo_2['@rot']),     
               'fov_horiz': fov_horiz,     
               'fov_vert': fov_vert,     
               }

        # initialise table
        ttt = Table()        
        #Update the table metadata
        ttt.meta.update(meta_dic)   
   
        #create time and main data arrays
        # Datetime is ISO 8601 UTC format
        datetime_array = []
    
        # Azimuth are East of North, in degrees
        azimuth_array = []
        # Altitudes are geometric (not apparent) angles above the horizon, in degrees
        altitude_array = []

        # Magnitude
        mag_array = []
 
        # Right Ascension / Declination coordinates read from file
        ra_array  = []
        dec_array = []
        
        N_records = len(ufo_6)
        print(f'{N_records} records found')
        for i in range(N_records):
            obs=ufo_6[i]
            az   = float(obs['@az'])
            elev = float(obs['@ev'])
            ra   = float(obs['@ra'])
            dec = float(obs['@dec'])
            mag = float(obs['@mag'])
            obs_time = (int(obs['@fno']) - timestamp_frame)/frame_rate
            time_stamp = str(timestamp + timedelta(seconds=obs_time))  
            azimuth_array.append(az)
            altitude_array.append(elev)
            ra_array.append(ra)
            dec_array.append(dec)
            mag_array.append(mag)
            datetime_array.append(time_stamp)
        
        ## Populate the table with the data created to date
        # create columns
        ttt['datetime'] = datetime_array
        ttt['ra']  = ra_array  * u.degree
        ttt['dec'] = dec_array * u.degree 
        ttt['azimuth'] = ((np.array(azimuth_array) + 180)%360) * u.degree
        ttt['altitude'] = altitude_array * u.degree    
        ttt['mag'] = mag_array 
        ttt['x_image'] = 0.0 
        ttt['y_image'] = 0.0
        
        # add some uncertainties, to be compatible with DFN pipeline
        DEF_AZ_ERR = 2./60
        DEF_ALT_ERR = 2./60
        ttt['err_minus_altitude'] = DEF_ALT_ERR
        ttt['err_plus_altitude'] = DEF_ALT_ERR
        ttt['err_minus_azimuth'] = DEF_AZ_ERR
        ttt['err_plus_azimuth'] = DEF_AZ_ERR
        
        DEF_TIME_ERR = 3./1000
        ttt['time_err_minus'] = DEF_TIME_ERR
        ttt['time_err_plus'] = DEF_TIME_ERR

        
        # if RA/Dec values are all zeros, then convert alt-az to Ra-Dec
        if not np.any(ttt['ra']):
            print('WARNING: RA/Dec values are all zeros, converting Alt-Az to RA-Dec, ignoring existing RA-Dec values')
            from DFN_to_GFE_conversion import Alt_Az_to_RA_Dec
            ra, dec = Alt_Az_to_RA_Dec(azimuth_array,
                                       altitude_array,
                                       datetime_array,
                                       obs_latitude,
                                        obs_longitude,
                                        obs_elevation)
            ttt['ra']  = ra  * u.degree
            ttt['dec'] = dec * u.degree 
            
    
        # now add ttt to the array of tables
        ttt_list.append(ttt)
    
    
    return(ttt_list, meteor_count);



def main(fname):
    # input is UFOAnalyzer.  
    print("UFO format being read")
    with open(fname) as fd:
        _obs_dic=xmltodict.parse(fd.read())
    ttt_list, meteor_count = ufo_to_std(_obs_dic)
    if meteor_count > 1:
        print('WARNING: more than one meteor in file, outputing just the first one')
    ttt = ttt_list[0]
    out_name = fname.replace('.xml','.ecsv')
    
    ttt.write(out_name, format='ascii.ecsv', delimiter=',', overwrite=False)
    print(f'table has been written to {out_name}')

if __name__ == "__main__":
    import argparse
    # parse arguments
    parser = argparse.ArgumentParser(description='Convert AMOS (UFO Analyzer) astrometry ECSV to GFE ECSV standard')
    parser.add_argument("-i", "--ifile", type=str, required=True, help="input filename")
    
    args = parser.parse_args()
    
    main(args.ifile)

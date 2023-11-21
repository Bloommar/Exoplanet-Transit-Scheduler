#importing libraries and tools

import urllib
import json
import pandas as pd
from astropy.time import Time
from astropy.coordinates import solar_system_ephemeris,  get_body_barycentric, EarthLocation
from astropy import units as u
from astropy.coordinates import SkyCoord, AltAz
import astropy.units as u
from astropy.coordinates import get_sun
import numpy as np

!pip install astroplan

from astroplan import Observer, FixedTarget
from astropy.coordinates import EarthLocation

#downloading the data from exoclock
exoclock_planets= pd.read_json("https://www.exoclock.space/database/planets_json")

exoclock_planets.to_csv("all_data.csv",header=True, index=True)


#creating data frame
df=exoclock_planets
df = pd.read_csv("all_data.csv", index_col=0, header=1)

#initialising the telescope characteristics 
min_aperture= 8.0
max_aperture=45.0
r_mag_min=8.0
telescope_latitude=-30.526309016377613
telescope_longitude=-70.85329602458853

#the extra time for before and after transit start in JD, here is one hour in JD
extra_time=0.04167000018



def sexagesimal_to_decimal(sexagesimal_str):



  #
    #Converts sexagesimal coordinates to decimal coordinates.

   # Parameters:
   # - sexagesimal_str (str): A string in the format 'degrees:minutes:seconds'.

   # Returns:
   # - str: A string representing the decimal coordinates with 5 decimal places.
   # 
   
    degrees, minutes, seconds = map(float, sexagesimal_str.split(':'))

    decimal_declination = degrees + minutes/60 + seconds/3600

    return '{:.5f}'.format(decimal_declination)


row_to_convert = df.iloc[56]
decimal_declinations = row_to_convert.apply(sexagesimal_to_decimal)

df.iloc[56] = decimal_declinations

df.to_csv("all_data.csv", header=1, index=True)


def sexagesimal_to_decimal(sexagesimal_str):


  #
    #Converts sexagesimal coordinates to decimal coordinates.

   # Parameters:
   # - sexagesimal_str (str): A string in the format 'degrees:minutes:seconds'.

   # Returns:
   # - str: A string representing the decimal coordinates with 5 decimal places.
   # 
   
    degrees, minutes, seconds = map(float, sexagesimal_str.split(':'))

    decimal_ra = degrees + minutes/60 + seconds/3600

    return '{:.5f}'.format(decimal_ra)


row_to_convert = df.iloc[55]
decimal_ra = row_to_convert.apply(sexagesimal_to_decimal)

df.iloc[55] = decimal_ra

df.to_csv("all_data.csv", header=1, index=True)


#transposing the dataframe
transposed_df = df.T
transposed_df.dtypes

#making columns float numbers instead of objects
columns_to_convert = ['depth_r_mmag','duration_hours', 'ephem_period','ra_j2000','dec_j2000','r_mag','v_mag','min_telescope_inches','ephem_mid_time','t0_bjd_tdb','period_days']  # Replace with the column names you want to convert

transposed_df[columns_to_convert] = transposed_df[columns_to_convert].astype(float)

transposed_df.min_telescope_inches.dtype


#first criteria needed depending on your telescope and needs
new_df= transposed_df[
    ((transposed_df["priority"] == "high") | (transposed_df["priority"] == "medium") | (transposed_df["priority"] == "alert")) &
    (transposed_df["min_telescope_inches"] >= min_aperture) & (transposed_df["min_telescope_inches"] <= max_aperture) & (transposed_df["r_mag"] >= r_mag_min)
]
new_df=transposed_df

new_df.to_csv("all_data.csv", header=1, index=True)

#new column order , keeping what we need
new_column_order = ['priority', 'current_oc_min','depth_r_mmag','star','duration_hours', 'ephem_period', 'ephem_period_units','ra_j2000','dec_j2000','r_mag','v_mag','min_telescope_inches','ephem_mid_time','t0_bjd_tdb','period_days']
new_df_2 = new_df[new_column_order]


#put when do you want to start observations 
date_time_start = "2023-11-20T12:00:00"

time_start = Time(date_time_start, format='isot', scale='utc')

with solar_system_ephemeris.set('builtin'):
    tdb_start = time_start.tdb.jd


    #put when you want to finish observations
date_time_end = "2023-11-21T10:00:00"
time_end = Time(date_time_end, format='isot', scale='utc')

with solar_system_ephemeris.set('builtin'):
    tdb_end = time_end.tdb.jd



#calculating the total observation period we want
total_period_of_observations=tdb_end-tdb_start


#calculating the current epoch
new_df_2['current_epoch']=(tdb_start-new_df_2.t0_bjd_tdb)/new_df_2.ephem_period


#calculating the next eclipses from current eppoch until the time we want to stop observations
def calc_next_eclipses(row, tdb_end_1):
    output_list = []
    m = 0
    while True:
        m = m + 1
        next_eclipse = row['t0_bjd_tdb'] + row['ephem_period'] * (int(row['current_epoch']) + m)
        if next_eclipse > tdb_end_1:
            break
        output_list.append(next_eclipse)
    return output_list
new_df_2['next_eclipse_times'] = new_df_2.apply(lambda row: calc_next_eclipses(row, tdb_end), axis=1)
new_df_2 = new_df_2[new_df_2['next_eclipse_times'].apply(lambda x: bool(x))]
new_df_2


# Define your observer's location using EarthLocation
observer_location = EarthLocation.from_geodetic(-70.85329602458853, -30.526309016377613, 1700.0)
observer = Observer(location=observer_location, name="My Observatory", timezone="UTC")

observable_targets_set = set()

for index, row in new_df_2.iterrows():
    # Round the coordinates to a certain precision
    ra_rounded = round(row['ra_j2000'], 5)
    dec_rounded = round(row['dec_j2000'], 5)

    # Extract target coordinates from the DataFrame
    target_coordinates = SkyCoord(ra=ra_rounded * u.deg, dec=dec_rounded * u.deg, frame='icrs')
    
    # Specify the time of observation
    observation_time = Time(row['t0_bjd_tdb'], format='jd', scale='tdb')

    # Check if the target is observable at the specified time
    is_observable = observer.target_is_up(observation_time, target_coordinates)

    if is_observable:
        # Add the target name to the set
        observable_targets_set.add(row['star'])

# Filter the DataFrame based on the set of unique observable targets
new_df_2 = new_df_2[new_df_2['star'].isin(observable_targets_set)]



#extra time in JD , 1 hour long so we can add/substract it form the timestamp
extra_time=0.04167000018
#calculating the duration days from the duration hours
new_df_2['duration_days']=new_df_2.duration_hours/24



#Calculating the start and end of each eclipse

def start_eclipse(row, extra_time):
    if not row['next_eclipse_times']:
        return None
    start_of_transit = [ecl - extra_time - 0.5 * row['duration_days'] for ecl in row['next_eclipse_times']]
    return start_of_transit

new_df_2['start_of_transit'] = new_df_2.apply(start_eclipse, extra_time=extra_time, axis=1)

new_df_2 = new_df_2.dropna(subset=['start_of_transit'])

if 'name' in new_df_2.columns:
    new_df_2.set_index('name', header=1, index=True)


def end_eclipse(row, extra_time):
    if not row['next_eclipse_times']:
        return None
    end_of_transit = [ecl + extra_time + 0.5 * row['duration_days'] for ecl in row['next_eclipse_times']]
    return end_of_transit

new_df_2['end_of_transit'] = new_df_2.apply(end_eclipse, extra_time=extra_time, axis=1)

new_df_2 = new_df_2.dropna(subset=['end_of_transit'])

if 'name' in new_df_2.columns:
    new_df_2.set_index('name', header=1, index=True)



#calculating the start altitude of the star based on your telescope coordinates
altitude_list = []

for index, row in new_df_2.iterrows():
    time_jd = Time(row['start_of_transit'], format='jd')
    bear_mountain = EarthLocation(lat=-30.526309016377613*u.deg, lon=-70.85329602458853*u.deg, height=1700*u.m)

    altitude = SkyCoord(ra=row['ra_j2000'], dec=row['dec_j2000'], unit=(u.hourangle, u.deg))
    altitude = SkyCoord(ra=row['ra_j2000'], dec=row['dec_j2000'], unit=(u.hourangle, u.deg)) \
        .transform_to(AltAz(obstime=time_jd, location=bear_mountain)).alt


    altitude_float = altitude.value
    altitude_list.append(altitude_float)

# Assign the list of altitudes to the DataFrame
new_df_2['start_altitude'] = altitude_list


#calculating the end altitude of the star based on your telescope coordinates

altitude_list_2 = []

for index, row in new_df_2.iterrows():
    time_jd = Time(row['end_of_transit'], format='jd')
    bear_mountain = EarthLocation(lat=-30.526309016377613*u.deg, lon=-70.85329602458853*u.deg, height=1700*u.m)

    altitude = SkyCoord(ra=row['ra_j2000'], dec=row['dec_j2000'], unit=(u.hourangle, u.deg))
    altitude = SkyCoord(ra=row['ra_j2000'], dec=row['dec_j2000'], unit=(u.hourangle, u.deg)) \
        .transform_to(AltAz(obstime=time_jd, location=bear_mountain)).alt


    altitude_float = altitude.value
    altitude_list_2.append(altitude_float)

# Assign the list of altitudes to the DataFrame
new_df_2['end_altitude'] = altitude_list_2


#getting the sun coordinates 
observatory_location = EarthLocation(lat = -30.526309016377613 * u.deg, lon = - 70.85329602458853 * u.deg, height=1700*u.m)

time = Time.now()


sun_coords = get_sun(time)

# Transform the Sun's coordinates to the AltAz frame for your observatory
sun_altaz = sun_coords.transform_to(AltAz(obstime=time, location=observatory_location))

# Access the altitude of the Sun
sun_altitude = sun_altaz.alt

#we want the twillight to be under -10 degrees 20 minutes before and after the transit
extra_time_for_twillight=0.01389000006  #extra time in JD we will calculate the sun altitude


#calculation where the sun will be before and after the transit
def twilight_before_transit(row, extra_time):
    # Initialize an empty list to store sun altitudes
    sun_altitudes_before = []

    # Iterate over each start_of_transit time in the list
    for transit_time_jd in row['start_of_transit']:
        # Convert the start time of transit to Time object
        time_jd = Time(transit_time_jd, format='jd')

        # Define the observatory location (replace with your observatory's coordinates)
        observatory_location = EarthLocation(lat=-30.526309016377613*u.deg, lon=-70.85329602458853*u.deg, height=1700*u.m)

        # Calculate the Sun's coordinates at the specified time
        sun_coords = get_sun(time_jd)

        # Transform the Sun's coordinates to the AltAz frame for your observatory
        sun_altaz = sun_coords.transform_to(AltAz(obstime=time_jd, location=observatory_location))

        # Access the altitude of the Sun and append to the list
        sun_altitudes_before.append(sun_altaz.alt.value)

    return sun_altitudes_before

new_df_2['sun_altitude_before'] = new_df_2.apply(twilight_before_transit, extra_time=extra_time, axis=1)




def twilight_after_transit(row, extra_time):
    # Initialize an empty list to store sun altitudes
    sun_altitudes_after = []

    # Iterate over each end_of_transit time in the list
    for transit_time_jd in row['end_of_transit']:
        # Convert the end time of transit to Time object
        time_jd = Time(transit_time_jd, format='jd')

        # Define the observatory location (replace with your observatory's coordinates)
        observatory_location = EarthLocation(lat=-30.526309016377613*u.deg, lon=-70.85329602458853*u.deg, height=1700*u.m)

        # Calculate the Sun's coordinates at the specified time
        sun_coords = get_sun(time_jd)

        # Transform the Sun's coordinates to the AltAz frame for your observatory
        sun_altaz = sun_coords.transform_to(AltAz(obstime=time_jd, location=observatory_location))

        # Access the altitude of the Sun and append to the list
        sun_altitudes_after.append(sun_altaz.alt.value)

    return sun_altitudes_after

new_df_2['sun_altitude_after'] = new_df_2.apply(twilight_after_transit, extra_time=extra_time, axis=1)



#in another code we have fitted the data we have [ra, dec,flux,mag] in an exponential function and the parameters are A,B for the 9 mag stars and C,B for the 11 mag
A=366615682.0727282
B=-0.8126434038858498
C=220555409.81953186
D=-0.899609043634542

data2 = pd.read_csv('/content/9mag_recalculated.csv')
data2.columns
df2 = pd.DataFrame(data2)

data3 = pd.read_csv('/content/11mag_recalculated.csv')
data3.columns
df3 = pd.DataFrame(data3)


#we will calculate the exposure time for 40000 counts 


def exposure_time(new_df_2, df2, df3, A, B, C, D):
    # Assuming r_mag is a column in new_df_2
    r_mag = new_df_2['r_mag']

    # Define the exponential functions
    def exponential_function(x, a, b):
        return a * np.exp(b * x)

    def exponential_function2(x, c, d):
        return c * np.exp(d * x)

    # Calculate exposure time based on conditions
    new_df_2["exp_time"] = np.where(r_mag <= 11.0, (10 * 40000) / exponential_function(r_mag, A, B),
                                    (60 * 40000) / exponential_function2(r_mag, C, D))

    # Return the DataFrame with the added exposure_time column
    return new_df_2

# Apply the function
new_df_2 = exposure_time(new_df_2, df2, df3, A, B, C, D)



#we extract the lists we have created 

def extract_lists(row, columns_with_lists):
    extracted_data = {}
    for col in columns_with_lists:
        extracted_data[col] = row[col]
    return pd.Series(extracted_data)

# Specify the columns in 'new_df_2' that have lists and those you want to keep
columns_with_lists = ['ra_j2000',"current_oc_min",'dec_j2000','exp_time','star','next_eclipse_times', 'start_of_transit', 'end_of_transit', 'start_altitude', 'end_altitude','r_mag' ,	'sun_altitude_before',	'sun_altitude_after']

# Apply the function to create 'new_df_3'
new_df_3 = new_df_2.apply(extract_lists, columns_with_lists=columns_with_lists, axis=1)



#we create a new dictionary with data from the existing dfs by itterating each row, each row in df is a new row 
#also we keep the columns we need

data = []


for index, row in new_df_3.iterrows():
    exploded_start_altitude = pd.Series(row['start_altitude']).explode() if isinstance(row['start_altitude'], list) else pd.Series([row['start_altitude']])
    exploded_end_altitude = pd.Series(row['end_altitude']).explode() if isinstance(row['end_altitude'], list) else pd.Series([row['end_altitude']])

    max_length = max(len(row["start_of_transit"]), len(exploded_start_altitude), len(row["end_of_transit"]))

    for i in range(max_length):
        data.append(
            {
                "name": row.name,
                "T_0": row['next_eclipse_times'][i] if isinstance(row['next_eclipse_times'], list) and len(row['next_eclipse_times']) > i else None,
                "start": row['start_of_transit'][i] if isinstance(row['start_of_transit'], list) and len(row['start_of_transit']) > i else None,
                "end": row['end_of_transit'][i] if isinstance(row['end_of_transit'], list) and len(row['end_of_transit']) > i else None,
                "start_altitude": row["start_altitude"][i] ,
                "end_altitude": row["end_altitude"][i],
                "sun_altitudes_before": row['sun_altitude_before'][i] if isinstance(row['sun_altitude_before'], list) and len(row['sun_altitude_before']) > i else None,
                "sun_altitudes_after": row['sun_altitude_after'][i] if isinstance(row['sun_altitude_after'], list) and len(row['sun_altitude_after']) > i else None,
                "ra_j2000": row.ra_j2000,
                "dec_j2000" : row.dec_j2000,
                "r_mag" : row.r_mag,
                "exp_time": row.exp_time,
                "current_oc_min" : row.current_oc_min
            }
        )

#putting data the df

new_df_3=pd.DataFrame.from_dict(data)


#filtering the data by the sun altitude 
#we want it under -10 degrees which is over nautical twillight but the sun light won't bother our our observations

new_df_3=new_df_3[
    ((new_df_3['sun_altitudes_before'] <= -10.0) & (new_df_3['sun_altitudes_after']<=-10.0))
]

#introducing a function that will convert the times to utc, 'cause afterwards we want to filter it 
def bjd_to_utc(bjd_value):
    bjd_time = Time(bjd_value, format='jd', scale='tdb')
    return bjd_time.utc.datetime

new_df_3['T_0_utc'] = new_df_3['T_0'].apply(bjd_to_utc)
new_df_3['start_utc'] = new_df_3['start'].apply(bjd_to_utc)
new_df_3['end_utc'] = new_df_3['end'].apply(bjd_to_utc)


#now we change the time to strftime so it's more usable

new_df_3['start_utc'] = pd.to_datetime(new_df_3['start_utc']).dt.strftime('%Y-%m-%d %H:%M:%S')

new_df_3['end_utc'] = pd.to_datetime(new_df_3['end_utc']).dt.strftime('%Y-%m-%d %H:%M:%S')




# Convert 'start_utc' and 'end_utc' columns to datetime objects
new_df_3['start_utc'] = pd.to_datetime(new_df_3['start_utc'])
new_df_3['end_utc'] = pd.to_datetime(new_df_3['end_utc'])

# Convert 'current_oc_min' to numeric and update 'start_utc' based on the condition
update_start_condition = pd.to_numeric(new_df_3['current_oc_min'], errors='coerce').abs() <= 15
new_df_3.loc[update_start_condition, 'start_utc'] += pd.to_timedelta(abs(pd.to_numeric(new_df_3.loc[update_start_condition, 'current_oc_min'], errors='coerce')), unit='m')

# Convert 'current_oc_min' to numeric and update 'end_utc' based on the condition
update_end_condition = pd.to_numeric(new_df_3['current_oc_min'], errors='coerce').abs() >= 15
new_df_3.loc[update_end_condition, 'end_utc'] += pd.to_timedelta(abs(pd.to_numeric(new_df_3.loc[update_end_condition, 'current_oc_min'], errors='coerce')), unit='m')

# Round 'start_utc' and 'end_utc' to the nearest 15 minutes and convert to string
new_df_3['obs_start_time'] = new_df_3['start_utc'].round('15min').astype(str)
new_df_3['obs_end_time'] = new_df_3['end_utc'].round('15min').astype(str)

# Split the strings in 'obs_start_time' and 'obs_end_time' to extract date and time components
new_df_3[['real_obs_date', 'obs_start_time']] = new_df_3['obs_start_time'].str.split(' ', expand=True)
new_df_3[['end_date_delete', 'obs_end_time']] = new_df_3['obs_end_time'].str.split(' ', expand=True)



new_df_3['start_utc'] = pd.to_datetime(new_df_3['start_utc'])

# Subtract one day from 'start_utc' if the hour is between 0 and 11, else keep the original value
new_df_3['obs_night'] = new_df_3['start_utc'].apply(lambda row_dt: row_dt - pd.Timedelta(days=1) if 0 <= row_dt.hour < 12 else row_dt)

# Convert 'obs_night' to string and split it into separate columns
new_df_3[['obs_night', 'obs_time_delete']] = new_df_3['obs_night'].astype(str).str.split(' ', expand=True)


#creating new columns for the final file for our project

new_df_3['FILTER']='R'
new_df_3['Focus']=11850
new_df_3['Binning']=2
new_df_3['GAIN']=0

new_df_3['project']='ExoClock UZ'

new_column_order = ['obs_night','obs_start_time', 'obs_end_time', 'name', 'project','ra_j2000','dec_j2000','FILTER','exp_time','Focus','Binning','GAIN']
new_df_3 = new_df_3[new_column_order]


#sorting the data by day and time
new_df_3 = new_df_3.sort_values(by=['obs_night', 'obs_start_time'], ascending=[True, True])
#storing them in the final file
new_df_3.to_csv("all_data_final.csv", header=1, index=True)

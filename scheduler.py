

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

URL = "https://www.exoclock.space/database/planets_json"
MIN_APERTURE= 8.0 # proposed minimum aperture
MAX_APERTURE=45.0 # proposed maximum aperture
R_MAG_MIN=8.0 # minimum R filter magnitude
telescope_LATITUDE=-30.526309016377613 # telescope latitude
telescope_longitude=-70.85329602458853 # telescope longitude
PRIORITY_criteria = ('high', 'alert', 'medium')  #criteria for the kind of priority objects you want to observe
columns_to_convert = ['depth_r_mmag','duration_hours', 'current_oc_min','ephem_period','ra_j2000','dec_j2000','r_mag','v_mag','min_telescope_inches','ephem_mid_time','t0_bjd_tdb','period_days']  # Replace with the column names you want to convert


date_time_start = "2023-11-25T12:00:00" # date and time you want to start your observations
date_time_end = "2023-11-30T10:00:00" # date and time you want to end your observations



observer_location = EarthLocation.from_geodetic(-70.85329602458853, -30.526309016377613, 1700.0) #exact location of the telescope
observatory_location = EarthLocation(lat = -30.526309016377613 * u.deg, lon = - 70.85329602458853 * u.deg, height=1700*u.m)

observer = Observer(location=observer_location, name="My Observatory", timezone="UTC") # use your preffered timezone
extra_time=0.04167000018 # extra time in JD , 1 hour long so we can add/substract it form the timestamp
extra_time_for_twillight=0.01389000006 # extra time in JD we will calculate the sun altitude

#in another code we have fitted the data we have [ra, dec,flux,mag] in an exponential function and the parameters are A,B for the 9 mag stars and C,B for the 11 mag
A = 366615682.0727282
B = -0.8126434038858498
C=220555409.81953186
D=-0.899609043634542

#columns that we want to keep and make them float
columns_to_convert = ['depth_r_mmag','duration_hours', 'ephem_period','ra_j2000','dec_j2000','r_mag','v_mag','min_telescope_inches','ephem_mid_time','t0_bjd_tdb','period_days']  # Replace with the column names you want to convert
new_column_order = ['priority', 'current_oc_min','depth_r_mmag','star','duration_hours', 'ephem_period', 'ephem_period_units','ra_j2000','dec_j2000','r_mag','v_mag','min_telescope_inches','ephem_mid_time','t0_bjd_tdb','period_days']

# Specify the columns in 'new_df' that have lists and those you want to keep
columns_with_lists = ['ra_j2000',"current_oc_min",'dec_j2000','exp_time','star','next_eclipse_times', 'start_of_transit', 'end_of_transit', 'start_altitude', 'end_altitude','r_mag' ,	'sun_altitude_before',	'sun_altitude_after']

time_start = Time(date_time_start, format='isot', scale='utc')  # defining the time in UTC
time_end = Time(date_time_end, format='isot', scale='utc')
with solar_system_ephemeris.set('builtin'):
    tdb_start = time_start.tdb.jd
    tdb_end = time_end.tdb.jd

def get_data(url):
  data = pd.read_json(url)
  data = data.T.reset_index(drop=True)

  return data



def sexagesimal_to_decimal(sexagesimal_str):
    #Converts sexagesimal dec to decimal coordinates.

   # Parameters:
   # - sexagesimal_str (str): A string in the format 'degrees:minutes:seconds'.

   # Returns:
   # - str: A string representing the decimal coordinates with 7 decimal places.


    degrees, minutes, seconds = map(float, sexagesimal_str.split(':'))

    decimal_coordinates = degrees + minutes/60 + seconds/3600

    return '{:.7f}'.format(decimal_coordinates)

def calc_epoch(new_df, date_time_start, date_time_end):
    time_start = Time(date_time_start, format='isot', scale='utc')  # defining the time in UTC
    time_end = Time(date_time_end, format='isot', scale='utc')

    with solar_system_ephemeris.set('builtin'):
        tdb_start = time_start.tdb.jd
        tdb_end = time_end.tdb.jd

    # calculating the total observation period we want
    total_period_of_observations = tdb_end - tdb_start

    # calculating the current epoch
    new_df['current_epoch'] = (tdb_start - new_df['t0_bjd_tdb'].astype(float)) / new_df['ephem_period'].astype(float)

    return new_df

#calculating the next eclipses from current eppoch until the time we want to stop observations
def calc_next_eclipses(row, tdb_end):
    output_list = []
    m = 0
    while True:
        m = m + 1
        next_eclipse = (row['t0_bjd_tdb'] + row['ephem_period'] * (int(row['current_epoch']) + m))
        print(f"t0_bjd_tdb: {row['t0_bjd_tdb']}, ephem_period: {row['ephem_period']}, current_epoch: {row['current_epoch']}, next_eclipse: {next_eclipse}")
        print(tdb_end)
        if next_eclipse > float(tdb_end):
            break
        output_list.append(float(next_eclipse))
    return output_list

def start_eclipse(row, extra_time):
    if not row['next_eclipse_times']:
        return None
    start_of_transit = [ecl - extra_time - 0.5 * row['duration_days'] for ecl in row['next_eclipse_times']]
    return start_of_transit

def end_eclipse(row, extra_time):
    if not row['next_eclipse_times']:
        return None
    end_of_transit = [ecl + extra_time + 0.5 * row['duration_days'] for ecl in row['next_eclipse_times']]
    return end_of_transit

def calc_altitude(row, time_column,location_telescope):
   time_jd = Time(row[time_column], format='jd')

   altitude = SkyCoord(ra=row['ra_j2000'], dec=row['dec_j2000'], unit=(u.hourangle, u.deg)) \
        .transform_to(AltAz(obstime=time_jd, location=location_telescope)).alt
   altitude_float = altitude.value
   return altitude_float

def twilight_according_transit(row, extra_time_for_twillight, time_column, observatory_location):
     #Initialize an empty list to store sun altitudes
    sun_altitudes = []

     #Iterate over each start_of_transit time in the list
    for time_value in row[time_column]:
         #Convert the start time of transit to Time object
        time_jd = Time(time_value, format='jd')

        # Calculate the Sun's coordinates at the specified time
        sun_coords = get_sun(time_jd)

        # Transform the Sun's coordinates to the AltAz frame for your observatory
        sun_altaz = sun_coords.transform_to(AltAz(obstime=time_jd, location=observatory_location))

        # Access the altitude of the Sun and append to the list
        sun_altitudes.append(sun_altaz.alt.value)

    return sun_altitudes

# Define the exponential function for exposure time calculation

def exponential_function(x, a, b):
  return a * np.exp(b * x)

def extract_lists(row, columns_with_lists):
    extracted_data = {}
    for col in columns_with_lists:
        extracted_data[col] = row[col]
    return pd.Series(extracted_data)

def create_new_dataframe(new_df):
    data = []

    for index, row in new_df.iterrows():
        exploded_start_altitude = pd.Series(row['start_altitude']).explode() if isinstance(row['start_altitude'], list) else pd.Series([row['start_altitude']])
        exploded_end_altitude = pd.Series(row['end_altitude']).explode() if isinstance(row['end_altitude'], list) else pd.Series([row['end_altitude']])

        max_length = max(len(row["start_of_transit"]), len(exploded_start_altitude), len(row["end_of_transit"]))

        for i in range(max_length):
            data.append(
                {
                    "name": row["name"],
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

    new_df = pd.DataFrame.from_dict(data)
    return new_df

def bjd_to_utc(bjd_value):
    bjd_time = Time(bjd_value, format='jd', scale='tdb')
    return bjd_time.utc.datetime













df = get_data(URL) # download the data

# apply the transformation from sexadecimal to decimal coordinates
df['ra_j2000'] = df.ra_j2000.apply(sexagesimal_to_decimal)
df['dec_j2000'] = df.dec_j2000.apply(sexagesimal_to_decimal)

#making the selected columns float so we can use them
df[columns_to_convert] = df[columns_to_convert].astype(float)

# first criteria needed depending on your telescope and needs
#df['priority'] = df['priority'].map(lambda x: x if x in PRIORITY_criteria else None)

new_df= df[
    ((df["priority"] == "high") | (df["priority"] == "medium") | (df["priority"] == "alert")) &
    (df["min_telescope_inches"] >= MIN_APERTURE) & (df["min_telescope_inches"] <= MAX_APERTURE) & (df["r_mag"] >= R_MAG_MIN)
]

# calculating currenting epoch of each object
new_df = calc_epoch(new_df, date_time_start, date_time_end)

# applying the fuction to calculate next eclipses and drop those objects that
# have no eclipses in the wanted observation period
new_df['next_eclipse_times'] = new_df.apply(lambda row: calc_next_eclipses(row, tdb_end), axis=1)
new_df = new_df[new_df['next_eclipse_times'].apply(lambda x: bool(x))]

#here we filter the stars that are observable by the position of the telescope

observable_targets_set = set()

for index, row in new_df.iterrows():
    # Extract target coordinates from the DataFrame
    target_coordinates = SkyCoord(ra=row['ra_j2000'] * u.deg, dec=row['dec_j2000'] * u.deg, frame='icrs')

    # Iterate over each next eclipse time
    for next_eclipse_time in row['next_eclipse_times']:
        # Specify the time of observation
        observation_time = Time(next_eclipse_time, format='jd', scale='tdb')

        # Check if the target is observable at the specified time
        is_observable = observer.target_is_up(observation_time, target_coordinates)

        if is_observable:
            # Add the target name to the set
            observable_targets_set.add(row['star'])

# Filter the DataFrame based on the set of unique observable targets
new_df = new_df[new_df['star'].isin(observable_targets_set)]

#calculating the duration days from the duration hours
new_df['duration_days']=new_df.duration_hours/24

#applying the start and end of the transit function
new_df['start_of_transit'] = new_df.apply(start_eclipse, extra_time=extra_time, axis=1)
new_df = new_df.dropna(subset=['start_of_transit'])

new_df['end_of_transit'] = new_df.apply(end_eclipse, extra_time=extra_time, axis=1)
new_df = new_df.dropna(subset=['end_of_transit'])

# calculating the start and end of each transit

new_df['start_altitude'] = new_df.apply(calc_altitude, args=('start_of_transit', observatory_location), axis=1)
new_df['end_altitude'] = new_df.apply(calc_altitude, args=('end_of_transit', observatory_location), axis=1)

# calculating before and after the sun altitude so we can filter out the no
# visible transist due to twilight
# we are putting 30 minutes before and after the transit
new_df['sun_altitude_before'] = new_df.apply(
    twilight_according_transit, args=(extra_time_for_twillight, 'start_of_transit', observatory_location) , axis=1)

new_df['sun_altitude_after'] = new_df.apply(
    twilight_according_transit, args=(extra_time_for_twillight, 'end_of_transit', observatory_location) , axis=1)

# using the 9mmag exponential fit we calculate the exposure time
new_df["exp_time"] = new_df['r_mag'].apply(lambda x: (10 * 40000) / exponential_function(x, A, B))

#creating a new dataframe with the columns we have put in the function
new_df = create_new_dataframe(new_df)

# the criteria that define if the altitude of the sun is acceptable
# and the altitude of the stars are visible

new_df = new_df[
    (
        (new_df['sun_altitudes_before'] <= -10.0) &
        (new_df['sun_altitudes_after'] <= -10.0) &
        (new_df['start_altitude'] >= 30.0) &
        (new_df['end_altitude'] >= 30.0)
    )
]

#calculating the time in utc and making it into string in our desired format

new_df['T_0_utc'] = new_df['T_0'].apply(bjd_to_utc)
new_df['start_utc'] = new_df['start'].apply(bjd_to_utc)
new_df['end_utc'] = new_df['end'].apply(bjd_to_utc)

new_df['start_utc'] = pd.to_datetime(new_df['start_utc']).dt.strftime('%Y-%m-%d %H:%M:%S')

new_df['end_utc'] = pd.to_datetime(new_df['end_utc']).dt.strftime('%Y-%m-%d %H:%M:%S')

#now let's make the time to datetime number
new_df['start_utc'] = pd.to_datetime(new_df['start_utc'])

#we want O-C to be between a specific range (15 > O-C > -15),
# if not to recalculate the exp_time ACCORDING TO THAT DIFFERENCE

condition = new_df['current_oc_min'] <= -15

# Update 'start_utc' based on the condition
new_df.loc[condition, 'start_utc'] += pd.to_timedelta(abs(new_df.loc[condition, 'current_oc_min']), unit='m')

#ROUNDING THE TIME TO EVERY 15 min to be more neat
new_df['obs_start_time'] = new_df['start_utc'].dt.round('15min')
new_df['obs_start_time'] = new_df['obs_start_time'].astype(str)

new_df[['real_obs_date', 'obs_start_time']] = new_df['obs_start_time'].str.split(' ', expand=True)


new_df['end_utc'] = pd.to_datetime(new_df['end_utc'])

condition = new_df['current_oc_min'] >= 15

# Update 'end_utc' based on the condition
new_df.loc[condition, 'end_utc'] += pd.to_timedelta(abs(new_df.loc[condition, 'current_oc_min']), unit='m')


new_df['obs_end_time'] = new_df['end_utc'].dt.round('15min')
new_df['obs_end_time'] = new_df['obs_end_time'].astype(str)

# splitting time so we can keep only the end time and not the date
new_df[['end_date_delete', 'obs_end_time']] = new_df['obs_end_time'].str.split(' ', expand=True)

#
new_df['start_utc'] = pd.to_datetime(new_df['start_utc'])

# if the transit is happening between 0am to 12pm substract a day from the date
# and naming it observation night
new_df['obs_night'] = new_df['start_utc'].apply(lambda row_dt: row_dt - pd.Timedelta(days=1) if 0 <= row_dt.hour < 12 else row_dt)

#making the observational time into string so we can slpit it and
# delete the time
new_df['obs_night'] = new_df['obs_night'].astype(str)
new_df[['obs_night', 'obs_time_delete']] = new_df['obs_night'].str.split(' ', expand=True)

# new columns to create a needed table for the scheduling
new_df['FILTER']='R'
new_df['Focus']=11850
new_df['Binning']=2
new_df['GAIN']=0
new_df['project']='ExoClock UZ'


# putting it in the order we want
new_column_order = ['obs_night','obs_start_time', 'obs_end_time', 'name', 'project','ra_j2000','dec_j2000','FILTER','exp_time','Focus','Binning','GAIN']
new_df = new_df[new_column_order]
# and sorting the transit based on the time and day
new_df = new_df.sort_values(by=['obs_night', 'obs_start_time'], ascending=[True, True])

# saving all the data
new_df.to_csv("all_data_schedule.csv", header=1, index=True)

new_df

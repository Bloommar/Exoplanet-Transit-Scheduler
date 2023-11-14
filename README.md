#Exoplanet Transit Scheduler 
###test

This python code has the purpose to create a schedule file that contains all the exoplanets transits that are observable by your telescope in your timeframe.

First it downloads all the data needed automatically from exoClock server. 

You have to put your telescope constraints:
-minimum and maximum aperture of the object-telescope observability
-the minimum magnitude of the R filter your telescope can observe and 
-your telescopes latitude

The ExoClock server downloads the coordinates to sexagecimal, and we are converting them to decimal for more clarity and usability for filtering the 56 that has the declination and 55 that has the right accession

Because, it’s hard to filter in that way the data, we transpose the df and we convert the needed columns to float numbers. 

CRITERIA

We put the first criteria based on our telescope constraints and we order the columns as wanted.

We calculate transit start and end of the 1st transit of each object given the datas last transit mid time. 

Converting the time we want to start  and end our observations to Julian date.

Of course, we have to calculate the current epoch using our start observation night, the period of the transit  and the integer of the current epoch. This is why:

The current epoch we calculate above should be a float number since it will be between 2 transits; int(current_epoch) and int(current_epoch)+1. The current epoch can be, of course, the current time, if the data from the server are up to date. 

The next eclipses/transits are calculated as the following and are stored in a cell as a list 
 
    next_eclipse=mid_time+ Period*{int(current_epoch)+m_times}
Where m_times increases by one every time the calculation is happening


We calculate the duration of the period of hours we want to have observations and we divided by 24 to convert in days. 

Now, with the duration of observation days we calculate the start and end of its transit for each list. We will need that to see the altitude of the object before and after each transit.

We do not need the objects in our dated-frame that don’t have a transit so we drop them.

At this point we use the altitude calculation of satrapy to calculate the altitude of each object in each transit and store it into a list. 

We need to see where is the sun in those transits , so we calculate it and then filter it if it’s over -12, which means over nautical twilight. 



**The magic happens for this point on**

We extract the list in rows by creating a dictionary. And finally we can filter it!!!

First, we filter based on the sun altitude as stated before. 

We convert times in UTC+0 and round them up to every 15’. 

We add the new columns that we have to include for the scheduling. 

We calculate the exposure time based on the data fitting we have previously done in another code, which takes into account our own data, and if we calculate it it flags up “1” so we know which one did we calculate. 

Order the columns as needed and FINALLY we sort the data by date!


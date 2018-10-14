function year_doy_sod_vector  = NRLMSISE_time_fromJD(jd)
% Return the year, day of year, and second of day in UT give a julian day
mytime= datetime(jd,'ConvertFrom','juliandate');
year= mytime.Year;
sod= seconds(timeofday(mytime));
doy= day(mytime,'dayofyear');
year_doy_sod_vector= [year doy sod];
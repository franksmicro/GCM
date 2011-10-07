/* This is the name of the data file we will create. */
#define FILE_NAME "gcm.nc"
     
/* We are writing 4D data, a 10 x 18 x 72 lvl-lat-lon grid */
#define NDIMS 4

#define UNITS "units"

/* Size of each dim */
#define START_TIME 0
#define TIME_STRIDE 1440
#define NTIME 144000
#define TIME_NAME "time"
#define TIME_UNITS "days"
#define TIME_DT 60

#define NLVL 10
#define LVL_NAME "level"
#define LVL_UNITS "layernum"

#define NLAT 36
#define LAT_NAME "latitude"
#define LAT_UNITS "degrees"

#define NLON 72
#define LON_NAME "longitude"
#define LON_UNITS "degrees"

#define UNITS "units"
#define LONG_NAME "long name"

/* Attributes for vars */
#define HEIGHT_NAME "height"
#define HEIGHT_UNITS "meters"
#define TEMP_NAME "temperature"
#define TEMP_UNITS "kelvin"
#define THETA_NAME "potential temperature"
#define THETA_UNITS "kelvin"
#define U_NAME "u_velocity"
#define U_UNITS "meters per second"
#define V_NAME "v_velocity"
#define V_UNITS "meters per second"
#define X_NAME "vorticity"
#define X_UNITS "per second"
#define A_NAME "vorticity advection"
#define A_UNITS "per second"
#define W_NAME "w_velocity"
#define W_UNITS "meters per second"


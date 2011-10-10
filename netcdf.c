/*
	netcdf.c - Create netcdf records for gcm
	Written by Frank Kienast in August, 2011
*/

#include <stdio.h>
#include <string.h>
#include <netcdf.h>

#include "griddefs.h"
     
/* Handle errors by printing an error message and exiting with a
 * non-zero status. */
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); return 2;}
     
int dimids[NDIMS];
int ncid, time_dimid, lvl_dimid, lat_dimid, lon_dimid;
int height_varid, temp_varid, u_varid, v_varid, x_varid, a_varid, w_varid;
int theta_varid;
     
/* The start and count arrays will tell the netCDF library where to
   write our data. */
size_t start[NDIMS], count[NDIMS];


int init_netcdf()
{
     
	/* These program variables hold the latitudes and longitudes. */
	float times[NTIME/TIME_STRIDE+1];
	float lats[NLAT], lons[NLON];
	long lvls[NLVL];

	/* Loop indexes. */
	int tim, lvl, lat, lon;

	/* Error handling. */
	int retval;

	/* Create index data for header */
	for(tim = 0; tim <= NTIME; tim += TIME_STRIDE)
		times[tim/TIME_STRIDE] = tim/TIME_STRIDE;
	for(lvl = 0; lvl < (NLVL); lvl++)
		lvls[lvl] = lvl;
	double latinc = 180. / NLAT;
	double loninc = 360. / NLON;
	for (lat = 0; lat < NLAT; lat++)
		lats[lat] = -90 + latinc * lat;
	for (lon = 0; lon < NLON; lon++)
		lons[lon] = loninc * lon;
     
        /* Create the file. */
	if ((retval = nc_create(FILE_NAME, NC_CLOBBER, &ncid)))
		ERR(retval);
 
	/* Define the dimensions. */
	if ((retval = nc_def_dim(ncid, TIME_NAME, NTIME/TIME_STRIDE+1, &time_dimid)))
		ERR(retval);
	if ((retval = nc_def_dim(ncid, LVL_NAME, (NLVL), &lvl_dimid)))
		ERR(retval);
	if ((retval = nc_def_dim(ncid, LAT_NAME, NLAT, &lat_dimid)))
		ERR(retval);
	if ((retval = nc_def_dim(ncid, LON_NAME, NLON, &lon_dimid)))
		ERR(retval);
     
	/* Define the coordinate variables. */
	if ((retval = nc_def_var(ncid, TIME_NAME, NC_FLOAT, 1, &time_dimid,
		&time_dimid)))
		ERR(retval);
	if ((retval = nc_def_var(ncid, LVL_NAME, NC_LONG, 1, &lvl_dimid,
		&lvl_dimid)))
		ERR(retval);
	if ((retval = nc_def_var(ncid, LAT_NAME, NC_FLOAT, 1, &lat_dimid,
		&lat_dimid)))
		ERR(retval);
	if ((retval = nc_def_var(ncid, LON_NAME, NC_FLOAT, 1, &lon_dimid,
		&lon_dimid)))
		ERR(retval);
 
	/* Assign units attributes to coordinate variables. */
	if ((retval = nc_put_att_text(ncid, time_dimid, UNITS,
		strlen(TIME_UNITS), TIME_UNITS)))
		ERR(retval);
	if ((retval = nc_put_att_text(ncid, lvl_dimid, UNITS,
		strlen(LVL_UNITS), LVL_UNITS)))
		ERR(retval);
	if ((retval = nc_put_att_text(ncid, lat_dimid, UNITS,
		strlen(LAT_UNITS), LAT_UNITS)))
		ERR(retval);
	if ((retval = nc_put_att_text(ncid, lon_dimid, UNITS,
		strlen(LON_UNITS), LON_UNITS)))
		ERR(retval);
 
	/* The dimids array is used to pass the dimids of the dimensions of
	the netCDF variables. Both of the netCDF variables we are
	creating share the same four dimensions. In C, the
	unlimited dimension must come first on the list of dimids. */
        dimids[0] = time_dimid;
        dimids[1] = lvl_dimid;
        dimids[2] = lat_dimid;
        dimids[3] = lon_dimid;
     
	/* Define the netCDF variables for height, temperature, u, v, w */
	if ((retval = nc_def_var(ncid, HEIGHT_NAME, NC_FLOAT, NDIMS,
		dimids, &height_varid)))
		ERR(retval);
	if ((retval = nc_def_var(ncid, TEMP_NAME, NC_FLOAT, NDIMS,
		dimids, &temp_varid)))
		ERR(retval);
	/*
	if ((retval = nc_def_var(ncid, THETA_NAME, NC_FLOAT, NDIMS,
		dimids, &theta_varid)))
		ERR(retval);
	*/
	if ((retval = nc_def_var(ncid, U_NAME, NC_FLOAT, NDIMS,
		dimids, &u_varid)))
		ERR(retval);
	if ((retval = nc_def_var(ncid, V_NAME, NC_FLOAT, NDIMS,
		dimids, &v_varid)))
		ERR(retval);
	if ((retval = nc_def_var(ncid, X_NAME, NC_FLOAT, NDIMS,
		dimids, &x_varid)))
		ERR(retval);
	/*
	if ((retval = nc_def_var(ncid, A_NAME, NC_FLOAT, NDIMS,
		dimids, &a_varid)))
		ERR(retval);
	*/
	if ((retval = nc_def_var(ncid, W_NAME, NC_FLOAT, NDIMS,
		dimids, &w_varid)))
		ERR(retval);

	/* Assign units attributes to the netCDF variables. */
	if ((retval = nc_put_att_text(ncid, height_varid, UNITS,
		strlen(HEIGHT_UNITS), HEIGHT_UNITS)))
		ERR(retval);
	if ((retval = nc_put_att_text(ncid, temp_varid, UNITS,
		strlen(TEMP_UNITS), TEMP_UNITS)))
		ERR(retval);
	/*
	if ((retval = nc_put_att_text(ncid, temp_varid, UNITS,
		strlen(THETA_UNITS), THETA_UNITS)))
		ERR(retval);
	*/
	if ((retval = nc_put_att_text(ncid, u_varid, UNITS,
		strlen(U_UNITS), U_UNITS)))
		ERR(retval);
	if ((retval = nc_put_att_text(ncid, v_varid, UNITS,
		strlen(V_UNITS), V_UNITS)))
		ERR(retval);
	if ((retval = nc_put_att_text(ncid, x_varid, UNITS,
		strlen(X_UNITS), X_UNITS)))
		ERR(retval);
	/*
	if ((retval = nc_put_att_text(ncid, a_varid, UNITS,
		strlen(A_UNITS), A_UNITS)))
		ERR(retval);
	*/
	if ((retval = nc_put_att_text(ncid, w_varid, UNITS,
		strlen(W_UNITS), W_UNITS)))
		ERR(retval);

	/* End define mode. */
	if ((retval = nc_enddef(ncid)))
		ERR(retval);
	/* Write the coordinate variable data. */
	if ((retval = nc_put_var_float(ncid, time_dimid, &times[0])))
		ERR(retval);
	if ((retval = nc_put_var_long(ncid, lvl_dimid, &lvls[0])))
		ERR(retval);
	if ((retval = nc_put_var_float(ncid, lat_dimid, &lats[0])))
		ERR(retval);
	if ((retval = nc_put_var_float(ncid, lon_dimid, &lons[0])))
		ERR(retval);

	/* These settings tell netcdf to write one timestep of data. (The
	setting of start[0] inside the loop below tells netCDF which
	timestep to write.) */
	count[0] = 1;
	count[1] = (NLVL);
	count[2] = NLAT;
	count[3] = NLON;
	start[1] = 0;
	start[2] = 0;
	start[3] = 0;


	return 0;
}

int write_height(int recno,double dvals[NLVL][NLAT][NLON])
{
	int lvl,lat,lon;
	float fvals[NLVL][NLAT][NLON];
	for(lvl = 0; lvl < (NLVL); lvl++)
		for(lat = 0; lat < NLAT; lat++)
			for(lon = 0; lon < NLON; lon++)
				fvals[lvl][lat][lon] = (float) dvals[lvl][lat][lon];
	start[0] = recno;
	return nc_put_vara_float(ncid, height_varid, start, count, &fvals[0][0][0]);
}

int write_temperature(int recno,double dvals[NLVL][NLAT][NLON])
{
	int lvl,lat,lon;
	float fvals[NLVL][NLAT][NLON];
	for(lvl = 0; lvl < (NLVL); lvl++)
		for(lat = 0; lat < NLAT; lat++)
			for(lon = 0; lon < NLON; lon++)
				fvals[lvl][lat][lon] = (float) dvals[lvl][lat][lon];
	start[0] = recno;
	return nc_put_vara_float(ncid, temp_varid, start, count, &fvals[0][0][0]);

}

/*
int write_theta(int recno,double dvals[NLVL][NLAT][NLON])
{
	int lvl,lat,lon;
	float fvals[NLVL][NLAT][NLON];
	for(lvl = 0; lvl < (NLVL); lvl++)
		for(lat = 0; lat < NLAT; lat++)
			for(lon = 0; lon < NLON; lon++)
				fvals[lvl][lat][lon] = (float) dvals[lvl][lat][lon];
	start[0] = recno;
	return nc_put_vara_float(ncid, theta_varid, start, count, &fvals[0][0][0]);

}
*/

int write_u(int recno,double dvals[NLVL][NLAT][NLON])
{
	start[0] = recno;
	int lvl,lat,lon;
	float fvals[NLVL][NLAT][NLON];
	for(lvl = 0; lvl < (NLVL); lvl++)
		for(lat = 0; lat < NLAT; lat++)
			for(lon = 0; lon < NLON; lon++)
				fvals[lvl][lat][lon] = (float) dvals[lvl][lat][lon];
	return nc_put_vara_float(ncid, u_varid, start, count, &fvals[0][0][0]);
}

int write_v(int recno,double dvals[NLVL][NLAT][NLON])
{
	int lvl,lat,lon;
	float fvals[NLVL][NLAT][NLON];
	for(lvl = 0; lvl < (NLVL); lvl++)
		for(lat = 0; lat < NLAT; lat++)
			for(lon = 0; lon < NLON; lon++)
				fvals[lvl][lat][lon] = (float) dvals[lvl][lat][lon];
	start[0] = recno;
	return nc_put_vara_float(ncid, v_varid, start, count, &fvals[0][0][0]);
}

int write_x(int recno,double dvals[NLVL][NLAT][NLON])
{
	int lvl,lat,lon;
	float fvals[NLVL][NLAT][NLON];
	for(lvl = 0; lvl < (NLVL); lvl++)
		for(lat = 0; lat < NLAT; lat++)
			for(lon = 0; lon < NLON; lon++)
				fvals[lvl][lat][lon] = (float) dvals[lvl][lat][lon];
	start[0] = recno;
	return nc_put_vara_float(ncid, x_varid, start, count, &fvals[0][0][0]);
}

/*
int write_a(int recno,double dvals[NLVL][NLAT][NLON])
{
	int lvl,lat,lon;
	float fvals[NLVL][NLAT][NLON];
	for(lvl = 0; lvl < (NLVL); lvl++)
		for(lat = 0; lat < NLAT; lat++)
			for(lon = 0; lon < NLON; lon++)
				fvals[lvl][lat][lon] = (float) dvals[lvl][lat][lon];
	start[0] = recno;
	return nc_put_vara_float(ncid, a_varid, start, count, &fvals[0][0][0]);
}
*/

int write_w(int recno,double dvals[NLVL][NLAT][NLON])
{
	int lvl,lat,lon;
	float fvals[NLVL][NLAT][NLON];
	for(lvl = 0; lvl < (NLVL); lvl++)
		for(lat = 0; lat < NLAT; lat++)
			for(lon = 0; lon < NLON; lon++)
				fvals[lvl][lat][lon] = (float) dvals[lvl][lat][lon];
	start[0] = recno;
	return nc_put_vara_float(ncid, w_varid, start, count, &fvals[0][0][0]);
}

int close_netcdf()
{
	return nc_close(ncid);
}


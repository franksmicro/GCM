/*
gcm.c - Attempt a "toy" global climate model
Written by Frank Kienast in 2011
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <netcdf.h>
#include <fftw3.h>

#include "griddefs.h"

#define PI 3.1415927

#define ERR(s) {printf("Error: %s\n", s); return 2;}

double U[NLVL][NLAT][NLON];
double UPrev[NLVL][NLAT][NLON];
double UNext[NLVL][NLAT][NLON];
double V[NLVL][NLAT][NLON];
double VPrev[NLVL][NLAT][NLON];
double VNext[NLVL][NLAT][NLON];
double X[NLVL][NLAT][NLON];
double XPrev[NLVL][NLAT][NLON];
double XNext[NLVL][NLAT][NLON];
double A[NLVL][NLAT][NLON];
double APrev[NLVL][NLAT][NLON];
double ANext[NLVL][NLAT][NLON];
double W[NLVL][NLAT][NLON];
double WPrev[NLVL][NLAT][NLON];
double WNext[NLVL][NLAT][NLON];
double T[NLVL][NLAT][NLON];
double TPrev[NLVL][NLAT][NLON];
double TNext[NLVL][NLAT][NLON];
double Theta[NLVL][NLAT][NLON];
double ThetaPrev[NLVL][NLAT][NLON];
double ThetaNext[NLVL][NLAT][NLON];
double Z[NLVL][NLAT][NLON];
double ZPrev[NLVL][NLAT][NLON];
double ZNext[NLVL][NLAT][NLON];

double EarthC[NLAT][NLON];
double EarthT[NLAT][NLON];


fftw_complex fftarr[NLAT][NLON];
fftw_plan p1,p2;

int main(int argc, char *argv[])
{
	long tFinal = NTIME;
	long t;
	int lev,lon,lat;

	void initVars();
	void doRadiation();
	void calcHeights();
	void AdvectX();
	void calcU();
	void calcV();
	void calcX();
	void calcA();
	void calcW();
	void calcT();
	void calcTheta();
	void calcEarth();
	int writeResults();

	if(init_netcdf())
		ERR("init_netcdf");
		
	initVars();
	writeResults(0);		
	for(t = 1; t <= tFinal; t++)
	{
		printf("t = %ld\n",t);
		doRadiation(t);	
		calcHeights();
		calcTheta();
		calcU();
		calcV();
		calcX();
		calcW();
		calcT(t);
		//calcEarth();

		for(lev = 0; lev < NLVL; lev++)
		for(lat = 0; lat < NLAT; lat++)
		for(lon = 0; lon < NLON; lon++)
		{
			UPrev[lev][lat][lon] = U[lev][lat][lon];
			U[lev][lat][lon] = UNext[lev][lat][lon];

			VPrev[lev][lat][lon] = V[lev][lat][lon];
			V[lev][lat][lon] = VNext[lev][lat][lon];

			XPrev[lev][lat][lon] = X[lev][lat][lon];
			X[lev][lat][lon] = XNext[lev][lat][lon];

			APrev[lev][lat][lon] = A[lev][lat][lon];
			A[lev][lat][lon] = ANext[lev][lat][lon];

			WPrev[lev][lat][lon] = W[lev][lat][lon];
			W[lev][lat][lon] = WNext[lev][lat][lon];

			TPrev[lev][lat][lon] = T[lev][lat][lon];
			T[lev][lat][lon] = TNext[lev][lat][lon];

			ThetaPrev[lev][lat][lon] = Theta[lev][lat][lon];
			Theta[lev][lat][lon] = ThetaNext[lev][lat][lon];

			ZPrev[lev][lat][lon] = Z[lev][lat][lon];
			Z[lev][lat][lon] = ZNext[lev][lat][lon];

		}

		if((t % TIME_STRIDE) == 0)
			if(writeResults(t/TIME_STRIDE))
				return 1;
	}
	if(close_netcdf())
		ERR("close_netcdf");

	return 0;
}

void initVars()
{
	int lev,lat,lon;

	for(lev = 0; lev < NLVL; lev++)
		for(lat = 0; lat < NLAT; lat++)
			for(lon = 0; lon < NLON; lon++)
			{
				switch(lev)
				{
					case 0: 
						T[lev][lat][lon] = 270;
						break;
					case 1:
						T[lev][lat][lon] = 250;
						break;
					case 2:
						T[lev][lat][lon] = 240;
						break;
					case 3:
						T[lev][lat][lon] = 230;
						break;
					case 4:
						T[lev][lat][lon] = 220;
						break;
					case 5:
						T[lev][lat][lon] = 210;
						break;
					case 6:
						T[lev][lat][lon] = 200;
						break;
					case 7:
						T[lev][lat][lon] = 190;
						break;
					case 8:
						T[lev][lat][lon] = 180;
						break;
					case 9:
						T[lev][lat][lon] = 160;
						break;
				}
				TPrev[lev][lat][lon] = T[lev][lat][lon];
				TNext[lev][lat][lon] = T[lev][lat][lon];
				U[lev][lat][lon] = 0;
				UPrev[lev][lat][lon] = 0;
				V[lev][lat][lon] = 0;
				VPrev[lev][lat][lon] = 0;
				X[lev][lat][lon] = 0;
				XPrev[lev][lat][lon] = 0;
				A[lev][lat][lon] = 0;
				APrev[lev][lat][lon] = 0;
				W[lev][lat][lon] = 0;
				WPrev[lev][lat][lon] = 0;
				WNext[lev][lat][lon] = 0;
				Z[lev][lat][lon] = 0;
				ZPrev[lev][lat][lon] = 0;
			}

	for(lat = 0; lat < NLAT; lat++)
	for(lon = 0; lon < NLON; lon++)
	{
		double latinc = 180. / NLAT;
		double loninc = 360. / NLON;
		double longitude = loninc * lon + loninc / 2.;
		double latitude = -90. + latinc * lat + latinc / 2.;
		EarthC[lat][lon] = 1000E3;
		EarthT[lat][lon] = T[0][lat][lon];

		/*
		//Water lon 135 to 225 lat 20 to 60
		if(longitude >= 135 && longitude <= 225
		&& latitude >= 20 && latitude <= 60)
			EarthC[lat][lon] *= 50;
		*/
	}

	p1 = fftw_plan_dft_2d(NLAT,NLON,
		&fftarr[0][0],&fftarr[0][0],
		FFTW_FORWARD, FFTW_ESTIMATE);
	p2 = fftw_plan_dft_2d(NLAT,NLON,
		&fftarr[0][0],&fftarr[0][0],
		FFTW_BACKWARD, FFTW_ESTIMATE);

}
	
int writeResults(long t)
{
	if(write_height(t,Z))
		ERR("write_height");
	if(write_temperature(t,T))
		ERR("write_temperature");
	if(write_u(t,U))
		ERR("write_u");
	/*
	if(write_theta(t,Theta))
		ERR("write_theta");
	*/
	if(write_v(t,V))
		ERR("write_v");
	if(write_x(t,X))
		ERR("write_x");
	/*
	if(write_a(t,A))
		ERR("write_a");
	*/
	if(write_w(t,W))
		ERR("write_w");
}	

void doRadiation(long t)
{
	int lev = 0;
	int lat, lon;
	double dt = TIME_DT;
	void GrayAtmosphere();

	long daylen = 86400 / TIME_DT;
	double degpert = 360. / daylen;
	double rotations = (t * degpert)/ 360.;
	double fraction = rotations - ((long)rotations);
	double sunlongitude = -fraction * 360.;

	long yearlen = 31557600 / TIME_DT;
	double deg2pert = 360. / yearlen;
	double revolutions = t * deg2pert / 365.25;
	double fraction2 = revolutions - ((long)revolutions);
	double sunlatitude = 22.5 * sin(2. * PI * fraction2);

	double latitude, longitude;
	double c;

	double calcRadIn();
	double calcRadOut();
		
	for(lat = 0; lat < NLAT; lat++)
	for(lon = 0; lon < NLON; lon++)
	{
		double latinc = 180. / NLAT;
		double loninc = 360. / NLON;
		double longitude = loninc * lon + loninc / 2.;
		double latitude = -90. + latinc * lat + latinc / 2.;

		/*
		c = 1010 * 1000; // Joules per cubic meter * thickness of layer

		TNext[lev][lat][lon] = T[lev][lat][lon] +  dt * 
		(calcRadIn(t,latitude,longitude,sunlatitude,sunlongitude) -
		calcRadOut(T[lev][lat][lon])) / c;
		*/
		double sunmagnitude = calcRadIn(t,latitude,longitude,
					sunlatitude,sunlongitude);
		GrayAtmosphere(lat,lon,sunmagnitude);
	}
	
}
	
double calcRadIn(long t,double latitude, double longitude, double sunlatitude, double sunlongitude)
{
	//Solar radiation in Watts
	
	double SUNMAG = 1400; // W/m^2
	double latrads, lonrads, sunlonrads,sunlatrads;
	double radIn;
	
	latrads = PI * latitude / 180;
	lonrads = PI * longitude / 180.;
	sunlonrads = PI * sunlongitude / 180.;
	sunlatrads = PI * sunlatitude / 180.;
		
	//If either of these is negaive, want 0, not + if both negative
	radIn =  SUNMAG * cos(latrads)
	* cos(lonrads - sunlonrads);
	if(radIn < 0) radIn = 0; 
	radIn *= cos(latrads - sunlatrads);
	if (radIn < 0) radIn = 0;

	return radIn;
}
	
/*
double calcRadOut(double T)
{
	//Radiation out in Watts
		
	return 5.67E-8 * T * T * T * T;
}
*/

void GrayAtmosphere(int lat, int lon, double sunmag)
{
	double tau = 1.8;
	double albedo = 0.7;
	double sigma = 5.67E-8;	
	double g = 9.82;
	double csurf = 1E6;
	double cp = 710;
	double p0 = 100000;

	double fu[NLVL+1]; //Upward flux
	double fd[NLVL+1]; //Downward flux
	double fudep[NLVL]; //Deposited upward flux
	double fddep[NLVL]; //Deposited downward flux
	double ftdep[NLVL]; //Deposited total flux

	int lev;
	double T1;

	T1 = T[0][lat][lon];
	fu[0] = sigma * (1-albedo) * T1 * T1 * T1 * T1;

	for(lev = 1; lev < NLVL+1; lev++)
	{
		T1 = T[lev-1][lat][lon];
		fu[lev] = fu[lev-1] + 
		(-0.5 * fu[lev-1] + sigma * T1 * T1 * T1 * T1 * (tau/NLVL)) / 
		(1 + 0.5 * tau/NLVL);
		
		fudep[lev-1] = fu[lev-1] - fu[lev];
	}

	for(lev = (NLVL-1); lev >= 0; lev--)
	{
		T1 = T[lev][lat][lon];
		fd[lev] = fd[lev+1] +
		(-0.5 * fd[lev+1] + sigma * T1 * T1 * T1 * T1 * (tau/NLVL)) /
		(1 + 0.5 * tau/NLVL);
		
		fddep[lev] = fd[lev+1] - fd[lev];
	}

	T1 = T[0][lat][lon];
	double fsurf = fd[0] + sunmag - sigma * T1 * T1 * T1 * T1;
	TNext[0][lat][lon] = T[0][lat][lon] + TIME_DT * fsurf / csurf;

	for(lev = 1; lev < NLVL; lev++)
	{
		/*
		printf("lev=%d,lat=%d,lon=%d,fudep=%f,fddep=%f\n",
			lev,lat,lon,fudep[lev],fddep[lev]);
		*/
		double ftot = fudep[lev] + fddep[lev];
		TNext[lev][lat][lon] = T[lev][lat][lon] 
		+ TIME_DT * ftot * g * NLVL / (cp * p0);
	}
}
	
void calcHeights()
{
	//Calculates the total height in meters at top of each layer k at point (i,j)
		
	double pressure, totheight;
	int lev,lat,lon;
		
	double calcZ();

	for(lat = 0; lat < NLAT; lat++)
	for(lon = 0; lon < NLON; lon++)
	{
		totheight = 0;
		for(lev = 0; lev < NLVL; lev++)
		{
			pressure = 100. - 50. / NLVL - (100. / NLVL) * lev;
			totheight += calcZ(pressure,TNext[lev][lat][lon]);
			ZNext[lev][lat][lon] = totheight;
		}
	}
}
	
double calcZ(double pressure, double temperature)
{
	//Thickness of pressure layer in meters (distance up to next pressure layer)
		
	double R = 287.053;
	double g = 9.80;	
	double p, T;
		
	p = pressure;	
	T = temperature;
	
	return (R * T) / (p * g) * 10; //10 is dp
}

double CalcDfDx(int lev, int lat, int lon, double var[NLVL][NLAT][NLON],
		int sign)
{
	double latinc = 180. / NLAT;
	double latitude = -90. + latinc * lat + latinc / 2.;
	double dx = 40E6 * cos(latitude * PI / 180) / NLON;
	double df;

	int lonprev = lon - 1; if(lonprev < 0) lonprev += NLON;
	int lonnext = lon + 1; if(lonnext >= NLON) lonnext -= NLON;

	
	if(sign > 0)
		df = var[lev][lat][lon] - var[lev][lat][lonprev];
	else
		df = var[lev][lat][lonnext] - var[lev][lat][lon];


	return df / dx;
}

double CalcDfDy(int lev, int lat, int lon, double var[NLVL][NLAT][NLON],
		int sign)
{
	double dy = 20E6 / NLAT;
	double df;

	if(sign > 0)
	{
		if(lat == 0)
			df = 0;
		else
			df = var[lev][lat][lon] - var[lev][lat-1][lon];
	}
	else 
	{
		if(lat == (NLAT-1))
			df = 0;
		else
			df = var[lev][lat+1][lon] - var[lev][lat][lon];
	}

	return df/dy;
}

double CalcDfDz(int lev, int lat, int lon, 
		double var[NLVL][NLAT][NLON], 
		double zvar[NLVL][NLAT][NLON],
		int sign)
{
	double dz,df;

	if(lev == 0)
		dz = zvar[lev][lat][lon];
	else
		dz = zvar[lev][lat][lon] - zvar[lev-1][lat][lon];

	if(sign > 0)
	{
		if(lev == 0)
			df = 0;
		else
			df = var[lev][lat][lon] - var[lev-1][lat][lon];
	}
	else
	{
		if(lev == (NLVL-1))
			df = 0;
		else
			df = var[lev+1][lat][lon] - var[lev][lat][lon];
	}

	if(dz <= 0)
		return 0;
	else
		return df/dz;
}

double sign(double val)
{
	if(val < 0) 
		return -1.;
	else
		return 1;
}

void AdvectX(double var[NLVL][NLAT][NLON],
		double uvar[NLVL][NLAT][NLON],
		double outvar[NLVL][NLAT][NLON])
{
	int lev, lat, lon;
	double var1[NLVL][NLAT][NLON];
	double var2[NLVL][NLAT][NLON];
	long dt = TIME_DT;

	for(lev = 0; lev < NLVL; lev++)
		for(lat = 0; lat < NLAT; lat++)
			for(lon = 0; lon < NLON; lon++)
			{
				double u = uvar[lev][lat][lon];
				var1[lev][lat][lon] = var[lev][lat][lon]
				- (1./3.) * u * dt * 
				CalcDfDx(lev,lat,lon,var,sign(u));
			}
				
	for(lev = 0; lev < NLVL; lev++)
		for(lat = 0; lat < NLAT; lat++)
			for(lon = 0; lon < NLON; lon++)
			{
				double u = uvar[lev][lat][lon];
				var2[lev][lat][lon] = var[lev][lat][lon]
				- (1./2.) * u * dt * 
				CalcDfDx(lev,lat,lon,var1,sign(u));
			}

	for(lev = 0; lev < NLVL; lev++)
		for(lat = 0; lat < NLAT; lat++)
			for(lon = 0; lon < NLON; lon++)
			{
				double u = uvar[lev][lat][lon];
				outvar[lev][lat][lon] = var[lev][lat][lon]
				- 1. * u * dt * 
				CalcDfDx(lev,lat,lon,var2,sign(u));
			}
}

void AdvectY(double var[NLVL][NLAT][NLON],
		double vvar[NLVL][NLAT][NLON],
		double outvar[NLVL][NLAT][NLON])
{
	int lev, lat, lon;
	double var1[NLVL][NLAT][NLON];
	double var2[NLVL][NLAT][NLON];
	long dt = TIME_DT;

	for(lev = 0; lev < NLVL; lev++)
		for(lat = 0; lat < NLAT; lat++)
			for(lon = 0; lon < NLON; lon++)
			{
				double v = vvar[lev][lat][lon];
				var1[lev][lat][lon] = var[lev][lat][lon]
				- (1./3.) * v * dt * 
				CalcDfDy(lev,lat,lon,var,sign(v));
			}
				
	for(lev = 0; lev < NLVL; lev++)
		for(lat = 0; lat < NLAT; lat++)
			for(lon = 0; lon < NLON; lon++)
			{
				double v = vvar[lev][lat][lon];
				var2[lev][lat][lon] = var[lev][lat][lon]
				- (1./2.) * v * dt * 
				CalcDfDy(lev,lat,lon,var1,sign(v));
			}

	for(lev = 0; lev < NLVL; lev++)
		for(lat = 0; lat < NLAT; lat++)
			for(lon = 0; lon < NLON; lon++)
			{
				double v = vvar[lev][lat][lon];
				outvar[lev][lat][lon] = var[lev][lat][lon]
				- 1. * v * dt * 
				CalcDfDy(lev,lat,lon,var2,sign(v));
			}
}

void AdvectZ(double var[NLVL][NLAT][NLON],
		double wvar[NLVL][NLAT][NLON],
		double outvar[NLVL][NLAT][NLON])
{
	int lev, lat, lon;
	double var1[NLVL][NLAT][NLON];
	double var2[NLVL][NLAT][NLON];
	long dt = TIME_DT;

	for(lev = 0; lev < NLVL; lev++)
		for(lat = 0; lat < NLAT; lat++)
			for(lon = 0; lon < NLON; lon++)
			{
				double w = wvar[lev][lat][lon];
				var1[lev][lat][lon] = var[lev][lat][lon]
				- (1./3.) * w * dt * 
				CalcDfDz(lev,lat,lon,var,Z,sign(w));
			}
				
	for(lev = 0; lev < NLVL; lev++)
		for(lat = 0; lat < NLAT; lat++)
			for(lon = 0; lon < NLON; lon++)
			{
				double w = wvar[lev][lat][lon];
				var2[lev][lat][lon] = var[lev][lat][lon]
				- (1./2.) * w * dt * 
				CalcDfDz(lev,lat,lon,var1,Z,sign(w));
			}

	for(lev = 0; lev < NLVL; lev++)
		for(lat = 0; lat < NLAT; lat++)
			for(lon = 0; lon < NLON; lon++)
			{
				double w = wvar[lev][lat][lon];
				outvar[lev][lat][lon] = var[lev][lat][lon]
				- 1. * w * dt * 
				CalcDfDz(lev,lat,lon,var2,Z,sign(w));
			}
}

void FFTFilter(double var[NLVL][NLAT][NLON])
{
	int lev, lat, lon;

	int lat1 = NLAT / 4;
	int lat2 = 3 * NLAT / 4;
	int lon1 = NLON / 4;
	int lon2 = 3 * NLON / 4;

	for(lev = 0; lev < NLVL; lev++)
	{
		//Store values for FFT
		for(lat = 0; lat < NLAT; lat++)
		for(lon = 0; lon < NLON; lon++)
		{
			fftarr[lat][lon][0] = var[lev][lat][lon];
			fftarr[lat][lon][1] = 0.;
		}

		//Forward FFT
		fftw_execute(p1);

		//Filter logic
		for(lat = 0; lat < NLAT; lat++)
		for(lon = 0; lon < NLON; lon++)
		{
			if((lat > lat1 && lat < lat2) 
			|| (lon > lon1 && lon < lon2))
			{
				fftarr[lat][lon][0] = 0.;
				fftarr[lat][lon][1] = 0.;
			}
		}
			
		//Backward FFT
		fftw_execute(p2);

		//Store values back 
		for(lat = 0; lat < NLAT; lat++)
		for(lon = 0; lon < NLON; lon++)
		{
			var[lev][lat][lon] = fftarr[lat][lon][0] /
			((double) NLAT * (double) NLON);
		}
	}
}


void calcU()
{
	int lev,lat,lon;
	double UX[NLVL][NLAT][NLON];
	double UY[NLVL][NLAT][NLON];

	double g = 9.80;
	double ohmega = 1.458E-4;
	double dt = TIME_DT; //Seconds
	double latinc = 180. / NLAT;

	AdvectX(U,U,UX); 
	AdvectY(UX,V,UY);
	AdvectZ(UY,W,UNext);

	for(lev = 0; lev < NLVL; lev++)
	for(lat = 0; lat < NLAT; lat++)
	for(lon = 0; lon < NLON; lon++)
	{
		double a = (lev == 0 ? 2E-8: 0);
		double latitude = -90. + latinc * lat + latinc / 2.;
		double f = 2. * ohmega * sin(PI * latitude / 180.);
		UNext[lev][lat][lon] += f * V[lev][lat][lon] * dt;
		UNext[lev][lat][lon] -= g * 
		CalcDfDx(lev,lat,lon,Z,sign(UNext[lev][lat][lon])) * dt;
		/*
		UNext[lev][lat][lon] -= a * sign(UNext[lev][lat][lon]) *
			UNext[lev][lat][lon] * UNext[lev][lat][lon] * dt;
		*/
		//UNext[lev][lat][lon] -= a * UNext[lev][lat][lon];
		//Enforce zero velocity at poles
		if(abs(latitude + 90.) < 5. 
		|| abs(latitude - 90.) < 5.)
			UNext[lev][lat][lon] = 0;
	} 
	//FFTFilter(UNext);

}

void calcV()
{
	int lev,lat,lon;
	double VX[NLVL][NLAT][NLON];
	double VY[NLVL][NLAT][NLON];

	double g = 9.80;
	double ohmega = 1.458E-4;
	double dt = TIME_DT; //Seconds
	double latinc = 180. / NLAT;

	AdvectX(V,U,VX); 
	AdvectY(VX,V,VY);
	AdvectZ(VY,W,VNext);

	for(lev = 0; lev < NLVL; lev++)
	for(lat = 0; lat < NLAT; lat++)
	for(lon = 0; lon < NLON; lon++)
	{
		double a = (lev == 0 ? 2E-8: 0);
		double latitude = -90. + latinc * lat + latinc / 2.;
		double f = 2. * ohmega * sin(PI * latitude / 180.);
		VNext[lev][lat][lon] -= f * U[lev][lat][lon] * dt;
		VNext[lev][lat][lon] -= g * 
		CalcDfDy(lev,lat,lon,Z,sign(VNext[lev][lat][lon])) * dt;
		/*
		VNext[lev][lat][lon] -= a * sign(VNext[lev][lat][lon]) *
			VNext[lev][lat][lon] * VNext[lev][lat][lon] * dt;
		*/
		//VNext[lev][lat][lon] -= a * VNext[lev][lat][lon];

		//Enforce zero velocity at poles
		if(abs(latitude + 90.) < 5. 
		|| abs(latitude - 90.) < 5.)
			VNext[lev][lat][lon] = 0;
	}

	//FFTFilter(VNext);

}
	
void calcTheta()
{
	//Calc height at middle of level

	int lev,lat,lon;
	double dz,dt;
	double R = 2870.53;
	double p = 100000 - 1./(2 * NLVL) -1./NLVL * lev;
	double c = 1010;

	for(lev = 0; lev < NLVL; lev++)
	for(lat = 0; lat < NLAT; lat++)
	for(lon = 0; lon < NLON; lon++)
	{
		if(lev == 0)
			dz = Z[lev][lat][lon] / 2.;
		else
			dz = (Z[lev][lat][lon] + Z[lev-1][lat][lon]) / 2.;


		ThetaNext[lev][lat][lon] = 
		TNext[lev][lat][lon] 
		+ 10. * (dz / 1000.);
	}
}

void calcX(int lev, int lat, int lon)
{
	// Vorticity calculation 

	double dvdx, dudy;

	for(lev = 0; lev < NLVL; lev++)
	for(lat = 1; lat < (NLAT-1); lat++) //Don't calc at poles
	for(lon = 0; lon < NLON; lon++)
	{
		dvdx = (CalcDfDx(lev,lat,lon,VNext,1) +
			CalcDfDx(lev,lat,lon,VNext,-1)) / 2.;
		dudy = (CalcDfDy(lev,lat,lon,UNext,1) +
			CalcDfDy(lev,lat,lon,UNext,-1)) / 2.;

		XNext[lev][lat][lon] = dvdx - dudy;
	}
}
	
/*
void calcA(int lev, int lat, int lon)
{
	//Vorticity advection

	double ohmega = 1.458E-4;
	double REarth = 6357E3;
	double latinc = 180. / NLAT;
	double latitude = -90. + latinc * lat + latinc / 2.;
	double f2 = 2. * ohmega * cos(PI * latitude / 180.) / REarth;

	double da;
	double dt = TIME_DT; //Seconds
	da = 	- calcAdvectX(lev,lat,lon,XNext)
		- calcAdvectY(lev,lat,lon,XNext);
		+ f2 * V[lev][lat][lon];
	ANext[lev][lat][lon] = da * dt;
}
*/

void calcW()
{
	int lev, lat, lon;
	double dt = TIME_DT;

	for(lev = 0; lev < NLVL; lev++)
	for(lat = 0; lat < NLAT; lat++)
	for(lon = 0; lon < NLON; lon++)
	{
		if(ZPrev[lev][lat][lon] != 0)
			WNext[lev][lat][lon] = 
			(Z[lev][lat][lon] 
			- ZPrev[lev][lat][lon]) / dt;
		else
			WNext[lev][lat][lon] = 0;
	}
}

void calcT(int iteration)
{
	int lev,lat,lon;
	double TX[NLVL][NLAT][NLON],TY[NLVL][NLAT][NLON],TZ[NLVL][NLAT][NLON];
	double WDummy[NLVL][NLAT][NLON];
	//So don't overwrite TOut during calculations
	//We use TNext because we have already updated TNext with radiation
	double TIn[NLVL][NLAT][NLON]; 
	double dt = TIME_DT;

	for(lev = 0; lev < NLVL; lev++)	
	for(lat = 0; lat < NLAT; lat++)
	for(lon = 0; lon < NLON; lon++)
	{
		TIn[lev][lat][lon] = TNext[lev][lat][lon];
	}

	AdvectX(TIn,U,TX);
	AdvectY(TX,V,TNext);
	//We can't advect T in Z because we must use theta, not T, so do below
	
	for(lev = 0; lev < NLVL; lev++)	
	for(lat = 0; lat < NLAT; lat++)
	for(lon = 0; lon < NLON; lon++)
	{
		TNext[lev][lat][lon] 
		-= CalcDfDz(lev,lat,lon,TNext,Z,sign(W[lev][lat][lon])) 
		* W[lev][lat][lon];
		TNext[lev][lat][lon]  
		-= W[lev][lat][lon] * 10. / 1000 * dt; //10K per km
		
		//Assume some sort of mixing or radiative transfer between
		//levels with time constant close to 2 day between levels
		/*
		double alpha = dt / 172800.;
		double deltat;
		if(lev == 0)
			deltat = ThetaNext[lev+1][lat][lon] 
			- ThetaNext[lev][lat][lon];
		else if(lev == (NLVL-1))
			deltat = ThetaNext[lev-1][lat][lon] 
			- ThetaNext[lev][lat][lon];
		else
			deltat = 
			((ThetaNext[lev+1][lat][lon] 
			- Theta[lev][lat][lon])
			+ (ThetaNext[lev-1][lat][lon] 
			- ThetaNext[lev][lat][lon])) / 2.;
		if(iteration > 5)
			TNext[lev][lat][lon] += alpha * deltat;
		*/
	}
}

void calcEarth()
{
	int lat, lon;
	double c1,c2,t1,t2,tf;

	for(lat = 0; lat < NLAT; lat++)
	for(lon = 0; lon < NLON; lon++)
	{
		c1 = 1010. * 1000.; 
		c2 = EarthC[lat][lon];
		t1 = TNext[0][lat][lon];
		t2 = EarthT[lat][lon];
		tf = (c1 * t1 + c2 * t2) / (c1 + c2);
		EarthT[lat][lon] = tf;
		TNext[0][lat][lon] = tf;
	}
}

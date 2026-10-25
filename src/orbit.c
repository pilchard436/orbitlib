#include "orbitlib_orbit.h"
#include "orbitlib_celestial.h"
#include <math.h>
#include <stdio.h>

struct Orbit constr_orbit_from_elements(double semimajor_axis, double eccentricity, double inclination, double raan, double arg_of_peri, double true_anomaly, struct Body *cb) {
	struct Orbit new_orbit;
	new_orbit.cb = cb;
	new_orbit.a = semimajor_axis;
	// avoid nan when converting to osv and back
	new_orbit.e = eccentricity != 0 ? eccentricity : eccentricity+1e-12;
	// avoid nan when converting to osv and back
	new_orbit.i = inclination != 0 ? inclination : inclination + 1e-12;
	new_orbit.raan = raan;
	new_orbit.arg_peri = arg_of_peri;
	new_orbit.ta = true_anomaly;
	
	return new_orbit;
}

struct Orbit constr_orbit_from_apsides(double apsis1, double apsis2, double inclination, struct Body *cb) {
	struct Orbit new_orbit;
	new_orbit.cb = cb;
	double ap, pe;
	if(apsis1 > apsis2) {
		ap = apsis1;
		pe = apsis2;
	} else {
		ap = apsis2;
		pe = apsis1;
	}
	new_orbit.a = (ap+pe)/2;
	new_orbit.i = inclination;
	new_orbit.e = (ap-pe)/(ap+pe);
	
	// not given, therefore zero
	new_orbit.raan = 0;
	new_orbit.arg_peri = 0;
	new_orbit.ta = 0;
	
	return new_orbit;
}

struct Orbit constr_orbit_from_osv(Vector3 r, Vector3 v, struct Body *cb) {
	double r_mag = mag_vec3(r);
	double v_mag = mag_vec3(v);
	double v_r = dot_vec3(v,r) / r_mag;
	double mu = cb->mu;
	
	double a = 1 / (2/r_mag - pow(v_mag,2)/mu);
	Vector3 h = cross_vec3(r,v);
	Vector3 e = scale_vec3(add_vec3(cross_vec3(v,h), scale_vec3(r, -mu/r_mag)), 1/mu);
	double e_mag = mag_vec3(e);
	
	Vector3 k = {0,0,1};
	Vector3 n_vec = cross_vec3(k, h);
	Vector3 n_norm = norm_vec3(n_vec);
	double RAAN, i, arg_peri;
	
	if(mag_vec3(n_vec) != 0) {
		RAAN = n_norm.y >= 0 ? acos(n_norm.x) : 2 * M_PI - acos(n_norm.x); // if n_norm.y is negative: raan > 180°
		i = acos(dot_vec3(k, norm_vec3(h)));
		double dp = dot_vec3(n_norm, e) / e_mag; if(dp > 1) dp = 1; if(dp < -1) dp = -1;	// if inside cos greater than 1 -> nan
		arg_peri = e.z >= 0 ? acos(dp) : 2 * M_PI - acos(dp);  // if r.z is positive: w > 180°
	} else {
		RAAN = 0;
		i = dot_vec3(k, norm_vec3(h)) > 0 ? 0 : M_PI;
		arg_peri = cross_vec3(r,v).z * e.y > 0 ? acos(e.x/e_mag) : 2*M_PI - acos(e.x/e_mag);
	}
	double dp = dot_vec3(e,r) / (e_mag*r_mag); if(dp > 1) dp = 1; if(dp < -1) dp = -1;	// if inside cos greater than 1 -> nan
	double ta = v_r >= 0 ? acos(dp) : 2*M_PI - acos(dp);
	
	struct Orbit new_orbit = constr_orbit_from_elements(a, e_mag, i, RAAN, arg_peri, ta, cb);
	
	return new_orbit;
}

Vector2 calc_vel_vec2(double r_mag, double v_mag, double true_anomaly, double flight_path_angle) {
	Vector2 r_norm = {cos(true_anomaly), sin(true_anomaly)};
	Vector2 r = scale_vec2(r_norm, r_mag);
	Vector2 n = {-r.y, r.x};
	Vector2 v = rotate_vec2(n, flight_path_angle);
	v = scale_vec2(v, v_mag/ mag_vec2(v));
	return v;
}

double calc_true_anomaly_from_mean_anomaly(struct Orbit orbit, double mean_anomaly) {
	// Solve Kepler's equation
	double ecc_anomaly = mean_anomaly; // Initial guess
	double delta;
	do {
		delta = (ecc_anomaly - orbit.e * sin(ecc_anomaly) - mean_anomaly) / (1 - orbit.e * cos(ecc_anomaly));
		ecc_anomaly -= delta;
	} while (fabs(delta) > 1e-6);
	
	// True anomaly from eccentric anomaly and eccentricity
	return 2 * atan(sqrt((1 + orbit.e) / (1 - orbit.e)) * tan(ecc_anomaly / 2));
}

Vector2 calc_orbital_speed_2d(double r_mag, double v_mag, double true_anomaly, double flight_path_angle) {
	Vector2 r_norm = {cos(true_anomaly), sin(true_anomaly)};
	Vector2 r = scale_vec2(r_norm, r_mag);
	Vector2 n = {-r.y, r.x};
	Vector2 v = rotate_vec2(n, flight_path_angle);
	v = scale_vec2(v, v_mag/ mag_vec2(v));
	
	return v;
}

double calc_orbit_flight_path_angle(double eccentricity, double true_anomaly) {
	return atan(eccentricity*sin(true_anomaly)/(1 + eccentricity*cos(true_anomaly)));
}

double calc_orbit_time_since_periapsis(Orbit orbit) {
	double n = sqrt(orbit.cb->mu / pow(fabs(orbit.a),3));
	double t;
	if(orbit.e < 1) {
		double E = 2*atan(sqrt((1 - orbit.e)/(1 + orbit.e))*tan(orbit.ta/2));
		t = (E - orbit.e*sin(E))/n;
		if(t < 0) {
			double T = 2*M_PI/n;
			t += T;
		}
	} else {
		double F = acosh((orbit.e + cos(orbit.ta)) / (1 + orbit.e * cos(orbit.ta)));
		t = (orbit.e * sinh(F) - F) / n;
		if(orbit.ta > M_PI) t *= -1;
	}
	return t;
}

double calc_orbital_period(Orbit orbit) {
	double n = sqrt(orbit.cb->mu / pow(fabs(orbit.a),3));
	if(orbit.e < 1) return 2*M_PI/n;
	else return INFINITY;
}

double calc_orbital_speed(struct Orbit orbit, double r) {
	double v2 = orbit.cb->mu * (2/r - 1/orbit.a);
	return sqrt(v2);
}

double calc_orbit_apoapsis(struct Orbit orbit) {
	return orbit.a*(1+orbit.e) - orbit.cb->radius;
}

double calc_orbit_periapsis(struct Orbit orbit) {
	return orbit.a*(1-orbit.e) - orbit.cb->radius;
}

Vector3 heliocentric_rot(Vector2 v, double raan, double argp, double incl) {
	double sin_raan = sin(raan);
	double sin_argp = sin(argp);
	double sin_incl = sin(incl);
	double cos_raan = cos(raan);
	double cos_argp = cos(argp);
	double cos_incl = cos(incl);
	
	double Q[3][3] = {
			{-sin_raan*cos_incl*sin_argp + cos_raan*cos_argp, -sin_raan*cos_incl*cos_argp - cos_raan*sin_argp,  sin_raan*sin_incl},
			{ cos_raan*cos_incl*sin_argp + sin_raan*cos_argp,  cos_raan*cos_incl*cos_argp - sin_raan*sin_argp, -cos_raan*sin_incl},
			{ sin_incl*sin_argp,                               sin_incl*cos_argp,                               cos_incl}};
	
	double v_vec[3] = {v.x, v.y, 0};
	double result[3] = {0,0,0};
	
	for(int i = 0; i < 3; i++) {
		for(int j = 0; j < 3; j++) {
			result[i] += Q[i][j]*v_vec[j];
		}
	}
	Vector3 result_v = {result[0], result[1], result[2]};
	return result_v;
}

OSV osv_from_orbit(Orbit orbit) {
	double flight_path_angle = calc_orbit_flight_path_angle(orbit.e, orbit.ta);
	double r_mag = orbit.a*(1-pow(orbit.e,2)) / (1+orbit.e*cos(orbit.ta));
	double v_mag = sqrt(orbit.cb->mu*(2/r_mag - 1/orbit.a));
	Vector2 r_2d = {cos(orbit.ta) * r_mag, sin(orbit.ta) * r_mag};
	Vector2 v_2d = calc_orbital_speed_2d(r_mag, v_mag, orbit.ta, flight_path_angle);
	
	Vector3 r = heliocentric_rot(r_2d, orbit.raan, orbit.arg_peri, orbit.i);
	Vector3 v = heliocentric_rot(v_2d, orbit.raan, orbit.arg_peri, orbit.i);
	
	OSV osv = {r, v};
	return osv;
}

OSV osv_from_elements(Orbit orbit, double epoch) {
	double dt = (epoch - orbit.cb->system->ut0) * (24 * 60 * 60);
	double e = orbit.e;
	if(e >= 0 && e < 1 && orbit.a > 0) {
		// Fast path: stay in eccentric anomaly and skip the true anomaly round trip.
		double mu = orbit.cb->mu;
		double n = sqrt(mu / (orbit.a * orbit.a * orbit.a));
		double sq_p = sqrt(1 + e), sq_m = sqrt(1 - e);
		double E0 = 2 * atan2(sq_m * sin(orbit.ta / 2), sq_p * cos(orbit.ta / 2));
		double M = fmod(E0 - e * sin(E0) + n * dt, 2 * M_PI);
		if(M < 0) M += 2 * M_PI;
		double E = e < 0.8 ? M : M_PI;
		int solved = 0;
		double sinE = 0, cosE = 1;
		for(int iteration = 0; iteration < 30; iteration++) {
			sinE = sin(E); cosE = cos(E);
			double next_E = E - (E - e * sinE - M) / (1 - e * cosE);
			if(!isfinite(next_E)) break;
			double change = fabs(next_E - E);
			E = next_E;
			if(change < 1e-14) { solved = 1; break; }
		}
		if(solved) {
			sinE = sin(E); cosE = cos(E);
			double a = orbit.a;
			double r_mag = a * (1 - e * cosE);
			double sq_mu_a = sqrt(mu * a) / r_mag;
			double sq_1e2 = sq_p * sq_m;
			Vector2 r_2d = {a * (cosE - e), a * sq_1e2 * sinE};
			Vector2 v_2d = {-sq_mu_a * sinE, sq_mu_a * sq_1e2 * cosE};
			OSV osv = {heliocentric_rot(r_2d, orbit.raan, orbit.arg_peri, orbit.i),
					heliocentric_rot(v_2d, orbit.raan, orbit.arg_peri, orbit.i)};
			return osv;
		}
	}
	orbit = propagate_orbit_time(orbit, dt);
	return osv_from_orbit(orbit);
}


Orbit propagate_orbit_time(Orbit orbit, double dt) {
	double ta = orbit.ta;
	double e = orbit.e;
	double a = orbit.a;
	double mu = orbit.cb->mu;
	double t = calc_orbit_time_since_periapsis(orbit);
	double target_t = t + dt;
	double T = calc_orbital_period(orbit);
	
	double n = sqrt(mu / pow(fabs(a),3));
	
	if(e >= 0 && e < 1 && isfinite(target_t) && isfinite(n)) {
		double reduced_t = target_t;
		while(reduced_t > T) reduced_t -= T;
		while(reduced_t < 0) reduced_t += T;
		double M = n * reduced_t;
		double E = e < 0.8 ? M : M_PI;
		int solved = 0;
		for(int iteration = 0; iteration < 30; iteration++) {
			double f = E - e * sin(E) - M;
			double next_E = E - f / (1 - e * cos(E));
			if(!isfinite(next_E)) break;
			double change = fabs(next_E - E);
			E = next_E;
			if(change < 1e-14) { solved = 1; break; }
		}
		if(solved) {
			orbit.ta = pi_norm(2 * atan2(sqrt(1 + e) * sin(E / 2), sqrt(1 - e) * cos(E / 2)));
			return orbit;
		}
	}
		double step = deg2rad(5);
	// if dt is basically 0, only add step, as this gets subtracted after the loop (not going inside loop)
	if(e<1) ta += fabs(t-target_t) > 1 ? dt/T * M_PI*2 : step;
	
	ta = pi_norm(ta);
	while(target_t > T && e < 1) target_t -= T;
	while(target_t < 0 && e < 1) target_t += T;
	
	int c = 0;
	
	while(fabs(t-target_t) > 1) {
		c++;
		// prevent endless loops (floating point imprecision can lead to not changing values for very small steps)
		if(c == 500) break;
		
		ta = pi_norm(ta);
		if(e < 1) {
			double E = acos((e + cos(ta))/(1 + e*cos(ta)));
			t = (E - e*sin(E))/n;
			if(ta > M_PI) t = T-t;
		} else {
			//printf("[%f %f %f %f]\n", t/(24*60*60), target_t/(24*60*60), (target_t-t)/(24*60*60), rad2deg(theta));
			double F = acosh((e + cos(ta))/(1 + e*cos(ta)));
			t = (e*sinh(F) - F)/n;
			if(ta > M_PI) t *= -1;
			if(isnan(t)) {
				step /= 2;
				ta -= step;
				t = 100;	// to not exit the loop;
				continue;
			}
		}
		
		// check in which half t is with respect to target_t (forwards or backwards from target_t) and move it closer
		if(target_t < T/2 || e > 1) {
			if(t > target_t && (t < target_t+T/2  || e > 1)) {
				if (step > 0) step *= -1.0 / 4;
			} else {
				if (step < 0) step *= -1.0 / 4;
			}
		} else {
			if(t < target_t && t > target_t-T/2) {
				if (step < 0) step *= -1.0 / 4;
			} else {
				if (step > 0) step *= -1.0 / 4;
			}
		}
		ta += step;
	}
	ta -= step; // reset theta1 from last change inside the loop
	orbit.ta = ta;
	return orbit;
}

OSV propagate_osv_time(OSV osv, Body *cb, double dt) {
	Orbit orbit = constr_orbit_from_osv(osv.r, osv.v, cb);
	orbit = propagate_orbit_time(orbit, dt);
	return osv_from_orbit(orbit);
}

OSV propagate_osv_ta(OSV osv, Body *cb, double delta_ta) {
	Orbit orbit = constr_orbit_from_osv(osv.r, osv.v, cb);
	orbit.ta = pi_norm(orbit.ta+delta_ta);
	return osv_from_orbit(orbit);
}



// Printing info #######################################################

void print_orbit_info(struct Orbit orbit) {
	struct Body *body = orbit.cb;
	printf("\n______________________\nORBIT:\n\n");
	printf("Orbiting: \t\t%s\n", body->name);
	printf("Apoapsis:\t\t%g km\n", (calc_orbit_apoapsis(orbit))/1000);
	printf("Periapsis:\t\t%g km\n", (calc_orbit_periapsis(orbit))/1000);
	printf("Semi-major axis:\t%g km\n", orbit.a /1000);
	printf("Inclination:\t\t%g°\n", rad2deg(orbit.i));
	printf("Eccentricity:\t\t%g\n", orbit.e);
	printf("RAAN:\t\t\t\t%g°\n", rad2deg(orbit.raan));
	printf("Arg of Periapsis:\t%g°\n", rad2deg(orbit.arg_peri));
	printf("Orbital Period:\t\t%gs\n", calc_orbital_period(orbit));
	printf("______________________\n\n");
}
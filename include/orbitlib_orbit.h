#ifndef ORBITLIB_ORBITLIB_ORBIT_H
#define ORBITLIB_ORBITLIB_ORBIT_H

#include "geometrylib.h"

typedef struct Body Body;

/*
 * ------------------------------------
 * Orbit & State Vector
 * ------------------------------------
 */

/**
 * @brief Represents an orbital path around a central body using classical orbital elements
 */
typedef struct Orbit {
	struct Body * cb;	/**< Central body being orbited */
	double e;           /**< Eccentricity of the orbit */
	double a;           /**< Semi-major axis [m] */
	double i; 			/**< Inclination [radians] */
	double raan;        /**< Right ascension of ascending node [radians] */
	double arg_peri; 	/**< Argument of periapsis [radians] */
	double ta;   		/**< True anomaly [radians] */
} Orbit;

/**
 * @brief Represents an orbital state vector (position and velocity)
 */
typedef struct OSV {
	Vector3 r; /**< Position vector [m] */
	Vector3 v; /**< Velocity vector [m/s] */
} OSV;


/*
 * ------------------------------------
 * Orbit Construction
 * ------------------------------------
 */

/**
 * @brief Constructs an orbit using classical orbital elements
 *
 * @param semimajor_axis Semi-major axis [m]
 * @param eccentricity Orbit eccentricity
 * @param inclination Inclination [radians]
 * @param raan Right ascension of ascending node [radians]
 * @param arg_of_peri Argument of periapsis [radians]
 * @param true_anomaly True anomaly [radians]
 * @param cb Pointer to the central body
 * @return Constructed Orbit struct
 */
struct Orbit constr_orbit_from_elements(double semimajor_axis, double eccentricity, double inclination, double raan, double arg_of_peri, double true_anomaly, struct Body *cb);

/**
 * @brief Constructs an orbit using the two apsides and inclination
 *
 * @param apsis1 Radius of the first apsis [m]
 * @param apsis2 Radius of the second apsis [m]
 * @param inclination Inclination [radians]
 * @param cb Pointer to the central body
 * @return Constructed Orbit struct
 */
struct Orbit constr_orbit_from_apsides(double apsis1, double apsis2, double inclination, Body *cb);

/**
 * @brief Constructs an orbit from an orbital state vector (position and velocity)
 *
 * @param r Position vector [m]
 * @param v Velocity vector [m/s]
 * @param cb Pointer to the central body
 * @return Constructed Orbit struct
 */
struct Orbit constr_orbit_from_osv(Vector3 r, Vector3 v, Body *cb);


/*
 * ------------------------------------
 * Orbit Calculations
 * ------------------------------------
 */

/**
 * @brief Calculates a 2D velocity vector from magnitudes and angles
 *
 * @param r_mag Magnitude of the position vector [m]
 * @param v_mag Magnitude of the velocity vector [m/s]
 * @param true_anomaly True anomaly [radians]
 * @param flight_path_angle Flight path angle [radians]
 * @return Velocity vector in 2D (Vector2)
 */
Vector2 calc_vel_vec2(double r_mag, double v_mag, double true_anomaly, double flight_path_angle);

/**
 * @brief Calculates the true anomaly from the mean anomaly for a given orbit
 *
 * @param orbit Orbit struct containing orbital elements
 * @param mean_anomaly Mean anomaly [radians]
 * @return True anomaly [radians]
 */
double calc_true_anomaly_from_mean_anomaly(struct Orbit orbit, double mean_anomaly);

/**
 * @brief Calculates the flight path angle of an orbit at a given true anomaly
 *
 * @param eccentricity Orbit eccentricity
 * @param true_anomaly True anomaly [radians]
 * @return Flight path angle [radians]
 */
double calc_orbit_flight_path_angle(double eccentricity, double true_anomaly);

/**
 * @brief Calculates time elapsed since periapsis passage for the orbit
 *
 * @param orbit Orbit struct containing orbital elements
 * @return Time since periapsis [s]
 */
double calc_orbit_time_since_periapsis(Orbit orbit);

/**
 * @brief Calculates the orbital period for the orbit
 *
 * @param orbit Orbit struct containing orbital elements
 * @return Orbital period [s]
 */
double calc_orbital_period(Orbit orbit);

/**
 * @brief Calculates the orbital speed at a given radius
 *
 * @param orbit Orbit struct containing orbital elements
 * @param r Distance from central body [m]
 * @return Orbital speed at radius r [m/s]
 */
double calc_orbital_speed(struct Orbit orbit, double r);

/**
 * @brief Calculates the apoapsis distance (farthest point) of the orbit
 *
 * @param orbit Orbit struct containing orbital elements
 * @return Apoapsis radius [m]
 */
double calc_orbit_apoapsis(struct Orbit orbit);

/**
 * @brief Calculates the periapsis distance (closest point) of the orbit
 *
 * @param orbit Orbit struct containing orbital elements
 * @return Periapsis radius [m]
 */
double calc_orbit_periapsis(struct Orbit orbit);


/*
 * ------------------------------------
 * Orbital State Vector (OSV) Construction
 * ------------------------------------
 */

/**
 * @brief Rotates a 2D vector into heliocentric coordinates using orbital angles
 *
 * @param v 2D vector to rotate
 * @param raan Right ascension of ascending node [radians]
 * @param argp Argument of periapsis [radians]
 * @param incl Inclination [radians]
 * @return Rotated 3D vector in heliocentric coordinates (Vector3)
 */
Vector3 heliocentric_rot(Vector2 v, double raan, double argp, double incl);

/**
 * @brief Constructs an orbital state vector (position and velocity) from an Orbit struct
 *
 * @param orbit Orbit struct containing orbital elements
 * @return Orbital state vector (OSV) with position and velocity
 */
OSV osv_from_orbit(Orbit orbit);

/**
 * @brief Constructs an orbital state vector (position and velocity) from orbital elements and epoch
 *
 * @param orbit Orbit struct containing orbital elements
 * @param epoch Time at which to compute the OSV (Julian date)
 * @return Orbital state vector (OSV) with position and velocity
 */
OSV osv_from_elements(Orbit orbit, double epoch);

/**
 * @brief Constructs an orbital state vector dt seconds after the state described by the orbital elements
 *
 * Elliptic orbits use a cached eccentric-anomaly propagation; other orbits fall back to
 * propagate_orbit_time. osv_from_elements uses this with dt measured from the system epoch.
 *
 * @param orbit Orbit struct containing orbital elements
 * @param dt Time after the elements' epoch [s]
 * @return Orbital state vector (OSV) with position and velocity
 */
OSV osv_from_elements_after(const Orbit *orbit, double dt);

/*
 * ------------------------------------
 * Orbit Propagation
 * ------------------------------------
 */

/**
 * @brief Propagates an orbit forward in time by a given duration
 *
 * @param orbit Initial orbit
 * @param dt Time step to propagate [s]
 * @return Orbit propagated by dt seconds
 */
Orbit propagate_orbit_time(Orbit orbit, double dt);

/**
 * @brief Propagates an orbital state vector forward in time by a given duration
 *
 * @param osv Initial orbital state vector
 * @param cb Central body of the orbit
 * @param dt Time step to propagate [s]
 * @return Orbital state vector propagated by dt seconds
 */
OSV propagate_osv_time(OSV osv, Body *cb, double dt);

/**
 * @brief Propagates an orbital state vector by a change in true anomaly
 *
 * @param osv Initial orbital state vector
 * @param cb Central body of the orbit
 * @param delta_ta Change in true anomaly [radians]
 * @return Orbital state vector propagated by delta_ta radians
 */
OSV propagate_osv_ta(OSV osv, Body *cb, double delta_ta);

/*
 * ------------------------------------
 * Debug / Print
 * ------------------------------------
 */

/**
 * @brief Prints detailed information about an orbit to the console
 *
 * @param orbit Orbit struct containing orbital elements
 */
void print_orbit_info(struct Orbit orbit);

#endif //ORBITLIB_ORBITLIB_ORBIT_H

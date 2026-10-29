#ifndef ORBITLIB_ORBITLIB_EPHEMERIS_H
#define ORBITLIB_ORBITLIB_EPHEMERIS_H

#include "geometrylib.h"
#include "orbitlib_celestial.h"
#include "orbitlib_datetime.h"

/**
 * @brief Represents the ephemeral data consisting of epoch, position, and velocity components
 */
typedef struct Ephem {
	double epoch;	/**< Epoch associated with the ephemeral data (Julian Date) */
	Vector3 r;		/**< Position Vector */
	Vector3 v;		/**< Velocity Vector */
} Ephem;

/**
 * @brief Prints the date, position and velocity vector of the given ephemeris
 *
 * @param ephem The given ephemeris
 */
void print_ephem(struct Ephem ephem);


/**
 * @brief Retrieves the ephemeral data of requested body for requested time (from JPL's Horizon API or from file)
 *
 * @param body The body for which the ephemerides should be stored
 * @param central_body The central body of the body's system
 */
void get_body_ephems(Body *body, Datetime min_date, Datetime max_date, Datetime time_step, const char *ephem_directory);

/**
 * @brief Returns the ephemeris entry whose epoch is nearest to the given epoch
 *
 * The list must be sorted ascending by epoch. Ties go to the later entry.
 */
Ephem get_closest_ephem(Ephem *ephem_list, int num_ephems, double epoch);

/**
 * @brief Interpolates an ephemeris list to get the state vector at a given epoch
 *
 * Uses the provided list of ephemeris entries to compute the OSV (position and velocity)
 * of a body at the specified epoch, relative to the given central body.
 *
 * @param ephem_list Array of ephemeris entries
 * @param num_ephems Number of entries in the ephemeris list
 * @param epoch Time at which to compute the state (Julian Date)
 * @param cb Pointer to the central body
 * @return OSV (position and velocity) of the body at the given epoch
 */
OSV osv_from_ephem(Ephem *ephem_list, int num_ephems, double epoch, Body *cb);

#endif //ORBITLIB_ORBITLIB_EPHEMERIS_H

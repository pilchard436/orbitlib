#ifndef ORBITLIB_ORBITLIB_TRANSFER_H
#define ORBITLIB_ORBITLIB_TRANSFER_H

#include "orbitlib_celestial.h"
#include "orbitlib_orbit.h"


/*
 * ------------------------------------
 * Transfer and Trajectory Types
 * ------------------------------------
 */

/**
 * @brief Represents a Hohmann transfer maneuver with duration and delta-vs
 */
typedef struct Hohmann {
	double dur;      /**< Duration of transfer [s] */
	double dv_dep;   /**< Delta-v at departure [m/s] */
	double dv_arr;   /**< Delta-v at arrival [m/s] */
} Hohmann;

/**
 * @brief Enumeration of transfer types between circularization and capture orbits
 */
enum Transfer_Type {
	capcap,   /**< Capture to capture */
	circcap,  /**< Circular to capture */
	capcirc,  /**< Capture to circular */
	circcirc, /**< Circular to circular */
	capfb,    /**< Capture with flyby */
	circfb    /**< Circular with flyby */
};

/**
 * @brief Enumeration of possible Lambert solver result statuses
 */
enum LAMBERT_SOLVER_SUCCESS {
	LAMBERT_SUCCESS,         /**< Solver converged successfully */
	LAMBERT_IMPRECISION,     /**< Solution found but with reduced precision */
	LAMBERT_MAX_ITERATIONS,  /**< Maximum iterations exceeded */
	LAMBERT_FAIL_NAN,        /**< Result contains NaN */
	LAMBERT_FAIL_ECC         /**< Invalid eccentricity */
};

enum LambertSolverMethod {
	LAMBERT_METHOD_LEGACY,
	LAMBERT_METHOD_IZZO
};

void set_lambert_solver_method(enum LambertSolverMethod method);
enum LambertSolverMethod get_lambert_solver_method();

/**
 * @brief Lambert solution for 2D orbits including true anomalies and success status
 */
typedef struct Lambert2 {
	Orbit orbit;                         /**< Orbit resulting from Lambert solution */
	double true_anomaly0;                /**< Initial true anomaly [radians] */
	double true_anomaly1;                /**< Final true anomaly [radians] */
	enum LAMBERT_SOLVER_SUCCESS success; /**< Status of Lambert solver */
} Lambert2;

/**
 * @brief Lambert solution for 3D trajectories including position and velocity vectors
 */
typedef struct Lambert3 {
	Vector3 r0, v0;                     /**< Initial position and velocity vectors */
	Vector3 r1, v1;                     /**< Final position and velocity vectors */
	enum LAMBERT_SOLVER_SUCCESS success; /**< Status of Lambert solver */
} Lambert3;

typedef int (*LambertCancelCheck)(void *context);

typedef struct PoweredFlybyResult {
	double periapsis_radius;
	double burn_delta_v;
	Vector3 burn_delta_v_vector;
	Vector3 periapsis_direction;
	Vector3 incoming_tangent_direction;
	Vector3 outgoing_tangent_direction;
	double soi_radius;
	double incoming_soi_time;
	double outgoing_soi_time;
	OSV incoming_soi_state;
	OSV outgoing_soi_state;
} PoweredFlybyResult;

typedef struct DepartureBurn3D {
	Vector3 burn_delta_v;
	Vector3 periapsis_direction;
	Vector3 tangential_direction;
	double burn_magnitude;
} DepartureBurn3D;

/**
 * @brief Enumeration of hyperbolic transfer orbit types
 */
enum HyperbolaType {
	HYP_DEPARTURE, /**< Hyperbola for departure maneuver */
	HYP_ARRIVAL,   /**< Hyperbola for arrival maneuver */
	HYP_FLYBY      /**< Hyperbola for flyby maneuver */
};

/**
 * @brief Parameters describing incoming or outgoing leg of a hyperbola maneuver
 */
typedef struct HyperbolaLegParams {
	double decl;        /**< Declination angle [radians] */
	double bplane_angle; /**< B-plane angle [radians] */
	double bvazi;       /**< B-vector azimuth [radians] */
} HyperbolaLegParams;

/**
 * @brief Parameters for a hyperbola maneuver including periapsis, C3 energy, and legs (declination, b-plane angle and bvazi)
 */
typedef struct HyperbolaParams {
	enum HyperbolaType type;            /**< Type of hyperbola maneuver */
	double rp;                         /**< Periapsis radius [m] */
	double c3_energy;                  /**< Characteristic energy (C3) [m²/s²] */
	struct HyperbolaLegParams incoming; /**< Incoming leg parameters (ignored if type == HYP_DEPARTURE) */
	struct HyperbolaLegParams outgoing; /**< Outgoing leg parameters (ignored if type == HYP_ARRIVAL) */
} HyperbolaParams;


/*
 * ------------------------------------
 * Maneuver Delta-V Calculations
 * ------------------------------------
 */

/**
 * @brief Calculates delta-v required for an apsis maneuver changing from initial to new apsis at a static apsis
 *
 * @param static_apsis The apsis that remains unchanged during the maneuver [m]
 * @param initial_apsis The current apsis radius [m]
 * @param new_apsis The desired new apsis radius [m]
 * @param cb Pointer to the central body influencing the orbit
 * @return Required delta-v for the maneuver [m/s]
 */
double calc_apsis_maneuver_dv(double static_apsis, double initial_apsis, double new_apsis, struct Body *cb);

/**
 * @brief Calculates delta-v needed to insert into a circular orbit from hyperbolic approach at periapsis radius
 *
 * @param body Pointer to the central body
 * @param rp Radius of periapsis of the orbit [m]
 * @param vinf Hyperbolic excess speed (speed at infinity) [m/s]
 * @return Delta-v required for circularization [m/s]
 */
double dv_circ(struct Body *body, double rp, double vinf);

/**
 * @brief Calculates delta-v needed to capture into orbit from hyperbolic approach at periapsis radius
 *
 * @param body Pointer to the central body
 * @param rp Radius of periapsis of the orbit [m]
 * @param vinf Hyperbolic excess speed (speed at infinity) [m/s]
 * @return Delta-v required for capture [m/s]
 */
double dv_capture(struct Body *body, double rp, double vinf);
int calculate_free_plane_departure_burn(Vector3 v_inf, Body *body, double periapsis_radius,
		DepartureBurn3D *result);


/*
 * ------------------------------------
 * Transfer Orbit Calculations
 * ------------------------------------
 */

/**
 * @brief Calculates parameters for a Hohmann transfer between two circular orbits around a central body
 *
 * @param r0 Radius of the initial circular orbit [m]
 * @param r1 Radius of the target circular orbit [m]
 * @param cb Pointer to the central body
 * @return Struct containing duration and delta-v values for the transfer
 */
Hohmann calc_hohmann_transfer(double r0, double r1, Body *cb);

/**
 * @brief Computes a 2D Lambert solution for an orbit transfer between two radii with a change in true anomaly over a target time
 *
 * @param r0 Radius of starting position [m]
 * @param r1 Radius of ending position [m]
 * @param delta_ta Change in true anomaly between positions [radians]
 * @param target_dt Desired transfer time [s]
 * @param cb Pointer to the central body
 * @return Lambert2 struct containing orbit and solver status
 */
Lambert2 calc_lambert2(double r0, double r1, double delta_ta, double target_dt, Body *cb);

/**
 * @brief Computes a 3D Lambert solution for an orbit transfer between two position vectors over a target time
 *
 * @param r0 Initial position vector [m]
 * @param r1 Final position vector [m]
 * @param target_dt Desired transfer time [s]
 * @param cb Pointer to the central body
 * @return Lambert3 struct containing position, velocity vectors and solver status
 */
Lambert3 calc_lambert3(Vector3 r0, Vector3 r1, double target_dt, Body *cb);
Lambert3 calc_lambert3_selected(Vector3 r0, Vector3 r1, double target_dt, Body *cb);
Lambert3 calc_lambert3_selected_cancellable(Vector3 r0, Vector3 r1, double target_dt, Body *cb,
		LambertCancelCheck cancel_check, void *cancel_context);

/**
 * @brief Computes a zero-revolution Lambert solution using Izzo's formulation
 *
 * @param r0 Initial position vector [m]
 * @param r1 Final position vector [m]
 * @param target_dt Desired transfer time [s]
 * @param cb Pointer to the central body
 * @return Lambert3 containing the endpoint velocities and solver status
 */
Lambert3 calc_lambert3_izzo(Vector3 r0, Vector3 r1, double target_dt, Body *cb);


/*
 * ------------------------------------
 * Flyby and Hyperbolic Orbit Calculations
 * ------------------------------------
 */

/**
 * @brief Computes periapsis radius of a flyby hyperbola based on incoming, outgoing, and central body velocity vectors,
 * all given in the same inertial reference frame
 *
 * @param v_arr Velocity vector at arrival [m/s]
 * @param v_dep Velocity vector at departure [m/s]
 * @param v_body Velocity vector of the central body [m/s]
 * @param body Pointer to the central body
 * @return Periapsis radius of the flyby hyperbola [m]
 */
double get_flyby_periapsis(Vector3 v_arr, Vector3 v_dep, Vector3 v_body, Body *body);
int calculate_powered_flyby(Vector3 v_arr, Vector3 v_dep, Vector3 v_body, Body *body,
		double minimum_periapsis_radius, PoweredFlybyResult *result);

/**
 * @brief Calculates inclination angle of a flyby based on incoming, outgoing, and central body velocity vectors,
 * all given in the same inertial reference frame
 *
 * @param v_arr Velocity vector at arrival [m/s]
 * @param v_dep Velocity vector at departure [m/s]
 * @param v_body Velocity vector of the central body [m/s]
 * @return Flyby inclination angle [radians]
 */
double get_flyby_inclination(Vector3 v_arr, Vector3 v_dep, Vector3 v_body, Plane3 equatorial_plane);

/**
 * @brief Calculates parameters of a hyperbolic orbit (departure, arrival, or flyby) given velocity vectors and periapsis radius.
 * Incoming and outgoing velocity vectors and central body velocity vector are all in the same inertial reference frame.
 * The incoming vector is ignored if type == HYP_DEPARTURE; the outgoing vector is ignored if type == HYP_ARRIVAL;
 * The radius of periapsis is ignored if type == HYP_FLYBY
 *
 * @param v_arr Velocity vector at arrival [m/s]
 * @param v_dep Velocity vector at departure [m/s]
 * @param v_body Velocity vector of the central body [m/s]
 * @param body Pointer to the central body
 * @param rp Radius of periapsis of the hyperbolic orbit [m]
 * @param type Type of hyperbolic orbit (departure, arrival, or flyby)
 * @return Struct containing hyperbolic orbit parameters such as periapsis radius, C3 energy, and leg parameters
 */
HyperbolaParams get_hyperbola_params(Vector3 v_arr, Vector3 v_dep, Vector3 v_body, struct Body *body, double rp, enum HyperbolaType type);

/**
 * @brief Determines if a flyby is viable by comparing the difference in excess speed between incoming and outgoing velocity vectors,
 * and the radius of periapsis (not being inside body or body's atmosphere)
 * All velocity vectors are in the same inertial reference frame.
 *
 * @param v_arr Velocity vector at arrival [m/s]
 * @param v_dep Velocity vector at departure [m/s]
 * @param v_body Velocity vector of the central body [m/s]
 * @param body Pointer to the central body
 * @param precision Allowed difference in excess speed (vinf) for the flyby to be considered viable [m/s]
 * @return true if the flyby maneuver is viable, false otherwise
 */
bool is_flyby_viable(Vector3 v_arr, Vector3 v_dep, Vector3 v_body, Body *body, double precision);


#endif //ORBITLIB_ORBITLIB_TRANSFER_H

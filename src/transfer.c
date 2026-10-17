#include "orbitlib_transfer.h"
#include "geometrylib.h"
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

static int is_finite_vector(Vector3 value) {
	return isfinite(value.x) && isfinite(value.y) && isfinite(value.z);
}

#define IZZO_TIME_OF_FLIGHT_TOLERANCE 1e-12

double calc_apsis_maneuver_dv(double static_apsis, double initial_apsis, double new_apsis, struct Body *cb) {
	struct Orbit initial_orbit;
	struct Orbit new_orbit;
	
	initial_orbit = constr_orbit_from_apsides(static_apsis, initial_apsis, 0, cb);
	new_orbit = constr_orbit_from_apsides(static_apsis, new_apsis, 0, cb);
	
	double v0 = calc_orbital_speed(initial_orbit, static_apsis);
	double v1 = calc_orbital_speed(new_orbit, static_apsis);
	
	return fabs(v1-v0);
}

double dv_circ(struct Body *body, double rp, double vinf) {
	return sqrt(2 * body->mu / rp + vinf * vinf) - sqrt(body->mu / rp);
}

double dv_capture(struct Body *body, double rp, double vinf) {
	return sqrt(2 * body->mu / rp + vinf * vinf) - sqrt(2 * body->mu / rp);
}

int calculate_free_plane_departure_burn(Vector3 v_inf, Body *body, double periapsis_radius,
		DepartureBurn3D *result) {
	if(body == NULL || result == NULL || !(body->mu > 0) || !(periapsis_radius > 0)) return 0;
	double v_inf_magnitude = mag_vec3(v_inf);
	if(!(v_inf_magnitude > 0) || !isfinite(v_inf_magnitude)) return 0;
	Vector3 asymptote_direction = scale_vec3(v_inf, 1 / v_inf_magnitude);
	Vector3 reference_axis = fabs(asymptote_direction.z) < 0.9 ? vec3(0, 0, 1) : vec3(0, 1, 0);
	Vector3 plane_normal = norm_vec3(cross_vec3(asymptote_direction, reference_axis));
	double radial_component = 1 / (1 + periapsis_radius * v_inf_magnitude * v_inf_magnitude / body->mu);
	double tangent_component = sqrt(fmax(0, 1 - radial_component * radial_component));
	Vector3 periapsis_direction = add_vec3(scale_vec3(asymptote_direction, radial_component),
		scale_vec3(plane_normal, tangent_component));
	Vector3 parking_tangent = norm_vec3(subtract_vec3(asymptote_direction,
		scale_vec3(periapsis_direction, radial_component)));
	double hyperbolic_periapsis_speed = sqrt(v_inf_magnitude * v_inf_magnitude + 2 * body->mu / periapsis_radius);
	double parking_orbit_speed = sqrt(body->mu / periapsis_radius);
	double burn_magnitude = hyperbolic_periapsis_speed - parking_orbit_speed;
	*result = (DepartureBurn3D) {
		.burn_delta_v = scale_vec3(parking_tangent, burn_magnitude),
		.periapsis_direction = periapsis_direction,
		.tangential_direction = parking_tangent,
		.burn_magnitude = fabs(burn_magnitude)
	};
	return isfinite(result->burn_magnitude) && isfinite(result->burn_delta_v.x) &&
		isfinite(result->burn_delta_v.y) && isfinite(result->burn_delta_v.z);
}

Hohmann calc_hohmann_transfer(double r0, double r1, struct Body *cb) {
	double sma_pow_3 = pow(((r0 + r1) / 2),3);
	double dur = M_PI * sqrt(sma_pow_3 / cb->mu);
	double dv_dep = calc_apsis_maneuver_dv(r0, r0, r1, cb);
	double dv_arr = calc_apsis_maneuver_dv(r1, r0, r1, cb);
	return (Hohmann) {dur, dv_dep, dv_arr};
}

// dt --> 0
double departure_true_anomaly_at_min_dt(double r0, double r1, double delta_ta) {
	delta_ta = pi_norm(delta_ta);
	double max_arr_ta;
	if(delta_ta < M_PI) {
		double r = r1/r0;
		double r1r2 = sqrt(1 + r*r - 2 * r * cos(delta_ta));
		double beta = acos((1 + r1r2 * r1r2 - r * r) / (2 * r1r2));
		double alpha = M_PI / 2 - beta;
		max_arr_ta = 2*M_PI-alpha;
	} else {
		max_arr_ta = 2*M_PI-delta_ta/2;
	}
	return pi_norm(max_arr_ta);
}

// dt --> inf
double departure_true_anomaly_at_max_dt(double r0, double r1, double delta_ta) {
	double min_arr_ta;
	Vector2 p1 = {r0, 0};
	Vector2 p2 = {cos(delta_ta)*r1, sin(delta_ta)*r1};
	double pxr = (p2.x-p1.x)/(r1 - r0);
	double pyr = (p2.y-p1.y)/(r1 - r0);
	double p = (2*pxr)/(pxr*pxr + pyr*pyr);
	double q = (1-pyr*pyr)/(pxr*pxr + pyr*pyr);
	
	// IF NAN IN 2D TRANSFER CALC, LOOK HERE. SQRT SHOULD NOT GET NEGATIVE
	double inside_sqrt = p*p/4 - q;
	if(inside_sqrt < 0) inside_sqrt *= -1;
	
	double mx1 = -p/2 - sqrt(inside_sqrt);
	double mx2 = -p/2 + sqrt(inside_sqrt);
	
	Vector2 m_n[] = {
			{mx1, +sqrt(1-mx1*mx1)},
			{mx1, -sqrt(1-mx1*mx1)},
			{mx2, +sqrt(1-mx2*mx2)},
			{mx2, -sqrt(1-mx2*mx2)}
	};
	
	double m_ml_dot[4];
	
	// dot product of m and ml should be 0, due to imprecision not accurately 0 -> find smallest of the pairs
	// pair (0,1) and (2,3) (see above; two times two solutions of square root)
	for(int i = 0; i < 4; i++) m_ml_dot[i] = m_n[i].x*m_n[i].x+pxr*m_n[i].x+m_n[i].y*m_n[i].y+pyr*m_n[i].y;
	if(fabs(m_ml_dot[1]) < fabs(m_ml_dot[0])) m_n[0] = m_n[1];
	m_n[1] = fabs(m_ml_dot[2]) < fabs(m_ml_dot[3]) ? m_n[2] : m_n[3];
	
	for(int i = 0; i < 2; i++) {
		min_arr_ta = angle_ccw_vec2_vec2(p1, m_n[i]);
		double min_theta2 = angle_ccw_vec2_vec2(p2, m_n[i]);
		if(min_arr_ta > min_theta2) {
			break;
		}
	}
	return pi_norm(min_arr_ta);
}

static double izzo_hypergeometric(double z) {
	double term = 1;
	double sum = 1;
	for(int j = 0; j < 100; j++) {
		term *= ((double)j + 3) / ((double)j + 2.5) * z;
		sum += term;
		if(fabs(term) < 1e-13) break;
	}
	return sum;
}

static double izzo_time_of_flight(double x, double lambda) {
	double lambda2 = lambda * lambda;
	double eccentricity_parameter = x * x - 1;
	double rho = fabs(eccentricity_parameter);
	double z = sqrt(1 + lambda2 * eccentricity_parameter);
	if(fabs(x - 1) < 0.01) {
		double eta = z - lambda * x;
		double series_argument = 0.5 * (1 - lambda - x * eta);
		double q = 4.0 / 3.0 * izzo_hypergeometric(series_argument);
		return (eta * eta * eta * q + 4 * lambda * eta) / 2;
	}

	double g = x * z - lambda * eccentricity_parameter;
	if(eccentricity_parameter < 0) {
		double angle = acos(fmax(-1.0, fmin(1.0, g)));
		return (x - lambda * z - angle / sqrt(rho)) / eccentricity_parameter;
	}
	double y = sqrt(rho);
	double f = y * (z - lambda * x);
	return (x - lambda * z - log(f + g) / y) / eccentricity_parameter;
}

static int izzo_time_of_flight_derivatives(double x, double lambda, double time_of_flight,
		double *first, double *second, double *third) {
	double one_minus_x2 = 1 - x * x;
	if(fabs(one_minus_x2) < 1e-5) {
		double h = 1e-4;
		double f_minus_two = izzo_time_of_flight(x - 2 * h, lambda);
		double f_minus_one = izzo_time_of_flight(x - h, lambda);
		double f_plus_one = izzo_time_of_flight(x + h, lambda);
		double f_plus_two = izzo_time_of_flight(x + 2 * h, lambda);
		double f_center = time_of_flight;
		*first = (f_minus_two - 8 * f_minus_one + 8 * f_plus_one - f_plus_two) / (12 * h);
		*second = (-f_plus_two + 16 * f_plus_one - 30 * f_center + 16 * f_minus_one - f_minus_two) / (12 * h * h);
		*third = (f_plus_two - 2 * f_plus_one + 2 * f_minus_one - f_minus_two) / (2 * h * h * h);
		return isfinite(*first) && isfinite(*second) && isfinite(*third) && fabs(*first) > DBL_MIN;
	}

	double lambda2 = lambda * lambda;
	double lambda3 = lambda2 * lambda;
	double y = sqrt(1 - lambda2 * one_minus_x2);
	double y2 = y * y;
	double y3 = y2 * y;
	*first = (3 * time_of_flight * x - 2 + 2 * lambda3 * x / y) / one_minus_x2;
	*second = (3 * time_of_flight + 5 * x * *first + 2 * (1 - lambda2) * lambda3 / y3) / one_minus_x2;
	*third = (7 * x * *second + 8 * *first -
		6 * (1 - lambda2) * lambda2 * lambda3 * x / (y3 * y2)) / one_minus_x2;
	return isfinite(*first) && isfinite(*second) && isfinite(*third) && fabs(*first) > DBL_MIN;
}

Lambert3 calc_lambert3_cancellable(Vector3 r0, Vector3 r1, double target_dt, Body *cb,
		LambertCancelCheck cancel_check, void *cancel_context) {
	Lambert3 failure = {.r0 = r0, .r1 = r1, .success = LAMBERT_FAIL_NAN};
	if(cancel_check != NULL && cancel_check(cancel_context)) return failure;
	if(cb == NULL || !(cb->mu > 0) || !(target_dt > 0) || !isfinite(target_dt)) return failure;
	double r0_magnitude = mag_vec3(r0);
	double r1_magnitude = mag_vec3(r1);
	if(!(r0_magnitude > 0) || !(r1_magnitude > 0)) return failure;

	Vector3 direction = cross_vec3(r0, r1);
	double direction_magnitude = mag_vec3(direction);
	if(!(direction_magnitude > 0) || !isfinite(direction_magnitude)) return failure;
	direction = scale_vec3(direction, 1 / direction_magnitude);

	Vector3 radial0 = scale_vec3(r0, 1 / r0_magnitude);
	Vector3 radial1 = scale_vec3(r1, 1 / r1_magnitude);
	Vector3 tangent0 = cross_vec3(direction, radial0);
	Vector3 tangent1 = cross_vec3(direction, radial1);
	double chord = mag_vec3(subtract_vec3(r1, r0));
	if(!(chord > 0)) return failure;
	double semiperimeter = (r0_magnitude + r1_magnitude + chord) / 2;
	double lambda_magnitude = sqrt(fmax(0, 1 - chord / semiperimeter));
	double lambda = direction.z < 0 ? -lambda_magnitude : lambda_magnitude;
	if(direction.z < 0) {
		tangent0 = scale_vec3(tangent0, -1);
		tangent1 = scale_vec3(tangent1, -1);
	}

	double lambda2 = lambda * lambda;
	double lambda3 = lambda2 * lambda;
	double dimensionless_time = sqrt(2 * cb->mu / (semiperimeter * semiperimeter * semiperimeter)) * target_dt;
	if(!(dimensionless_time > 0) || !isfinite(dimensionless_time)) return failure;
	double parabolic_time = acos(lambda) + lambda * sqrt(1 - lambda2);
	double short_time = 2.0 / 3.0 * (1 - lambda3);
	double x;
	if(dimensionless_time >= parabolic_time) {
		double difference = dimensionless_time - parabolic_time;
		x = -difference / (difference + 4);
	} else if(dimensionless_time <= short_time) {
		double denominator = 0.4 * (1 - lambda2 * lambda3) * dimensionless_time;
		if(fabs(denominator) < DBL_MIN) return failure;
		x = short_time * (short_time - dimensionless_time) / denominator + 1;
	} else {
		x = pow(dimensionless_time / parabolic_time,
			0.69314718055994530942 / log(short_time / parabolic_time)) - 1;
	}
	if(!isfinite(x) || x <= -1) return failure;

	int converged = 0;
	for(int iteration = 0; iteration < 15; iteration++) {
		if(cancel_check != NULL && cancel_check(cancel_context)) return failure;
		double computed_time = izzo_time_of_flight(x, lambda);
		double delta = computed_time - dimensionless_time;
			if(!isfinite(delta)) return failure;
		if(fabs(delta) <= IZZO_TIME_OF_FLIGHT_TOLERANCE * fmax(1, dimensionless_time)) {
			converged = 1;
			break;
		}
		double first, second, third;
		if(!izzo_time_of_flight_derivatives(x, lambda, computed_time, &first, &second, &third))
			return failure;
		double first2 = first * first;
		double denominator = first * (first2 - delta * second) + third * delta * delta / 6;
		if(!isfinite(denominator) || fabs(denominator) < DBL_MIN)
			return failure;
		double next_x = x - delta * (first2 - delta * second / 2) / denominator;
		if(!isfinite(next_x) || next_x <= -1)
			return failure;
		x = next_x;
	}
	if(!converged) {
		Lambert3 not_converged = {.r0 = r0, .r1 = r1, .success = LAMBERT_MAX_ITERATIONS};
		return not_converged;
	}

	double y = sqrt(1 - lambda2 * (1 - x * x));
	double rho = (r0_magnitude - r1_magnitude) / chord;
	double sigma = sqrt(fmax(0, 1 - rho * rho));
	double gamma = sqrt(cb->mu * semiperimeter / 2);
	double radial_velocity0 = gamma * ((lambda * y - x) - rho * (lambda * y + x)) / r0_magnitude;
	double radial_velocity1 = -gamma * ((lambda * y - x) + rho * (lambda * y + x)) / r1_magnitude;
	double transverse_velocity = gamma * sigma * (y + lambda * x);
	Vector3 v0 = add_vec3(scale_vec3(radial0, radial_velocity0),
		scale_vec3(tangent0, transverse_velocity / r0_magnitude));
	Vector3 v1 = add_vec3(scale_vec3(radial1, radial_velocity1),
		scale_vec3(tangent1, transverse_velocity / r1_magnitude));
	if(!isfinite(v0.x) || !isfinite(v0.y) || !isfinite(v0.z) ||
		!isfinite(v1.x) || !isfinite(v1.y) || !isfinite(v1.z))
		return failure;
	if(cancel_check != NULL && cancel_check(cancel_context)) return failure;
	return (Lambert3){r0, v0, r1, v1, LAMBERT_SUCCESS};
}

Lambert3 calc_lambert3(Vector3 r0, Vector3 r1, double target_dt, Body *cb) {
	return calc_lambert3_cancellable(r0, r1, target_dt, cb, NULL, NULL);
}

// Elliptic time of flight with whole revolutions, x in (-1, 1).
static double izzo_time_of_flight_revs(double x, double lambda, int revolutions) {
	double one_minus_x2 = 1 - x * x;
	double y = sqrt(1 - lambda * lambda * one_minus_x2);
	double angle = acos(fmax(-1.0, fmin(1.0, x * y + lambda * one_minus_x2)));
	return ((angle + revolutions * M_PI) / sqrt(one_minus_x2) - x + lambda * y) / one_minus_x2;
}

static int izzo_revs_derivatives(double x, double lambda, double time_of_flight,
		double *first, double *second, double *third) {
	double one_minus_x2 = 1 - x * x;
	if(fabs(one_minus_x2) < 1e-5) return 0;
	double lambda2 = lambda * lambda;
	double lambda3 = lambda2 * lambda;
	double y = sqrt(1 - lambda2 * one_minus_x2);
	double y2 = y * y;
	double y3 = y2 * y;
	*first = (3 * time_of_flight * x - 2 + 2 * lambda3 * x / y) / one_minus_x2;
	*second = (3 * time_of_flight + 5 * x * *first + 2 * (1 - lambda2) * lambda3 / y3) / one_minus_x2;
	*third = (7 * x * *second + 8 * *first -
		6 * (1 - lambda2) * lambda2 * lambda3 * x / (y3 * y2)) / one_minus_x2;
	return isfinite(*first) && isfinite(*second) && isfinite(*third);
}

// The fastest transfer with this many revolutions sits at the minimum of the time of flight.
static int izzo_minimum_time(double lambda, int revolutions, double *x_minimum, double *time_minimum) {
	double x = 0;
	for(int iteration = 0; iteration < 60; iteration++) {
		double time = izzo_time_of_flight_revs(x, lambda, revolutions);
		double first, second, third;
		if(!izzo_revs_derivatives(x, lambda, time, &first, &second, &third)) return 0;
		double denominator = 2 * second * second - first * third;
		if(!isfinite(denominator) || fabs(denominator) < DBL_MIN) return 0;
		double next = fmax(-0.999, fmin(0.999, x - 2 * first * second / denominator));
		if(!isfinite(next)) return 0;
		if(fabs(next - x) < 1e-13) { x = next; break; }
		x = next;
	}
	*x_minimum = x;
	*time_minimum = izzo_time_of_flight_revs(x, lambda, revolutions);
	return isfinite(*time_minimum);
}

Lambert3 calc_lambert3_revs_cancellable(Vector3 r0, Vector3 r1, double target_dt, Body *cb,
		int revolutions, int right_branch, LambertCancelCheck cancel_check, void *cancel_context) {
	if(revolutions <= 0)
		return calc_lambert3_cancellable(r0, r1, target_dt, cb, cancel_check, cancel_context);
	Lambert3 failure = {.r0 = r0, .r1 = r1, .success = LAMBERT_FAIL_NAN};
	if(cancel_check != NULL && cancel_check(cancel_context)) return failure;
	if(cb == NULL || !(cb->mu > 0) || !(target_dt > 0) || !isfinite(target_dt)) return failure;
	double r0_magnitude = mag_vec3(r0);
	double r1_magnitude = mag_vec3(r1);
	if(!(r0_magnitude > 0) || !(r1_magnitude > 0)) return failure;

	Vector3 direction = cross_vec3(r0, r1);
	double direction_magnitude = mag_vec3(direction);
	if(!(direction_magnitude > 0) || !isfinite(direction_magnitude)) return failure;
	direction = scale_vec3(direction, 1 / direction_magnitude);

	Vector3 radial0 = scale_vec3(r0, 1 / r0_magnitude);
	Vector3 radial1 = scale_vec3(r1, 1 / r1_magnitude);
	Vector3 tangent0 = cross_vec3(direction, radial0);
	Vector3 tangent1 = cross_vec3(direction, radial1);
	double chord = mag_vec3(subtract_vec3(r1, r0));
	if(!(chord > 0)) return failure;
	double semiperimeter = (r0_magnitude + r1_magnitude + chord) / 2;
	double lambda_magnitude = sqrt(fmax(0, 1 - chord / semiperimeter));
	double lambda = direction.z < 0 ? -lambda_magnitude : lambda_magnitude;
	if(direction.z < 0) {
		tangent0 = scale_vec3(tangent0, -1);
		tangent1 = scale_vec3(tangent1, -1);
	}
	double lambda2 = lambda * lambda;
	double dimensionless_time = sqrt(2 * cb->mu / (semiperimeter * semiperimeter * semiperimeter)) * target_dt;
	if(!(dimensionless_time > 0) || !isfinite(dimensionless_time)) return failure;

	double x_minimum, time_minimum;
	if(!izzo_minimum_time(lambda, revolutions, &x_minimum, &time_minimum) ||
		dimensionless_time < time_minimum) return failure;
	double x;
	if(!right_branch) {
		double a = pow((revolutions * M_PI + M_PI) / (8 * dimensionless_time), 2.0 / 3.0);
		x = (a - 1) / (a + 1);
	} else {
		double b = pow(8 * dimensionless_time / (revolutions * M_PI), 2.0 / 3.0);
		x = (b - 1) / (b + 1);
	}
	x = fmax(-0.999, fmin(0.999, x));

	int converged = 0;
	for(int iteration = 0; iteration < 60; iteration++) {
		if(cancel_check != NULL && cancel_check(cancel_context)) return failure;
		double computed_time = izzo_time_of_flight_revs(x, lambda, revolutions);
		double delta = computed_time - dimensionless_time;
		if(!isfinite(delta)) return failure;
		if(fabs(delta) <= 1e-11 * fmax(1, dimensionless_time)) {
			converged = 1;
			break;
		}
		double first, second, third;
		if(!izzo_revs_derivatives(x, lambda, computed_time, &first, &second, &third)) return failure;
		double first2 = first * first;
		double denominator = first * (first2 - delta * second) + third * delta * delta / 6;
		if(!isfinite(denominator) || fabs(denominator) < DBL_MIN) return failure;
		double next_x = x - delta * (first2 - delta * second / 2) / denominator;
		if(!isfinite(next_x)) return failure;
		if(next_x <= -0.9999) next_x = (x - 1) / 2;
		if(next_x >= 0.9999) next_x = (x + 1) / 2;
		x = next_x;
	}
	if(!converged) return failure;
	if((!right_branch && x > x_minimum + 1e-9) || (right_branch && x < x_minimum - 1e-9)) return failure;

	double y = sqrt(1 - lambda2 * (1 - x * x));
	double rho = (r0_magnitude - r1_magnitude) / chord;
	double sigma = sqrt(fmax(0, 1 - rho * rho));
	double gamma = sqrt(cb->mu * semiperimeter / 2);
	double radial_velocity0 = gamma * ((lambda * y - x) - rho * (lambda * y + x)) / r0_magnitude;
	double radial_velocity1 = -gamma * ((lambda * y - x) + rho * (lambda * y + x)) / r1_magnitude;
	double transverse_velocity = gamma * sigma * (y + lambda * x);
	Vector3 v0 = add_vec3(scale_vec3(radial0, radial_velocity0),
		scale_vec3(tangent0, transverse_velocity / r0_magnitude));
	Vector3 v1 = add_vec3(scale_vec3(radial1, radial_velocity1),
		scale_vec3(tangent1, transverse_velocity / r1_magnitude));
	if(!isfinite(v0.x) || !isfinite(v0.y) || !isfinite(v0.z) ||
		!isfinite(v1.x) || !isfinite(v1.y) || !isfinite(v1.z)) return failure;
	return (Lambert3){r0, v0, r1, v1, LAMBERT_SUCCESS};
}

Lambert3 calc_lambert3_revs(Vector3 r0, Vector3 r1, double target_dt, Body *cb,
		int revolutions, int right_branch) {
	return calc_lambert3_revs_cancellable(r0, r1, target_dt, cb, revolutions, right_branch, NULL, NULL);
}

double get_flyby_periapsis(Vector3 v_arr, Vector3 v_dep, Vector3 v_body, Body *body) {
	Vector3 v1 = subtract_vec3(v_arr, v_body);
	Vector3 v2 = subtract_vec3(v_dep, v_body);
	double beta = (M_PI - angle_vec3_vec3(v1, v2))/2;
	return (1 / cos(beta) - 1) * (body->mu / (pow(mag_vec3(v1), 2)));
}

// only_side is -1 or +1 to fix the side of the plane of the two excess velocities, 0 for the cheaper one.
static int powered_flyby_at_radius_side(Vector3 incoming, Vector3 outgoing, double incoming_speed,
		double outgoing_speed, double mu, double periapsis_radius, int only_side, PoweredFlybyResult *result) {
	Vector3 incoming_direction = scale_vec3(incoming, 1 / incoming_speed);
	Vector3 outgoing_direction = scale_vec3(outgoing, 1 / outgoing_speed);
	double direction_dot = fmax(-1, fmin(1, dot_vec3(incoming_direction, outgoing_direction)));
	double denominator = 1 - direction_dot * direction_dot;
	if(denominator < 1e-12) return 0;
	double incoming_sine = 1 / (1 + periapsis_radius * incoming_speed * incoming_speed / mu);
	double outgoing_sine = 1 / (1 + periapsis_radius * outgoing_speed * outgoing_speed / mu);
	double incoming_cosine = sqrt(fmax(0, 1 - incoming_sine * incoming_sine));
	double outgoing_cosine = sqrt(fmax(0, 1 - outgoing_sine * outgoing_sine));
	double incoming_scale = (incoming_sine + direction_dot * outgoing_sine) / denominator;
	double outgoing_scale = (-outgoing_sine - direction_dot * incoming_sine) / denominator;
	Vector3 periapsis_planar = add_vec3(scale_vec3(incoming_direction, incoming_scale),
		scale_vec3(outgoing_direction, outgoing_scale));
	double planar_norm_squared = dot_vec3(periapsis_planar, periapsis_planar);
	if(planar_norm_squared > 1 + 1e-10) return 0;
	Vector3 normal = cross_vec3(incoming_direction, outgoing_direction);
	double normal_magnitude = mag_vec3(normal);
	if(!(normal_magnitude > 0)) return 0;
	normal = scale_vec3(normal, 1 / normal_magnitude);
	double normal_component = sqrt(fmax(0, 1 - planar_norm_squared));
	double incoming_periapsis_speed = sqrt(incoming_speed * incoming_speed + 2 * mu / periapsis_radius);
	double outgoing_periapsis_speed = sqrt(outgoing_speed * outgoing_speed + 2 * mu / periapsis_radius);
	PoweredFlybyResult best = {.periapsis_radius = periapsis_radius, .burn_delta_v = DBL_MAX};
	for(int side = -1; side <= 1; side += 2) {
		if(only_side != 0 && side != only_side) continue;
		Vector3 periapsis_direction = add_vec3(periapsis_planar,
			scale_vec3(normal, side * normal_component));
		Vector3 incoming_tangent = norm_vec3(subtract_vec3(incoming_direction,
			scale_vec3(periapsis_direction, incoming_sine)));
		Vector3 outgoing_tangent = norm_vec3(add_vec3(outgoing_direction,
			scale_vec3(periapsis_direction, outgoing_sine)));
		Vector3 incoming_velocity = scale_vec3(incoming_tangent, incoming_periapsis_speed);
		Vector3 outgoing_velocity = scale_vec3(outgoing_tangent, outgoing_periapsis_speed);
		Vector3 burn = subtract_vec3(outgoing_velocity, incoming_velocity);
		double burn_magnitude = mag_vec3(burn);
		if(burn_magnitude < best.burn_delta_v) {
			best.periapsis_direction = periapsis_direction;
			best.incoming_tangent_direction = incoming_tangent;
			best.outgoing_tangent_direction = outgoing_tangent;
			best.burn_delta_v_vector = burn;
			best.burn_delta_v = burn_magnitude;
		}
	}
	if(!isfinite(best.burn_delta_v)) return 0;
	*result = best;
	return 1;
}

static int powered_flyby_at_radius(Vector3 incoming, Vector3 outgoing, double incoming_speed,
		double outgoing_speed, double mu, double periapsis_radius, PoweredFlybyResult *result) {
	return powered_flyby_at_radius_side(incoming, outgoing, incoming_speed, outgoing_speed, mu,
		periapsis_radius, 0, result);
}

static int powered_flyby_soi_state(Body *body, Vector3 periapsis_direction,
		Vector3 tangent_direction, double periapsis_radius, double excess_speed,
		double soi_radius, int incoming, double *travel_time, OSV *state) {
	double eccentricity = 1 + periapsis_radius * excess_speed * excess_speed / body->mu;
	double semi_latus_rectum = periapsis_radius * (1 + eccentricity);
	double cosine_anomaly = (semi_latus_rectum / soi_radius - 1) / eccentricity;
	if(!(eccentricity > 1) || cosine_anomaly < -1 || cosine_anomaly > 1) return 0;
	double anomaly_magnitude = acos(cosine_anomaly);
	double anomaly = incoming ? -anomaly_magnitude : anomaly_magnitude;
	double cosine = cos(anomaly_magnitude);
	double hyperbolic_anomaly = acosh((eccentricity + cosine) / (1 + eccentricity * cosine));
	double semi_major_axis_magnitude = body->mu / (excess_speed * excess_speed);
	double mean_motion = sqrt(body->mu / pow(semi_major_axis_magnitude, 3));
	*travel_time = (eccentricity * sinh(hyperbolic_anomaly) - hyperbolic_anomaly) / mean_motion;
	double velocity_scale = sqrt(body->mu / semi_latus_rectum);
	state->r = add_vec3(scale_vec3(periapsis_direction, soi_radius * cos(anomaly)),
		scale_vec3(tangent_direction, soi_radius * sin(anomaly)));
	state->v = add_vec3(scale_vec3(periapsis_direction, -velocity_scale * sin(anomaly)),
		scale_vec3(tangent_direction, velocity_scale * (eccentricity + cosine)));
	return isfinite(*travel_time) && *travel_time > 0 &&
		is_finite_vector(state->r) && is_finite_vector(state->v);
}

static int populate_powered_flyby_soi(Body *body, double incoming_speed, double outgoing_speed,
		PoweredFlybyResult *result) {
	Body *central_body = body->orbit.cb;
	if(central_body == NULL || !(central_body->mu > 0) || !(body->orbit.a != 0)) return 0;
	double soi_radius = fabs(body->orbit.a) * pow(body->mu / central_body->mu, 0.4);
	if(!isfinite(soi_radius) || soi_radius <= result->periapsis_radius) return 0;
	if(!powered_flyby_soi_state(body, result->periapsis_direction,
		result->incoming_tangent_direction, result->periapsis_radius, incoming_speed,
		soi_radius, 1, &result->incoming_soi_time, &result->incoming_soi_state) ||
		!powered_flyby_soi_state(body, result->periapsis_direction,
			result->outgoing_tangent_direction, result->periapsis_radius, outgoing_speed,
			soi_radius, 0, &result->outgoing_soi_time, &result->outgoing_soi_state)) return 0;
	result->soi_radius = soi_radius;
	return 1;
}

int calculate_powered_flyby(Vector3 v_arr, Vector3 v_dep, Vector3 v_body, Body *body,
		double minimum_periapsis_radius, PoweredFlybyResult *result) {
	if(body == NULL || result == NULL || !(body->mu > 0) || !(minimum_periapsis_radius > 0)) return 0;
	Vector3 incoming = subtract_vec3(v_arr, v_body);
	Vector3 outgoing = subtract_vec3(v_dep, v_body);
	double incoming_speed = mag_vec3(incoming);
	double outgoing_speed = mag_vec3(outgoing);
	if(!(incoming_speed > 0) || !(outgoing_speed > 0)) return 0;
	PoweredFlybyResult shortcut;
	int shortcut_valid = 0;
	if(fabs(incoming_speed - outgoing_speed) < 1) {
		// Single-periapsis shortcut; it can fail or overstate the burn when the speeds differ slightly,
		// so it is only returned directly when the burn is negligible.
		double unpowered_periapsis = get_flyby_periapsis(v_arr, v_dep, v_body, body);
		if(unpowered_periapsis < minimum_periapsis_radius) return 0;
		shortcut_valid = powered_flyby_at_radius(incoming, outgoing, incoming_speed, outgoing_speed,
			body->mu, unpowered_periapsis, &shortcut) &&
			populate_powered_flyby_soi(body, incoming_speed, outgoing_speed, &shortcut);
		if(shortcut_valid && shortcut.burn_delta_v < 2.0) {
			*result = shortcut;
			return 1;
		}
	}

	double maximum_periapsis = minimum_periapsis_radius * 1e5;
	if(body->orbit.cb != NULL && body->orbit.cb->mu > 0 && fabs(body->orbit.a) > minimum_periapsis_radius) {
		double sphere_of_influence = fabs(body->orbit.a) * pow(body->mu / body->orbit.cb->mu, 0.4);
		if(isfinite(sphere_of_influence) && sphere_of_influence > minimum_periapsis_radius)
			maximum_periapsis = fmin(maximum_periapsis, 0.95 * sphere_of_influence);
	}
	if(!(maximum_periapsis > minimum_periapsis_radius)) return 0;

	int found_solution = 0;
	const int sample_count = 64;
	double log_minimum = log(minimum_periapsis_radius);
	double log_range = log(maximum_periapsis / minimum_periapsis_radius);
	double lower = 0, upper = 0;
	PoweredFlybyResult best = {.burn_delta_v = DBL_MAX};
	for(int sample = 0; sample < sample_count; sample++) {
		double periapsis_radius = exp(log_minimum + log_range * sample / (sample_count - 1));
		PoweredFlybyResult candidate;
		if(powered_flyby_at_radius(incoming, outgoing, incoming_speed, outgoing_speed,
			body->mu, periapsis_radius, &candidate) && candidate.burn_delta_v < best.burn_delta_v) {
			best = candidate;
			lower = log_minimum + log_range * fmax(0, sample - 1) / (sample_count - 1);
			upper = log_minimum + log_range * fmin(sample_count - 1, sample + 1) / (sample_count - 1);
			found_solution = 1;
		}
	}
	// The coarse grid can miss a narrow feasible window near the turn-matching periapsis.
	double local_lower = 0, local_upper = 0, local_best_burn = DBL_MAX;
	int local_found = 0;
	double turn_periapsis = get_flyby_periapsis(v_arr, v_dep, v_body, body);
	if(isfinite(turn_periapsis) && turn_periapsis > 0 && best.burn_delta_v > 5.0) {
		const int local_count = 49;
		double local_span = log(2.0);
		for(int sample = 0; sample < local_count; sample++) {
			double log_radius = log(turn_periapsis) + local_span * (2.0 * sample / (local_count - 1) - 1);
			double periapsis_radius = exp(log_radius);
			if(periapsis_radius < minimum_periapsis_radius || periapsis_radius > maximum_periapsis) continue;
			PoweredFlybyResult candidate;
			if(powered_flyby_at_radius(incoming, outgoing, incoming_speed, outgoing_speed,
				body->mu, periapsis_radius, &candidate) && candidate.burn_delta_v < local_best_burn) {
				local_best_burn = candidate.burn_delta_v;
				double step = 2 * local_span / (local_count - 1);
				local_lower = log_radius - step;
				local_upper = log_radius + step;
				if(candidate.burn_delta_v < best.burn_delta_v) best = candidate;
				found_solution = 1;
				local_found = 1;
			}
		}
	}
	if(!found_solution) {
		if(!shortcut_valid) return 0;
		*result = shortcut;
		return 1;
	}

	for(int bracket = 0; bracket < 2; bracket++) {
		if(bracket == 1) {
			if(!local_found) break;
			lower = fmax(local_lower, log_minimum);
			upper = fmin(local_upper, log_minimum + log_range);
		}
		for(int iteration = 0; iteration < 32; iteration++) {
			double first_log = lower + (upper - lower) * 0.3819660112501051;
			double second_log = lower + (upper - lower) * 0.6180339887498949;
			PoweredFlybyResult first, second;
			int first_valid = powered_flyby_at_radius(incoming, outgoing, incoming_speed, outgoing_speed,
				body->mu, exp(first_log), &first);
			int second_valid = powered_flyby_at_radius(incoming, outgoing, incoming_speed, outgoing_speed,
				body->mu, exp(second_log), &second);
			double first_cost = first_valid ? first.burn_delta_v : DBL_MAX;
			double second_cost = second_valid ? second.burn_delta_v : DBL_MAX;
			if(first_cost <= second_cost) upper = second_log;
			else lower = first_log;
		}
		PoweredFlybyResult refined;
		if(powered_flyby_at_radius(incoming, outgoing, incoming_speed, outgoing_speed,
			body->mu, exp((lower + upper) / 2), &refined) && refined.burn_delta_v < best.burn_delta_v)
			best = refined;
	}
	int best_valid = populate_powered_flyby_soi(body, incoming_speed, outgoing_speed, &best);
	if(shortcut_valid && (!best_valid || shortcut.burn_delta_v <= best.burn_delta_v)) {
		*result = shortcut;
		return 1;
	}
	if(!best_valid) return 0;
	*result = best;
	return 1;
}

int calculate_powered_flyby_at_periapsis(Vector3 v_arr, Vector3 v_dep, Vector3 v_body, Body *body,
		double periapsis_radius, int side, PoweredFlybyResult *result) {
	if(body == NULL || result == NULL || !(body->mu > 0) || !(periapsis_radius > 0) || (side != -1 && side != 1)) return 0;
	Vector3 incoming = subtract_vec3(v_arr, v_body);
	Vector3 outgoing = subtract_vec3(v_dep, v_body);
	double incoming_speed = mag_vec3(incoming);
	double outgoing_speed = mag_vec3(outgoing);
	if(!(incoming_speed > 0) || !(outgoing_speed > 0)) return 0;
	return powered_flyby_at_radius_side(incoming, outgoing, incoming_speed, outgoing_speed,
		body->mu, periapsis_radius, side, result) &&
		populate_powered_flyby_soi(body, incoming_speed, outgoing_speed, result);
}

double get_flyby_inclination(Vector3 v_arr, Vector3 v_dep, Vector3 v_body, Plane3 equatorial_plane) {
	Vector3 v1 = subtract_vec3(v_arr, v_body);
	Vector3 v2 = subtract_vec3(v_dep, v_body);
	Plane3 hyperbola_plane = constr_plane3(vec3(0, 0, 0), v1, v2);
	return angle_plane3_plane3(equatorial_plane, hyperbola_plane);
}

HyperbolaLegParams get_dep_hyperbola_params(Vector3 v_inf) {
	Plane3 xy = constr_plane3(vec3(0,0,0), vec3(1,0,0), vec3(0,1,0));
	Plane3 xz = constr_plane3(vec3(0,0,0), vec3(1,0,0), vec3(0,0,1));
	double decl = angle_plane3_vec3(xy, v_inf);
	decl = v_inf.z < 0 ? -fabs(decl) : fabs(decl);
	double bplane_angle = -angle_plane3_vec3(xz, v_inf);
	if(cross_vec3(v_inf,vec3(0,1,0)).z < 0) bplane_angle = M_PI-bplane_angle;
	bplane_angle = pi_norm(bplane_angle);
	// bvazi only relevant with given inclination (not yet implemented here, but implemented for fly-by)
	HyperbolaLegParams hyp_params = {decl, bplane_angle, M_PI/2};
	return hyp_params;
}

HyperbolaParams get_hyperbola_params(Vector3 v_arr, Vector3 v_dep, Vector3 v_body, struct Body *body, double h_pe, enum HyperbolaType type) {
	Vector3 vinf_arr = subtract_vec3(v_arr, v_body);
	Vector3 vinf_dep = subtract_vec3(v_dep, v_body);
	HyperbolaParams hyperbola_params = {.type = type};
	
	if(type != HYP_DEPARTURE) {
		hyperbola_params.incoming = get_dep_hyperbola_params(vinf_arr);
		// invert hyperbola direction (departure to arrival hyperbola)
		hyperbola_params.incoming.decl *= -1;
		hyperbola_params.incoming.bplane_angle = pi_norm(M_PI + hyperbola_params.incoming.bplane_angle);
	}
	if(type != HYP_ARRIVAL) hyperbola_params.outgoing = get_dep_hyperbola_params(vinf_dep);
	if(type == HYP_FLYBY) {
		Vector3 N = cross_vec3(vinf_arr, vinf_dep);
		Vector3 B_arr = cross_vec3(vinf_arr, N);
		Vector3 B_dep = cross_vec3(vinf_dep, scale_vec3(N,-1));
		
		// BVAZI (Azimuth of B-vector) calculated from south
		hyperbola_params.incoming.bvazi = angle_vec3_vec3(vec3(0,0,-1),B_arr);
		hyperbola_params.outgoing.bvazi = angle_vec3_vec3(vec3(0,0,-1),B_dep);
		
		// if retrograde orbit
		if(N.z < 0) {
			hyperbola_params.incoming.bvazi *= -1;
			hyperbola_params.outgoing.bvazi *= -1;
		}
	}
	
	double rp = type == HYP_FLYBY ? get_flyby_periapsis(v_arr, v_dep, v_body, body) : h_pe+body->radius;
	double c3_energy = sq_mag_vec3(type != HYP_DEPARTURE ? vinf_arr : vinf_dep);
	
	hyperbola_params.rp = rp;
	hyperbola_params.c3_energy = c3_energy;
	
	return hyperbola_params;
}

bool is_flyby_viable(Vector3 v_arr, Vector3 v_dep, Vector3 v_body, Body *body, double precision) {
	double rp = get_flyby_periapsis(v_arr, v_dep, v_body, body);
	if(rp < body->radius+body->atmo_alt) return false;
	if(fabs(mag_vec3(subtract_vec3(v_dep, v_body)) - mag_vec3(subtract_vec3(v_arr, v_body))) > precision) return false;
	return true;
}

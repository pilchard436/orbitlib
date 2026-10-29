#include "orbitlib_ephemeris.h"
#include "orbitlib_fileio.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

static Ephem *closest_ephem_entry(Ephem *ephem, int num_ephems, double epoch);


void print_ephem(struct Ephem ephem) {
	printf("Date: %f  (", ephem.epoch);
	print_date(convert_JD_date(ephem.epoch, DATE_ISO), 0);
	printf(")\nx: %g m,   y: %g m,   z: %g m\n"
		   "vx: %g m/s,   vy: %g m/s,   vz: %g m/s\n\n",
		   ephem.r.x, ephem.r.y, ephem.r.z, ephem.v.x, ephem.v.y, ephem.v.z);
}

void get_ephem_data_filepath(int id, char *filepath, const char *ephem_directory) {
	sprintf(filepath, "%s/%d.ephem", ephem_directory, id);
}

int is_ephem_available(int body_code, const char *ephem_directory) {
	char filepath[50];
	get_ephem_data_filepath(body_code, filepath, ephem_directory);
	FILE *file = fopen(filepath, "r");  // Try to open file in read mode
	if (file) {
		fclose(file);  // Close file if it was opened
		return 1;      // File exists
	}
	return 0;          // File does not exist
}

void get_body_ephems(Body *body, Datetime min_date, Datetime max_date, Datetime time_step, const char *ephem_directory) {
	if(body->orbit.cb == NULL) return;
	if(!directory_exists(ephem_directory)) create_directory(ephem_directory);
	
	char filepath[50];
	get_ephem_data_filepath(body->id, filepath, ephem_directory);
	
	if(!is_ephem_available(body->id, ephem_directory)) {
		char d0_s[32];
		char d1_s[32];
		date_to_string(min_date, d0_s, 1);
		date_to_string(max_date, d1_s, 1);
		// Construct the URL with your API key and parameters
		
		char timestep_s[8];
		if(time_step.y > 0) {
			sprintf(timestep_s, "%d y", time_step.y);
		} else if(time_step.m > 0) {
			sprintf(timestep_s, "%d mo", time_step.m);
		} else if(time_step.d > 0) {
			sprintf(timestep_s, "%d d", time_step.d);
		} else return;
		
		char body_id[24];
		if(body->id >= 20000000) sprintf(body_id, "DES=%d", body->id);
		else sprintf(body_id, "%d", body->id);
		
		char url[256];
		sprintf(url, "https://ssd.jpl.nasa.gov/api/horizons.api?"
					 "format=text&"
					 "COMMAND='%s'&"
					 "OBJ_DATA='YES'&"
					 "MAKE_EPHEM='YES'&"
					 "EPHEM_TYPE='VECTORS'&"
					 "CENTER='500@%d'&"
					 "START_TIME='%s'&"
					 "STOP_TIME='%s'&"
					 "STEP_SIZE='%s'&"
					 "VEC_TABLE='2'", body_id, body->orbit.cb->id, d0_s, d1_s, timestep_s);
		
		download_file(url, filepath);
	}
	
	FILE *file;
	char line[256];  // Assuming lines are no longer than 255 characters
	
	file = fopen(filepath, "r");
	
	if(file == NULL) {
		perror("Unable to open file");
		return;
	}
	
	// Read lines from the file until the end is reached
	while(fgets(line, sizeof(line), file) != NULL) {
		line[strcspn(line, "\n")] = '\0';
		if(strcmp(line, "$$SOE") == 0) {
			break; // Exit the loop when "$$SOE" is encountered
		}
	}
	
	int max_num_ephems = 12;
	
	if(body->ephem != NULL) free(body->ephem);
	body->ephem = calloc(max_num_ephems, sizeof(Ephem));
	body->num_ephems = 0;
	
	fgets(line, sizeof(line), file);
	line[strcspn(line, "\n")] = '\0';
	
	while(strcmp(line, "$$EOE") != 0) {
		char *endptr;
		double date = strtod(line, &endptr);
		
		fgets(line, sizeof(line), file);
		double x, y, z;
		sscanf(line, " X =%lf Y =%lf Z =%lf", &x, &y, &z);
		fgets(line, sizeof(line), file);
		double vx, vy, vz;
		sscanf(line, " VX=%lf VY=%lf VZ=%lf", &vx, &vy, &vz);
		
		if(body->num_ephems == max_num_ephems) {
			max_num_ephems *= 2;
			Ephem *temp = realloc(body->ephem, max_num_ephems*sizeof(Ephem));
			if(temp != NULL) body->ephem = temp;
		}
		
		body->ephem[body->num_ephems].epoch = date;
		body->ephem[body->num_ephems].r.x = x*1e3;
		body->ephem[body->num_ephems].r.y = y*1e3;
		body->ephem[body->num_ephems].r.z = z*1e3;
		body->ephem[body->num_ephems].v.x = vx*1e3;
		body->ephem[body->num_ephems].v.y = vy*1e3;
		body->ephem[body->num_ephems].v.z = vz*1e3;
		body->num_ephems++;
		
		fgets(line, sizeof(line), file);
		line[strcspn(line, "\n")] = '\0';
	}
	fclose(file);
}

Ephem get_closest_ephem(Ephem *ephem, int num_ephems, double epoch) {
	return *closest_ephem_entry(ephem, num_ephems, epoch);
}

// The entry whose epoch is nearest; the list is sorted ascending by epoch.
static Ephem *closest_ephem_entry(Ephem *ephem, int num_ephems, double epoch) {
	int low = 0, high = num_ephems;	// first index with ephem[index].epoch > epoch
	while(low < high) {
		int middle = low + (high - low) / 2;
		if(epoch < ephem[middle].epoch) high = middle;
		else low = middle + 1;
	}
	if(low == num_ephems) return &ephem[num_ephems - 1];
	if(low == 0) return &ephem[0];
	return fabs(ephem[low - 1].epoch - epoch) < fabs(ephem[low].epoch - epoch) ? &ephem[low - 1] : &ephem[low];
}

// Elements of recently used ephemeris entries, so repeated propagation from the same entry
// skips the state-to-elements conversion. Keyed by entry address and contents; per thread.
#define EPHEM_ELEMENT_CACHE_SIZE 16
typedef struct EphemElementCache {
	const Ephem *entry;
	Ephem contents;
	struct Body *cb;
	Orbit orbit;
} EphemElementCache;

static _Thread_local EphemElementCache ephem_element_cache[EPHEM_ELEMENT_CACHE_SIZE];
static _Thread_local int ephem_element_cache_next;

static int same_ephem(const Ephem *a, const Ephem *b) {
	return a->epoch == b->epoch && a->r.x == b->r.x && a->r.y == b->r.y && a->r.z == b->r.z &&
		a->v.x == b->v.x && a->v.y == b->v.y && a->v.z == b->v.z;
}

OSV osv_from_ephem(Ephem *ephem_list, int num_ephems, double epoch, struct Body *cb) {
	const Ephem *entry = closest_ephem_entry(ephem_list, num_ephems, epoch);
	double dt = (epoch - entry->epoch) * (24 * 60 * 60);
	for(int k = 0; k < EPHEM_ELEMENT_CACHE_SIZE; k++) {
		if(ephem_element_cache[k].entry == entry && ephem_element_cache[k].cb == cb &&
			same_ephem(&ephem_element_cache[k].contents, entry))
			return osv_from_elements_after(&ephem_element_cache[k].orbit, dt);
	}
	EphemElementCache *cache = &ephem_element_cache[ephem_element_cache_next];
	ephem_element_cache_next = (ephem_element_cache_next + 1) % EPHEM_ELEMENT_CACHE_SIZE;
	cache->entry = entry;
	cache->contents = *entry;
	cache->cb = cb;
	cache->orbit = constr_orbit_from_osv(entry->r, entry->v, cb);
	return osv_from_elements_after(&cache->orbit, dt);
}
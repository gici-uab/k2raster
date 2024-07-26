#define __STDC_LIMIT_MACROS 1

//#include <linux/limits.h>
#include <libgen.h>
//#include <stdio.h>
#include <string.h>
//#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <math.h>
#include <arpa/inet.h>
#include <sys/time.h>
#include <bits/stdc++.h>
#define USE_LIMITED_LEVEL 0
#if USE_LIMITED_LEVEL
	#include "./limited_levels/dacs.h"
#else
	#include "./no_restrictions/dacs.h"
#endif
//#include "./no_restrictions/basics.h"
#include "utils.h"
#include "utils/HashTable.h"
#include <vector>
#include <cstdlib>
#include <cstdio>
#include <climits>

#define compile_with_DAC	(1)
#define NL					(putchar('\n'))
#define mask31 				0x0000001F
#if (!compile_with_DAC)
	#define W      			32
#endif
#define TREE_DEPTH			((int)((double)log10(nx)/(double)log10(k_val)) + 1)

#if !defined(CHAR_BIT)
#define CHAR_BIT	    	8
#endif

struct voc {
	union {
		int32_t   *submatrix;
		int32_t   submatrix_value;
	};
	uint32_t  weight;
};

typedef struct _maxmin {
	union {
		int16_t  max_s16;
		uint16_t max_u16;
	};
	union {
		int16_t  min_s16;
		uint16_t min_u16;
	};
} maxmin;


typedef struct {
	union {
		int8_t		*s8;
		uint8_t		*u8;
		int16_t		*s16;
		uint16_t	*u16;
		int32_t		*s32;
		uint32_t	*u32;
	};
} vector_type;

uint16_t	k_val;

/////////////////////////
uint 		*pmax_2d;
uint 		*pmin_2d;
int32_t 	**vmax_2d;
int32_t 	**vmin_2d;
int32_t 	*lmax_2d;
int32_t 	*lmin_2d;
uint8_t 	**t3_2d;
uint 		**T_2d;
maxmin 		*rMaxMin_2d;
/////////////////////////

/*********************************************************/
int64_t	    get_file_size(const char *fname);
int32_t		file_exists(const char* fname);
int32_t  	combine_bits_into_int(uint *input_vector, int input_size, uint *output_vector);
int32_t		all_zeros(int32_t *v, uint32_t v_len);
int32_t  	get_bit(uint *T, int size, int bit_pos); // bit_pos 0...t_size-1
uint32_t 	get_rank(const uint *list, uint pos);
int32_t  	get_base_name(char *input_path_name, char *output_base_name);
void 		print_size(const char *prog_dir, const char *file_name, char heuristic, uint32_t nz);
void        print_vector(int16_t *v, int b1, int b2, int r1, int r2, int c1, int c2);
double      compute_entropy(int16_t *vect, uint32_t vect_len, uint32_t band);
void        remove_data(const char * dir, const char *file_name, uint16_t k_val, uint32_t nz);
/*********************************************************/

int32_t get_cell(uint *T, FTRep *ftrep, int32_t T_size, uint32_t n, uint32_t r, uint32_t c, int32_t z, int32_t max_val)
{
	z = get_rank(T, z) * k_val * k_val;
	
	z = z + r / (n/k_val) * k_val + c / (n/k_val);
	
	int32_t val = accessFT(ftrep, z);
	
	max_val -= val;
	
	if (z >= T_size || bitget(T, z) == 0)
		return max_val;
	else 
		return get_cell(T, ftrep, T_size, n/k_val, r%(n/k_val), c%(n/k_val), z, max_val);
}

/*
void get_window(uint32_t n,
				uint32_t b,
                uint32_t r1,
				uint32_t r2,
				uint32_t c1,
				uint32_t c2,
				int32_t  z,
				int16_t  max_val)
{
	uint32_t  i = 0;
	uint32_t  j = 0;

	uint32_t r1_new      = 0;
	uint32_t r2_new      = 0;
	uint32_t c1_new      = 0;
	uint32_t c2_new      = 0;
	int32_t  z_new       = 0;
	int16_t  max_val_new = 0;

	z = get_rank(T[b], z) * k_val * k_val;

	for (j = r1/(n/k_val); j <= r2/(n/k_val); j++) {
		if (j == r1/(n/k_val))
			r1_new = r1 % (n/k_val);
		else
			r1_new = 0;
		if (j == r2/(n/k_val))
			r2_new = r2 % (n/k_val);
		else
			r2_new = (n/k_val) - 1;
		for (i = c1/(n/k_val); i <= c2/(n/k_val); i++) {
			if (i == c1/(n/k_val))
				c1_new = c1 % (n/k_val);
			else
				c1_new = 0;
			if (i == c2/(n/k_val))
				c2_new = c2 % (n/k_val);
			else
				c2_new = (n/k_val) - 1;
			z_new = z + k_val*j + i;
			max_val_new = max_val - accessFT(dacMax[b], z_new+1);
			int bit = get_bit(T[b], T_size, z_new);
			if (z_new >= T_size || bit == 0) {
				// printf("%d %d times\n",
					   // max_val_new,
					   // ((r2_new-r1_new)+1)*((c2_new-c1_new)+1));
			}
			else {
				// puts("*******");
				// puts("");
				// printf("j: %d\n", j);
				// printf("i: %d\n", i);
				// printf("n/k_val: %d\n", n/k_val);
				// printf("r1_new: %d\n", r1_new);
				// printf("r2_new: %d\n", r2_new);
				// printf("c1_new: %d\n", c1_new);
				// printf("c2_new: %d\n", c2_new);
				// printf("z_new: %d\n", z_new);
				// printf("max_val_new: %d\n", max_val_new);
				// puts("*******");
				get_window(n/k_val, b, r1_new, r2_new, c1_new, c2_new, z_new, max_val_new);
			}
		}
	}
}

void search_values(uint32_t n,
                   uint32_t b,
				   uint32_t r1,
				   uint32_t r2,
				   uint32_t c1,
				   uint32_t c2,
				   int16_t  vb,
				   int16_t  ve,
				   int16_t  max_val,
				   int16_t  min_val,
				   int32_t  z)
{
	uint32_t i          = 0;
	uint32_t j          = 0;
	uint32_t r1_new      = 0;
	uint32_t r2_new      = 0;
	uint32_t c1_new      = 0;
	uint32_t c2_new      = 0;
	int32_t  z_new       = 0;
	int16_t  max_val_new = 0;
	int16_t  min_val_new = 0;

	z = get_rank(T[b], z) * k_val * k_val;
	for (j = r1/(n/k_val); j <= r2/(n/k_val); j++) {
		if (j == r1/(n/k_val))
			r1_new = r1 % (n/k_val);
		else
			r1_new = 0;
		if (j == r2/(n/k_val))
			r2_new = r2 % (n/k_val);
		else
			r2_new = (n/k_val) - 1;

		for (i = c1/(n/k_val); i <= c2/(n/k_val); i++) {
			if (i == c1/(n/k_val))
				c1_new = c1 % (n/k_val);
			else
				c1_new = 0;
			if (i == c2/(n/k_val))
				c2_new = c2 % (n/k_val);
			else
				c2_new = (n/k_val) - 1;
			z_new = z + k_val*j + i;
			max_val_new = max_val - accessFT(dacMax[b], z_new+1);

			int bit = get_bit(T[b], T_size, z_new);
			if (z_new >= T_size || bit == 0) { // leaf
				min_val_new = max_val_new;
				if (min_val_new >= vb && max_val_new <= ve) {
					printf("LEAF: %d %d %d %d max_val_new=%d row=%d column=%d.\n",
					    r1, r2, c1, c2, max_val_new, j, i);
				}
				printf("%d %d times\n",
					   max_val_new,
					   ((r2_new-r1_new)+1)*((c2_new-c1_new)+1));
			}
			else { // internal node
				min_val_new = min_val + accessFT(dacMin[b], (get_rank(T[b], z))+1);
				if (min_val_new >= vb && max_val_new <= ve) {
					//printf("INTERNAL: min_val_new=%d row=%d column=%d.\n", min_val_new, j, i);
				}
				if (min_val_new > ve || max_val_new < vb) {
					//printf("INTERNAL: no cells meet the condition in this branch.\n");
				}
				if (min_val_new < ve || max_val_new > vb) {
					search_values(n/k_val, b, r1_new, r2_new, c1_new, c2_new, vb, ve, max_val_new, min_val_new, z_new);
				}
			}
		}
	}
}
*/

//#define img_vector[x][y] (*(img_vector+(y*512+x)))

maxmin build_u16(uint16_t *vec, uint32_t nx_new, uint32_t n, uint32_t l, uint32_t r, uint32_t c)
{
	maxmin ret;
	 // std::numeric_limits<std::int16_t>::max(); // INT16_MAX;
	uint16_t	maxval = 0;
	uint16_t	minval = UINT16_MAX;

	 // std::numeric_limits<std::int16_t>::min(); //INT16_MIN;
	uint32_t	i  = 0;
	uint32_t	j  = 0;
	int32_t	x      = 0;
	int32_t y      = 0;

	for (j = 0; j < k_val; ++j) {
		for (i = 0; i < k_val; ++i) {
			if (k_val == n) { // last level
				y = r+j;
				x = c+i;
				if (minval > *(vec + (y * nx_new + x)))
					minval = *(vec + (y * nx_new + x));
				if (maxval < *(vec + (y * nx_new + x)))
					maxval = *(vec + (y * nx_new + x));
				vmax_2d[l][pmax_2d[l]] = *(vec + (y * nx_new + x));
				pmax_2d[l]++;
			}
			else { // internal node
				maxmin child = build_u16(vec, nx_new, n/k_val, l+1, r+j*(n/k_val), c+i*(n/k_val));
				vmax_2d[l][pmax_2d[l]] = child.max_u16;
				vmin_2d[l][pmin_2d[l]] = child.min_u16;
				pmin_2d[l]++;

				if (child.max_u16 != child.min_u16)
					t3_2d[l][pmax_2d[l]] = 1;

				pmax_2d[l]++;

				if (minval > child.min_u16)
					minval = child.min_u16;
				if (maxval < child.max_u16)
					maxval = child.max_u16;
			}
		}
	}

	if (minval == maxval) {
		pmax_2d[l] -= k_val*k_val;
		pmin_2d[l] -= k_val*k_val;
	}

	ret.max_u16 = maxval;
	ret.min_u16 = minval;

	return ret;
}

maxmin build_s16(int16_t *vec, uint32_t nx_new, uint32_t n, uint32_t l, uint32_t r, uint32_t c)
{
	maxmin ret;
	 // std::numeric_limits<std::int16_t>::max(); // INT16_MAX;
	int16_t	maxval = INT16_MIN;
	int16_t	minval = INT16_MAX;

	 // std::numeric_limits<std::int16_t>::min(); //INT16_MIN;
	uint32_t	i = 0;
	uint32_t	j = 0;
	int32_t	    x = 0;
	int32_t     y = 0;

	for (j = 0; j < k_val; ++j) {
		for (i = 0; i < k_val; ++i) {
			if (k_val == n) { // last level
				y = r+j;
				x = c+i;
				if (minval > *(vec + (y * nx_new + x)))
					minval = *(vec + (y * nx_new + x));
				if (maxval < *(vec + (y * nx_new + x)))
					maxval = *(vec + (y * nx_new + x));
				vmax_2d[l][pmax_2d[l]] = *(vec + (y * nx_new + x));
				pmax_2d[l]++;
			}
			else { // internal node
				maxmin child = build_s16(vec, nx_new, n/k_val, l+1, r+j*(n/k_val), c+i*(n/k_val));
				vmax_2d[l][pmax_2d[l]] = child.max_s16;
				vmin_2d[l][pmin_2d[l]] = child.min_s16;
				pmin_2d[l]++;

				if (child.max_s16 != child.min_s16)
					t3_2d[l][pmax_2d[l]] = 1;

				pmax_2d[l]++;

				if (minval > child.min_s16)
					minval = child.min_s16;
				if (maxval < child.max_s16)
					maxval = child.max_s16;
			}
		}
	}

	if (minval == maxval) {
		pmax_2d[l] -= k_val*k_val;
		pmin_2d[l] -= k_val*k_val;
	}

	ret.max_s16 = maxval;
	ret.min_s16 = minval;

	return ret;
}

int initialize(int8_t td, uint32_t *t_size, int32_t *T_size, uint32_t nz)
{
	uint32_t i    = 0;
	uint32_t j    = 0;
	uint32_t k    = 0;
	int32_t  num  = 1;
	uint32_t num_elements = 0;
	
	if (td <= 0) {
		fprintf(stderr, "tree depth cannot be 0 or less than 0.");
		return -1;
	}
	
	//int16_t td;

	//td = floating_point_correction(log(nx), log((double)k_val)) + 1;
	//td = (uint8_t)(ceil(log((double)nx) / log((double)k_val)));

	//*tree_depth = td;
	printf("Tree Depth: %d (Tree Height: %d)\n", td, td-1);

	for (i = 0; i < (uint32_t)td; ++i) {
		if (i != 0 && i != (uint32_t)(td-1))
			*T_size += num;
		num_elements += num;
		num *= k_val * k_val;
	}

	// int num2 = 1;
	// int num_elements2 = 0;
	// for (i = 0; i < td; ++i) {
		// num_elements2 += num2;
		// num2 *= k_val * k_val * k_val;
	// }

	//printf("num_elements2: %d\n", num_elements2);

	if ((T_2d = (uint **)calloc(nz, sizeof(uint *))) == NULL)
		return -12;

	printf("num_elements: %d\n", num_elements);
	printf("T_size: %d\n", *T_size);
	*t_size = (int)ceil((double)(*T_size) / (double)(sizeof(uint32_t) * 8));
	printf("t_size: %d\n", *t_size);
	
	for (k = 0; k < nz; ++k) {
		if ((T_2d[k] = (uint *)calloc(*t_size, sizeof(uint))) == NULL)
			return -24;		
	}

	if ((rMaxMin_2d = (maxmin *)calloc(nz, sizeof(maxmin))) == NULL)
		return -10;
	if ((vmax_2d = (int32_t **)calloc(td, sizeof(int32_t *))) == NULL)
		return -13;
	if ((vmin_2d = (int32_t **)calloc(td, sizeof(int32_t *))) == NULL)
		return -14;
	if ((t3_2d = (uint8_t **)calloc(td, sizeof(uint8_t *))) == NULL)
		return -17;
	if ((pmax_2d = (uint *)calloc(td, sizeof(uint))) == NULL)
		return -18;
	if ((pmin_2d = (uint *)calloc(td, sizeof(uint))) == NULL)
		return -19;
	if ((lmax_2d = (int32_t *)calloc(num_elements-1, sizeof(int32_t))) == NULL) // excludes the element at level 0
		return -22;
	if ((lmin_2d = (int32_t *)calloc(num_elements-1, sizeof(int32_t))) == NULL) // excludes the element at level 0
		return -23;

	num = 1;
	for (j = 0; j < (uint32_t)td; ++j) {
		if ((vmax_2d[j] = (int32_t *)calloc(num, sizeof(int32_t))) == NULL)
			return -25;
		if ((vmin_2d[j] = (int32_t *)calloc(num, sizeof(int32_t))) == NULL)
			return -26;
		if ((t3_2d[j] = (uint8_t *)calloc(num, sizeof(uint8_t))) == NULL)
			return -29;
		num *= k_val * k_val;
	}

	return 0;	
}

void free_mem_2d(int16_t tree_depth, uint32_t nz)
{
	uint32_t j;
	uint32_t k;
	
	for (k = 0; k < nz; ++k) {
		free(T_2d[k]);
	}
	
	for (j = 0; j < (uint32_t)tree_depth; ++j) {
		free(vmax_2d[j]);
		free(vmin_2d[j]);
		free(t3_2d[j]);
	}

	free(rMaxMin_2d);
	free(vmax_2d);
	free(vmin_2d);
	free(t3_2d);
	free(pmax_2d);
	free(pmin_2d);
	free(lmax_2d);
	free(lmin_2d);
	free(T_2d);
}

/*
int16_t ***create_matrix_int16_t(int32_t nx, int32_t ny, int32_t nz)
{
	int16_t ***m;
	int32_t	j, k;

	if ((m = (int16_t ***)calloc(nz, sizeof(int16_t **))) == NULL) {
		return NULL;
	} else {
		for (k = 0; k < nz; ++k) {
			if ((m[k] = (int16_t **)calloc(ny, sizeof(int16_t *))) == NULL)
				return NULL;
			for (j = 0; j < ny; ++j) {
				if ((m[k][j] = (int16_t *)calloc(nx, sizeof(int16_t))) == NULL)
					return NULL;
			}
		}
	}

	return m;
}

int free_matrix_int16_t(int16_t ***matrix, int32_t ny, int32_t nz)
{
	
	if (matrix == NULL)
		 return 1;
	 
	int32_t j;
	int32_t k;

	for (k = 0; k < nz; ++k) {
		for (j = 0; j < ny; ++j) {
			free((matrix)[k][j]);
		}
		free((matrix)[k]);
	}
	free(matrix);

	matrix = NULL;

	return 0;
}
*/

/*
int load_vector(char *command_line_file_name, char *prog_file_name)
{
	int32_t i;
	int32_t j;
	int32_t k;
	int32_t nx_new;
	int32_t ny_new;
	int32_t nz_new;
	char	buf[PATH_MAX] = {0};
	int16_t ***matrix_temp = NULL;
	int16_t *vector_old = NULL;
	int32_t	nx_t;
	int32_t ny_t;
	int32_t nz_t;

	if (*command_line_file_name != '\0')
		strcpy(buf, command_line_file_name);
	else
		strcpy(buf, FILENAME);

	if (file_exists(buf) < 0) {
		fprintf(stderr, "\nFile does not exist. ");
		return -1;
	}

	if (get_base_name(buf, prog_file_name) < 0) {
		fprintf(stderr, "\nCannot get file name. ");
		return -1;
	}

	if ((vector_old = (int16_t *)calloc(nx * ny * nz, sizeof(int16_t))) == NULL) {
		fprintf(stderr, "\nError allocating memory. ");
		return -2;
	}

	FILE *fp;
	fp = fopen(buf, "rb");
	size_t orig_image_size = nx * ny * nz;
	if (orig_image_size != fread(vector_old, sizeof(int16_t), orig_image_size, fp)) {
		fprintf(stderr, "\nError reading file. ");
		return -3;
	}
	fclose(fp);

	nx_t = nx;
	ny_t = ny;
	nz_t = nz;

	//double x_log = floating_point_correction2(log10((double)(nx_t)), log10((double)k_val));
	//double y_log = floating_point_correction2(log10((double)(ny_t)), log10((double)k_val));
	double x_log = ceil(log((double)(nx_t)) / log((double)k_val));
	double y_log = ceil(log((double)(ny_t)) / log((double)k_val));
	
	printf("x_log: %0.19f\n", x_log);
	printf("y_log: %0.19f\n", y_log);

	ny_new = nx_new = pow((double)k_val, max(x_log, y_log));
	nz_new = nz_t;

	printf("nx_new: %d\n", nx_new);
	printf("ny_new: %d\n", ny_new);
	printf("nz_new: %d\n", nz_new);

	if ((img_vector = (int16_t *)calloc(nx_new * ny_new * nz_new, sizeof(int16_t))) == NULL) {
		fprintf(stderr, "\nError allocating memory. ");
		return -2;
	}

	for (k = 0; k < nz_t; ++k) {
		for (j = 0; j < ny_t; ++j) {
			for (i = 0; i < nx_t; ++i) {
				img_vector[k*ny_new*nx_new + j*nx_new + i] = (int16_t)ntohs(vector_old[k*ny_t*nx_t + j*nx_t + i]);
				//img_vector[k*ny_new*nx_new + j*nx_new + i] = matrix_temp[k][j][i];
			}
		}
	}

	nx_old = nx_t;
	ny_old = ny_t;
	nz_old = nz_t;

	nx = nx_new;
	ny = ny_new;
	nz = nz_new;


	free(vector_old);
	free_matrix_int16_t(matrix_temp, ny_t, nz_t);

	return (EXIT_SUCCESS);
}
*/

int find_width(int max_num)
{
	return ((int)log10(max_num)) + 1;
}

void remove_data(const char * dir, const char *file_name, uint16_t k_val, uint32_t nz)
{
	uint32_t k;
	char     buf[PATH_MAX] = {0};
	
	for (k = 0; k < nz; k++) {
		sprintf(buf, "%s/data/%s_k%d_Max_%04d.dat", dir, file_name, k_val, k);
		remove(buf);
		sprintf(buf, "%s/data/%s_k%d_Min_%04d.dat", dir, file_name, k_val, k);
		remove(buf);
		sprintf(buf, "%s/data/%s_k%d_voc_%04d.dat", dir, file_name, k_val, k);
		remove(buf);
		sprintf(buf, "%s/data/%s_k%d_plain_%04d.dat", dir, file_name, k_val, k);
		remove(buf);
		sprintf(buf, "%s/data/%s_k%d_data.dat", dir, file_name, k_val);
		remove(buf);
		sprintf(buf, "%s/data/%s_k%d_encoded_values_%04d.dat", dir, file_name, k_val, k);
		remove(buf);
		sprintf(buf, "%s/data/%s_k%d_is_in_voc_%04d.dat", dir, file_name, k_val, k);
		remove(buf);
	}
}

void print_size(const char *prog_dir, const char *file_name, char heuristic, uint32_t nz)
{
	char buf[PATH_MAX] = {0};
	uint64_t size1 = 0;
	uint64_t size2 = 0;
	uint64_t size3 = 0;
	uint64_t size4 = 0;
	uint64_t size5 = 0;
	uint64_t size6 = 0;
	uint64_t size7 = 0;
	uint32_t k;

	printf("==== SIZE (k=%d) ====\n", k_val);

	sprintf(buf, "%s/data/%s_k%d_data.dat", prog_dir, file_name, k_val);
	size1 = get_file_size(buf);

	printf("T+rMax+rMin:   %10ld Bytes\n", size1);

	// Max
	for (k = 0; k < nz; ++k) {
		sprintf(buf, "%s/data/%s_k%d_Max_%04d.dat", prog_dir, file_name, k_val, k);
		size2 += get_file_size(buf);
	}

	printf("Max:           %10ld Bytes\n", size2);

	for (k = 0; k < nz; ++k) {
		sprintf(buf, "%s/data/%s_k%d_Min_%04d.dat", prog_dir, file_name, k_val, k);
		size3 += get_file_size(buf);
	}

	printf("Min:           %10ld Bytes\n", size3);
	
	if (heuristic) {
		for (k = 0; k < nz; ++k) {
			sprintf(buf, "%s/data/%s_k%d_voc_%04d.dat", prog_dir, file_name, k_val, k);
			size4 += get_file_size(buf);
		}

		printf("Voc:           %10ld Bytes\n", size4);
		
		for (k = 0; k < nz; ++k) {
			sprintf(buf, "%s/data/%s_k%d_plain_%04d.dat", prog_dir, file_name, k_val, k);
			size5 += get_file_size(buf);
		}

		printf("Plain          %10ld Bytes\n", size5);
		
		for (k = 0; k < nz; ++k) {
			sprintf(buf, "%s/data/%s_k%d_encoded_values_%04d.dat", prog_dir, file_name, k_val, k);
			size6 += get_file_size(buf);
		}

		printf("Encoded_words: %10ld Bytes\n", size6);
		
		for (k = 0; k < nz; ++k) {
			sprintf(buf, "%s/data/%s_k%d_is_in_voc_%04d.dat", prog_dir, file_name, k_val, k);
			size7 += get_file_size(buf);
		}

		printf("Is_in_voc:     %10ld Bytes\n", size7);
	}
	printf("==================================================\n");
	if (heuristic)
		printf("Total Size:    %10ld Bytes (%2.3f MB) [k=%d] [File: %s]\n", size1+size2+size3+size4+size5+size6+size7, (double)(size1+size2+size3+size4+size5+size6+size7)/(double)1024/(double)1024, k_val, file_name);
	else
		printf("Total Size:    %10ld Bytes (%2.3f MB) [k=%d] [File: %s]\n", size1+size2+size3, (double)(size1+size2+size3)/(double)1024/(double)1024, k_val, file_name);
}

int64_t get_file_size(const char *fname)
{
	int64_t fsize = -1L;
	FILE *fp = fopen(fname, "rb");

    if (fp) {
		fseek(fp, 0, SEEK_END);
		fsize = ftell(fp);
		fclose(fp);
	}
	return fsize;
}

int file_exists(const char* fname)
{
	return access(fname, F_OK) != -1 ? 0 : -1;
}

int combine_bits_into_int(uint *input_vector, int input_size, uint *output_vector)
{
	int i;
	int counter = 0;

	for (i = 0; i < input_size; ++i) {
		output_vector[counter/W] |= (input_vector[i] << (counter%W));
		counter++;
	}

	return counter;
}

int32_t all_zeros(int32_t *v, uint32_t v_len)
{
	if (v_len == 0)
		return 1;

	int32_t i;

	for (i = 0; i < (int32_t)v_len; i++) {
		if (v[i] != 0) {
			return 0;
		}
	}

	return 1;
}

int get_bit(uint *T, int size, int bit_pos) // bit_pos 0...t_size-1
{
	//printf("\nget_bit: bit_pos: %d\n", bit_pos);
	if (bit_pos >= size)
		return -1;

	return (*(T + (bit_pos/32)) >> (bit_pos%32)) & 0x0001;
}

uint get_rank(const uint *list, uint pos)
{
	uint count = 0;
	uint i;

	if (pos + 1 == 0)
		return 0;
	++pos;

	for (i = 0; i < pos/W; i++) {
		count += __builtin_popcount(*(list+i));
	}
	count += __builtin_popcount(*(list+pos/32) & ((1 << (pos & mask31))-1));

	return count;
}

int get_base_name(char *input_path_name, char *output_base_name)
{
	char* local_file = input_path_name;

	//char* ts1 = strdup(local_file);
	char* ts2 = strdup(local_file);

	//char* dir = dirname(ts1);
	char* fname = basename(ts2);
	if (fname == NULL)
		return -1;
	strcpy(output_base_name, fname);

	printf("File name: %s\n", output_base_name);
	
	free(ts2);
	return 0;
}

double compute_entropy(int16_t *vect, uint32_t vect_len, uint32_t band)
{
	double ret = 0.0;
	int32_t hist[1024*64] = {0};
	int32_t histlen = makehist_16(vect, hist, vect_len);

	ret = entropy(hist, histlen, vect_len);
	
    //printf("Entropy value for band %d = %f\n", band, ret);
	
	return ret;
}


void fill_parent_2d(uint level, uint8_t **t, int32_t *v_selected, int32_t **v_all, uint *p)
{
	int32_t count = 0;
	uint    i = 0;
	for (i = 0; i < pmax_2d[level]; i++) {
		if (t[level][i] == 1) {
			v_selected[count] = v_all[level][i];
			count++;
		}
	}
}

int test_get_cell(const char *command_line_file_name, const char *prog_dir, const char *file_name, uint32_t k_val, char data_type, uint32_t nx_old, uint32_t ny_old, uint32_t nx)
{
	uint32_t	**T_local;
	FTRep     	**ftrep_local;
	FILE		*fp1;
	uint32_t	nz_local;
	uint32_t    T_size_local;
	uint32_t	t_size_local;
	maxmin		*rMaxMin_local;
	uint32_t	i;
	//uint32_t	j;
	uint32_t    k;
	char		buf[PATH_MAX] = {0};
	size_t      result;
	vector_type vector_original = {{0}};

	if (data_type == 2) {
		if ((vector_original.u16 = (uint16_t *)calloc(ny_old * nx_old, sizeof(uint16_t))) == NULL) {
			fprintf(stderr, "\nError allocating memory. ");
			return (EXIT_FAILURE);
		}		
	}
	else if (data_type == 3) {
		if ((vector_original.s16 = (int16_t *)calloc(ny_old * nx_old, sizeof(int16_t))) == NULL) {
			fprintf(stderr, "\nError allocating memory. ");
			return (EXIT_FAILURE);
		}
	}	
	
	sprintf(buf, "%s/data/%s_k%d_data.dat", prog_dir, file_name, k_val);
	fp1 = fopen(buf, "rb");
	result = fread(&nz_local, sizeof(uint32_t), 1, fp1); if (result != 1) {fputs ("Reading error", stderr); return (EXIT_FAILURE);}
	result = fread(&T_size_local, sizeof(uint32_t), 1, fp1); if (result != 1) {fputs ("Reading error", stderr); return (EXIT_FAILURE);}

	t_size_local = (int)ceil((double)(T_size_local) / (double)(sizeof(uint32_t) * 8));

	T_local = (uint32_t **)calloc(nz_local, sizeof(uint32_t *));
	for (k = 0; k < nz_local; ++k)
		T_local[k] = (uint32_t *)calloc(t_size_local, sizeof(uint32_t));
	rMaxMin_local = (maxmin *)calloc(nz_local, sizeof(maxmin));	
		
	ftrep_local = (FTRep **)calloc(nz_local, sizeof(FTRep *));
	
	for (k = 0; k < nz_local; ++k) {
		result = fread(T_local[k], sizeof(uint32_t), t_size_local, fp1); if (result != t_size_local) {fputs ("Reading error", stderr); return (EXIT_FAILURE);}
		if (data_type == 2) {
			result = fread(&rMaxMin_local[k].max_u16, sizeof(uint16_t), 1, fp1); if (result != 1) {fputs ("Reading error", stderr); return (EXIT_FAILURE);}
			result = fread(&rMaxMin_local[k].min_u16, sizeof(uint16_t), 1, fp1); if (result != 1) {fputs ("Reading error", stderr); return (EXIT_FAILURE);}
		}
		else if (data_type == 3) {
			result = fread(&rMaxMin_local[k].max_s16, sizeof(int16_t), 1, fp1); if (result != 1) {fputs ("Reading error", stderr); return (EXIT_FAILURE);}
			result = fread(&rMaxMin_local[k].min_s16, sizeof(int16_t), 1, fp1);	if (result != 1) {fputs ("Reading error", stderr); return (EXIT_FAILURE);}
		}
		sprintf(buf, "%s/data/%s_k%d_Max_%04d.dat", prog_dir, file_name, k_val, k);
		
		ftrep_local[k] = loadFT(buf);		
	}
	fclose(fp1);
	
///////
	
	//size_t	orig_raster_size = ny_old * nx_old;	
	if (data_type == 2) {
		struct timeval time_begin;
		struct timeval time_end;
		double         elapsed_time;
		clock_t        t = clock();
		uint16_t       x_val;
		uint16_t       y_val;
		uint16_t       z_val;

		gettimeofday(&time_begin, NULL);
		//printf("Verifying the first 50 iterations of random access in img_vector[%d][%d][%d]...\n", nx, ny, nz);
		//printf("%d\t%d\t%d\n", nx_old, ny_old, nz);
		//printf("[ x ][ y ][ z ]\n");
		srand(time(0));
		for (i = 0; i < 100000; ++i) {
			x_val = rand() % nx_old;
			y_val = rand() % ny_old;
			z_val = rand() % nz_local;
			if (rMaxMin_local[z_val].max_u16 == 0)
				continue;

			int32_t val = get_cell(T_local[z_val], ftrep_local[z_val], T_size_local, nx, y_val, x_val, -1, rMaxMin_local[z_val].max_u16);
			
			(void) val;
			
			// if (i < 100) {
				// //printf("rMaxMin_local[%d].max_u16 = %d\n", z_val, rMaxMin_local[z_val].max_u16);
				// fp2 = fopen(command_line_file_name, "rb");
				// fseek(fp2, orig_raster_size * z_val * sizeof(uint16_t), SEEK_SET);
				// if (orig_raster_size != fread(vector_original.u16, sizeof(uint16_t), orig_raster_size, fp2)) {
					// fprintf(stderr, "\nError reading file. ");
					// return (EXIT_FAILURE);
				// }
				// for (j = 0; j < orig_raster_size; ++j)
					// vector_original.u16[j] = ntohs(vector_original.u16[j]);
				// printf("%d\t%d\t%d\t%s\t%d\t%d\t%d\n", i, val, vector_original.u16[y_val * nx_old + x_val], 
					// val == vector_original.u16[y_val * nx_old + x_val] ? "same" : "different", x_val, y_val, z_val);
				// fclose(fp2);
			// }
		}
		gettimeofday(&time_end, NULL);
		elapsed_time = (time_end.tv_sec - time_begin.tv_sec) +
					   (time_end.tv_usec - time_begin.tv_usec) * 1.0e-6;

		t = clock() - t;
		printf("Elapsed time : %.3f (s)\n", elapsed_time);
		printf("Time taken for get_cell() to run: %8.5f milliseconds\n", ((double)t/CLOCKS_PER_SEC) * 1000);

		// FILE *fp3 = fopen(command_line_file_name, "rb");
		// for (k = 0; k < nz_local; ++k) {
			// fread(vector_original.u16, sizeof(uint16_t), orig_raster_size, fp3);
			// for (j = 0; j < orig_raster_size; ++j)
				// vector_original.u16[j] = ntohs(vector_original.u16[j]);			
			// int found = 0; 
			// for (j = 0; j < ny_old; ++j)
				// for (i = 0; i < nx_old; ++i) {
					// int32_t val = get_cell(T_local, ftrep_local[k], T_size_local, nx, j, i, -1, rMaxMin_local[k].max_u16);
					// if (val != vector_original.u16[j * nx_old + i]) {
						// //printf("Not the same:\t%d (original)\t%d (restored)\n", vector_original.u16[j * nx_old + i], val);
						// found = 1;
					// }
				// }
			// if (found) 
				// printf("Elements are not the same in band %d\n", k + 1);
			// else
				// printf("Elements are OK in band %d\t", k + 1);
			// fflush(stdout);
		// }
		// fclose(fp3);

	}
	else if (data_type == 3) {
		struct timeval time_begin;
		struct timeval time_end;
		double         elapsed_time;
		clock_t        t = clock();
		uint16_t       x_val;
		uint16_t       y_val;
		uint16_t       z_val;

		gettimeofday(&time_begin, NULL);
		srand(time(0));
		for (i = 0; i < 100000; ++i) {
			x_val = rand() % nx_old;
			y_val = rand() % ny_old;
			z_val = rand() % nz_local;
			if (rMaxMin_local[z_val].max_s16 == 0)
				continue;

			int32_t val = get_cell(T_local[z_val], ftrep_local[z_val], T_size_local, nx, y_val, x_val, -1, rMaxMin_local[z_val].max_s16);
			(void) val;
			// if (i < 100) {
				// fp2 = fopen(command_line_file_name, "rb");
				// fseek(fp2, orig_raster_size * z_val * sizeof(uint16_t), SEEK_SET);
				// if (orig_raster_size != fread(vector_original.s16, sizeof(uint16_t), orig_raster_size, fp2)) {
					// fprintf(stderr, "\nError reading file. ");
					// return (EXIT_FAILURE);
				// }
				// for (j = 0; j < orig_raster_size; ++j)
					// vector_original.s16[j] = ntohs(vector_original.s16[j]);
				// printf("%d\t%d\t%d\t%s\t%d\t%d\t%d\n", i, val, vector_original.s16[y_val * nx_old + x_val], val == vector_original.s16[y_val * nx_old + x_val] ? "same" : "different", x_val, y_val, z_val);
				// fclose(fp2);
			// }
		}
		gettimeofday(&time_end, NULL);
		elapsed_time = (time_end.tv_sec - time_begin.tv_sec) +
					   (time_end.tv_usec - time_begin.tv_usec) * 1.0e-6;

		t = clock() - t;
		printf("Elapsed time : %.3f (ms)\n", elapsed_time * 1000);
		printf("Time taken for get_cell() to run: %8.5f milliseconds\n", ((double)t/CLOCKS_PER_SEC) * 1000);

		// FILE *fp3 = fopen(command_line_file_name, "rb");
		// for (k = 0; k < nz_local; ++k) {
			// fread(vector_original.s16, sizeof(uint16_t), orig_raster_size, fp3);
			// for (j = 0; j < orig_raster_size; ++j)
				// vector_original.s16[j] = ntohs(vector_original.s16[j]);			
			// int found = 0; 
			// for (j = 0; j < ny_old; ++j)
				// for (i = 0; i < nx_old; ++i) {
					// int32_t val = get_cell(T_local, ftrep_local[k], T_size_local, nx, j, i, -1, rMaxMin_local[k].max_s16);
					// if (val != vector_original.s16[j * nx_old + i]) {
						// //printf("Not the same:\t%d (original)\t%d (restored)\n", vector_original.u16[j * nx_old + i], val);
						// found = 1;
					// }
				// }
			// if (found) 
				// printf("Elements are not the same in band %d\n", k + 1);
			// else
				// printf("Elements are OK in band %d\t", k + 1);
			// fflush(stdout);
		// }
		// fclose(fp3);
	}	
	
	
///////
	for (k = 0; k < nz_local; ++k) {
		free(T_local[k]);
		destroyFT(ftrep_local[k]);
	}
	
//printf("A\n");
	
	free(ftrep_local);
	free(T_local);
	free(rMaxMin_local);
//printf("B\n");
	if (data_type == 2) {
		free(vector_original.u16);
	}
	else if (data_type == 3) {
		free(vector_original.s16);
	}
//printf("C\n");
	return (EXIT_SUCCESS);
}

int verify_saved_data(const char *command_line_file_name, const char *prog_dir, const char *file_name, uint32_t k_val, char data_type, uint32_t nx_old, uint32_t ny_old, uint32_t nx)
{
	uint32_t	*T_local;
	FTRep     	*ftrep_local;
	FILE		*fp1;
	FILE		*fp2;
	uint32_t	nz_local;
	uint32_t    T_size_local;
	uint32_t	t_size_local;
	maxmin		*rMaxMin_local;
	uint32_t	i;
	uint32_t	j;
	uint32_t	k;
	char		buf[PATH_MAX] = {0};
	size_t      result;
	vector_type vector_original = {{0}};
	
	if (data_type == 2) {
		if ((vector_original.u16 = (uint16_t *)calloc(ny_old * nx_old, sizeof(uint16_t))) == NULL) {
			fprintf(stderr, "\nError allocating memory. ");
			return (EXIT_FAILURE);
		}		
	}
	else if (data_type == 3) {
		if ((vector_original.s16 = (int16_t *)calloc(ny_old * nx_old, sizeof(int16_t))) == NULL) {
			fprintf(stderr, "\nError allocating memory. ");
			return (EXIT_FAILURE);
		}
	}	
	
	sprintf(buf, "%s/data/%s_k%d_data.dat", prog_dir, file_name, k_val);
	fp1 = fopen(buf, "rb");
	result = fread(&nz_local, sizeof(uint32_t), 1, fp1); if (result != 1) {fputs ("Reading error", stderr); return (EXIT_FAILURE);}
	result = fread(&T_size_local, sizeof(uint32_t), 1, fp1); if (result != 1) {fputs ("Reading error", stderr); return (EXIT_FAILURE);}

	t_size_local = (int)ceil((double)(T_size_local) / (double)(sizeof(uint32_t) * 8));

	T_local = (uint32_t *)calloc(t_size_local, sizeof(uint32_t));
	rMaxMin_local = (maxmin *)calloc(nz_local, sizeof(maxmin));	
		
	fp2 = fopen(command_line_file_name, "rb");
	size_t	orig_raster_size = ny_old * nx_old; 	
	for (k = 0; k < nz_local; ++k) {
		result = fread(T_local, sizeof(uint32_t), t_size_local, fp1); if (result != t_size_local) {fputs ("Reading error", stderr); return (EXIT_FAILURE);}
		if (data_type == 2) {
			result = fread(&rMaxMin_local[k].max_u16, sizeof(uint16_t), 1, fp1); if (result != 1) {fputs ("Reading error", stderr); return (EXIT_FAILURE);}
			result = fread(&rMaxMin_local[k].min_u16, sizeof(uint16_t), 1, fp1); if (result != 1) {fputs ("Reading error", stderr); return (EXIT_FAILURE);}
		}
		else if (data_type == 3) {
			result = fread(&rMaxMin_local[k].max_s16, sizeof(int16_t), 1, fp1); if (result != 1) {fputs ("Reading error", stderr); return (EXIT_FAILURE);}
			result = fread(&rMaxMin_local[k].min_s16, sizeof(int16_t), 1, fp1);	if (result != 1) {fputs ("Reading error", stderr); return (EXIT_FAILURE);}
		}
		sprintf(buf, "%s/data/%s_k%d_Max_%04d.dat", prog_dir, file_name, k_val, k);
		ftrep_local = loadFT(buf);		
		if (data_type == 2) {
			if (orig_raster_size != fread(vector_original.u16, sizeof(uint16_t), orig_raster_size, fp2)) {
				fprintf(stderr, "\nError reading file. ");
				return (EXIT_FAILURE);
			}
			for (i = 0; i < orig_raster_size; ++i)
				vector_original.u16[i] = ntohs(vector_original.u16[i]);

			int found = 0; 
			for (j = 0; j < ny_old; ++j)
				for (i = 0; i < nx_old; ++i) {
					int32_t val = get_cell(T_local, ftrep_local, T_size_local, nx, j, i, -1, rMaxMin_local[k].max_u16);
					if (val != vector_original.u16[j * nx_old + i]) {
						printf("Not the same:\t%d (original)\t%d (restored)\n", vector_original.u16[j * nx_old + i], val);
						found = 1;
					}
				}
			if (found) 
				printf("Elements are not the same in band %d\n", k + 1);
			else
				printf("Elements are OK in band %d\t", k + 1);
			fflush(stdout);
		}
		else if (data_type == 3) {
			if (orig_raster_size != fread(vector_original.s16, sizeof(int16_t), orig_raster_size, fp2)) {
				fprintf(stderr, "\nError reading file. ");
				return (EXIT_FAILURE);
			}
			for (i = 0; i < orig_raster_size; ++i)
				vector_original.s16[i] = ntohs(vector_original.s16[i]);
			
			int found = 0; 
			for (j = 0; j < ny_old; ++j)
				for (i = 0; i < nx_old; ++i) {
					int32_t val = get_cell(T_local, ftrep_local, T_size_local, nx, j, i, -1, rMaxMin_local[k].max_s16);
					if (val != vector_original.s16[j * nx_old + i]) {
						printf("Not the same:\t%d (original)\t%d (restored)\n", vector_original.s16[j * nx_old + i], val);
						found = 1;
					}
				}
			if (found) 
				printf("Elements are not the same in band %d\n", k + 1);
			else
				printf("Elements are OK in band %d\t", k + 1);
			fflush(stdout);
		}
		destroyFT(ftrep_local);
	}
	fclose(fp1);
	fclose(fp2);
	
	free(T_local);
	free(rMaxMin_local);
	
	if (data_type == 2) {
		free(vector_original.u16);
	}
	else if (data_type == 3) {
		free(vector_original.s16);
	}

	return (EXIT_SUCCESS);
}

void print_vector(int16_t *v, int nx, int ny, int b1, int b2, int r1, int r2, int c1, int c2)
{

	if (b2-b1+1 >= 256 || r2-r1+1 >= 256 || c2-c1+1 >= 256) {
		printf("Can only print a matrix with a band or column or row size of less than 256.\n");
		return;
	}

	int32_t i;
	int32_t j;
	int32_t k;

	for (k = b1; k <= b2; ++k) {
		for (j = r1; j <= r2; ++j) {
			for (i = c1; i <= c2; ++i)
				printf("%d\t", (int16_t)ntohs(v[k*ny*nx + j*nx + i]));
			printf("\n");
		}
		printf("************************************\n");
	}
}
/*
int test_get_cell(int32_t nx_old, int32_t ny_old, int32_t nz_old)
{
	clock_t        t;
    struct timeval time_begin;
	struct timeval time_end;
	double         elapsed_time;
	int            i;
	uint16_t       x_val;
	uint16_t       y_val;
	uint16_t       z_val;
	
	srand(time(0));

	printf("*********** %d\t%d\n", nx, ny);

	size_t offset = ny*nx;
	t = clock();

	gettimeofday(&time_begin, NULL);
	printf("Verifying the first 50 iterations of random access in img_vector[%d][%d][%d]...\n", nx, ny, nz);
	printf("%d\t%d\t%d\n", nx_old, ny_old, nz_old);
	printf("[ x ][ y ][ z ]\n");
	for (i = 0; i < 100000; ++i) {
		x_val = rand() % nx_old;
		y_val = rand() % ny_old;
		z_val = rand() % nz_old;
		if (rMaxMin_2d[z_val].max_s16 == 0)
			continue;

		get_cell(nx, z_val, y_val, x_val, -1, rMaxMin_2d[z_val].max_s16);
		//if (iter < 50)
			//printf("[%3d][%3d][%3d]: %3d - %3d\n", x_val, y_val, z_val, get_cell(nx, z_val, y_val, x_val, -1, rMaxMin[z_val].max), orig_value);
		//iter++;
	}
	gettimeofday(&time_end, NULL);
	elapsed_time = (time_end.tv_sec - time_begin.tv_sec) +
				   (time_end.tv_usec - time_begin.tv_usec) * 1.0e-6;

	t = clock() - t;
	printf("Elapsed time : %.3f (s)\n", elapsed_time);
	printf("Time taken for get_cell() to run: %8.5f seconds\n", ((double)t/CLOCKS_PER_SEC));
	
	// repeat without timing
	printf("Repeat without timing...\n");
	srand(time(0));
	t = clock();
	int iter = 0;
	for (i = 0; i < 100000; ++i) {
		x_val = rand() % nx_old;
		y_val = rand() % ny_old;
		z_val = rand() % nz_old;
		if (rMaxMin_2d[z_val].max_s16 == 0)
			continue;

		int16_t curr_value = get_cell(nx, z_val, y_val, x_val, -1, rMaxMin_2d[z_val].max_s16);
		int16_t orig_value = img_vector[z_val * offset + y_val * nx + x_val];
		if (iter < 50)
			printf("%2d: [%3d][%3d][%4d]:\t%5d\t%5d\t(%s)\n", iter+1, x_val, y_val, z_val, get_cell(nx, z_val, y_val, x_val, -1, rMaxMin_2d[z_val].max_s16), orig_value, curr_value == orig_value ? "Same" : "Different" );
		else if (curr_value != orig_value)
			printf("PIXEL NOT THE SAME:  [%3d][%3d][%4d]:\t%5d\t%5d\n", x_val, y_val, z_val, get_cell(nx, z_val, y_val, x_val, -1, rMaxMin_2d[z_val].max_s16), orig_value);
		// if (iter < 20 || curr_value != orig_value)
			// printf("[%3d][%3d][%3d]: %3d - %3d\n", x_val, y_val, z_val, get_cell(nx, z_val, y_val, x_val, -1, rMaxMin[z_val].max), orig_value);
		iter++;
	}
	
	return (EXIT_SUCCESS);
}
*/

int get_dir(char *dir)
{
	char pBuf[PATH_MAX];
	size_t len = sizeof(pBuf);
	char szTmp[32];
	sprintf(szTmp, "/proc/%d/exe", getpid());
	int bytes = MIN((int)readlink(szTmp, pBuf, len), (int)(len - 1));
	if (bytes >= 0)
		pBuf[bytes] = '\0';
	
	
	printf("%s\n", pBuf);
	
	char* ts1 = strdup(pBuf);
	//char* ts2 = strdup(pBuf);

	strcpy(dir, dirname(ts1));
	//char* filename = basename(ts2);	
	
	printf("%s\n", dir);
	//printf("%s\n", filename);	
	
	if (ts1 != NULL)
		free(ts1);
	
	return (EXIT_SUCCESS);
}

int convert_to_es_format(char *string, int len)
{
	int i;
	
	for (i = 0; i < len; ++i) {
		if (string[i] == '.')
			string[i] = ',';
		else if (string[i] == ',')
			string[i] = '.';
	}
	
	return 1;
}

void save_entropy_bbp_to_file(int k, int16_t *img_vector, uint64_t fsMax, uint64_t fsMin, uint t_size, uint32_t nx, uint32_t ny, FILE *fp)
{
	double entropy_val = compute_entropy(img_vector, nx * ny, k);
	double k2r_bbp_val = ((double)((fsMax + fsMin + sizeof(uint) * t_size + sizeof(int16_t) * 2) * 8) / (double)(nx * ny)); // includes max and min of root

	char number_format[512] = {0};
	sprintf(number_format, "%f\t%f\n", entropy_val, k2r_bbp_val);
	convert_to_es_format(number_format, strlen(number_format));
	fprintf(fp, "%s", number_format);	
}

int cmpfunc(const void *a, const void *b)
{
	uint32_t x = ((struct voc *) a)->weight;
	uint32_t y = ((struct voc *) b)->weight;		
	
    if (x > y) {
        return -1;
    } else if (x < y) {
        return 1;
    } else {
        return 0;
    }
}

bool sort_by_count(const struct voc &a, const struct voc &b)
{
	if (a.weight > b.weight)
		return true;
	else if (a.weight < b.weight)
		return false;
	
	return false;
}

struct voc *compute_submatrices_freq(int32_t *vmax_last, uint32_t vmax_last_count, uint32_t &vocab_count, int band, uint32_t submatrix_size, double &entropy_val)
{
	uint32_t i = 0;
	uint32_t j = 0;
	uint32_t submatrix_count = vmax_last_count / submatrix_size;
	
	map<list<int32_t>, uint32_t> vocab_map;
			
	for (i = 0; i < submatrix_count; ++i) {
		list<int32_t> ls;
		for (j = 0; j < submatrix_size; ++j) {
			ls.push_back(vmax_last[i * submatrix_size + j]);
		}
		vocab_map[ls]++;
	}
	
	vocab_count = vocab_map.size();
	
	//printf("%d\tvocabulary_count: %d\tvmax_last_count: %d\n", band + 1, vocab_count, submatrix_count);
	
	struct voc *vocabulary = (struct voc *) calloc(vocab_count, sizeof(struct voc));
	
	memset(vocabulary, 0, sizeof(struct voc) * (vocab_count));
	int counter1 = 0;
	for (auto it = vocab_map.begin(); it != vocab_map.end(); ++it) {
		vocabulary[counter1].submatrix = (int32_t *)calloc(submatrix_size, sizeof(int32_t));
		int counter2 = 0;
		//for (auto lit = it->first.begin(); lit != it->first.end(); lit++) {
		for (auto x : it->first) {
			vocabulary[counter1].submatrix[counter2] = x;
			counter2++;
		}
		vocabulary[counter1].weight = it->second;
		counter1++;
	}
		
	int32_t v[vocab_count];
	for (i = 0; i < vocab_count; ++i)
		v[i] = vocabulary[i].weight;
	entropy_val = entropy(v, vocab_count, submatrix_count);
	//printf("Entropy - submatrix: %f\n", entropy_val);		
	
	return vocabulary;
}

struct voc *compute_values_freq(int32_t *vmax_last, uint32_t vmax_last_count, uint32_t &vocab_count, int band, double &entropy_val)
{
	uint32_t i = 0;
	
	map<int32_t, uint32_t> vocab_map;
			
	for (i = 0; i < vmax_last_count; ++i) {
		vocab_map[vmax_last[i]]++;
	}			
	
	vocab_count = vocab_map.size();
	
	//printf("%d\tvocabulary_count: %d\tvmax_last_count: %d\n", band + 1, vocab_count, vmax_last_count);
	
	struct voc *vocabulary = (struct voc *) calloc(vocab_count, sizeof(struct voc));
	
	memset(vocabulary, 0, sizeof(struct voc) * (vocab_count));
	int counter1 = 0;
	for (auto it = vocab_map.begin(); it != vocab_map.end(); ++it) {
		vocabulary[counter1].submatrix_value = it->first;
		vocabulary[counter1].weight = it->second;
		counter1++;
	}
	
	// double entropy_val = compute_entropy(vmax_last, vmax_last_count, band);
	// printf("Entropy: %f\n", entropy_val);
	
	int32_t v[vocab_count];
	for (i = 0; i < vocab_count; ++i)
		v[i] = vocabulary[i].weight;
	entropy_val = entropy(v, vocab_count, vmax_last_count);
	//printf("Entropy: %f\n", entropy_val);
	
	return vocabulary;
} 


uint16_t find_max(uint32_t *lmax_2d, uint32_t lmax_len_2d)
{
	uint32_t i;
	uint16_t max_val = 0;
	
	for (i = 0; i < lmax_len_2d; ++i)
		if (max_val < lmax_2d[i])
			max_val = lmax_2d[i];
		
	return max_val;
}

uint8_t count_bits(uint32_t n) 
{ 
	uint8_t count = 0; 
	while (n) { 
        count++; 
        n >>= 1; 
    } 
    return count; 
} 

struct k_val_data {
	uint8_t  index;
	uint32_t size;
} k_val_data;

bool sort_k_val(const struct k_val_data &a, const struct k_val_data &b)
{
	if (a.size < b.size)
		return true;
	else if (a.size > b.size)
		return false;
	
	return false;
}

int16_t find_best_k_val(uint32_t nx, uint32_t ny)
{
	//int16_t ret = 2;
	uint32_t i;
	
	struct k_val_data d[19];
	
	for (i = 2; i <= 20; ++i) {
		double x_log = (int)ceil(log((double)(nx)) / log((double)i));
		double y_log = (int)ceil(log((double)(ny)) / log((double)i));
		d[i-2].index = i;
		d[i-2].size = pow((double)i, (double)max(x_log, y_log));
	}
	
	std::sort(d, d+19, sort_k_val);
	
	return d[0].index;
}

int main(int argc, char **argv)
{
	printf("%d %d %d \n", 1 << 1, 1 << 2, 1 <<3);
	
	uint32_t  	i;
	uint32_t  	j;
	uint32_t  	k;
	int  		T_counter = 0;
	int			err_code;
	int32_t     *vmax_parent_2d;
	int32_t     *vmin_parent_2d;
	
	char        prog_dir[PATH_MAX] = {0};
	char		buf[PATH_MAX]      = {0};
	char    	command_line_file_name[PATH_MAX] = {0};
	char		prog_file_name[PATH_MAX] = {0}; // File name of the image used in the program, doesn't include path.
	int8_t		tree_depth = 0;                 // Tree Depth = Tree Height + 1
	int32_t     nx_new;
	int32_t     ny_new;
	uint32_t    nx_old;
	uint32_t    ny_old;
	uint32_t    t_size = 0;
	int32_t		T_size = 0;
	uint32_t	lmax_len_2d;
	uint32_t	lmin_len_2d;
	vector_type vector_old = {{0}};
	vector_type vector_new = {{0}};
	FTRep 		*dacMax_2d = NULL;
	FTRep 		*dacMin_2d = NULL;
	uint32_t	nx = 1024;
	uint32_t 	ny = 1024;
	uint32_t 	nz = 1;

	// user options
	char        data_type    = 2; // 2 = unsigned 16-bit; 3 = signed 16-bit.
	char        heuristic    = 0; // 0 = original; 1 = heuristic.
//	char        show_entropy = 0; // 0 = not showing entropy; 1 = showing entropy;
	
	k_val = 2; // default is 2

	if (argc < 7) {
		fprintf(stderr, "Required parameters missing.\n");
		fprintf(stderr, "Optimal k values - \n");
		fprintf(stderr, "airs:\t\t\t6\n");
		fprintf(stderr, "aviris_calibrated:\t6\n");
		fprintf(stderr, "aviris_uncalibrated:\t9\n");
		fprintf(stderr, "crism:\t\t\t5, 6\n");
		fprintf(stderr, "hyperion_calibrated:\t8\n");
		fprintf(stderr, "hyperion_uncalibrated:\t8\n");
		return (EXIT_FAILURE);
	}
	strcpy(command_line_file_name, argv[1]);
	nz        = atoi(argv[2]); // get  first command line parameter nz
	ny        = atoi(argv[3]); // get second command line parameter ny
	nx        = atoi(argv[4]); // get  third command line parameter nx
	data_type = atoi(argv[5]); // get  third command line parameter data_type
	k_val     = atoi(argv[6]); // get k value from command line
	
	//fprintf(stderr, "k = %d\n", k_val);
	
	printf("Number of bands:   nz = %d\n", nz);
	printf("Number of rows:    ny = %d\n", ny);
	printf("Number of columns: nx = %d\n", nx);
	
	printf("Best k value: %d\n", find_best_k_val(nx, ny));
	
	get_dir(prog_dir);

	//////////////////////////// BEGINNING OF PART 1

	if (data_type == 2) {
		if ((vector_old.u16 = (uint16_t *)calloc(ny * nx, sizeof(uint16_t))) == NULL) {
			fprintf(stderr, "\nError allocating memory. ");
			return -2;
		}		
	}
	else if (data_type == 3) {
		if ((vector_old.s16 = (int16_t *)calloc(ny * nx, sizeof(int16_t))) == NULL) {
			fprintf(stderr, "\nError allocating memory. ");
			return -2;
		}
	}

	double x_log = (int)ceil(log((double)(nx)) / log((double)k_val));
	double y_log = (int)ceil(log((double)(ny)) / log((double)k_val));
	
	printf("1st step x_log: %0.19f\n", x_log);
	printf("1st step y_log: %0.19f\n", y_log);

	tree_depth = (uint8_t)max(x_log, y_log) + 1;

	ny_new = nx_new = pow((double)k_val, tree_depth-1);

	printf("nx_new: %d\n", nx_new);
	printf("ny_new: %d\n", ny_new);

	if (data_type == 2) {
		if ((vector_new.u16 = (uint16_t *)calloc(ny_new * nx_new, sizeof(uint16_t))) == NULL) {
			fprintf(stderr, "\nError allocating memory. ");
			return (EXIT_FAILURE);
		}		
	}
	else if (data_type == 3) {
		if ((vector_new.s16 = (int16_t *)calloc(ny_new * nx_new, sizeof(int16_t))) == NULL) {
			fprintf(stderr, "\nError allocating memory. ");
			return (EXIT_FAILURE);
		}
	}

	nx_old = nx;
	ny_old = ny;

	nx = nx_new;
	ny = ny_new;

	printf("main() 1 - Initialize()\n");
	if ((err_code = initialize(tree_depth, &t_size, &T_size, nz)) < 0) {
		fprintf(stderr, "Memory allocation error. Error Code: %d\n", err_code);
		return (EXIT_FAILURE);
	}
	
	strcpy(buf, command_line_file_name);

	if (file_exists(buf) < 0) {
		fprintf(stderr, "\nFile does not exist. ");
		return -1;
	}

	if (get_base_name(buf, prog_file_name) < 0) {
		fprintf(stderr, "\nCannot get file name. ");
		return -1;
	}

	FILE *fp1 = fopen(buf, "rb");
	sprintf(buf, "%s/data/%s_k%d_data.dat", prog_dir, prog_file_name, k_val);
	FILE *fp2 = fopen(buf, "wb");
	fwrite(&nz, sizeof(uint32_t), 1, fp2);
	fwrite(&T_size, sizeof(uint32_t), 1, fp2);
	//sprintf(buf, "./entropy-bitrate.txt");
	//FILE *fp4 = fopen(buf, "w");
	
	// sprintf(buf, "%s/data/%s_k%d_max_size.dat", prog_dir, prog_file_name, k_val);
	// FILE *fp100 = fopen(buf, "w");
	
	// sprintf(buf, "%s/data/%s_k%d_min_size.dat", prog_dir, prog_file_name, k_val);
	// FILE *fp101 = fopen(buf, "w");
		
	
	

	// start timing
	clock_t        t;
	struct timeval time_begin;
	struct timeval time_end;
	double         elapsed_time;
	
	// variables for showing entropy
	//uint64_t       fsMax = 0;
	//uint64_t       fsMin = 0;

	t = clock();
	gettimeofday(&time_begin, NULL);
	
	uint how_many = 0;
	
	size_t		orig_raster_size = ny_old * nx_old; 
	uint32_t	total_code_word_size = 0;
	for (k = 0; k < nz; k++) {
		if (k % 10 == 0) {
			printf("\nk=%d\n", k);
			fflush(stdout);
		}
		if (data_type == 2) {
			if (orig_raster_size != fread(vector_old.u16, sizeof(uint16_t), orig_raster_size, fp1)) {
				fprintf(stderr, "\nError reading file. ");
				return (EXIT_FAILURE);
			}
		}
		else if (data_type == 3) {
			if (orig_raster_size != fread(vector_old.s16, sizeof(int16_t), orig_raster_size, fp1)) {
				fprintf(stderr, "\nError reading file. ");
				return (EXIT_FAILURE);
			}
		}

		for (j = 0; j < ny_old; ++j) {
			for (i = 0; i < nx_old; ++i) {
				if (data_type == 2) {
					
					vector_new.u16[j * nx + i] = (uint16_t)ntohs(vector_old.u16[j * nx_old + i]);
				}
				else if (data_type == 3) {
					vector_new.s16[j * nx + i] = (int16_t)ntohs(vector_old.s16[j * nx_old + i]);
				}
			}
		}

		uint32_t contador = 1;
		for (j = 0; j < (uint32_t)tree_depth; ++j) {
			pmax_2d[j] = 0;
			pmin_2d[j] = 0;
			for (i = 0; i < contador; ++i) {
				vmax_2d[j][i] = 0;
				vmin_2d[j][i] = 0;
				t3_2d[j][i] = 0; 
			}
			contador *= k_val * k_val;
		}
	
		if (data_type == 2) {
			rMaxMin_2d[k] = build_u16(vector_new.u16, nx, nx, 1, 0, 0);
		}
		else if (data_type == 3) {
			rMaxMin_2d[k] = build_s16(vector_new.s16, nx, nx, 1, 0, 0);
		}
		//printf("rMax: %d  rMin: %d\n", rMaxMin_2d[k].max, rMaxMin_2d[k].min);
		
		// Stuffing the T3 array elements into T
		T_counter = 0;
		for (j = 1; j < (uint32_t)(tree_depth-1); ++j) {
			for (i = 0; i < pmax_2d[j]; i++) {
				T_2d[k][T_counter/W] |= (t3_2d[j][i] << (T_counter%W));
				T_counter++;
			}
		}

		if (data_type == 2) {
			vmax_2d[0][0] = rMaxMin_2d[k].max_u16;
			vmin_2d[0][0] = rMaxMin_2d[k].min_u16;
		}
		else if (data_type == 3) {
			vmax_2d[0][0] = rMaxMin_2d[k].max_s16;
			vmin_2d[0][0] = rMaxMin_2d[k].min_s16;
		}
		pmax_2d[0] = 1;
		pmin_2d[0] = 1;
		t3_2d[0][0] = 1;
		
		int lmax_total = 0;
		int lmin_total = 0;
		
		int32_t  size_vmax_last = pow(k_val * k_val, tree_depth-1);
		int32_t  vmax_last[size_vmax_last];
		uint32_t vmax_last_count = 0;
		// // printf("k2r_plain - Concatenate vmax into lmax\n");
	    // maximum
		lmax_len_2d = 0;
		for (j = 0; j < (uint32_t)(tree_depth-1); ++j) {
			vmax_parent_2d = (int32_t *)calloc(pow(k_val*k_val, j), sizeof(int32_t));
			fill_parent_2d(j, t3_2d, vmax_parent_2d, vmax_2d, pmax_2d);
			for (i = 0; i < pmax_2d[j+1]; i++) {
				if (heuristic) {
					if (j == (uint32_t)(tree_depth - 2)) {
						vmax_last[vmax_last_count] = vmax_parent_2d[i/(k_val*k_val)] - vmax_2d[j+1][i];
						vmax_last_count++;
					} else {
						lmax_2d[lmax_len_2d] = vmax_parent_2d[i/(k_val*k_val)] - vmax_2d[j+1][i];
						lmax_len_2d++;
					}
				} else {
					lmax_2d[lmax_len_2d] = vmax_parent_2d[i/(k_val*k_val)] - vmax_2d[j+1][i];
					lmax_len_2d++;
				}
			}
			free(vmax_parent_2d);
		}

		lmax_total += lmax_len_2d;

		// // printf("k2r_plain - Concatenate vmin into lmin\n");
		// minimum
		lmin_len_2d = 0;
		for (j = 0; j < (uint32_t)(tree_depth-2); j++) {
			vmin_parent_2d = (int32_t *)calloc(pow(k_val * k_val, j), sizeof(int32_t));
			fill_parent_2d(j, t3_2d, vmin_parent_2d, vmin_2d, pmin_2d);
			for (i = 0; i < pmin_2d[j+1]; i++) {
				if (t3_2d[j+1][i] == 1) {
					lmin_2d[lmin_len_2d] = vmin_2d[j+1][i] - vmin_parent_2d[i/(k_val*k_val)];
					lmin_len_2d++;
				}
			}
			free(vmin_parent_2d);
		}
		lmin_total += lmin_len_2d;
		
		// // printf("Saving all the files\n");
		sprintf(buf, "%s/data/%s_k%d_Max_%04d.dat", prog_dir, prog_file_name, k_val, k);
		if (all_zeros((int32_t *)lmax_2d, lmax_len_2d)) {
			printf("ALL ZEROS MAX: %d\n", k);
			FILE * fp3 = fopen(buf, "w");
			fclose(fp3);
		} else {
			#if USE_LIMITED_LEVEL
				dacMax_2d = createFT((uint32_t *)lmax_2d, lmax_len_2d, 100);
			#else
				dacMax_2d = createFT((uint32_t *)lmax_2d, lmax_len_2d);
			#endif
			saveFT(dacMax_2d, buf);
			//fsMax = get_file_size(buf);
			destroyFT(dacMax_2d);
			
			//uint16_t max_val_1 = find_max((uint32_t *)lmax_2d, lmax_len_2d);
			//if (k < 100)
				//printf("max_val_1 = %d - %d - %d\n", max_val_1, count_bits(max_val_1), lmax_len_2d);
			// sprintf(buf2, "%s/data/nodacs_%s_k%d_Max_%04d.dat", prog_dir, prog_file_name, k_val, k);
			// FILE *fp6 = fopen(buf2, "w");
			// fwrite(lmax_2d, sizeof(uint32_t), lmax_len_2d, fp6);
			// fclose(fp6);
			// fwrite(&lmax_len_2d, sizeof(uint32_t), 1, fp100);
		}

		sprintf(buf, "%s/data/%s_k%d_Min_%04d.dat", prog_dir, prog_file_name, k_val, k);
		if (all_zeros((int32_t *)lmin_2d, lmin_len_2d)) {
			printf("ALL ZEROS MIN: %d\n", k);
			FILE * fp4 = fopen(buf, "w");
			fclose(fp4);
		} else {
			#if USE_LIMITED_LEVEL
				dacMin_2d = createFT((uint32_t *)lmin_2d, lmin_len_2d, 100);
			#else
				dacMin_2d = createFT((uint32_t *)lmin_2d, lmin_len_2d);
			#endif
			saveFT(dacMin_2d, buf);
			//fsMin = get_file_size(buf);
			destroyFT(dacMin_2d);
			
			//uint16_t max_val_2 = find_max((uint32_t *)lmin_2d, lmin_len_2d);
			//if (k < 100)
				//printf("max_val_2 = %d - %d - %d\n", max_val_2, count_bits(max_val_2), lmin_len_2d);
			// sprintf(buf2, "%s/data/nodacs_%s_k%d_Min_%04d.dat", prog_dir, prog_file_name, k_val, k);
			// FILE *fp6 = fopen(buf2, "w");
			// fwrite(lmin_2d, sizeof(uint32_t), lmin_len_2d, fp6);
			// fclose(fp6);
			// fwrite(&lmin_len_2d, sizeof(uint32_t), 1, fp101);
		}

		fwrite(T_2d[k], sizeof(uint32_t), t_size, fp2);
		if (data_type == 2) {
			fwrite(&rMaxMin_2d[k].max_u16, sizeof(uint16_t), 1, fp2);
			fwrite(&rMaxMin_2d[k].min_u16, sizeof(uint16_t), 1, fp2);
		}
		else if (data_type == 3) {
			fwrite(&rMaxMin_2d[k].max_s16, sizeof(int16_t), 1, fp2);
			fwrite(&rMaxMin_2d[k].min_s16, sizeof(int16_t), 1, fp2);
		}

		// if (show_entropy) {
			// if (data_type == 2)
				// save_entropy_bbp_to_file(k, (int16_t *)vector_old.u16, fsMax, fsMin, t_size, nx_old, ny_old, fp4);
			// else
				// save_entropy_bbp_to_file(k, vector_old.s16, fsMax, fsMin, t_size, nx_old, ny_old, fp4);
		// }
		
		if (heuristic) {
			// variables for heuristic
			list<struct voc> final_voc;
			list<bool>       in_voc;
			list<uint32_t>   encoded_code_words;
			list<int32_t>    plain_value;
			uint32_t         vocab_count = 0;
			struct voc       *vocabulary = NULL;
			double           entropy_voc = 0.0;

			///// 1st step heuristic
			vocabulary = compute_submatrices_freq(vmax_last, vmax_last_count, vocab_count, k, k_val * k_val, entropy_voc);
			
			// printf("Before sorting...\n");
			// for (i = 0; i < (int32_t) vocab_count; ++i)
			// {
				// printf("%d\t%d\t%d\t%d\t%d\t%d\n", i, vocabulary[i].submatrix[0], vocabulary[i].submatrix[1],
					// vocabulary[i].submatrix[2], vocabulary[i].submatrix[3], vocabulary[i].weight);
			// }
			
			std::sort(vocabulary, vocabulary+vocab_count, sort_by_count);

			// printf("After sorting...\n");
			// for (i = 0; i < (int32_t) vocab_count; ++i)
			// {
				// printf("%d\t%d\t%d\t%d\t%d\t%d\n", i, vocabulary[i].submatrix[0], vocabulary[i].submatrix[1],
					// vocabulary[i].submatrix[2], vocabulary[i].submatrix[3], vocabulary[i].weight);
			// }		
			
			///// 2nd step heuristic
			uint32_t vocab_count2 = 0;
			double entropy_value = 0.0;
			struct voc *vocabulary2 = compute_values_freq(vmax_last, vmax_last_count, vocab_count2, k, entropy_value);
			// printf("Before sorting...\n");
			// for (i = 0; i < (int32_t) vocab_count2; ++i) {
				// printf("%d\t%d\t%d\n", i+1, vocabulary2[i].submatrix_value, vocabulary2[i].weight);
			// }
			
			std::sort(vocabulary2, vocabulary2 + vocab_count2, sort_by_count);

			// printf("After sorting...\n");
			// for (i = 0; i < (int32_t) vocab_count2; ++i) {
				// printf("%d\t%d\t%d\n", i+1, vocabulary2[i].submatrix_value, vocabulary2[i].weight);
			// }		

			////// 4th step heuristic
			uint32_t code_word = 0;
			uint32_t max_compacted_values = 0; 
			uint32_t max_no_compacted_values = 0;

			double    pPlain;
			double    pVoc;
			uint32_t  p;
			for (i = 0; i < vocab_count; ++i) {
				p = vocabulary[i].weight;
				pPlain = entropy_value * p * (k_val * k_val);
				pVoc   = entropy_voc * p + (k_val * k_val) * sizeof(int16_t) * 8;
				
				if (pVoc <= pPlain) {
					in_voc.push_back(1);
					encoded_code_words.push_back(code_word);
					struct voc temp;
					temp.submatrix = (int32_t *)calloc(k_val * k_val, sizeof (int32_t));
					//temp.index = vocabulary[i].index = code_word;
					memcpy(temp.submatrix, vocabulary[i].submatrix, k_val * k_val * sizeof(int32_t));
					temp.weight = vocabulary[i].weight;
					final_voc.push_back(temp);
					code_word++;
					max_compacted_values += p;
					//free(temp.submatrix);
				} else {
					in_voc.push_back(0);
					for (j = 0; j < (uint32_t)(k_val * k_val); ++j) {
						plain_value.push_back((vocabulary[i].submatrix[j]));
					}
					max_no_compacted_values += p;
				}
			}
			
			// if (final_voc.size() > 0) {
				// for (auto x: final_voc) {
					// printf("Result: \n");
					// for (j = 0; j < k_val * k_val; ++j)
						// printf("%d\t", x.submatrix[j]);
					// printf("\n");
				// }
			// }
			total_code_word_size += code_word;

			// vocabulary
			sprintf(buf, "%s/data/%s_k%d_voc_%0d.dat", prog_dir, prog_file_name, k_val, k);
			if (final_voc.size() > 0) {
				
				int32_t tmp1[final_voc.size() * k_val * k_val];
				int c1 = 0;

				for (auto it1 : final_voc) {
					for (j = 0; j < (uint32_t)(k_val * k_val); ++j) {
						printf("%d\n", it1.submatrix[j]);
						tmp1[c1++] = it1.submatrix[j];
					}
				}
				#if USE_LIMITED_LEVEL
					FTRep *vocFT = createFT((uint32_t *) tmp1, c1, 100);
				#else
					FTRep *vocFT = createFT((uint32_t *) tmp1, c1);
				#endif
				saveFT(vocFT, buf);
				destroyFT(vocFT);		
			} else {
				FILE * fp4 = fopen(buf, "w");
				fclose(fp4);			
			}
			// plain
			printf("k = %d\n", k);
			sprintf(buf, "%s/data/%s_k%d_plain_%04d.dat", prog_dir, prog_file_name, k_val, k);
			if (plain_value.size() > 0) {
				uint32_t tmp2[plain_value.size()];
				int c2 = 0;
				for (auto it : plain_value) {
					//if (plain_value.size() == 347064 && it < 0)
						//printf("%d\t", it);
					if (it < 0) {
						how_many++;
						printf("k = %d\n", k);
					}
					tmp2[c2++] = it;
				}
				if (all_zeros((int32_t *)tmp2, c2)) {
					FILE * fp4 = fopen(buf, "w");
					fclose(fp4);
				} else {
					// FILE *fp = fopen(buf, "w");
					// fwrite(tmp2, sizeof(uint32_t), c2, fp);
					// fclose(fp);
					#if USE_LIMITED_LEVEL
					FTRep *plainFT = createFT(tmp2, c2, 100);
					#else
					FTRep *plainFT = createFT(tmp2, c2);
					#endif
					if (plainFT) {
						saveFT(plainFT, buf);
						destroyFT(plainFT);
					}
				}
			} else {
				FILE * fp4 = fopen(buf, "w");
				fclose(fp4);				
			}
			// encoded values
			sprintf(buf, "%s/data/%s_k%d_encoded_values_%04d.dat", prog_dir, prog_file_name, k_val, k);
			FILE *fp5 = fopen(buf, "w");
			if (encoded_code_words.size() > 0) {
				for (auto x : encoded_code_words) {
					fwrite(&x, sizeof(uint32_t), 1, fp5);
				}
			}
			fclose(fp5);
			
			// is it in vocabulary?
			sprintf(buf, "%s/data/%s_k%d_is_in_voc_%04d.dat", prog_dir, prog_file_name, k_val, k);
			FILE *fp6 = fopen(buf, "w");
			if (in_voc.size() > 0) {
				uint32_t nelements = ceil(in_voc.size() / W);
				uint32_t is_in_voc[nelements];
				int counter = 0;
				for (auto x : in_voc) {
					is_in_voc[counter/W] |= x << (counter % W);
					counter++;
				}

				fwrite(&is_in_voc, sizeof(uint32_t), nelements, fp6);
			}
			fclose(fp6);

			// final - free memory:	
			for (auto x : final_voc)
				if (x.submatrix)
					free(x.submatrix);
			for (i = 0; i < vocab_count; ++i) {
				if (vocabulary[i].submatrix)
					free(vocabulary[i].submatrix);
				vocabulary[i].weight = 0;
				//vocabulary[i].index = 0;
			}
		} // end heuristic

	} // end k
	printf("\ntotal_code_word %d\n", total_code_word_size);
	
	fclose(fp1);
	fclose(fp2);
	//fclose(fp4);
	
	//fclose(fp100);
	//fclose(fp101);
	
	t = clock() - t;
	gettimeofday(&time_end, NULL);
	elapsed_time = (time_end.tv_sec - time_begin.tv_sec) +
				   (time_end.tv_usec - time_begin.tv_usec) * 1.0e-6;

	printf("Elapsed time for building k2-raster (using gettimeofday()): %8.5f milliseconds\n", elapsed_time * 1000);
	printf("Elapsed time for building k2-raster (using clock()):        %8.5f milliseconds\n", ((double)t/CLOCKS_PER_SEC) * 1000);
	//end timing
	
	printf("HOW_MANY: %d\n", how_many);
	
	print_size(prog_dir, prog_file_name, heuristic, nz);

	free_mem_2d(tree_depth, nz);
	
	//////////////////////////// END OF PART 1
	
	// compare the saved data and the original data
	// printf("\nReload data and compare them with the original data...\n");
	// verify_saved_data(command_line_file_name, prog_dir, prog_file_name, k_val, data_type, nx_old, ny_old, nx);

	NL;
	
	// testing get_cell
	//printf("\nTesting GetCell()...\n");
	//test_get_cell(command_line_file_name, prog_dir, prog_file_name, k_val, data_type, nx_old, ny_old, nx);

	NL;

	if (data_type == 2) {
		free(vector_new.u16);
		free(vector_old.u16);
	}
	else if (data_type == 3) {
		free(vector_new.s16);
		free(vector_old.s16);
	}

	NL;

	printf("Now deleting data files to save space... Please wait\n");
	remove_data(prog_dir, prog_file_name, k_val, nz);
	printf("Done deleting.\n");

	NL;
	
	return (EXIT_SUCCESS);
} // END main()


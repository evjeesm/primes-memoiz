#ifndef _PRIMES_H_
#define _PRIMES_H_

#include "dynarr.h"

#include <stdbool.h>

#define REFLECT(name, code) \
    __attribute__((weak)) const char* name = #code; \
    code

/*
* Can be tweaked depending on file system limitation.
*/
#define MAX_FILE_SIZE 2199023255552

/*
* Cache control struct, maps single file block at a time.
* Utilizes sparce file storage, where a file can grow up to fs limit measured in TiB,
* but since it is mostly empty, zero pages won't be physically allocated.
*/
typedef struct cache
{
    int       fd;
    size_t    file_idx;
    size_t    page_offset; /* current page offset */
    char      *page;       /* mapped */
}
cache_t;

/*
* Each odd number in the cache takes 2 bits of storage.
* Which means that each cache page stores primerility for (3 * page_size) odd numbers.
*/
typedef enum cache_value
{
    UNDEFINED = 0,
    PRIME,
    NOT_PRIME,
    VALUES_TOTAL
}
cache_value_t;


REFLECT(pair_definition,
    typedef struct pair
    {
        size_t first;
        size_t second;
    }
    pair_t
);

/*
* Function that performs primerility test for a number.
* Complexity: sqrt(N)/3, where N -> number
*/
bool is_prime(size_t number);

/*
* Initialize/Deinitialize cache for primes memoization.
*/
void init_cache(void);
void fini_cache(void);

/*
* Finding prime using memoization and updating cache on the way.
*/
bool is_prime_cached(size_t number);

/*
* Factory function for vector that stores primes.
*/
dynarr_t *create_primes_array(void);

/*
* Function checks range of numbers from `begin` to `end` inclusive for primerility.
* Result stored in a vector `out` that has to be created prior this call.
*/
void get_primes_range(size_t begin, size_t end, dynarr_t **out);

/*
* Function calculates and returns lowest primitive root of a `prime` number.
* If prime has no primitive roots, zero value will be returned.
*/
size_t get_lowest_primitive_root(size_t prime);

/*
* Function calculates and returns all primitive roots of the `prime`,
* given lowest primitive root calculated earlier.
* Result will be collected in preallocated vector `out` in ascending order.
*/
void get_primitive_roots(size_t prime, size_t lowest_root, dynarr_t **out);

/*
* Returns middle element from vector of primitive roots.
*/
size_t get_medium_range_proot(dynarr_t *proots);

/*
* Function that calculates medium range primitive root for the `prime`,
* performing all necessary steps from start to end.
* If prime has no primitive roots, zero value will be returned.
*/
size_t calc_medium_range_proot(size_t prime);

/*
* Factory function that creates vector of size_t pairs.
*/
dynarr_t *create_pair_array(void);

/*
* Function calculates table of prime numbers and their medium range primitive roots,
* storing result in preallocated vector `out`.
*/
void calc_PMPR_table(size_t begin, size_t end, dynarr_t **out);

/*
* Function generates PMPR c header file that contais all necessary definitions
* so it can be included and used by other programs.
*/
void gen_PMPR_c_header(size_t begin, size_t end, const char *filename);


#endif/*_PRIMES_H_*/

#include "primes.h"

#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <assert.h>

#define FILENAME_PREFIX "primes.dat"
#define MAX_FILENAME_SIZE 19
#define MAX_HEADER_NAME_SIZE 32
#define FILE_CAPACITY (MAX_FILE_SIZE)
#define MAX_CONTENTS_LINE_SIZE 64

static void open_cache(cache_t *cache);
static void change_file(cache_t *cache, size_t file_idx);
static void open_page(cache_t *cache, size_t offset);
static void close_cache(cache_t *cache);
static cache_value_t check_prime(cache_t *cache, size_t number);
static void set_prime(cache_t *cache, size_t prime, cache_value_t value);

static void find_prime_factors(size_t prime, dynarr_t **out);
static void check_factor(size_t factor, size_t *number, dynarr_t **out);

/*
* Modular exponentiation algorithm.
* Allows calculations of modular division for very large numbers
* that can be only expressed in form of x^y.
*/
static ssize_t modpow(size_t base, size_t exp, size_t mod);

/*
* Euclidean recursive greatest commont divider algorithm.
*/
static size_t gcd(size_t a, size_t b);

/*
* Used for binary sorted insertion algorithm.
*/
static ssize_t cmp_asc(const void *value, const void *element, void *param);


static void generate_include_guard(const char *filename, char *out);
static void generate_table_contents(dynarr_t *table, dynarr_t **contents);

/*
* Page size of the system.
*/
static size_t s_page_size;

/*
* global prime Cache.
*/
static cache_t s_cache = {};


bool is_prime(size_t number)
{
    if (number == 1) return false;

    for (size_t factor = 2; (factor != number && factor <= 3); ++factor)
    {
        if (number % factor == 0)
        {
            return false;
        }
    }

    for (size_t k = 0, p = 4; p <= sqrt(number); ++k, p = 6 * k)
    {
        if (number % (p + 1) == 0
         || number % (p + 5) == 0)
        {
            return false;
        }
    }

    return true;
}


void init_cache(void)
{
    s_page_size = sysconf(_SC_PAGESIZE);

    open_cache(&s_cache);
}


void fini_cache(void)
{
    close_cache(&s_cache);
}


bool is_prime_cached(size_t number)
{
    if (number == 2) return true;
    if (number % 2 == 0) return false;
    switch (check_prime(&s_cache, number))
    {
        case UNDEFINED: {
            bool prime = is_prime(number);
            set_prime(&s_cache, number, prime ? PRIME : NOT_PRIME);
            return prime;
        }
        case PRIME:     return true;
        case NOT_PRIME: return false;
        default:        exit(EXIT_FAILURE);
    }
}


dynarr_t *create_primes_array(void)
{
    dynarr_t *array = dynarr_create(.element_size = sizeof(size_t));
    return array;
}


void get_primes_range(size_t begin, size_t end, dynarr_t **out)
{
    // dynarr_clear(*out);
    for (; begin <= end; ++begin)
    {
        if (is_prime_cached(begin))
        {
            dynarr_append(out, &begin);
        }
    }
}


size_t get_lowest_primitive_root(size_t prime)
{
    dynarr_t *factors = create_pair_array();
    find_prime_factors(prime, &factors);

    /* reduce */
    size_t total_count = 0;
    for (size_t i = 0; i < dynarr_size(factors); ++i)
    {
        total_count += ((pair_t*) dynarr_get(factors, i))->second;
    }

    size_t unique_count = dynarr_size(factors);
    size_t lowest_factor = 0; /* zero means no primitive roots */

    for (size_t k = 0; k < total_count; ++k)
    {
        size_t test = k + 2;
        bool check = true;
        for (size_t power = 0; check && power < unique_count; ++power)
        {
            size_t factor = ((pair_t*) dynarr_get(factors, power))->first;
            check &= (1 != modpow(test, (prime - 1) / factor, prime));
        }
        if (check)
        {
            lowest_factor = test;
            break;
        }
    }

    dynarr_destroy(factors);
    return lowest_factor;
}


void get_primitive_roots(size_t prime, size_t lowest_root, dynarr_t **out)
{
    size_t power = 2;

    dynarr_clear(*out);
    dynarr_append(out, &lowest_root);

    for (; power < prime; ++power)
    {
        if (1 == gcd(power, prime - 1))
        {
            size_t root = modpow(lowest_root, power, prime);
            dynarr_binary_insert(out, &root, cmp_asc, NULL, NULL);
        }
    }
}


size_t get_medium_range_proot(dynarr_t *proots)
{
    return *(size_t*)dynarr_get(proots, dynarr_size(proots) / 2);
}


size_t calc_medium_range_proot(size_t prime)
{
    size_t lowest_root = get_lowest_primitive_root(prime);
    if (!lowest_root)  return 0; /* zero means no primitive roots */

    dynarr_t *proots = create_primes_array();

    get_primitive_roots(prime, lowest_root, &proots);

    size_t mid_proot = get_medium_range_proot(proots);

    dynarr_destroy(proots);
    return mid_proot;
}


dynarr_t *create_pair_array(void)
{
    dynarr_t *array = dynarr_create(.element_size = sizeof(pair_t));
    return array;
}


void calc_PMPR_table(size_t begin, size_t end, dynarr_t **out)
{
    dynarr_t *primes = create_primes_array();
    get_primes_range(begin, end, &primes);

    size_t size = dynarr_size(primes);
    for (size_t i = 0; i < size; ++i)
    {
        size_t* prime = (size_t*) dynarr_get(primes, i);
        size_t mid_proot = calc_medium_range_proot(*prime);

        if (!mid_proot) continue; /* skip prime with no primitive roots */

        pair_t *pair = &(pair_t){
            .first = *prime,
            .second = mid_proot
        };

        dynarr_append(out, pair);
    }
    dynarr_destroy(primes);
}


void gen_PMPR_c_header(size_t begin, size_t end, const char *filename)
{
    if (filename == NULL) filename = "pmpr.h";
    assert(MAX_HEADER_NAME_SIZE >= strlen(filename));

    dynarr_t *PMPR_table = create_pair_array();
    calc_PMPR_table(begin, end, &PMPR_table);

    FILE *file = fopen(filename, "w");

    char include_guard[MAX_HEADER_NAME_SIZE + 2];
    generate_include_guard(filename, include_guard);

    dynarr_t *contents = dynarr_create(.element_size = 1);
    generate_table_contents(PMPR_table, &contents);

    fprintf(file,
        "#ifndef %s\n"
        "#define %s\n"
        "%s\n" /* pair definition */
        "pair_t pmpr_table[] = {\n"
        "%s" /* table contents */
        "};\n"
        "#endif/*%s*/",
        include_guard,
        include_guard,
        pair_definition,
        (char*)dynarr_first(contents),
        include_guard
        );
    
    fclose(file);
    dynarr_destroy(contents);
    dynarr_destroy(PMPR_table);
}


static ssize_t cmp_asc(const void *value, const void *element, void *param)
{
    (void)param;
    return (ssize_t)*(size_t*)value - *(size_t*)element;
}


static size_t gcd(size_t a, size_t b)
{
    size_t temp;
    while (b != 0)
    {
        temp = a % b;

        a = b;
        b = temp;
    }
    return a;
}


static ssize_t modpow(size_t base, size_t exp, size_t mod)
{
    ssize_t product;
    ssize_t pseq;
    product = 1;
    pseq = base % mod;
    while (exp > 0)
    {
        if (exp & 1)
        {
            product = (product * pseq) % mod;
        }
        pseq = (pseq * pseq) % mod;
        exp >>= 1;
    }
    return product;
}


static void open_cache(cache_t *cache)
{
    change_file(cache, 0);

    if (0 == lseek(cache->fd, 0, SEEK_END)) /* just created */
    {
        if (-1 == ftruncate(cache->fd, s_page_size)) exit(EXIT_FAILURE);
    }

    cache->page_offset = -1ul;
}


static void open_page(cache_t *cache, size_t offset)
{
    size_t page_offset = (offset - offset % s_page_size);

    /* already opened */
    if (cache->page_offset == page_offset) return;

    size_t file_offset = page_offset % FILE_CAPACITY;
    size_t file_idx = page_offset / FILE_CAPACITY;

    if (cache->file_idx != file_idx)
    {
        change_file(cache, file_idx);
    }

    /* need to extend file */
    if ((size_t)lseek(cache->fd, SEEK_END, 0) <= file_offset)
    {
        if (-1 == ftruncate(cache->fd, file_offset + s_page_size))
        {
            exit(EXIT_FAILURE);
        }
    }

    char *page = (char*) mmap(NULL,
        s_page_size,
        PROT_READ|PROT_WRITE,
        MAP_SHARED,
        cache->fd,
        file_offset
    );

    if (MAP_FAILED == page)
    {
        close(cache->fd);
        exit(EXIT_FAILURE);
    }

    if (cache->page) munmap(cache->page, s_page_size);

    cache->page = page;
    cache->page_offset = page_offset;
}


static void change_file(cache_t *cache, size_t file_idx)
{
    char filename[MAX_FILENAME_SIZE];
    if (MAX_FILENAME_SIZE < sprintf(filename, "%s.%zu", FILENAME_PREFIX, file_idx))
    {
        exit(EXIT_FAILURE);
    }

    int fd = open(filename, O_RDWR|O_CREAT, S_IRUSR|S_IWUSR);
    if (-1 == fd)
    {
        exit(EXIT_FAILURE);
    }

    if (cache->fd)
    {
        close(cache->fd);
    }

    cache->file_idx = file_idx;
    cache->fd = fd;
}


static void close_cache(cache_t *cache)
{
    close(cache->fd);
    if (cache->page) munmap(cache->page, s_page_size);
}


static cache_value_t check_prime(cache_t *cache, size_t number)
{
    size_t odd_idx = number / 2;
    size_t byte_offset = odd_idx / 4;
    size_t in_page_offset = byte_offset % s_page_size;
    size_t bit = (odd_idx % 4) * 2;

    open_page(cache, byte_offset);
    return (cache->page[in_page_offset] >> bit) & VALUES_TOTAL;
}


static void set_prime(cache_t *cache, size_t prime, cache_value_t value)
{
    size_t odd_idx = prime / 2;
    size_t byte_offset = odd_idx / 4;
    size_t in_page_offset = byte_offset % s_page_size;
    size_t bit = (odd_idx % 4) * 2;

    open_page(cache, byte_offset);
    cache->page[in_page_offset] |= (value << bit);
}


static void find_prime_factors(size_t prime, dynarr_t **out)
{
    size_t number = prime - 1;
    dynarr_clear(*out);

    while (number > 1)
    {
        for (size_t factor = 2; factor <= 3; ++factor)
        {
            check_factor(factor, &number, out);
        }

        for (size_t k = 0, p = 4; p < number; ++k, p = 6 * k)
        {
            if (number % (p + 1) == 0)
            {
                check_factor(p + 1, &number, out);
            }

            if (number % (p + 5) == 0)
            {
                check_factor(p + 5, &number, out);
            }
        }
    }
}


static void check_factor(size_t factor, size_t *number, dynarr_t **out)
{
    size_t power = 0;

    while (0 == *number % factor)
    {
        ++power;
        *number /= factor;
    }

    if (power >= 1)
    {
        pair_t *prime_factor = &(pair_t){
            .first = factor,
            .second = power
        };
        dynarr_append(out, prime_factor);
    }
}


static void generate_include_guard(const char *filename, char *out)
{
    size_t len = strlen(filename) + 1;

    out[0] = '_';
    strncpy(&out[1], filename, len);
    
    size_t i = 1;
    for (; i < len; ++i)
    {
        if (isalpha(out[i])) out[i] = toupper(out[i]);
        if (out[i] == '.') out[i] = '_';
    }

    out[i] = '_';
    out[i + 1] = '\0';
}


static void generate_table_contents(dynarr_t *table, dynarr_t **contents)
{
    size_t table_size = dynarr_size(table);
    for (size_t i = 0; i < table_size; ++i)
    {
        pair_t *pair = (pair_t*) dynarr_get(table, i);
        char buf[MAX_CONTENTS_LINE_SIZE] = {0};
        size_t len = sprintf(buf, "    {%zu, %zu},\n", pair->first, pair->second);
        for (size_t c = 0; c < len; ++c)
        {
            dynarr_append(contents, &buf[c]);
        }
    }
    dynarr_append(contents, TMP_REF(char, '\0'));
} 

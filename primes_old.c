#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/types.h>

#include "vector.h"

#define PRIME_TABLE_INITIAL_CAP 262144
#define MAX_PRIME_FACTORS 16
#define MAX_PRIMITIVE_ROOTS 4096

typedef struct
{
    size_t primes_num;
    size_t table[];
}
ptable_mapping_t;

typedef struct
{
    size_t cap;
    ptable_mapping_t *mapping;
}
ptable_t;

typedef struct
{
    size_t prime;
    size_t total_count;
    vector_t *factors; /* vector.size => unique_count */
}
pfactors_t;

static bool pfactors_create(pfactors_t *pfactors);
static void pfactors_destroy(pfactors_t *pfactors);

static void find_prime_factors(pfactors_t *pfactors);

static bool open_ptable(ptable_t *ptable, size_t initial_cap)
{
    /* create if not exist */
    int fd = open("primes.dat", O_RDWR|O_CREAT, S_IRUSR|S_IWUSR);
    off_t size = lseek(fd, 0, SEEK_END);
    if (size == 0)
    {
        size = (initial_cap + 1) * sizeof(size_t);
        ftruncate(fd, size);
    }
    lseek(fd, 0, SEEK_SET);
    ptable->mapping = (ptable_mapping_t*) mmap(NULL, size, PROT_READ|PROT_WRITE, MAP_SHARED, fd, 0);
    close(fd);
    if (MAP_FAILED == ptable->mapping)
    {
        return false;
    }

    ptable->cap = size - sizeof(size_t);
    return true;
}

static bool expand_ptable(ptable_t *ptable)
{
    int fd = open("primes.dat", O_RDWR);
    off_t offset = lseek(fd, 0, SEEK_END);
    size_t part_size = offset - sizeof(size_t);

    ftruncate(fd, offset + part_size);
 
    void *extension = mremap(ptable->mapping, offset, offset + part_size, 0);
    close(fd);

    if (MAP_FAILED == extension)
    {
        return false;
    }

    ptable->mapping = extension;
    ptable->cap += part_size;

    return true;
}

static void add_to_ptable(ptable_t *ptable, size_t prime)
{
    if (ptable->mapping->primes_num * sizeof(size_t) >= ptable->cap)
    {
        if (!expand_ptable(ptable))
        {
            exit(1);
        }
    }

    ptable->mapping->table[ ptable->mapping->primes_num++ ] = prime;
}

static void close_ptable(ptable_t *ptable)
{
    munmap(ptable->mapping, ptable->cap + sizeof(size_t));
}

bool ptable_find(size_t *prime_table, size_t number, size_t start, size_t end)
{
    if (start > end)
    {
        return false;
    }

    size_t middle = (start + end) / 2;
    if (prime_table[middle] < number)
    {
        return ptable_find(prime_table, number, middle + 1, end);
    }

    if (prime_table[middle] > number)
    {
        return ptable_find(prime_table, number, start, middle - 1);
    }

    return true; /* prime_table[middle] == number */
}

bool is_prime(ptable_t *ptable, size_t number)
{
    if (number == 2)
    {
        return true;
    }

    if (number == 1 || number % 2 == 0 || number % 3 == 0)
    {
        return false;
    }

    if (ptable_find(ptable->mapping->table, number, 0, ptable->mapping->primes_num))
    {
        return true;
    }

    for (size_t k = 0, p = 0; p <= sqrt(number); ++k)
    {
        p = 6 * k;

        if (p != 0 && number % (p + 1) == 0)
            return false;

        if (number % (p + 5) == 0)
            return false;
    }

    add_to_ptable(ptable, number);

    return true;
}

size_t find_prime(ptable_t *ptable, size_t minimum, size_t maximum)
{
    for (size_t number = minimum; number <= maximum; ++number)
    {
        if (is_prime(ptable, number))
        {
            return number;
        }
    }

    return 0;
}

bool is_primitive_root(size_t prime, size_t candidate)
{
    return false;
}

/* size_t gcd(size_t x, size_t y) */
/* { */
/*     if (x == 0 || y == 0)  return x + y; */
/*     if (x == y)  return x; */
/*     if (x > y)   return gcd(x-y, y); */
/*     return gcd(x, y-x); */
/* } */


size_t gcd(size_t a, size_t b)
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

ssize_t modpow(size_t base, size_t exp, size_t mod)
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

static size_t find_lowest_primitive_root(size_t prime)
{
    size_t unique_count = 0;
    size_t total_count = 0;
    pfactors_t pfactors = {};
    pfactors_create(&pfactors);
    pfactors.prime = prime;
    
    find_prime_factors(&pfactors);
    unique_count = vector_size(pfactors.factors);
    total_count = pfactors.total_count;
    
    /* for (size_t i = 0; i < unique_count; ++i) */
    /* { */
    /*     printf(" %zu,", *(size_t*)vector_get(pfactors.factors, i)); */
    /* } */
    /* printf("\n"); */

    size_t lowest_factor;
    for (size_t k = 0; k < total_count; ++k)
    {
        size_t test = k + 2;
        bool check = true;
        for (size_t i = 0; check && i < unique_count; ++i)
        {
            check &= (1 != modpow(test, (prime - 1) / *(size_t*)vector_get(pfactors.factors, i), prime));
        }

        if (check)
        {
            lowest_factor = test; /* lowest primitive root found */
        }
    }

    pfactors_destroy(&pfactors);
    return lowest_factor;
}

static void check_factor(pfactors_t *pfactors, size_t factor, size_t *number)
{
    size_t power = 0;
    while (0 == *number % factor)
    {
        ++power;
        *number /= factor;
    }

    if (power >= 1)
    {
        vector_append_back(&pfactors->factors, &factor);
        pfactors->total_count += power;
    }
}

static void find_prime_factors(pfactors_t *pfactors)
{
    size_t number = pfactors->prime - 1;

    while (number > 1)
    {
        for (size_t factor = 2; factor <= 3; ++factor)
        {
            check_factor(pfactors, factor, &number);
        }

        for (size_t k = 0, p = 0; p <= number; ++k)
        {
            p = 6 * k;

            if (p != 0 && number % (p + 1) == 0)
            {
                check_factor(pfactors, p + 1, &number);
            }

            if (number % (p + 5) == 0)
            {
                check_factor(pfactors, p + 5, &number);
            }
        }
    }
}

static size_t find_all_proots(size_t prime, size_t lowest_proot, size_t proots[MAX_PRIMITIVE_ROOTS])
{
    size_t count = 1;
    size_t power = 2;
    
    proots[0] = lowest_proot;

    for (; power < prime; ++power)
    {
        if (gcd(power, prime - 1) == 1)
        {
            proots[count++] = modpow(lowest_proot, power, prime);
        }
    }

    return count;
}

static int compare(const void *a, const void *b)
{
    return *(size_t*)a - *(size_t*)b;
}

static bool pfactors_create(pfactors_t *pfactors)
{
    vector_create(pfactors->factors, .esize = sizeof(size_t));
    if (pfactors->factors) return true;
    return false;
}

static void pfactors_destroy(pfactors_t *pfactors)
{
    vector_destroy(pfactors->factors);
}

int main(void)
{
    size_t proots[MAX_PRIMITIVE_ROOTS];
    size_t prime = 761;//9214703;
    size_t lowest_proot = find_lowest_primitive_root(prime);
    size_t proots_num = find_all_proots(prime, lowest_proot, proots);

    qsort(proots, proots_num, sizeof(size_t), compare);
    
    size_t medium_proot = proots[proots_num/2];
    /*  */
    /* for (size_t i = 0; i < proots_num; ++i) */
    /* { */
    /*     printf("%zu, ", proots[i]); */
    /* } */
    
    /* printf("medium_proot: %zu\n", medium_proot); */

    size_t state = medium_proot;
    for (size_t n = 0; n < prime - 1; ++n)
    {
        printf("%zu\n", (state=(medium_proot * state) % prime));
    }

    return 0;
}

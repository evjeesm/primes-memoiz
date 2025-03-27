#include "primes.h"
#include "vector.h"
#include <stdio.h>

int main(void)
{
    init_cache();
#if 0
    printf("is %zu a prime? %d\n", 1, is_prime(1));
    printf("is %zu a prime? %d\n", 2, is_prime(2));
    printf("is %zu a prime? %d\n", 3, is_prime(3));
    printf("is %zu a prime? %d\n", 5, is_prime(5));
    printf("is %zu a prime? %d\n", 6, is_prime(6));
    printf("is %zu a prime? %d\n", 7, is_prime(7));
    printf("is %zu a prime? %d\n", 8, is_prime(8));
    printf("is %zu a prime? %d\n", 9, is_prime(9));
    printf("is %zu a prime? %d\n", 13, is_prime(13));
    printf("is %zu a prime? %d\n", 53, is_prime(53));
    printf("is %zu a prime? %d\n", 761, is_prime(761));
#endif

#if 0
    printf("is %zu a prime? %d\n", 761, is_prime_cached(761));
    printf("is %zu a prime? %d\n", 4294967296, is_prime_cached(4294967296));
    printf("is %zu a prime? %d\n", 17179869184, is_prime_cached(17179869184));
    printf("is %zu a prime? %d\n", 17179869185, is_prime_cached(17179869185));
    printf("is %zu a prime? %d\n", 4398046507007, is_prime_cached(4398046507007));
    printf("is %zu a prime? %d\n", 4398046507009, is_prime_cached(4398046507009));
    printf("is %zu a prime? %d\n", 348394224319, is_prime_cached(348394224319));
    printf("is %zu a prime? %d\n", 1471800996949, is_prime_cached(1471800996949));
    printf("is %zu a prime? %d\n", 14084568907963, is_prime_cached(14084568907963));
    printf("is %zu a prime? %d\n", 16289967626711, is_prime_cached(16289967626711));
    printf("is %zu a prime? %d\n", 17592186044417, is_prime_cached(17592186044417));
    printf("is %zu a prime? %d\n", 5587434143095588297, is_prime_cached(5587434143095588297));
    printf("is %zu a prime? %d\n", 11580786069677071957ul, is_prime_cached(11580786069677071957ul));

#endif

#if 1
    dynarr_t *primes = create_primes_array();
    get_primes_range(   900000000000,    900000000100, &primes);
    get_primes_range( 21000000000101,  21000000000150, &primes);
    get_primes_range(110000000001001, 110000000001101, &primes);
    get_primes_range(300000000001001, 300000000001101, &primes);
    get_primes_range(610000000001001, 610000000001101, &primes);
    for (size_t i = 0; i < dynarr_size(primes); ++i)
    {
        printf("%zu, ", *(size_t*) dynarr_get(primes, i));
    }
    puts("\n");

    vector_destroy(primes);
#endif

#if 0

    size_t root = get_lowest_primitive_root(761);
    printf("lowest primitive root of 761 is %zu\n", root);

    vector_t *roots = create_primes_vector();
    get_primitive_roots(761, root, &roots);

    for (int i = 0; i < vector_size(roots); ++i)
    {
        printf("%zu, ", *(size_t*) vector_get(roots, i));
    }
    puts("\n");

    size_t mproot = get_medium_range_proot(roots);
    printf("medium range primitive root for 761 is %zu\n", mproot);

    vector_destroy(roots);
#endif

#if 0
    
    vector_t *PMPR_table = create_pair_vector();
    calc_PMPR_table(761, 931, &PMPR_table);

    for (int i = 0; i < vector_size(PMPR_table); ++i)
    {
        pair_t *pair = vector_get(PMPR_table, i);
        printf("%zu - %zu\n", pair->first, pair->second);
    }

    vector_destroy(PMPR_table);

#endif

#if 0
    gen_PMPR_c_header(101, 1023, "test_pmpr_1.h");
#endif

    fini_cache();
    return 0;
}

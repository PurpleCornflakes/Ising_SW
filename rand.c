#undef  NDEBUG                       /* #define will turn off assertion code */
#define NR_RANDOM_NUMBER              /* #define use "Numerical Recipes" ran2 */

#ifdef NR_RANDOM_NUMBER   /* then use "Numerical Recipes" random number ran2 */

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-13
#define RNMX (1.0-EPS)

static long int idum = 1;

double drand64(void)               /* this is "Numerical Recipes" ran2(idum) */
{ 
    int j;
    long k;
    static long idum2=123456789;
    static long iy=0;
    static long iv[NTAB];
    double temp;
    
    if( idum <= 0 )
    {
        if( ((-1)*idum ) < 1 )
            idum = 1;
        else
            idum  = (-1)*idum;
        idum2 = idum;
        for( j = NTAB+7 ; j>=0 ; j-- )
        {
            k = idum/IQ1;
            idum = IA1 * (idum-k*IQ1)- k*IR1;
            if( idum < 0 )
                idum += IM1;
            if( j< NTAB )
                iv[ j ] = idum;
        }
        iy = iv[0];
    }
    
    k = idum/IQ1;
    idum = IA1 * (idum-k*IQ1)- k*IR1;
    if( idum < 0 )
        idum += IM1;
    k = idum2/IQ2;
    idum2 = IA2 * (idum2-k*IQ2)- k*IR2;
    if( idum2 < 0 )
        idum2 += IM2;
    j = iy/NDIV;
    iy = iv[ j ]-idum2;
    iv[ j ] = idum;
    if( iy < 1 )
        iy += IMM1;
    if( (temp = AM*iy ) > RNMX )
        return RNMX;
    else
        return temp; 
}  

void srand64(unsigned long seed)
{
   // assert(sizeof(long int) == 4);
   assert(sizeof(long long int) == 8);
   assert(sizeof(long double) > sizeof(double));
   idum = -seed;
   if(idum > 0) idum = - idum;
   drand64();                       /* call with negative idum to initialize */
   // printf("NR ran2, initial seed x = %d\n", (int) seed);
}

#endif
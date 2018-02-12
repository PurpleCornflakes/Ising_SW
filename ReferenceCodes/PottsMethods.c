/* Do Metropolis/Swendsen-Wang/Sweeny/Gliozzi/Wolff dynamics

   D dimensions, general Q Potts model.

   Compute decorrelation time tau from the definition
   var(A_M) = tau var(A_1)/M, where A_M = (1/M) sum( A_i), A is energy

   Speed: 0.9 (700MHz Pentium 3), 0.6 (600MHz alpha), 0.4 (1.5G Athlon) 
          * 10^-6 sec per SW spin flip

   14 June 2002
   by Wang Jian-Sheng, wangjs@cz3.nus.edu.sg
*/

#undef  NDEBUG                       /* #define will turn off assertion code */
#undef MERSENNE_TWISTER                    /* random number, Mersenne or lcg */
#define NR_RANDOM_NUMBER              /* #define use "Numerical Recipes" ran2 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <assert.h>
                                                                   /* macros */
#define  ZMAX  10                                    /* max Z, good upto D=5 */
#define  LOOP  65536                                   /* number of ``runs'' */
#define  N_WRITE  1                                    /* how often to write */
#define  TK  16                                         /* buffer array size */
#define  MMAX  65536                        /* max M, must be 2^TK or larger */
#ifdef MERSENNE_TWISTER
#define srand64 sgenrand
#define drand64 genrand
#endif
                                                                  /* globals */
int Algo;                /* 0 metropolis, 1 SW, 2 sweeny, 3 gliozzi, 4 wolff */
int D;                                                          /* dimension */
int L;                                                /* lattice linear size */
int N;                                              /* total number of spins */
int Z;                                                  /* 2*D, coordination */
double Q;                                          /* number of Potts states */
int MCTOT;                                       /* total Monte Carlo sweeps */
int MCDIS;                              /* sweeps discarded in the beginning */
double Tr, P;                                                   /* Tc and Pc */
int *s;                                              /* Potts spin, 0 to Q-1 */
int *bond;                             /* bond[D*i+d] at site i, direction d */
int *list;                                       /* site i has label list[i] */
int Csize;                                                   /* Cluster size */
long double buff[MMAX];                        /* circular buffer for energy */
long double sum[TK+1];                                            /* sum E_i */
long double ave[TK+1];                                        /* <sum E_i/M> */
long double var[TK+1];                                    /* <(sum E_i/M)^2> */
                                                      /* function prototypes */
void init(void);
void sweeny(void);
int connect(int i, int j, int *bond);
void sw(void);
void metrop(int step);
void flip(int i, int s0, int q);
void wolff(int n);
int e_change(int i, int q);
long double energy(void);
void variance(long double e);
void neighbor(int i, int nn[ ]);
double drand64(void); 
void srand64(unsigned long seed);

int main(int argc, char **argv)
{
   long double u, C_av, C_av_all, E_av, E_av_all, e, eng, eng2, engk, ckT;
   long double tau[TK+1], err[TK+1], err1;
   long double ave_all[TK+1], var_all[TK+1]; 
   int mc, k, d, dt, loop, fac = 1;
   FILE *fin, *fout;
 
   fin = stdin;
   fout = stdout;
   if(argc >= 2) {
      fin = fopen(argv[1], "r");
      assert(fin != NULL);
   }

   fprintf(stderr, "enter Algo (0 Metropolis, 1 SW, 2 Sweeny, 3 Gliozzi, 4 Wolff):\n");
   fscanf(fin, "%d", &Algo);
   fprintf(stderr, "enter dimension D:\n");
   fscanf(fin, "%d", &D);
   fprintf(stderr, "enter Potts states Q:\n");
   fscanf(fin, "%lf", &Q);
   fprintf(stderr, "enter size L:\n");
   fscanf(fin, "%d", &L);
   fprintf(stderr, "enter MCTOT:\n");
   fscanf(fin, "%d", &MCTOT);
   fprintf(stderr, "enter MCDIS:\n");
   fscanf(fin, "%d", &MCDIS);
   if(argc >= 2) {
      fclose(fin);
   }
   assert(Q>0.0 && D>0 && L>0);
   Z = 2 * D;
   assert(ZMAX >= Z);               /* enlarge ZMAX if this assertion failed */
   N = 1;
   for(d = 0; d < D; ++d) {              /* N = L^D is total number of sites */
      N *= L;
   }
   
   loop = 0;
   for(k = 0; k <= TK; ++k) {
      tau[k] = 0.0;
      err[k] = 0.0;
      ave_all[k] = 0.0;
      var_all[k] = 0.0;
   }
   C_av = 0.0;
   C_av_all = 0.0;
   E_av = 0.0;
   E_av_all = 0.0;

   init();                 /* initialize s[], buff[], sum[], var[], A_ex etc */

Loop:

   for(mc = 1; mc <= MCTOT; ++mc) {
      if(Algo == 2 || Algo == 3) {
         sweeny();
      } else if (Algo == 1) {
         sw();
      } else if (Algo == 4) {
         Csize = 0;
         wolff(1);
         C_av += Csize;                              /* average cluster size */
      } else {
         metrop(N);
      }

      e = energy();
      variance(e);
 
      if(mc == MCDIS) {
         for(k = 0; k <= TK; ++k) {
            ave[k] = 0.0;
            var[k] = 0.0;                     /* discard accumulated numbers */
         }                                      /* for mc = 1, 2, ..., MCDIS */
         E_av = 0.0;
         C_av = 0.0;
      }
      if(mc > MCDIS) {                                      /* do statistics */
         E_av += e;
      }
   }
  
   ++loop;

   dt = 1;
   eng = E_av/(MCTOT-MCDIS);
   for(k = 0; k <= TK; ++k) {
      engk = ave[k]/(MCTOT-MCDIS);
      u = dt * (var[k]/(MCTOT-MCDIS)-engk*engk)/(var[0]/(MCTOT-MCDIS)-eng*eng);
      tau[k] += u;                            /* binned value for error only */
      err[k] += u*u;                            /* err[] is actually <tau^2> */
      dt *= 2;
      ave_all[k] += ave[k];                          /* sum up the loop runs */
      var_all[k] += var[k];      
   }
   E_av_all += E_av;
   C_av_all += C_av;

   if(loop % (fac*N_WRITE) == 0) {
      fac *= 2;

      if(argc == 3) {                /* open file to append only when write */
         fout = fopen(argv[2], "a");
         assert(fout != NULL);
      }
      switch (Algo) {
         case 0:
            fprintf(fout, "Single spin flip - ");
            break;
         case 1:
            fprintf(fout, "SW dynamics - ");
            break;
         case 2:
            fprintf(fout, "Sweeny-heat-bath dynamics - ");
            break;
         case 3:
            fprintf(fout, "Gliozzi dynamics - ");
            break;
         case 4:
            fprintf(fout, "Wolff dynamics - ");
            break;
         default:
            fprintf(stderr, "unknown algo!\n");
            exit(1);
      }
      fprintf(fout, "%dD Q=%g Potts, L=%d, %d MCS used, %d discard, %d runs\n", 
              D, Q, L, MCTOT-MCDIS, MCDIS, loop);
      if(Algo == 4) {
         fprintf(fout, "<C> = %Lg\n", C_av_all/(MCTOT-MCDIS)/loop);
      }
      eng = E_av_all/(MCTOT-MCDIS)/loop;
      if(Algo == 2 || Algo == 3) {
        ckT = ( (var_all[0]/(MCTOT-MCDIS)/loop-eng*eng) + eng*(1-P)/P)/N;
      } else {
        ckT = (var_all[0]/(MCTOT-MCDIS)/loop-eng*eng)/N;
      }
      fprintf(fout, "<E> = %Lg   cT^2 = %Lg\n", eng/N, ckT); 
      fprintf(fout, "    M      tau       err           <e_M^2>            <e_M>\n");
      dt = 1;
      for(k = 0; k <= TK; ++k) {
         engk = ave_all[k]/(MCTOT-MCDIS)/loop;
         eng2 = var_all[k]/(MCTOT-MCDIS)/loop;
         u = dt * (eng2-engk*engk) / (var_all[0]/(MCTOT-MCDIS)/loop-eng*eng);
         fprintf(fout, "%5d  %Lg", dt, u);
         if (loop > 1) { 
            err1 = sqrtl((err[k]/loop-(tau[k]/loop)*(tau[k]/loop))/(loop-1));
            fprintf(fout, " +/- %Lg   ", err1);
         } else {
            fprintf(fout, " +/- ?   ");
         } 
         fprintf(fout, "   %.15Lg    %.15Lg\n", eng2/N/N, engk/N);
         dt *= 2;
      }   
      fflush(fout);
      if(argc == 3) {
         fclose(fout);
      }
   }
   
   if(loop < LOOP) goto Loop;

   return 0;
}

                
void sweeny(void)                                /* Sweeny/Gliozzi algorithm */
{
   int mc, i, j, k, d;
   int nn[ZMAX];
   double r;

   for(mc = 0; mc < D*N; ++mc) {            /* go over each bonds on average */

      k = D * N * drand64();                        /* pick a bond at random */

      i = k / D;                                              /* site number */
      d = k % D;                                                /* direction */
      assert(i >= 0 && i < N);
      neighbor(i, nn);
      j = nn[2*d];
      
      r = drand64();
      if(Algo == 2) {                               /* Sweeny - PRB 27, 4445 */
         bond[k] = 0;
         if (connect(i,j, bond)) {
            bond[k] = r < P;
         } else {
            bond[k] =  r < P/( (1-P)*Q+P ); 
         }
      } else if (Algo == 3) {                   /* Gliozzi, cond-mat/0201285 */
         if( bond[k] || connect(i,j,bond) ) {
            bond[k] = r < P;
         } else {
            bond[k] = r < P/Q;
         }
      }
   }
}

long long int *mark;
int *stack;

                           /* Decide if site i and j are on the same cluster */
int connect(int i, int j, int *bond)
{
   int pt, jj, k, find, nn[ZMAX], b[ZMAX];
   static long long int inc = 0;

   if(inc == 0) { 
      for(k = 0; k < N; ++k) {
         assert(mark[k] == 0);
      }
   }

   ++inc;

   stack[0] = i;
   mark[i] = inc;
   pt = 1;                               /* points to the 1st available slot */
   find = 0;
   while(pt > 0) {

      assert(pt < N);
      k = stack[--pt];                                /* pop one from bottom */ 
      neighbor(k, nn);

      for(jj = 0; jj < D; ++jj) {
         b[2*jj] = bond[D*k + jj];
         b[2*jj+1] = bond[D*nn[2*jj+1]+jj];
      }

      for(jj = 0; jj < Z; ++jj) {
         if (b[jj] && mark[nn[jj]] != inc) {
            mark[nn[jj]] = inc;
            if (j == nn[jj]) {
               find = 1;
               goto Stop;
            } else {
               stack[pt++] = nn[jj];
            }
         }
      }
    }      
 
Stop:
    return find;
}

                 
void sw(void)     /* Swendsen-Wang algorithm, perform one Swendsen-Wang step */
{
   int i, ip, j, cnt, inc, a, b, min, max;
   int r, p, q;

   for(i = 0; i < N; ++i) {                      /* clear or initialize list */
      list[i] = i;             /* initially each site is a cluster by itself */
   }       /* later, site i has the same label as site list[i] (recursively) */
      
   for(i = 0; i < N; ++i) {                      /* set the bond with prob P */
      cnt = 0;
      r = i;
      p = 1 - L;
      q = 1;
Repeat:                                                           /* D times */ 
      ip = (r + 1) % L == 0 ? i + p : i + q;
      if(s[i] == s[ip] && drand64() < P) {      /* implement Hoshen-Kopelman */
         a = list[i];
         while (a > list[a]) {             /* run through until a == list[a] */
            a = list[a];
         } 
         b = list[ip];
         while (b > list[b]) {
            b = list[b];
         } 
         if (a > b) {                      /* find min and max of two labels */
            min = b;
            max = a;
         } else {
            min = a;
            max = b;
         }
         list[max] = min;
         list[i] = min;
      }
      ++cnt;
      if(cnt < D) {
         r = r/L;
         p *= L;
         q *= L;
         goto Repeat;
      }
   }

   inc = 0;           /* last sweep to make list pointing to the final label */
   for(i = 0; i < N; ++i) {
      if(i == list[i]) {
         s[i] = inc;                     /* spin value over-written by label */
         ++inc;                                        /* a new cluster find */
      } else {
        j = list[i];
        while (j > list[j]) {
           j = list[j];
        } 
        assert(j < i);
        s[i] = s[j];                                   /* one of old cluster */
      }
   }
   assert(inc <= N);
   for(j = 0; j < inc; ++j) {
      list[j] = Q*drand64();                    /* new spin for each cluster */
   }
   for(i = 0; i < N; ++i) {
      assert(s[i] < inc);
      s[i] = list[s[i]];                       /* old s[i] is cluster number */
   }
}


/* What is the energy change if site i is changed to state q from s[i].
*/
int e_change(int i, int q)
{
   int nn[ZMAX];
   int e, si, sj, j;

   neighbor(i, nn);
   si = s[i];
   e = 0;
   for(j = 0; j < Z; ++j) {
      sj = s[nn[j]];
      e += (si == sj) - (q == sj);
   }

   return e;
}


/* do Metropolis canonical Monte Carlo, temperarure Tr, size N are globals. */
void metrop(int step)
{
   int mc, i, q, de;

   for(mc = 0; mc < step; ++mc) {
      i = N * drand64();
      q = (Q-1) * drand64();
      if(q >= s[i]) {
         ++q;
      }
      de = e_change(i, q);
      if (de <= 0 || drand64() < exp(-de/Tr)) {           /* Metropolis rate */
         s[i] = q;
      }
   }
}

/* This function wolff() performs n cluster flips, equivalent to
   few Monte Carlo steps in standard single-spin-flip algorithms.
   It picks a seed site at random and calls the flip function to generate
   one cluster.  N is total number of spin (a global).  */

void wolff(int n)
{
   int i, k, q;

   assert(n == 1);
   for(k = 0; k < n; ++k) {
      i = drand64() * (double) N;
      q = (Q-1) * drand64();                      /* choose a new spin value */
      if(q >= s[i]) {
         ++q;
      }
      flip(i, s[i], q);
   }
}

/* Perform a Wolff single cluster flip. s[], P, and Z are passed globally.
   The first argument i of flip function is the site to be flipped, the
   second argument is the spin of the cluster before flipping.  q is the
   spin value to flip into.  Recursive version may cause stack overflow. */

void flip(int i, int s0, int q)
{
   int j, nn[ZMAX];

   ++Csize;                                     /* count the size of cluster */

   s[i] = q;                                    /* flip the spin immediately */
   neighbor(i, nn);                            /* find nearest neighbor of i */
   for(j = 0; j < Z; ++j)                       /* flip the neighbor if ...  */
      if(s0 == s[nn[j]] && drand64() < P)
         flip(nn[j], s0, q);
}

/* Neighbor returns in the array nn[ ] the neighbor sites of i.  The sites
are labelled sequentially, starting from 0.  It works for any hypercubic
lattice.  Z (=2*D) is the coordination number, L is linear size, all passed 
from globals. */

void neighbor(int i, int nn[ ])
{
   int j, r, p, q;

   r = i;
   p = 1 - L;
   q = 1;

   for(j = 0; j < Z; j += 2) {
      nn[j] = (r + 1) % L == 0 ? i + p : i + q;
      nn[j+1]     = r % L == 0 ? i - p : i - q;
      r = r/L;
      p *= L;
      q *= L;
   }
}


/* This function computes the energy of the current configuration s[].  
   The energy is defined as E = - sum ( delta (s[i], s[j]),
*/

long double energy(void)
{
   int b, k;
   int i, ip, ie, si, j;
   int r, p, q;

   if(Algo == 2 || Algo == 3) {
      b = 0;
      for(k = 0; k < D*N; ++k) {
         b += bond[k];
      }
      return -(long double) b/P; 
   } else {
      ie = 0;
      for(i = 0; i < N; ++i) {
         si = s[i];
         assert(si >= 0 && si < Q);
         r = i;
         p = 1 - L;
         q = 1;
         for(j = 0; j < D; ++j) {
            ip = (r + 1) % L == 0 ? i + p : i + q;
            ie += (si==s[ip]);
            r = r/L;
            p *= L;
            q *= L;
         }
      }
      return (long double) -ie;
   }
}


void variance(long double e)
{
   int dt, k, t;
   long double del;
   static int tpt = 0;

      /* current value is e, most recent value was at tpt (before increment) */
   dt = 1;                                  /* how many terms in partial sum */
   tpt = (tpt+1) % MMAX; 
   for(k = 0; k <= TK; ++k) {                         /* from dt = 1 to 2^TK */
      assert(dt <= MMAX);
      t = (tpt - dt + MMAX) % MMAX;              /* subtract dt, wrap around */
      sum[k] += (e - buff[t]);        /* sum_{t=t0+1}^{t0+dt} e_t for all t0 */
      dt *= 2; 
   }
   buff[tpt] = e;

   dt = 1;
   for(k = 0; k <= TK; ++k) {                                 /* k is lg2 dt */
      del = sum[k]/dt;                                       /* i.e. delta A */
      ave[k] += del;
      var[k] += del*del;
      dt *= 2;
   }
}

                         /* 2D Ising exact energy and specific heat per spin */
#if L==4
#define  EN_EXACT  -1.56562378763832
#define  CA_EXACT   0.78326682592891
#elif L==8
#define  EN_EXACT  -1.491589107439706
#define  CA_EXACT   1.145559239894408
#elif L==16
#define  EN_EXACT  -1.453064852813477
#define  CA_EXACT   1.498704959400026
#elif L==32
#define  EN_EXACT  -1.433658466146248
#define  CA_EXACT   1.846767590039559
#elif L==64
#define  EN_EXACT  -1.423938389833011
#define  CA_EXACT   2.192211393140571
#elif L==128
#define  EN_EXACT  -1.419076272084985
#define  CA_EXACT   2.536331335108602
#elif L==256
#define  EN_EXACT  -1.416644954196832 
#define  CA_EXACT   2.879786255202499
#elif L==512
#define  EN_EXACT  -1.4154292629050 
#define  CA_EXACT   3.22290795
#elif L==1024
#define  EN_EXACT  -1.41482141321652494
#define  CA_EXACT   3.565862873717
#endif


/* Do initialization of spins, coupling, mat[], etc */
void init(void)
{
   int i, k;

   if(Algo == 1) { 
      list = (int *) malloc(N*sizeof(int));
      assert(list != NULL);
   }
   if(Algo == 2 || Algo == 3) {
      bond = (int *) malloc(D*N*sizeof(int));
      mark = (long long int *) calloc(N, sizeof(long long int));
      stack = (int *) malloc(N*sizeof(int));
      assert(bond != NULL);
      assert(mark != NULL);
      assert(stack != NULL);
   } else {                                            /* for Algo = 0, 1, 4 */
      s = (int *) malloc(N*sizeof(int));
      assert(s != NULL);
   }

   if(D == 2) {
      Tr = (1.0/(log(sqrt((double)Q)+1.0)));
      P = (1.0-1.0/(sqrt((double)Q)+1.0));
   } else if (D == 3 && Q == 2) {
      Tr = 1.0/(0.221657*2.0);                         /* using Potts scaled */
      P = 1.0 - exp(-1.0/Tr);
   } else {
      Tr = D/2.0;
      P = 1.0 - exp(-1.0/Tr);
   }

   srand64(time(NULL));

   if(Algo == 2 || Algo == 3) {
      for(k = 0; k < D*N; ++k) {
         bond[k] = (drand64() < P);
      }
   } else {
      for(i = 0; i < N; ++i) {
         s[i] = Q * drand64();
      }
   }

   for(i = 0; i < MMAX; ++i) {
      buff[i] = 0;
   }
   for(k = 0; k <= TK; ++k) {
      sum[k] = 0.0;
      ave[k] = 0.0;
      var[k] = 0.0;
   }
}


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
   assert(sizeof(long int) == 4);
   assert(sizeof(long long int) == 8);
   assert(sizeof(long double) > sizeof(double));
   idum = -seed;
   if(idum > 0) idum = - idum;
   drand64();                       /* call with negative idum to initialize */
   printf("NR ran2, initial seed x = %d\n", (int) seed);
}

#else
#ifndef MERSENNE_TWISTER
/* 64-bit random number generator */

static unsigned long long int x = 1;

double drand64(void) 
{
   x = 6364136223846793005ll * x + (long long int) 1;
   return (double) x * 5.4210108624275218e-20; 
}

void srand64(unsigned long seed)
{
   assert(sizeof(long long int) == 8);
   x = seed;
   printf("drand 64-bit, initial seed x = %d\n", (int) x);
}
#endif
#endif
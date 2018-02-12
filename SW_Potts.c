#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include "rand.c"
// 2D Ising: z=0.35 with sw, z=2.215 with metropolis
// 3D Ising: z=0.75 with sw, z=2.0 with metropolis        

int D;                                                          /* dimension */
int L;                                                /* lattice linear size */
int N;                                              /* total number of spins */
int Q;                                          /* number of Potts states */
int MCTOT;                                       /* total Monte Carlo sweeps */
int MCDIS;                              /* sweeps discarded in the beginning */
double Tr, P;                                                   /* Tc and Pc */
int *s;                                              /* Potts spin, 0 to Q-1 */
int *bond;                             /* bond[D*i+d] at site i, direction d */
int *list;                                       /* site i has label list[i] */

void sw(void);
void init(void);
int energy(void);

int main(int argc, char **argv)
{
   double E_av;
   int e;
   int mc, d, ii;
   FILE *fin, *fout, *e_file;
 
   fin = stdin;
   fout = stdout;
   if(argc >= 2) {
      fin = fopen(argv[1], "r");
      assert(fin != NULL);
   }
   fprintf(stderr, "enter Potts states Q:\n");
   fscanf(fin, "%d", &Q);
   fprintf(stderr, "\nQ:%d\n",Q);
   fprintf(stderr, "enter dimension D:\n");
   fscanf(fin, "%d", &D);
   fprintf(stderr, "\nD:%d\n",D);
   fprintf(stderr, "enter size L:\n");
   fscanf(fin, "%d", &L);
   fprintf(stderr, "\nL:%d\n",L);
   fprintf(stderr, "enter MCTOT:\n");
   fscanf(fin, "%d", &MCTOT);
   fprintf(stderr, "\nMCTOT:%d\n",MCTOT);
   fprintf(stderr, "enter MCDIS:\n");
   fscanf(fin, "%d", &MCDIS);
   fprintf(stderr, "\nMCDIS:%d\n",MCDIS);
   assert(D>0 && L>0);
   if(argc >= 2) 
      fclose(fin);
   
   if(argc == 3) {                /* "r"read, "a"append, "w"write */
      fout = fopen(argv[2], "w");
      assert(fout != NULL);
   }

   N = 1;
   for(d = 0; d < D; ++d)            /* N = L^D is total number of sites */
      N *= L;

   init();                 /* initialize s[], buff[], sum[], var[], A_ex etc */
   E_av = 0.0;
   e_file = fopen("e_dist","w");
   assert(e_file != NULL);
   for(mc = 1; mc <= MCTOT; ++mc) {
      sw();
      e = energy();
      
      // print s[ii] to file
      if(mc%200 == 0){
         for(ii = 0; ii < N; ++ii)
            fprintf(fout, "%d",s[ii]); 
         fprintf(fout, "\n");
      } 

      //start record energy
      if(mc > MCDIS){                               /* do statistics */
         E_av += e;
         fprintf(e_file, "%d\n", e);
      }
      fprintf(stdout, "mc=%d", mc);
   }
   fclose(e_file);
   if(argc == 3) 
      fclose(fout);

   E_av /= (MCTOT-MCDIS);
   fprintf(stderr, "\nE_av = %g\n", E_av);
   return 0;
}

/* Do initialization of spins, coupling, mat[], etc */
void init(void)
{
   int i, k;

   list = (int *) malloc(N*sizeof(int));
   assert(list != NULL);

   s = (int *) malloc(N*sizeof(int));
   assert(s != NULL);

   if(D == 2) {
      Tr = (1.0/(log(sqrt((double)Q)+1.0))); //2.612 != Tc/2=2.269/2
      P = (1.0-1.0/(sqrt((double)Q)+1.0)); //0.58
   } else if (D == 3 && Q == 2) {
      Tr = 1.0/(0.221657*2.0);             //2.25 = Tc/2=4.5/2            /* using Potts scaled */
      P = 1.0 - exp(-1.0/Tr);
   } else {
      Tr = D/2.0;
      P = 1.0 - exp(-1.0/Tr);
   }

   srand64(time(NULL));
   
   for(i = 0; i < N; ++i) 
      s[i] = Q * drand64();
      
}


int energy(void)
{
   int b, k;
   int i, ip, ie, si, j;
   int r, p, q;


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
   return  -ie;
   
}

void sw(void)     /* Swendsen-Wang algorithm, perform one Swendsen-Wang step */
{
   int i, ip, j, cnt, inc, a, b, min, max;
   int r, p, q;

   for(i = 0; i < N; ++i)                   
      list[i] = i;             /* initially each site is a cluster by itself */  
      
   for(i = 0; i < N; ++i) {                      /* set the bond with prob P */
      r = i;
      p = 1 - L;
      q = 1;
      for(cnt = 0; cnt < D; ++cnt) {
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
         r = r/L;
         p *= L;
         q *= L;      
      }
   }


   inc = 0;           /* last sweep to make list pointing to the final label */
   for(i = 0; i < N; ++i) {
      if(i == list[i]) {
         s[i] = inc;                     /* spin value over-written by label */
         ++inc;                                        /* a new cluster find */
      } 
      else {
        j = list[i];
        while (j > list[j])
           j = list[j];
        assert(j < i);
        s[i] = s[j];                                   /* one of old cluster */
      }
   }

   assert(inc <= N);
   for(j = 0; j < inc; ++j) 
      list[j] = Q*drand64();                    /* new spin for each cluster */
   for(i = 0; i < N; ++i) {
      assert(s[i] < inc);
      s[i] = list[s[i]];                       /* old s[i] is cluster number */
   }

}


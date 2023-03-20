/* genotypecall
   Ziheng Yang, 6 August 2019
   This calculates the genotyping error rate for both homozygotes and heterozygotes as a
   function of read depth and single-read error rate e.  Genotype is inferred using ML based
   on reads for one sample (one individual), and the true single-read error rate e is assumed
   to be known and given.  Error rates are assumed to be the same among reads, independent of
   the true nucleotides.
   The model is one of simple binomial sampling.  See equations (4.12), (4.13) and (4.14) in
   the mathematical notes of Li (2011), for a more complex model.

   Li, H., 2011 A statistical framework for SNP calling, mutation discovery, association mapping
   and population genetical parameter estimation from sequencing data. Bioinformatics 27: 2987-2993.

   genotype 11:  1 is the true base, 0 is error.  Data is k, the number of 1's on p.20.

   cl -Ox genotypecall.c tools.c
   icc -o genotypecall -O3 genotypecall.c tools.c -lm
   genotypecall
*/

#include "paml.h"
#define Maruki_Lynch 0

int call_genotype(double lnL[3], double* dlnL);
void GTerror_Exact(void);
void GTerror_MC(void);
void bias_GT_Exact(double GQstar);

int call_genotype(double lnL[3], double* dlnL)
{
   int index[3], space[3], GT;

   indexing(lnL, 3, index, 1, space);
   GT = index[0];
   *dlnL = lnL[index[0]] - lnL[index[1]];
   if (*dlnL < 1e-6) printf("\nDo we have a tie?\n");

   if (Maruki_Lynch && GT == 1 && *dlnL < 2.7055 / 2) {  /* 1.9207 */
      GT = index[1];
      /* printf("\aconverting heterozygote to homozygote! "); */
   }
   return (GT);
}

void GTerror_MC(void)
{
   int nr = 1000000, ir, n = 5, k, ie, iGT, GTe = 0;
   double e, lne, ln1e, lnL[3], dlnL, GTerror[2];

   printf("simulating %d replicates\n", nr);
   printf("\ninput per-read per-base error rate? ");
   scanf("%lf", &e);
   printf("\nper-read per-base error rate e = %.6g\n", e);
   SetSeed(-1, 0);
   lne = log(e); ln1e = log(1 - e);
   printf("\n%s\t%s\t%s\n\n", "reads", "homo_error", "hetero_error");
   for (n = 2; n <= 30; n++) {
      GTerror[0] = GTerror[1] = 0;
      lnL[1] = -n * log(2.0);
      for (iGT = 0; iGT < 2; iGT++) {  /* iGT0: homozygote (11), iGT1: heterozygote (01) */
         for (ir = 0; ir < nr; ir++) {
            if (iGT == 0) k = rndBinomial(n, 1 - e);  /* iGT0: homozygote (11) */
            else          k = rndBinomial(n, 0.5);    /* iGT1: heterozygote (01) */
            lnL[0] = (n - k) * ln1e + k * lne;  /* given GT = 00, data is (k, n - k)  */
            lnL[2] = k * ln1e + (n - k) * lne;  /* given GT = 11, data is (k, n - k)  */
            GTe = call_genotype(lnL, &dlnL);
            if ((iGT == 0 && GTe != 2) || (iGT == 1 && GTe != 1))
               GTerror[iGT] ++;
         }
      }
      printf("%d\t%.6f\t%.6f\n", n, GTerror[0] / nr, GTerror[1] / nr);
   }
}

void GTerror_Exact(void)
{
   int maxnreads = 30, n, k, GTe;
   double e = 0.01, lne, ln1e, lnL[3], dlnL, pk[2], bnk, z, GTerror[2];
   char* GTstr[3] = { "00", "01", "11" };
   FILE* fout;
   char line[1024] = "genotypecallingerror_e0.005.txt";

   if (Maruki_Lynch)
      printf("\n*** This mimics Maruki & Lynch 2017, figure 1 and applies LRT cutoff of 1.9207 *** ");

   printf("calculating genotype calling error rate given read depth n and base-call error e\n");
   printf("input per-read per-base error rate? ");
   scanf("%lf", &e);
   printf("\nper-read per-base error rate e = %.6g\n", e);
   sprintf(line, "genotypecallingerror_e%g.txt", e);
   fout = (FILE*)fopen(line, "w");
   if (fout == NULL) zerror("file open error");

   lne = log(e); ln1e = log(1 - e);
   fprintf(fout, "\n%s\t%s\t%s\n\n", "reads", "homo_error", "hetero_error");
   for (n = 2; n <= maxnreads; n++) {
      printf(" n  k     lnL00    lnL01    lnL11   dlnL  GT  pk0(GT11) pk1(GT01)\n");
      GTerror[0] = GTerror[1] = 0;
      for (k = 0; k <= n; k++) {
         lnL[1] = -n * log(2.0);
         lnL[0] = (n - k) * ln1e + k * lne;
         lnL[2] = k * ln1e + (n - k) * lne;
         GTe = call_genotype(lnL, &dlnL);
         bnk = Binomial(n, k, &z);
         if (z) zerror("perhaps n is too large for this?");
         pk[0] = bnk * exp(lnL[2]);   /* prob(k|GT 11) */
         pk[1] = bnk * exp(lnL[1]);   /* prob(k|GT 01) */
         if (GTe != 2) GTerror[0] += pk[0];
         if (GTe != 1) GTerror[1] += pk[1];
         printf("%2d %2d %9.3f%9.3f%9.3f%7.2f %3s %9.5f%9.5f\n", n, k, lnL[0], lnL[1], lnL[2], dlnL, GTstr[GTe], pk[0], pk[1]);
      }
      printf("\nn = %2d  GTcallerror (for GT 11, 01) = %10.6f %10.6f\n\n", n, GTerror[0], GTerror[1]);
      fprintf(fout, "%d\t%.6f\t%.6f\n", n, GTerror[0], GTerror[1]);
   }
   fclose(fout);
}

void bias_GT_Exact(double GQstar)
{
   /* GQ=20 means log(100)=4.60517 */
   int maxnreads = 30, n, k, GTe;
   double e = 0.01, theta = 0.01, dlnLstar = GQstar / 10 * log(10), pGTe[3], sGTe;
   double lne, ln1e, lnL[3], dlnL, pk[2], bnk, z;

   printf("Estimating heterozygosity when GT is called using a GQ cutoff of %.1f\n", GQstar);
   printf("\nper-read per-base error rate e = %.6g\nheterozygosity = theta = %.6g\n", e, theta);
   printf("\n n   called genotypes (00 01 11)\n");

   lne = log(e); ln1e = log(1 - e);
   for (n = 2; n <= maxnreads; n++) {
      pGTe[0] = pGTe[1] = pGTe[2] = 0;
      for (k = 0; k <= n; k++) {
         lnL[1] = -n * log(2.0);              /* 1: 01 */
         lnL[0] = (n - k) * ln1e + k * lne;   /* 0: 00 */
         lnL[2] = k * ln1e + (n - k) * lne;   /* 2: 11 */
         GTe = call_genotype(lnL, &dlnL);
         bnk = Binomial(n, k, &z);
         if (z) zerror("perhaps n is too large for this?");
         pk[0] = bnk * exp(lnL[2]);   /* prob(k|GT 11) */
         pk[1] = bnk * exp(lnL[1]);   /* prob(k|GT 01) */
         if (dlnL < dlnLstar) continue;
         pGTe[GTe] += (1 - theta) * pk[0] + theta * pk[1];
      }
      sGTe = pGTe[0] + pGTe[1] + pGTe[2];
      printf("%2d %9.5f%9.5f%9.5f  h = %9.5f\n", n, pGTe[0], pGTe[1], pGTe[2], pGTe[1] / sGTe);
   }
}

void main(int argc, char *argv[])
{
   int option;
   if (argc < 2) {
      printf("Usage: \n\tgenotypecall 0: GTerror exact\n\tgenotypecall 1: GT error simulation\n\tgenotypecall 2: GT bias cutoff\n");
      exit(0);
   }
   option = atoi(argv[1]);
   if (option==0)         GTerror_Exact();
   else if (option == 1)  GTerror_MC();
   else if (option == 2)  bias_GT_Exact(20);
   exit(0);
}

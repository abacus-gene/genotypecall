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


int GTerror_Exact(void)
{
   int maxnreads = 30, n, k, iGT, GT;
   double e = 0.01, lne, ln1e, lnL11, lnL01, lnL00, p, psum[2], GTerror[2], z;
   char *GTstr[2] = { "homozygote 11", "heterozygote 01" }, *GTcalled[3] = { "00", "01", "11" };
   FILE *fout;
   char line[1024] = "genotypecallingerror_e0.005.txt";

   printf("\ninput per-read per-base error rate? ");
   scanf("%lf", &e);
   printf("\nper-read per-base error rate e = %.6g\n", e);
   sprintf(line, "genotypecallingerror_e%g.txt", e);
   fout = (FILE*)fopen(line, "w");
   if (fout == NULL) zerror("file open error");

   lne = log(e), ln1e = log(1 - e);
   fprintf(fout, "\n%s\t%s\t%s\n\n", "reads", "homo_error", "hetero_error");
   for (n = 2; n <= maxnreads; n++) {
      printf("\n%s****** %s ****** %21s ****** %s ****** \n", "", GTstr[0], "", GTstr[1]);
      printf(" n  k     prob    lnL11   lnL01   lnL00 GTcall    ");
      printf(" n  k     prob    lnL11   lnL01   lnL00 GTcall\n");
      GTerror[0] = GTerror[1] = 0;
      psum[0] = psum[1] = 0;
      for (k = 0; k <= n; k++) {
         for (iGT = 0; iGT < 2; iGT++) {    /* iGT0: homozygote(11), iGT1: heterozygote(01) */
            lnL01 = -n*log(2.0);
            lnL00 = (n - k)*ln1e + k*lne;
            lnL11 = k*ln1e + (n - k)*lne;
            if (iGT == 0) p = Binomial(n, k, &z)*exp(lnL11);
            else          p = Binomial(n, k, &z)*exp(lnL01);
            if (z) zerror("perhaps n is too large for this?");
            if (lnL11 > lnL00 && lnL11 > lnL01)      GT = 2;
            else if (lnL01 > lnL00 && lnL01 > lnL11) GT = 1;
            else if (lnL00 > lnL01 && lnL00 > lnL11) GT = 0;
            else
               printf("\nI see a tie?\n");

            if ((iGT == 0 && GT != 2) || (iGT == 1 && GT != 1))
               GTerror[iGT] += p;
            psum[iGT] += p;
            printf("%2d %2d %8.4f %8.3f%8.3f%8.3f %5s", n, k, p, lnL11, lnL01, lnL00, GTcalled[GT]);
            printf(iGT ? "\n" : "     ");
         }
      }
      if (fabs(1 - psum[0]) > 1e-9 || fabs(1 - psum[1]) > 1e-9)
         zerror("pk does not sum to 1.");
      printf("\n%10sn = %2d  GTcallerror = %9.7f %10s %2d  GTcallerror = %9.7f\n", "", n, GTerror[0], "n =", n, GTerror[1]);
      fprintf(fout, "%d\t%.7f\t%.7f\n", n, GTerror[0], GTerror[1]);
   }
   fclose(fout);
   exit(0);
}

int GTerror_MC(void)
{
   int nr = 1000000, ir, n = 5, k, ie, iGT, GT = 0;
   double e1x[3] = { 0.005, 0.01, 0.05 }, e, lne, ln1e;
   double lnL11, lnL01, lnL00, GTerror[2];
   FILE *fout;
   char line[1024] = "genotypecallingerror_e0.005.txt";

   SetSeed(-1, 0);
   for (ie = 0; ie < 3; ie++) {
      e = e1x[ie], lne = log(e), ln1e = log(1 - e);
      sprintf(line, "genotypecallingerror_e%.g.txt", e);
      fout = (FILE*)fopen(line, "w");
      if (fout == NULL) zerror("file open error");
      printf("\nper-read per-base error rate = %.7f\n", e);
      printf("\n%s\t%s\t%s\n\n", "reads", "homo_error", "hetero_error");
      fprintf(fout, "\n%s\t%s\t%s\n\n", "reads", "homo_error", "hetero_error");
      for (n = 2; n <= 30; n++) {
         GTerror[0] = GTerror[1] = 0;
         lnL01 = -n*log(2.0);
         for (iGT = 0; iGT < 2; iGT++) {  /* iGT0: homozygote (11), iGT1: heterozygote (01) */
            for (ir = 0; ir < nr; ir++) {
               if (iGT == 0) k = rndBinomial(n, 1 - e);  /* iGT0: homozygote (11) */
               else          k = rndBinomial(n, 0.5);    /* iGT1: heterozygote (01) */
               lnL00 = (n - k)*ln1e + k*lne;  /* given GT = 00, data is (k, n - k)  */
               lnL11 = k*ln1e + (n - k)*lne;  /* given GT = 11, data is (k, n - k)  */
               if (lnL11 > lnL00 && lnL11 > lnL01)      GT = 2;
               else if (lnL01 > lnL00 && lnL01 > lnL11) GT = 1;
               else if (lnL00 > lnL01 && lnL00 > lnL11) GT = 0;
               else                                     printf("\nI see a tie?\n");
               if (iGT == 0) {      /* iGT0: homozygote (11) */
                  if (GT != 2)
                     GTerror[iGT]++;
               }
               else {               /* iGT1: heterozygote (01) */
                  if (GT != 1)
                     GTerror[iGT]++;
               }
            }
         }
         printf("%d\t%.7f\t%.7f\n", n, GTerror[0] / nr, GTerror[1] / nr);
         fprintf(fout, "%d\t%.7f\t%.7f\n", n, GTerror[0] / nr, GTerror[1] / nr);
      }
      fclose(fout);
   }
   exit(0);
}

void main(void)
{
   GTerror_Exact();
   /* GTerror_MC(); */
}

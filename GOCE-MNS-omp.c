#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

//vystupny subor B a L a u do stvorcovej
//q - posledny stlpec suboru bude prava strana

#define N 8102 // pocet zdrojovych bodov
#define M 160002 // pocet meracich bodov

#define TOL 0.00000000000001
#define KONST 1000000.0

#define RE 6378137.0
#define E2 0.00669437999014

// prevod jednotiek
double rad(double degrees)
{
   return (degrees * M_PI) / 180.0;
}

// načítanie bodov zo súboru pre Sj  (zo súboru BL-3602.dat) //tieto budu vzdy tie vacsie
void loadSourceData(const char *filename, double *B, double *L, double *H, double *dg, double *f)
{
   FILE *file = fopen(filename, "r");
   if (file == NULL)
   {
       printf("Unable to open file\n");
       exit(1);
   }
   for (int i = 0; i < N; i++) // Povedzme M je počet zdrojových bodov
   {
       fscanf(file, "%lf %lf %lf %lf %lf", &B[i], &L[i], &H[i], &dg[i], &f[i]);
   }
   fclose(file);
}

// načítanie bodov merania pre Xi (zo súboru BL-8102.dat) //tieto budu vzdy tie mensie
void loadMeasurementData(const char *filename, double *B, double *L, double *RR, double *Tzz)
{
   FILE *file = fopen(filename, "r");
   if (file == NULL)
   {
       printf("Unable to open file\n");
       exit(1);
   }
   for (int i = 0; i < M; i++) // Povedzme N je počet bodov merania
   {
       fscanf(file, "%lf %lf %lf %lf %lf", &B[i], &L[i], &RR[i], &Tzz[i]);
   }
   fclose(file);
}

// dot product
double dot_product(const double *vec1, const double *vec2)
{
   double result = 0.0;
   for (int i = 0; i < N; i++)
   {
       result += vec1[i] * vec2[i];
   }
   return result;
}

// matica*vektor
void matrix_vector_mult(double **A, const double *vec, double *result)
{
   for (int i = 0; i < N; i++)
   {
       result[i] = 0.0;
       for (int j = 0; j < N; j++)
       {
           result[i] += A[i][j] * vec[j];
       }
   }
}

int main()
{

   const double GM = 398600.5;
   double R = 6371.0;
   double alt = 250.0;
   double dGM = ((2 * GM) / ((R + alt) * (R + alt) * (R + alt))); // vzorec
   // printf("%.10lf", dGM); //vyslo 0.0000027629

   // Dynamicky alokovane polia
   double(*coordinatesS)[3] = (double(*)[3])malloc(N * sizeof(double[3])); // S coordinates
   double(*coordinatesX)[3] = (double(*)[3])malloc(M * sizeof(double[3])); // X coordinates
   double(*coordinatesE)[3] = (double(*)[3])malloc(M * sizeof(double[3])); // E coordinates

   double *B_source = (double *)malloc(N * sizeof(double));
   double *L_source = (double *)malloc(N * sizeof(double));
   double *H_source = (double *)malloc(N * sizeof(double));
   double *dg_source = (double *)malloc(N * sizeof(double));
   double *f_source = (double *)malloc(N * sizeof(double));

   double *dGMarray_measurement = (double *)malloc(M * sizeof(double));

   double *B_measurement = (double *)malloc(M * sizeof(double));
   double *L_measurement = (double *)malloc(M * sizeof(double));
   double *RR = (double *)malloc(M * sizeof(double));
   double *Tzz = (double *)malloc(M * sizeof(double));
   //double *f_measurement = (double *)malloc(M * sizeof(double));
   double *u = calloc(M, sizeof(double));

//   omp_set_num_threads(16); //pôjde to na 4 jadrach

   // Alokácia pre maticu A a S
   double **A = (double **)malloc(M * sizeof(double *));
   for (int i = 0; i < M; i++)
   {
       A[i] = (double *)malloc(N * sizeof(double));
   }
   double **S = (double **)malloc(N * sizeof(double *));
   for (int i = 0; i < N; i++)
   {
       S[i] = (double *)malloc(N * sizeof(double));
   }

   /*if (coordinatesS == NULL || coordinatesX == NULL || coordinatesE == NULL ||
   A == NULL || S == NULL || B_source == NULL || L_source == NULL ||
   H_source == NULL || B_measurement == NULL || L_measurement == NULL ||
   H_measurement == NULL || dg == NULL || f == NULL || dGMarray == NULL ||
   b == NULL || x == NULL || r == NULL || r_hat == NULL || p == NULL ||
   s == NULL || t == NULL || v == NULL) {
   printf("Memory allocation failed\n");
   return 1;
}*/

   loadSourceData("BL-8102.dat", B_source, L_source, H_source, dg_source, f_source);
   loadMeasurementData("BLRTzz-160002-GOCE-filtered.dat", B_measurement, L_measurement, RR, Tzz);
   printf("Data loaded successfully.\n");

   // Nacitanie potrebnych dat

   for (int i = 0; i < M; i++) // pre source pointy
   {
       dGMarray_measurement[i] = dGM; // konstantna matica, co bude vlastne nasa parav strana, A*alpha = dGMarray
   }

       double BSrad, LSrad, BXrad, LXrad;
    n = RE / sqrt(1 - E2 * pow(sin(Brad), 2));

   // koordinaty
   /* for (int i = 0; i < N; i++)
    {
       BSrad = rad(B_source[i]);
       LSrad = rad(L_source[i]);

       // S coordinates
       coordinatesS[i][0] = (n + H[i]) * cos(BSrad) * cos(LSrad); // X
       coordinatesS[i][1] = (n + H[i]) * cos(BSrad) * sin(LSrad); // Y
       coordinatesS[i][2] = (n + H[i]) * sin(BSrad);              // Z

    }

   for (int i = 0; i < M; i++)
    {
       BXrad = rad(B_measurement[i]);
       LXrad = rad(L_measurement[i]);
       // X coordinates
       coordinatesX[i][0] = (R + alt) * cos(BXrad) * cos(LXrad); // X
       coordinatesX[i][1] = (R + alt) * cos(BXrad) * sin(LXrad); // Y
       coordinatesX[i][2] = (R + alt) * sin(BXrad);              // Z

       // E coordinates

       coordinatesE[i][0] = cos(BXrad) * cos(LXrad); // X
       coordinatesE[i][1] = cos(BXrad) * sin(LXrad); // Y
       coordinatesE[i][2] = sin(BXrad);              // Z
   } */
    for (int i = 0; i < N; i++)
        {
            double BSrad = rad(B_sourcce[i]);
            double LSrad = rad(L_source[i]);
            double n;
            n = RE / sqrt(1 - E2 * pow(sin(Brad), 2));

            // S coordinates
            coordinatesS[i][0] = (n + H_source[i]) * cos(BSrad) * cos(LSrad); // X
            coordinatesS[i][1] = (n + H_source[i]) * cos(BSrad) * sin(LSrad); // Y
            coordinatesS[i][2] = (n + H_source[i]) * sin(BSrad);             // Z
        }

    // measurement points koordinaty
    for (int i = 0; i < M; i++)
        {
            double BXrad = rad(B_measurement[i]);
            double LXrad = rad(L_measurement[i]);

            // X coordinates
            coordinatesS[i][0] = RR * cos(BXrad) * cos(LXrad); // X
            coordinatesS[i][1] = RR * cos(BXrad) * sin(LXrad); // Y
            coordinatesS[i][2] = RR * sin(BXrad);             // Z

            // E coordinates

        coordinatesE[i][0] = cos(BXrad) * cos(LXrad); // X
        coordinatesE[i][1] = cos(BXrad) * sin(LXrad); // Y
        coordinatesE[i][2] = sin(BXrad);         
        }
   
   printf("Prešli koordinaty, Continue...\n");

   // MATICA A
    double dx, dy, dz, rij, dotProduct;

    #pragma omp parallel for private (dx, dy, dz, rij, dotProduct)

    for (int i = 0; i < M; i++)
     {
       for (int j = 0; j < N; j++)
       {
           double p1[3] = {coordinatesX[i][0], coordinatesX[i][1], coordinatesX[i][2]};
           double p2[3] = {coordinatesS[j][0], coordinatesS[j][1], coordinatesS[j][2]};

           // komponenty
           dx = p1[0] - p2[0];
           dy = p1[1] - p2[1];
           dz = p1[2] - p2[2];
           rij = sqrt(dx * dx + dy * dy + dz * dz);

           // dot product
           dotProduct = dx * coordinatesE[i][0] + dy * coordinatesE[i][1] + dz * coordinatesE[i][2];

           // hodnoty A[i][j]
           A[i][j] = (1 / pow(rij, 3)) - ((3 * pow(dotProduct, 2)) / pow(rij, 5));
       }
     }
   printf("prešlo vytvorenie A, Continue...\n");

   // Overenie, prvych 5 prvkov
   for (int i = 0; i < 5; i++)
   {
       for (int j = 0; j < 5; j++)
       {
           printf("A[%d][%d] = %.15f\n", i, j, A[i][j]);
       }
   }
   printf("Vytlačili sa prvky A, Continue...\n");

   // Výpočet súčinu A^T * A (A_T * A)

#pragma omp parallel for
   for (int i = 0; i < N; i++) {
       for (int j = i; j < N; j++) { // začíname od i, aby sme vypĺňali len hornú časť
           S[i][j] = 0.0;
           for (int k = 0; k < M; k++) {
               S[i][j] += A[k][j] * A[k][i];
           }
       }
   }

#pragma omp parallel for
   for (int i = 1; i < N; i++) {
       for (int j = 0; j < i; j++) {
           S[i][j] = S[j][i];
       }
   }

   printf("Prešiel súčin AT*A, teda vznikla S, Continue...\n");

   for (int i = 0; i < 6; i++)
   {
       for (int j = i; j < 6; j++)
       {
           printf("S[%d][%d] = %.15f\n", i, j, KONST*S[i][j]);
       }
   }

   printf("Vytlačili sa prvé prvky S, , Continue...\n");

   double *b = (double *)malloc(N * sizeof(double));
   double *x = (double *)malloc(N * sizeof(double)); // Initial guess x(0)
   double *r = malloc(N * sizeof(double));
   double *r_hat = malloc(N * sizeof(double));
   double *p = malloc(N * sizeof(double));
   double *s = malloc(N * sizeof(double));
   double *t = malloc(N * sizeof(double));
   double *v = malloc(N * sizeof(double));
   double *at_dg = malloc(N * sizeof(double));

   // Perform A_T * dgM_source
   // A_T NxM dgM_source Mx1
   for (int k = 0; k < N; k++)
   {
       at_dg[k] = 0.0;
   for (int i = 0; i < M; i++)
   {
       
       at_dg[k] += A[i][k] * dGMarray_measurement[i];
     
   }
   }

   /* for (int i = 0; i < N; i++) {
   printf("b[%d] = %.10f\n", i, b[i]);
   } */

 
   printf("Vytvorila sa PS pre BCG, AT*dgM_source, Continue...\n");

   // A*x=b
   for (int i = 0; i < N; i++)
   {
       b[i] = - Tzz[i] * 0.0000000000001;             // do b vlozim hodnoty dGM, prava strana //resp AT*dgM(vektor velkosti M)
       x[i] = 0.0;                 // Initial guess x(0) = 0
       r[i] = r_hat[i] = b[i]; // r = r s vlnkou
       // PS bude teraz AT*dg
       //Tzz - 4. stlpec z noveho suboru bude prava strana
   }

   double alpha = 1.0, omega = 1.0, rho_new, beta, rho_old = dot_product(r_hat, r);

   int max_iter = 50; // pocet iteracii
   for (int iter = 1; iter <= max_iter; iter++)
   {
       // [UTF-8?]Ď(i-1) = [UTF-8?]r_hatáľ r(i-1)
       rho_new = dot_product(r_hat, r);
       if (rho_new == 0)
       {
           printf("BiCGSTAB failed: rho je 0.\n");
           break;
       }

       if (iter == 1)
       {
           for (int j = 0; j < N; j++)
           {
               p[j] = r[j]; // Initial p = r
           }
       }
       else
       {
           // [UTF-8?]Î˛(i-1) = [UTF-8?](Ď(i-1) / [UTF-8?]Ď(i-2)) * [UTF-8?](Îą(i-1) / [UTF-8?]Ď(i-1))
           beta = (rho_new / rho_old) * (alpha / omega);
           for (int j = 0; j < N; j++)
           {
               p[j] = r[j] + beta * (p[j] - omega * v[j]);
           }
       }

       // v(i) = A * p(i) //S*p(i)
       

       for (int i = 0; i < N; i++)
       {
           v[i] = 0.0;
           for (int j = 0; j < N; j++)
           {
               v[i] += S[i][j] * p[j];
           }
       }

       // [UTF-8?]Îą(i) = [UTF-8?]Ď(i-1) / [UTF-8?](r_hatáľ v(i))
       alpha = rho_new / dot_product(r_hat, v);

       // s = r(i-1) - [UTF-8?]Îą(i) * v(i)
       for (int j = 0; j < N; j++)
       {
           s[j] = r[j] - alpha * v[j];
       }

       // je norma prilis mala?
       double norm_s = sqrt(dot_product(s, s));
       printf("%d iteracia \t%.15lf\n",iter,norm_s);
       if (norm_s < TOL)
       {
           for (int j = 0; j < N; j++)
           {
               x[j] += alpha * p[j];
           }
           printf("BiCGSTAB konvergoval v %d iteracii.\n", iter);
           break;
       }

       // t = A * s //S÷s
       
       for (int i = 0; i < N; i++)
       {
           t[i] = 0.0;
           for (int j = 0; j < N; j++)
           {
               t[i] += S[i][j] * s[j];
           }
       }

       // [UTF-8?]Ď(i) = [UTF-8?](táľ s) / [UTF-8?](táľ t)
       omega = dot_product(t, s) / dot_product(t, t);

       // Update x(i) = x(i-1) + [UTF-8?]Îą(i) * p(i) + [UTF-8?]Ď(i) * s
       for (int j = 0; j < N; j++)
       {
           x[j] += alpha * p[j] + omega * s[j];
       }

       // r(i) = s - [UTF-8?]Ď(i) * t
       for (int j = 0; j < N; j++)
       {
           r[j] = s[j] - omega * t[j];
       }

       // Check convergence based on norm of r(i)
       double norm_r = sqrt(dot_product(r, r));

        printf("%d iteráciq:\t\t%.12lf\n", iter,norm_r);

       if (norm_r < TOL)
       {
           printf("BiCGSTAB konvergoval v %d iterácií.\n", iter);
           break;
       }

       // Update [UTF-8?]Ď(i-1) pre dalsiu iteraciu
       rho_old = rho_new;

       // ci je omega 0
       if (omega == 0)
       {
           printf("BiCGSTAB padol: omega je 0.\n");
           break;
       }
   }

   // vektor x, co je vlastne nasa hladana alpha
   for (int i = 0; i < 10; i++)
   {
       printf("x[%d] = %.15f\n", i, x[i]);
   }
    for(int i = 0; i < M; i++)
    {
        double BXrad = rad(B_measurement[i]);
        double LXrad = rad(L_measurement[i]);
        double n = RE / sqrt(1 - E2 * pow(sin(Brad), 2));
        // X coordinates
        coordinatesX[i][0] = (n + alt) * cos(BXrad) * cos(LXrad); // X
        coordinatesX[i][1] = (n + alt) * cos(BXrad) * sin(LXrad); // Y
        coordinatesX[i][2] = (n + alt) * sin(BXrad);              // Z
        

    }
    double G_ij;

    #pragma omp parallel for private (dx, dy, dz, rij, G_ij)
 
   for (int i = 0; i < M; i++)
   {
       for (int j = 0; j < N; j++)
       {
           double dx = coordinatesX[i][0] - coordinatesS[j][0];
           double dy = coordinatesX[i][1] - coordinatesS[j][1];
           double dz = coordinatesX[i][2] - coordinatesS[j][2];
           double rij = sqrt(dx * dx + dy * dy + dz * dz);

           if (rij != 0)
           {
               G_ij = 1.0 / rij;
               u[i] += x[j] * G_ij;
           }
       }
   }

   // overenie hodnot u
   for (int i = 0; i < 10; i++)
   {
       printf("u[%d] = %.10f\n", i, u[i]);
   }
    FILE *out1 = fopen("outBLx.dat", "w");
    FILE *out2 = fopen("outBLu.dat", "w");
    if (out1 == NULL)
    {
         printf("nepodarilo sa otvorit\n");
         return 1;
    }
    if (out2 == NULL)
    {
         printf("nepodarilo sa otvorit\n");
         return 1;
    }
    for (int i = 0; i < N; i++)
    {
         fprintf(out1, "%.2f %.5f %.3f\n", B_source[i], L_source[i], x[i]);
    }
    for (int i = 0; i < M; i++)
    {
         fprintf(out2, "%.2f %.5f %.3f\n", B_measurement[i], L_measurement[i], u[i]);
    }

    fclose(out1);
    fclose(out2);
 
   for (int i = 0; i < N; i++)
   {
       free(S[i]);
   }
   free(S);

 

   // uvolnenie pamate
   free(coordinatesS);
   free(coordinatesX);
   free(coordinatesE);
   free(at_dg);
   free(B_source);
   free(L_source);
   free(H_source);
   free(B_measurement);
   free(dGMarray_measurement);
   free(L_measurement);
   free(RR);
   free(Tzz);
   //free(f_source);
   free(dg_source);
   free(f_measurement);
   free(b);
   free(x);
   free(r);
   free(r_hat);
   free(p);
   free(u);
   for (int i = 0; i < M; i++)
   {
       free(A[i]);
   }
   free(A);
   free(s);
   free(t);
   free(v);
   return 0;
} 


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define N 1600002 // pocet bodov
#define TOL 1e-9
#define RE 6378137.0
#define E2 0.00669437999014

//     gcc-14 -fopenmp g_realnehod.c -o g_realnehod
//    ./g_realnehod
// prevod jednotiek
double rad(double degrees)
{
     return (degrees * M_PI) / 180.0;
}

// Nacitanie zo suboru
void loadData(const char *filename, double *B, double *L, double *H, double *dg, double *f)
{
     FILE *file = fopen(filename, "r");
     if (file == NULL)
     {
          printf("Unable to open file\n");
          exit(1);
     }
     for (int i = 0; i < N; i++)
     {
          fscanf(file, "%lf %lf %lf %lf %lf", &B[i], &L[i], &H[i], &dg[i], &f[i]);
     }
     fclose(file);
}

void loadData2(const char *filename, double *B2, double *L2, double *RR, double *Tzz)
{
     FILE *file = fopen(filename, "r");
     if (file == NULL)
     {
          printf("Unable to open file\n");
          exit(1);
     }
     for (int i = 0; i < N; i++)
     {
          fscanf(file, "%lf %lf %lf %lf %lf", &B2[i], &L2[i], &RR[i], &Tzz[i]);
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
#pragma omp parallel
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
     omp_set_num_threads(4);

     const double GM = 398600.5;
     double R = 6378000.0;
     double alt = 230000.0;
     double dGM = ((2 * GM) / ((R + alt) * (R + alt) * (R + alt))); // vzorec
     // printf("%.10lf", dGM); //vyslo 0.0000027629

     // Dynamicky alokovane polia
     double(*coordinatesS)[3] = (double(*)[3])malloc(N * sizeof(double[3])); // S coordinates
     double(*coordinatesX)[3] = (double(*)[3])malloc(N * sizeof(double[3])); // X coordinates
     double(*coordinatesE)[3] = (double(*)[3])malloc(N * sizeof(double[3])); // E coordinates
     double *B = (double *)malloc(N * sizeof(double));
     double *L = (double *)malloc(N * sizeof(double));
     double *H = (double *)malloc(N * sizeof(double));
     double *dg = (double *)malloc(N * sizeof(double));
     double *f = (double *)malloc(N * sizeof(double));
     double *B2 = (double *)malloc(N * sizeof(double));
     double *L2 = (double *)malloc(N * sizeof(double));
     double *RR = (double *)malloc(N * sizeof(double));
     double *Tzz = (double *)malloc(N * sizeof(double));

     // MaticA A
     double **A = (double **)malloc(N * sizeof(double *));
     for (int i = 0; i < N; i++)
     {
          A[i] = (double *)malloc(N * sizeof(double));
     }
     if (B == NULL || L == NULL || H == NULL || dg == NULL || f == NULL)
     {
          printf("Memory allocation failed\n");
          return 1;
     }

     loadData("/", B, L, H, dg, f);
     loadData2("/", B2, L2, RR, Tzz);

     // loadDataAlt("/Users/hannah/Desktop/GGGGG/gocee/BL-8102.dat", Alt);

     /* //
     for (int i = 0; i < 5; i++)
     {
         printf("B[%d] = %.2f, L[%d] = %.2f, H[%d] = %.2f, dg[%d] = %.2f, f[%d] = %.2f\n",
                i, B[i], i, L[i], i, H[i], i, dg[i], i, f[i]);
     }
  */
     // Nacitanie potrebnych dat

     double *dGMarray = (double *)malloc(N * sizeof(double));
     for (int i = 0; i < N; i++)
     {
          dGMarray[i] = dGM; // konstantna matica, co bude vlastne nasa parav strana, A*alpha = dGMarray
     }

     // koordinaty
     for (int i = 0; i < N; i++)
     {
          double Brad = rad(B[i]);
          double Lrad = rad(L[i]);
          double B2rad = rad(B2[i]);
          double L2rad = rad(L2[i]);
          double n;
          n = RE / sqrt(1 - E2 * pow(sin(Brad), 2));

          // S coordinates
          coordinatesS[i][0] = (n + H[i]) * cos(Brad) * cos(Lrad); // X
          coordinatesS[i][1] = (n + H[i]) * cos(Brad) * sin(Lrad); // Y
          coordinatesS[i][2] = (n + H[i]) * sin(Brad);             // Z

          // X coordinates
          coordinatesX[i][0] = RR[i] * cos(B2rad) * cos(L2rad); // X
          coordinatesX[i][1] = RR[i] * cos(B2rad) * sin(L2rad); // Y
          coordinatesX[i][2] = RR[i] * sin(B2rad);              // Z

          // E coordinates
          coordinatesE[i][0] = cos(Brad) * cos(Lrad); // X
          coordinatesE[i][1] = cos(Brad) * sin(Lrad); // Y
          coordinatesE[i][2] = sin(Brad);             // Z
     }

     // MATICA A
     double dx, dy, dz, rij; // Declare variables before the parallel region
#pragma omp parallel for private(dx, dy, dz, rij)
     for (int i = 0; i < N; i++)
     {
          for (int j = 0; j < N; j++)
          {
               double p1[3] = {coordinatesX[i][0], coordinatesX[i][1], coordinatesX[i][2]};
               double p2[3] = {coordinatesS[j][0], coordinatesS[j][1], coordinatesS[j][2]};
               double dx, dy, dz, rij;

               // komponenty
               dx = p1[0] - p2[0];
               dy = p1[1] - p2[1];
               dz = p1[2] - p2[2];
               rij = sqrt(dx * dx + dy * dy + dz * dz);

               // dot product
               double dotProduct = dx * coordinatesE[i][0] + dy * coordinatesE[i][1] + dz * coordinatesE[i][2];

               // hodnoty A[i][j]
               A[i][j] = (1 / pow(rij, 3)) - ((3 * pow(dotProduct, 2)) / pow(rij, 5));
          }
     }

     // Overenie, prvych 5 prvkov
     for (int i = 0; i < 5; i++)
     {
          for (int j = 0; j < 5; j++)
          {
               printf("A[%d][%d] = %.15f\n", i, j, A[i][j]);
          }
     }

     double *b = (double *)malloc(N * sizeof(double));
     double *x = (double *)malloc(N * sizeof(double)); // Initial guess x(0)
     double *r = malloc(N * sizeof(double));
     double *r_hat = malloc(N * sizeof(double));
     double *p = malloc(N * sizeof(double));
     double *s = malloc(N * sizeof(double));
     double *t = malloc(N * sizeof(double));
     double *v = malloc(N * sizeof(double));

     // A*x=b

     for (int i = 0; i < N; i++)
     {
          // poseldny stplec druheho suboru, opravit
          b[i] = -f[i] * 0.000000001; // do b vlozim realne hodnoty z vektora dg
          x[i] = 0.0;                 // Initial guess x(0) = 0
          r[i] = r_hat[i] = b[i];     // r = r s vlnkou
     }

     double alpha = 1.0, omega = 1.0, rho_new, beta, rho_old = dot_product(r_hat, r);

     int max_iter = 1000; // pocet iteracii
     for (int iter = 1; iter <= max_iter; iter++)
     {
          // ρ(i-1) = r_hatᵀ r(i-1)
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
               // β(i-1) = (ρ(i-1) / ρ(i-2)) * (α(i-1) / ω(i-1))
               beta = (rho_new / rho_old) * (alpha / omega);
               for (int j = 0; j < N; j++)
               {
                    p[j] = r[j] + beta * (p[j] - omega * v[j]);
               }
          }

          // v(i) = A * p(i)
          matrix_vector_mult(A, p, v);

          // α(i) = ρ(i-1) / (r_hatᵀ v(i))
          alpha = rho_new / dot_product(r_hat, v);

          // s = r(i-1) - α(i) * v(i)
          for (int j = 0; j < N; j++)
          {
               s[j] = r[j] - alpha * v[j];
          }

          // Check if norm(s) is small enough to stop
          double norm_s = sqrt(dot_product(s, s));
          if (norm_s < TOL)
          {
               for (int j = 0; j < N; j++)
               {
                    x[j] += alpha * p[j];
               }
               printf("BiCGSTAB konvergoval v %d iteracii.\n", iter);
               break;
          }

          // t = A * s
          matrix_vector_mult(A, s, t);

          // ω(i) = (tᵀ s) / (tᵀ t)
          omega = dot_product(t, s) / dot_product(t, t);

          // Update x(i) = x(i-1) + α(i) * p(i) + ω(i) * s
          for (int j = 0; j < N; j++)
          {
               x[j] += alpha * p[j] + omega * s[j];
          }

          // r(i) = s - ω(i) * t
          for (int j = 0; j < N; j++)
          {
               r[j] = s[j] - omega * t[j];
          }

          // Check convergence based on norm of r(i)
          double norm_r = sqrt(dot_product(r, r));
          if (norm_r < TOL)
          {
               printf("BiCGSTAB konvergoval v %d iterácií.\n", iter);
               break;
          }

          // Update ρ(i-1) pre dalsiu iteraciu
          rho_old = rho_new;

          // ci je omega 0
          if (omega == 0)
          {
               printf("BiCGSTAB padol: omega je 0.\n");
               break;
          }
     }

     // vektor x, co je vlastne nasa hladana alpha
     for (int i = 0; i < 5; i++)
     {
          printf("x[%d] = %.3f\n", i, x[i]);
     }
     double *u = calloc(N, sizeof(double)); // malloc, ale aj vynuluje vektor u, ktory bude konst

     for (int i = 0; i < N; i++)
     {
          for (int j = 0; j < N; j++)
          {
               double dx = coordinatesX[i][0] - coordinatesS[j][0];
               double dy = coordinatesX[i][1] - coordinatesS[j][1];
               double dz = coordinatesX[i][2] - coordinatesS[j][2];
               double rij = sqrt(dx * dx + dy * dy + dz * dz);

               if (rij != 0)
               {
                    double G_ij = 1.0 / rij;
                    u[i] += x[j] * G_ij;
               }
          }
     }

     // overenie hodnot u
     for (int i = 0; i < 5; i++)
     {
          printf("u[%d] = %.10f\n", i, u[i]);
     }

     FILE *output_file = fopen("output_u_BL.dat", "w");
     if (output_file == NULL)
     {
          printf("nepodarilo sa otvorit\n");
          return 1;
     }

     for (int i = 0; i < N; i++)
     {
          fprintf(output_file, "%.3f %.3f %.3f %.3f\n", B[i], L[i], u[i], x[i]);
     }

     fclose(output_file);

     // uvolnenie pamate
     free(coordinatesS);
     free(coordinatesX);
     free(coordinatesE);
     free(B);
     free(L);
     free(H);
     free(dg);
     free(f);
     free(dGMarray);
     free(b);
     free(x);
     free(r);
     free(r_hat);
     free(p);
     free(u);

     free(s);

     free(t);
     free(v);

     return 0;
}

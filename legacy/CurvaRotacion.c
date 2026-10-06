/*
 * Legacy educational C implementation of Metropolis–Hastings fitting
 * for a three-component galactic rotation curve.
 *
 * This is a cleaned-up descendant of the original homework code. The
 * modern, recommended pipeline lives in the Python package
 * ``galaxy_mcmc`` (see the repository README).
 *
 * Build:  cc CurvaRotacion.c -o CurvaRotacion.x -lm -std=c99
 * Run:    ./CurvaRotacion.x
 *
 * Expects ../data/RadialVelocities.dat (or RadialVelocities.dat in cwd).
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define MAX_LINE 128
#define MAX_POINTS 512
#define N_STEPS 20000
#define BURN_IN 5000

/* Fixed geometric scales [kpc] */
static const double B_b = 0.2497;
static const double B_d = 5.16;
static const double A_d = 0.3105;
static const double A_h = 64.3;
static const double SIGMA = 2.2; /* km/s observational noise */

static double R_exp[MAX_POINTS];
static double V_exp[MAX_POINTS];
static int n_data = 0;

static double urand(void) {
    return (double)rand() / ((double)RAND_MAX + 1.0);
}

/* Physically consistent circular velocity: sqrt(vb² + vd² + vh²) */
static double v_circ(double R, double Mb, double Md, double Mh) {
    double R2 = R * R;
    double vb2 = Mb * R2 / pow(R2 + B_b * B_b, 1.5);
    double vd2 = Md * R2 / pow(R2 + (B_d + A_d) * (B_d + A_d), 1.5);
    double vh2 = Mh / sqrt(R2 + A_h * A_h);
    double s = vb2 + vd2 + vh2;
    return (s > 0.0) ? sqrt(s) : 0.0;
}

/* Full-dataset Gaussian log-likelihood (up to an additive constant) */
static double log_likelihood(double Mb, double Md, double Mh) {
    double chi2 = 0.0;
    for (int i = 0; i < n_data; ++i) {
        double resid = (V_exp[i] - v_circ(R_exp[i], Mb, Md, Mh)) / SIGMA;
        chi2 += resid * resid;
    }
    return -0.5 * chi2;
}

static int load_data(const char *path) {
    FILE *fp = fopen(path, "r");
    if (!fp) {
        return 0;
    }
    char line[MAX_LINE];
    n_data = 0;
    while (fgets(line, MAX_LINE, fp) != NULL) {
        if (line[0] == '#' || line[0] == '\n') {
            continue;
        }
        double r, v;
        if (sscanf(line, "%lf %lf", &r, &v) == 2) {
            if (n_data >= MAX_POINTS) {
                fprintf(stderr, "Too many data points\n");
                fclose(fp);
                return 0;
            }
            R_exp[n_data] = r;
            V_exp[n_data] = v;
            n_data++;
        }
    }
    fclose(fp);
    return n_data > 0;
}

int main(void) {
    const char *candidates[] = {
        "../data/RadialVelocities.dat",
        "data/RadialVelocities.dat",
        "RadialVelocities.dat",
        NULL
    };
    int loaded = 0;
    for (int i = 0; candidates[i] != NULL; ++i) {
        if (load_data(candidates[i])) {
            printf("Loaded %d points from %s\n", n_data, candidates[i]);
            loaded = 1;
            break;
        }
    }
    if (!loaded) {
        fprintf(stderr, "Could not open RadialVelocities.dat\n");
        return 1;
    }

    srand((unsigned)time(NULL));

    /* Sample in log-mass space with a Gaussian random-walk proposal */
    double logMb = log(1.0);
    double logMd = log(1.4e4);
    double logMh = log(2.6e4);
    double lp = log_likelihood(exp(logMb), exp(logMd), exp(logMh));

    double sum_Mb = 0.0, sum_Md = 0.0, sum_Mh = 0.0;
    int accepted = 0;
    int kept = 0;

    FILE *chain = fopen("chain.dat", "w");
    if (!chain) {
        fprintf(stderr, "Cannot write chain.dat\n");
        return 1;
    }
    fprintf(chain, "# step Mb Md Mh logL\n");

    const double step_b = 0.8;
    const double step_d = 0.04;
    const double step_h = 0.04;

    for (int i = 0; i < N_STEPS; ++i) {
        double prop_b = logMb + step_b * (urand() * 2.0 - 1.0);
        double prop_d = logMd + step_d * (urand() * 2.0 - 1.0);
        double prop_h = logMh + step_h * (urand() * 2.0 - 1.0);

        double Mb = exp(prop_b);
        double Md = exp(prop_d);
        double Mh = exp(prop_h);

        /* Broad improper-flat prior on log M in a huge box */
        if (prop_b < log(1e-8) || prop_b > log(1e9) ||
            prop_d < log(1e-8) || prop_d > log(1e9) ||
            prop_h < log(1e-8) || prop_h > log(1e9)) {
            /* reject by leaving state unchanged */
        } else {
            double lp_prop = log_likelihood(Mb, Md, Mh);
            double log_alpha = lp_prop - lp;
            if (log(urand()) < log_alpha) {
                logMb = prop_b;
                logMd = prop_d;
                logMh = prop_h;
                lp = lp_prop;
                accepted++;
            }
        }

        if (i >= BURN_IN) {
            sum_Mb += exp(logMb);
            sum_Md += exp(logMd);
            sum_Mh += exp(logMh);
            kept++;
            fprintf(chain, "%d %.6e %.6e %.6e %.6f\n",
                    i, exp(logMb), exp(logMd), exp(logMh), lp);
        }
    }
    fclose(chain);

    double Mb_mean = sum_Mb / kept;
    double Md_mean = sum_Md / kept;
    double Mh_mean = sum_Mh / kept;

    printf("Acceptance rate: %.3f\n", (double)accepted / N_STEPS);
    printf("Posterior means (after burn-in):\n");
    printf("  Mb = %.6e\n", Mb_mean);
    printf("  Md = %.6e\n", Md_mean);
    printf("  Mh = %.6e\n", Mh_mean);

    FILE *fit = fopen("Ajuste.dat", "w");
    FILE *datos = fopen("Datos.dat", "w");
    if (!fit || !datos) {
        fprintf(stderr, "Cannot write output files\n");
        return 1;
    }
    for (int i = 0; i < n_data; ++i) {
        fprintf(datos, "%.8f %.8f\n", R_exp[i], V_exp[i]);
        fprintf(fit, "%.8f\n", v_circ(R_exp[i], Mb_mean, Md_mean, Mh_mean));
    }
    fclose(fit);
    fclose(datos);

    printf("Wrote Datos.dat, Ajuste.dat, chain.dat\n");
    return 0;
}

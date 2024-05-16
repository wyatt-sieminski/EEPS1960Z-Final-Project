#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>

#define TIME_STEP_ERUPTION 100.0
#define TIME_MAX 1e9

#define CONDUIT_HEIGHT 5000.0
#define CONDUIT_RADIUS 15.0
#define SHAPE_FACTOR 0.1
#define GRAVITY 9.81

#define MASS_FRACTION_H20 0.02
#define MASS_FRACTION_CRYSTAL 0.4

#define CRYSTAL_DENSITY 2600.0
#define MELT_DENSITY 2300.0
#define HOST_ROCK_DENSITY 2500.0

#define SIEVERTS_CONSTANT (4.0*pow(10.0, -6))
#define GAS_CONSTANT 462.0
#define SILICIC_MAGMA_TEMPERATURE 1073.15

#define BETA_R pow(10.0, 10.0)
#define DYNAMIC_VISCOSITY pow(10.0, 7.0)

struct eruption_parameters
{
    double t;
    double Qo;
    double rho;
    double V;
    double p;
    double n;
    double drho_dp;
    FILE* eruptions_ptr;
} typedef eruption_parameters;

// primary functions for moderating chamber processes over time
void chamber(void);
void eruption(eruption_parameters* e);

void forward_euler(eruption_parameters *e);


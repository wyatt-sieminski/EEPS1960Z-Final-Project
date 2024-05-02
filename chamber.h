#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>

#define TIME_STEP_ERUPTION 1000.0
#define TIME_STEP_DORMANT 1.0
#define TIME_MAX 3.154e8

#define CONDUIT_HEIGHT 5000.0
#define CONDUIT_RADIUS 15.0
#define MASS_FRACTION_H20 0.02
#define MASS_FRACTION_CRYSTAL 0.4
#define SHAPE_FACTOR 0.1
#define DYNAMIC_VISCOSITY 10e7
#define GRAVITY 9.81

#define CRYSTAL_DENSITY 2600.0
#define MELT_DENSITY 2300.0
#define HOST_ROCK_DENSITY 2500.0

#define SIEVERTS_CONSTANT 5.67e-8
#define GAS_CONSTANT 8.314
#define SILICIC_MAGMA_TEMPERATURE 800.0

#define BETA_R 10e10

struct eruption_parameters
{
    double t;
    double Qo;
    double rho;
    double V;
    double p;
    double n;
    double drho_dp;
    bool erupting;
    FILE* eruptions_ptr;
} typedef eruption_parameters;

// primary functions for moderating chamber processes over time
void chamber(void);


void eruption(eruption_parameters* e);


void dormant(eruption_parameters* e);


// function responsible for printing current values to appropriate output files
void log_values(eruption_parameters* e);


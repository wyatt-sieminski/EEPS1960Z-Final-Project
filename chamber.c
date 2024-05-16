#include <stdlib.h>
#include <stdio.h>

#include <stdbool.h>
#include <math.h>

#include "chamber.h"

int main(void)
{
    chamber();
}

void chamber(void)
{
    eruption_parameters *eruption_params = malloc(sizeof(eruption_parameters));

    double p_litho = HOST_ROCK_DENSITY * GRAVITY * CONDUIT_HEIGHT;

    // SET INITAL CONDITIONS HERE
    eruption_params->t = 0.0;
    eruption_params->V = 10.0e9;
    eruption_params->p = 10.0e6 + p_litho;

    // t = 0 values for n, rho, drho_dp, and Q0
    eruption_params->n = MASS_FRACTION_H20 - SIEVERTS_CONSTANT * sqrt(eruption_params->p) * (1.0 - MASS_FRACTION_CRYSTAL);

    // ensure n does not go below 0
    if (eruption_params->n < 0.0)
    {
        eruption_params->n = 0.0;
    }

    // t = 0 for rho
    eruption_params->rho = 1.0 / ((eruption_params->n / (eruption_params->p / (GAS_CONSTANT * SILICIC_MAGMA_TEMPERATURE))) + (1.0 - eruption_params->n) * ((MASS_FRACTION_CRYSTAL / CRYSTAL_DENSITY) + ((1.0 - MASS_FRACTION_CRYSTAL) / MELT_DENSITY)));

    // creating parts of the drho_dp equation
    double alpha = (MASS_FRACTION_CRYSTAL / CRYSTAL_DENSITY) + ((1.0 - MASS_FRACTION_CRYSTAL) / MELT_DENSITY);
    double numerator = eruption_params->n + ((1 - eruption_params->n) * alpha * (eruption_params->p / (GAS_CONSTANT * SILICIC_MAGMA_TEMPERATURE))) - ((MASS_FRACTION_H20 - eruption_params->n) / 2.0) * (1.0 - ((alpha * eruption_params->p) / (GAS_CONSTANT * SILICIC_MAGMA_TEMPERATURE))) + (1.0 - eruption_params->n) * (alpha / (GAS_CONSTANT * SILICIC_MAGMA_TEMPERATURE));
    double denominator = GAS_CONSTANT * SILICIC_MAGMA_TEMPERATURE * pow(eruption_params->n + (1.0 - eruption_params->n) * alpha * (eruption_params->p / (GAS_CONSTANT * SILICIC_MAGMA_TEMPERATURE)), 2.0);

    eruption_params->drho_dp = numerator / denominator;

    // calculations for Qo
    eruption_params->Qo = (eruption_params->p - p_litho) * ((eruption_params->rho * SHAPE_FACTOR * pow(M_PI * pow(CONDUIT_RADIUS, 2.0), 2.0)) / (CONDUIT_HEIGHT * DYNAMIC_VISCOSITY));

    // open file for writing model results
    eruption_params->eruptions_ptr = fopen("eruption_h20_testing.txt", "w");
 
    // print t = 0 values
    fprintf(eruption_params->eruptions_ptr, "%f, %f, %f, %f, %f, %f,\n", eruption_params->t, eruption_params->p, eruption_params->V, eruption_params->rho, eruption_params->n, eruption_params->drho_dp);

    eruption(eruption_params);

    free(eruption_params);
    fclose(eruption_params->eruptions_ptr);
}

void eruption(eruption_parameters *e)
{
    while (e->t < TIME_MAX)
    {
        forward_euler(e);

        e->t += TIME_STEP_ERUPTION;

        fprintf(e->eruptions_ptr, "%f, %f, %f, %f, %f, %f,\n", e->t, e->p, e->V, e->rho, e->n, e->Qo);
    }
}

void forward_euler(eruption_parameters *e)
{
    // calculate the pressure in the conduit
    double p_litho = HOST_ROCK_DENSITY * GRAVITY * CONDUIT_HEIGHT;

    // calculate the new mass flux
    e->Qo = (e->p - p_litho) * ((e->rho * SHAPE_FACTOR * pow(M_PI * pow(CONDUIT_RADIUS, 2.0), 2.0)) / (CONDUIT_HEIGHT * DYNAMIC_VISCOSITY));

    // calculate the new pressure
    double new_p = e->p - ((TIME_STEP_ERUPTION * e->Qo) / (e->rho * e->V)) * pow(((1.0 / BETA_R) + ((1.0 / e->rho) * e->drho_dp)), -1.0);

    // update the volume
    double new_V = e->V + e->V * ((new_p - e->p) / BETA_R);

    // update the mass_fraction_gas_bubbles
    double new_n = MASS_FRACTION_H20 - SIEVERTS_CONSTANT * sqrt(e->p) * (1.0 - MASS_FRACTION_CRYSTAL);

    if (new_n < 0.0)
    {
        new_n = 0.0;
    }

    // update the density
    double new_rho = pow((e->n / (e->p / (GAS_CONSTANT * SILICIC_MAGMA_TEMPERATURE))) + (1.0 - e->n) * ((MASS_FRACTION_CRYSTAL / CRYSTAL_DENSITY) + ((1.0 - MASS_FRACTION_CRYSTAL) / MELT_DENSITY)), -1.0);

    // creating parts of the drho_dp equation
    double alpha = (MASS_FRACTION_CRYSTAL / CRYSTAL_DENSITY) + ((1.0 - MASS_FRACTION_CRYSTAL) / MELT_DENSITY);
    double numerator = e->n + ((1 - e->n) * alpha * (e->p / (GAS_CONSTANT * SILICIC_MAGMA_TEMPERATURE))) - ((MASS_FRACTION_H20 - e->n) / 2.0) * (1.0 - ((alpha * e->p) / (GAS_CONSTANT * SILICIC_MAGMA_TEMPERATURE))) + ((1.0 - e->n) * (alpha / (GAS_CONSTANT * SILICIC_MAGMA_TEMPERATURE)));
    double denominator = GAS_CONSTANT * SILICIC_MAGMA_TEMPERATURE * pow(e->n + (1.0 - e->n) * alpha * (e->p / (GAS_CONSTANT * SILICIC_MAGMA_TEMPERATURE)), 2.0);

    // update the drho_dp
    e->drho_dp = numerator / denominator;

    e->p = new_p;
    e->V = new_V;
    e->n = new_n;
    e->rho = new_rho;
}

#include <stdlib.h>
#include <stdio.h>

#include <stdbool.h>
#include <math.h>

#include "chamber.h"
#include <Kernel/stdbool.h>

int main(void)
{
    chamber();
}

void chamber(void)
{
    eruption_parameters *eruption_params = malloc(sizeof(eruption_parameters));

    // initial conditions
    eruption_params->t = 0.0;
    eruption_params->V = 10.0e9;
    eruption_params->p = 1.32e8;
    eruption_params->erupting = true;

    // t = 0 values for n, rho, and drho_dp
    eruption_params->n = MASS_FRACTION_H20 - SIEVERTS_CONSTANT * sqrt(eruption_params->p) * (1.0 - MASS_FRACTION_CRYSTAL);

    eruption_params->rho = 1.0 / ((eruption_params->n / (eruption_params->p / (GAS_CONSTANT * SILICIC_MAGMA_TEMPERATURE))) + (1.0 - eruption_params->n) * ((MASS_FRACTION_CRYSTAL / CRYSTAL_DENSITY) + ((1.0 - MASS_FRACTION_CRYSTAL) / MELT_DENSITY)));

    double alpha = (MASS_FRACTION_CRYSTAL / CRYSTAL_DENSITY) + ((1.0 - MASS_FRACTION_CRYSTAL) / MELT_DENSITY);
    double numerator = eruption_params->n + ((1 - eruption_params->n) * alpha * (eruption_params->p / (GAS_CONSTANT * SILICIC_MAGMA_TEMPERATURE))) - ((MASS_FRACTION_H20 - eruption_params->n) / 2.0) * (1.0 - ((alpha * eruption_params->p) / (GAS_CONSTANT * SILICIC_MAGMA_TEMPERATURE))) + (1.0 - eruption_params->n) * (alpha / (GAS_CONSTANT * SILICIC_MAGMA_TEMPERATURE));
    double denominator = GAS_CONSTANT * SILICIC_MAGMA_TEMPERATURE * (eruption_params->n + pow(1.0 - eruption_params->n * alpha * (eruption_params->p / (GAS_CONSTANT * SILICIC_MAGMA_TEMPERATURE)), 2.0));

    eruption_params->drho_dp = numerator / denominator;

    // open file for writing model results
    eruption_params->eruptions_ptr = fopen("eruptions.txt", "w");
    fprintf(eruption_params->eruptions_ptr, "Time, Pressure, Volume, Density, Mass Fraction Gas Bubbles, Mass Erupted\n");

    while (eruption_params->t < TIME_MAX)
    {
        if (eruption_params->erupting)
        {
            eruption(eruption_params);
        }
        else
        {
            dormant(eruption_params);
        }
    }

    free(eruption_params);
    fclose(eruption_params->eruptions_ptr);
}

void forward_euler(eruption_parameters *e);

void eruption(eruption_parameters *e)
{
    double eruption_time_scale = TIME_MAX;
    double eruption_start_time = e->t;

    while (e->t < eruption_start_time + eruption_time_scale)
    {
        forward_euler(e);

        log_values(e);

        e->t += TIME_STEP_ERUPTION;
    }
}

void dormant(eruption_parameters *e)
{
    printf("%f", e->t);
    // TODO: implement
}

void log_values(eruption_parameters *e)
{
    fprintf(e->eruptions_ptr, "%f, %f, %f, %f, %f, %f\n", e->t, e->p, e->V, e->rho, e->n, e->Qo);
}

void forward_euler(eruption_parameters *e)
{
    // calculate the pressure in the conduit
    double p_litho = HOST_ROCK_DENSITY * GRAVITY * CONDUIT_HEIGHT;

    // calculate the mass flux
    e->Qo = (e->p - p_litho) * ((e->rho * SHAPE_FACTOR * pow(M_PI * pow(CONDUIT_RADIUS, 2.0), 2.0)) / (CONDUIT_HEIGHT * DYNAMIC_VISCOSITY));

    // calculate the new pressure
    double new_p = e->p - ((TIME_STEP_ERUPTION * e->Qo) / (e->rho * e->V)) * ((1.0 / BETA_R) + ((1.0 / e->rho) * e->drho_dp));

    // update the volume
    e->V += e->V*((new_p - e->p) / BETA_R);

    // update the mass_fraction_gas_bubbles
    e->n = MASS_FRACTION_H20 - SIEVERTS_CONSTANT * sqrt(e->p) * (1.0 - MASS_FRACTION_CRYSTAL);

    if (e->n < 0.0)
    {
        e->n = 0.0;
    }

    // update the density
    e->rho = pow(((e->n / (e->p / (GAS_CONSTANT * SILICIC_MAGMA_TEMPERATURE))) + (1.0 - e->n) * ((MASS_FRACTION_CRYSTAL / CRYSTAL_DENSITY) + ((1.0 - MASS_FRACTION_CRYSTAL) / MELT_DENSITY))), -1.0);

    // creating parts of the drho_dp equation
    double alpha = (MASS_FRACTION_CRYSTAL / CRYSTAL_DENSITY) + ((1.0 - MASS_FRACTION_CRYSTAL) / MELT_DENSITY);
    double numerator = e->n + ((1 - e->n) * alpha * (e->p / (GAS_CONSTANT * SILICIC_MAGMA_TEMPERATURE))) - ((MASS_FRACTION_H20 - e->n) / 2.0) * (1.0 - ((alpha * e->p) / (GAS_CONSTANT * SILICIC_MAGMA_TEMPERATURE))) + (1.0 - e->n) * (alpha / (GAS_CONSTANT * SILICIC_MAGMA_TEMPERATURE));
    double denominator = GAS_CONSTANT * SILICIC_MAGMA_TEMPERATURE * (e->n + pow(1.0 - e->n * alpha * (e->p / (GAS_CONSTANT * SILICIC_MAGMA_TEMPERATURE)), 2.0));

    // update the drho_dp
    e->drho_dp = numerator / denominator;

    e->p = new_p;
}

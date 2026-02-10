// copyright ############################### //
// This file is part of the Xtrack Package.  //
// Copyright (c) CERN, 2021.                 //
// ######################################### //

#ifndef XTRACK_ELENS_H
#define XTRACK_ELENS_H

#include "xtrack/headers/track.h"

GPUFUN
float T_Chebyshev(int n, float u){
  return cos(n * acos(u));
}

GPUFUN
float T_prime_Chebyshev(int n, float u){
  return n * sin(n * acos(u)) / sqrt(1 - pow((double)u, (double)2));
}

GPUFUN
void Elens_track_local_particle(ElensData el, LocalParticle* part0){

    double elens_length = ElensData_get_elens_length(el);

    if (LocalParticle_check_track_flag(part0, XS_FLAG_BACKTRACK)) {
        elens_length = -elens_length;
    }

    int const elens_on = ElensData_get_elens_on(el);

    double const inner_radius = ElensData_get_inner_radius(el);
    double const outer_radius = ElensData_get_outer_radius(el);
    double const offset_x = ElensData_get_offset_x(el);
    double const offset_y = ElensData_get_offset_y(el);
    double const current = ElensData_get_current(el);
    double const voltage = ElensData_get_voltage(el);
    double const residual_kick_x = ElensData_get_residual_kick_x(el);
    double const residual_kick_y = ElensData_get_residual_kick_y(el);

    int const polynomial_order = ElensData_get_polynomial_order(el);

    int const e_beam_dir = ElensData_get_e_beam_dir(el);

    GPUGLMEM double const* polynomial_coefficients =
                                ElensData_getp1_polynomial_coefficients(el, 0);

    int const chebyshev_max_order = ElensData_get_chebyshev_max_order(el);
    float const chebyshev_reference_radius = ElensData_get_chebyshev_reference_radius(el);

    GPUGLMEM double const* chebyshev_coefficients =
                                ElensData_getp1_chebyshev_coefficients(el, 0);

    START_PER_PARTICLE_BLOCK(part0, part);
        // electron mass
        double const EMASS  = 510998.928;
        // speed of light

        double x      = LocalParticle_get_x(part);
        double y      = LocalParticle_get_y(part);
        double xc     = x - offset_x;
        double yc     = y - offset_y;

        // delta
        // double delta  = LocalParticle_get_delta(part);
        // charge ratio: q/q0
        // double qratio = LocalParticle_get_charge_ratio(part);
        // chi = q/q0 * m0/m
        double const chi    = LocalParticle_get_chi(part);
        // reference particle momentum
        double const p0c    = LocalParticle_get_p0c(part);
        // particle momentum
        // double pc     = (1+delta)*(chi/qratio)*(p0c);
        // reference particle charge
        double const q0     = LocalParticle_get_q0(part);

        // rpp = P0/P
        double const rpp     = LocalParticle_get_rpp(part);


        // transverse radius
        double r      = sqrt(xc*xc + yc*yc);

        double rvv    = LocalParticle_get_rvv(part);
        double beta0  = LocalParticle_get_beta0(part);

        // # magnetic rigidity
        double const Brho0  = p0c/(q0*C_LIGHT);

        // # Electron properties
        // total electron energy
        double const etot_e       = voltage + EMASS;
        // // electron momentum
        double const p_e          = sqrt(etot_e*etot_e - EMASS*EMASS);
        // // relativistic beta of electron
        double const beta_e       = p_e/etot_e;
        //
        // // # relativistic beta  of protons
        double beta_p = rvv*beta0;

        // keep the formulas more compact
        double const r1 = inner_radius;
        double const r2 = outer_radius;

        // # geometric factor frr uniform distribution
        double frr = 0.;

        //
        if( r < r1 )
        {
          frr = 0.;
        }
        else if ( r > r2 )
        {
          frr = 1.;
        }
        else
        {
          // frr = ((r*r - r1*r1)/(r2*r2 - r1*r1));
          if (polynomial_order ==0)
          {
            frr = ((r*r - r1*r1)/(r2*r2 - r1*r1));
          }
          else
          {
            frr = 0;
            for(int i=0; i<(polynomial_order+1); ++i){
              frr += polynomial_coefficients[i]*(
                pow((double)r, (double)(i+2)) - pow((double)r1, (double)(i+2))
              ) / (i+2);
            }
          }

        }
        

        // # calculate the kick at r2 (maximum kick)
        double theta_max = ((1.0/(4.0*PI*EPSILON_0)));
        theta_max = theta_max*(2*elens_length*current);

        // theta_max depens on direction of e beam
        if (e_beam_dir < 0){
          theta_max = theta_max*(1+beta_e*beta_p);
        }
        else {
          theta_max = theta_max*(1-beta_e*beta_p);
        }

        theta_max = theta_max/(outer_radius*Brho0*beta_e*beta_p);
        theta_max = theta_max/(C_LIGHT*C_LIGHT);
        // theta max is now completed
        // theta_max = (-1)*theta_max/(rpp*chi);


        // now the actual kick the particle receives

        double theta_pxpy = 0.;
        double dpx = 0.;
        double dpy = 0.;
        //


        if ( r > r1 )
        {
          theta_pxpy = (-1)*frr*theta_max*(outer_radius/r)*(1/(rpp*chi));
          dpx        = xc*theta_pxpy/r;
          dpy        = yc*theta_pxpy/r;
        }
        else
        {
          // if the particle is not inside the e-beam, it will only
          // be subject to the residual kick
          if (chebyshev_max_order == 0) {
            dpx = residual_kick_x;
            dpy = residual_kick_y;
          }
          else {
            // Use Chebychev polynomials to evaluate non-linear residual kick 
            // https://lss.fnal.gov/archive/test-fn/0000/fermilab-fn-0972-apc.pdf

            double ux = x / chebyshev_reference_radius;
            double uy = y / chebyshev_reference_radius;

            int coeff_num = 0;
            for (int n=0; n<(chebyshev_max_order+1); ++n){
              for (int m=0; m<(n+1); ++m){
                dpx += chebyshev_coefficients[coeff_num] * T_prime_Chebyshev(m, ux) * T_Chebyshev(n-m, uy);
                dpy += chebyshev_coefficients[coeff_num] * T_Chebyshev(m, ux) * T_prime_Chebyshev(n-m, uy);
                coeff_num += 1;
              }
            }
            if (e_beam_dir < 0){
              dpx *= -(1+beta_e*beta_p) / (beta_p * C_LIGHT * Brho0 * chebyshev_reference_radius);
              dpy *= -(1+beta_e*beta_p) / (beta_p * C_LIGHT * Brho0 * chebyshev_reference_radius);
            }
            else {
              dpx *= -(1-beta_e*beta_p) / (beta_p * C_LIGHT * Brho0 * chebyshev_reference_radius);
              dpy *= -(1-beta_e*beta_p) / (beta_p * C_LIGHT * Brho0 * chebyshev_reference_radius);
            }
          }
        }

        dpx *= elens_on;
        dpy *= elens_on;

        LocalParticle_add_to_px(part, dpx );
        LocalParticle_add_to_py(part, dpy );

        // we can update the particle properties or add to the particle properties
        // LocalParticle_add_to_px(part, dpx);
        // LocalParticle_add_to_py(part, dpy);

        // LocalParticle_set_py(part, py_hat);
    END_PER_PARTICLE_BLOCK;
}

#endif
// copyright ############################### //
// This file is part of the Xtrack Package.  //
// Copyright (c) CERN, 2023.                 //
// ######################################### //

#ifndef XTRACK_TRACK_THICK_BEND_H
#define XTRACK_TRACK_THICK_BEND_H

#define POW2(X) ((X)*(X))
#define NONZERO(X) ((X) != 0.0)

/*gpufun*/
void track_thick_bend(
        LocalParticle* part,  // LocalParticle to track
        const double length,  // length of the element
        const double k,       // normal dipole strength
        const double h        // curvature
) {
    const double rvv = LocalParticle_get_rvv(part);

    // Particle coordinates
    const double x = LocalParticle_get_x(part);
    const double y = LocalParticle_get_y(part);
    const double px = LocalParticle_get_px(part);
    const double py = LocalParticle_get_py(part);
    const double ell = LocalParticle_get_zeta(part) / rvv;
    const double pt = LocalParticle_get_ptau(part);
    const double beta0 = LocalParticle_get_beta0(part);
    const double t = LocalParticle_get_zeta(part) / beta0;
    const double s = length;

    double new_x, new_px, new_y, new_ell;

    const double rho = 1 / h;
    const double ang = length / rho;
    const double sa = sin(ang);
    const double ca = cos(ang);

    const double pw2 = 1 + 2*pt/beta0 + POW2(pt) - POW2(py);
    const double pz = sqrt(pw2 - POW2(px));
    const double pzx = pz - k*(rho+x);
    const double npx = sa*pzx + ca*px;
    const double dpx = ca*pzx - sa*px;
    const double pzs = sqrt(pw2 - POW2(npx));
    const double _ptt = 1.0 / sqrt(pw2);
    const double dxs = (ang + asin(px*_ptt) - asin(npx*_ptt))/k;
    new_x = (pzs - dpx)/k - rho;
    new_px = npx;
    new_y = y + dxs*py;
    // new_t = 0;//t - dxs*(1/beta+pt) + (1-T)*(ld/beta)

    // Update Particles object
    LocalParticle_set_x(part, new_x);
    LocalParticle_set_px(part, new_px);
    LocalParticle_set_y(part, new_y);
    //LocalParticle_set_zeta(part, new_ell * rvv);
    LocalParticle_add_to_s(part, s);
}

#endif // XTRACK_TRACK_THICK_BEND_H
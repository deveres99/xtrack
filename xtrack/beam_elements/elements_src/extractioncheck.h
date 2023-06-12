// copyright ############################### //
// This file is part of the Xtrack Package.  //
// Copyright (c) CERN, 2021.                 //
// ######################################### //

#ifndef XTRACK_EXTRACTIONCHECK_H
#define XTRACK_EXTRACTIONCHECK_H


/*gpufun*/
void Extractioncheck_track_local_particle(ExtractioncheckData el, LocalParticle* part0){

    double const x_min = ExtractioncheckData_get_x_min(el);
    double const x_max = ExtractioncheckData_get_x_max(el);
    double const px_min = ExtractioncheckData_get_px_min(el);
    double const px_max = ExtractioncheckData_get_px_max(el);

    //start_per_particle_block (part0->part)

        double x = LocalParticle_get_x(part);
        double px = LocalParticle_get_px(part);

        if ((x_min <= x) && (x_max >= x) && (px_min <= px) && (px_max >= px))
        {
            LocalParticle_set_state(part, -1000);
        }

    //end_per_particle_block

}


#endif /* XTRACK_EXTRACTIONCHECK_H */

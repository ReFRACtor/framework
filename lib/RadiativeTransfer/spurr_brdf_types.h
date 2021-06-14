#ifndef SPURR_BRDF_TYPES_H
#define SPURR_BRDF_TYPES_H

namespace FullPhysics {

    /// These match numbers used internal to the Fortran code for Spurr based RTs
    enum SpurrBrdfType {
        LAMBERTIAN  = 1,
        ROSSTHIN    = 2,
        ROSSTHICK   = 3,
        LISPARSE    = 4,
        LIDENSE     = 5,
        HAPKE       = 6,
        ROUJEAN     = 7,
        RAHMAN      = 8,
        COXMUNK     = 9,
        BPDFSOIL    = 10,
        BPDFVEGN    = 11,
        BRDFNDVI    = 12,
        NEWCMGLINT  = 13,
        EMISSIVITY  = 100, // Not a real index from LIDORT family, but
			   // included to meet interface
	UNKNOWN = -999	   // Type that will trigger an error, if we
			   // have a ground type that doesn't map to
			   // the spurr BRDF
    };
 
}

#endif

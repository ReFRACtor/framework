module first_order_debug_output_m

use iso_c_binding

implicit none

contains

! Provides the same interface as fo_dtgeometry_master_m_fo_dtgeometry_master_wrap from first_order_interface.F90
! Modify first_order_interface.h to call this before calling that routine to enable debugging.

subroutine first_order_debug_output(FO_maxgeoms, &
                                    FO_maxlayers, &
                                    FO_maxfine, &
                                    FO_do_planpar, &
                                    FO_do_regular_ps, &
                                    FO_do_enhanced_ps, &
                                    FO_donadir, &
                                    FO_do_sleave, &
                                    FO_ngeoms, &
                                    FO_nlayers, &
                                    FO_nfinedivs, &
                                    FO_aclevel, &
                                    FO_reflec, &
                                    FO_slterm, &
                                    FO_extincs, &
                                    FO_deltaus, &
                                    FO_exactscat, &
                                    FO_flux, &
                                    FO_mu0, &
                                    FO_mu1, &
                                    FO_ncrit, &
                                    FO_xfine, &
                                    FO_wfine, &
                                    FO_csqfine, &
                                    FO_cotfine, &
                                    FO_raycon, &
                                    FO_cota, &
                                    FO_sunpaths, &
                                    FO_ntraverse, &
                                    FO_sunpaths_fine, &
                                    FO_ntraverse_fine, &
                                    FO_intensity_up, &
                                    FO_intensity_db, &
                                    FO_cumsource_up) bind(C)

    ! Arguments
    integer(c_int), intent(in) :: FO_maxgeoms
    integer(c_int), intent(in) :: FO_maxlayers
    integer(c_int), intent(in) :: FO_maxfine
    logical(c_bool), intent(in) :: FO_do_planpar
    logical(c_bool), intent(in) :: FO_do_regular_ps
    logical(c_bool), intent(in) :: FO_do_enhanced_ps
    logical(c_bool), dimension(FO_MAXGEOMS), intent(in) :: FO_donadir
    logical(c_bool), intent(in) :: FO_do_sleave
    integer(c_int), intent(in) :: FO_ngeoms
    integer(c_int), intent(in) :: FO_nlayers
    integer(c_int), dimension(FO_MAXLAYERS, FO_MAXGEOMS), intent(in) :: FO_nfinedivs
    integer(c_int), intent(in) :: FO_aclevel
    real(c_double), dimension(FO_MAXGEOMS), intent(in) :: FO_reflec
    real(c_double), dimension(FO_MAXGEOMS), intent(in) :: FO_slterm
    real(c_double), dimension(FO_MAXLAYERS), intent(in) :: FO_extincs
    real(c_double), dimension(FO_MAXLAYERS), intent(in) :: FO_deltaus
    real(c_double), dimension(FO_MAXLAYERS, FO_MAXGEOMS), intent(in) :: FO_exactscat
    real(c_double), intent(in) :: FO_flux
    real(c_double), dimension(FO_MAXGEOMS), intent(in) :: FO_mu0
    real(c_double), dimension(FO_MAXGEOMS), intent(in) :: FO_mu1
    integer(c_int), dimension(FO_MAXGEOMS), intent(in) :: FO_ncrit
    real(c_double), dimension(FO_MAXLAYERS, FO_MAXFINE, FO_MAXGEOMS), intent(in) :: FO_xfine
    real(c_double), dimension(FO_MAXLAYERS, FO_MAXFINE, FO_MAXGEOMS), intent(in) :: FO_wfine
    real(c_double), dimension(FO_MAXLAYERS, FO_MAXFINE, FO_MAXGEOMS), intent(in) :: FO_csqfine
    real(c_double), dimension(FO_MAXLAYERS, FO_MAXFINE, FO_MAXGEOMS), intent(in) :: FO_cotfine
    real(c_double), dimension(FO_MAXGEOMS), intent(in) :: FO_raycon
    real(c_double), dimension(0:FO_MAXLAYERS, FO_MAXGEOMS), intent(in) :: FO_cota
    real(c_double), dimension(0:FO_MAXLAYERS, FO_MAXLAYERS, FO_MAXGEOMS), intent(in) :: FO_sunpaths
    integer(c_int), dimension(0:FO_MAXLAYERS, FO_MAXGEOMS), intent(in) :: FO_ntraverse
    real(c_double), dimension(FO_MAXLAYERS, FO_MAXLAYERS, FO_MAXFINE, FO_MAXGEOMS), intent(in) :: FO_sunpaths_fine
    integer(c_int), dimension(FO_MAXLAYERS, FO_MAXFINE, FO_MAXGEOMS), intent(in) :: FO_ntraverse_fine
    real(c_double), dimension(FO_MAXGEOMS), intent(out) :: FO_intensity_up
    real(c_double), dimension(FO_MAXGEOMS), intent(out) :: FO_intensity_db
    real(c_double), dimension(0:FO_MAXLAYERS, FO_MAXGEOMS), intent(out) :: FO_cumsource_up

    ! Index variables
    integer :: lay, lay2, fine, geom

    ! nfine used equals maxfine in our usage of FO
    integer :: FO_nfineinput
    FO_nfineinput = FO_maxfine

    open(38,file='FO_WRITE_STD_INPUT.dbg',status='replace')

    WRITE(38,*)
    WRITE(38,'(A)') '-----------------------'
    WRITE(38,'(A)') ' Dimensions/Flags/Ctrl '
    WRITE(38,'(A)') '-----------------------'

    WRITE(38,*)
    WRITE(38,*)  'FO_maxgeoms            = ',FO_maxgeoms
    WRITE(38,*)  'FO_maxlayers           = ',FO_maxlayers
    WRITE(38,*)  'FO_maxfine             = ',FO_maxfine
    WRITE(38,*)  'FO_do_planpar          = ',FO_do_planpar
    WRITE(38,*)  'FO_do_regular_ps       = ',FO_do_regular_ps
    WRITE(38,*)  'FO_do_enhanced_ps      = ',FO_do_enhanced_ps
    WRITE(38,*)  'FO_doNadir             = ',FO_doNadir(1:FO_ngeoms)
    WRITE(38,*)  'FO_ngeoms              = ',FO_ngeoms
    WRITE(38,*)  'FO_nlayers             = ',FO_nlayers
    do lay = 1, FO_nlayers
       do geom = 1, FO_ngeoms
          WRITE(38,*)  ' LAY = ',LAY,' GEOM = ',GEOM,&
          ' FO_nfinedivs(LAY,GEOM) = ',&
            FO_nfinedivs(LAY,GEOM)
       enddo
    enddo

    WRITE(38,*)
    WRITE(38,'(A)') '-----------------------'
    WRITE(38,'(A)') '     Optical Inputs    '
    WRITE(38,'(A)') '-----------------------'
    WRITE(38,*)
    WRITE(38,*)  'FO_reflec              = ',FO_reflec(1:FO_ngeoms)
    WRITE(38,*)
    DO LAY=1,FO_NLAYERS
      WRITE(38,*)  'LAY = ',LAY,&
        ' FO_extincs(LAY)      = ',FO_extincs(LAY)
    END DO
    WRITE(38,*)
    DO LAY=1,FO_NLAYERS
      WRITE(38,*)  'LAY = ',LAY,&
        ' FO_deltaus(LAY)      = ',FO_deltaus(LAY)
    END DO
    WRITE(38,*)
    do lay = 1, FO_nlayers
       do geom = 1, FO_ngeoms
          WRITE(38,*)  ' LAY = ',LAY,' GEOM = ',GEOM,&
          ' FO_exactscat(LAY,GEOM) = ',&
            FO_exactscat(LAY,GEOM)
       enddo
    enddo
    WRITE(38,*)
    WRITE(38,*)  'FO_flux                = ',FO_flux


    WRITE(38,*)
    WRITE(38,'(A)') '-----------------------'
    WRITE(38,'(A)') '    Geometry Inputs    '
    WRITE(38,'(A)') '-----------------------'
    WRITE(38,*)
    WRITE(38,*)  'FO_Mu0            = ',FO_Mu0(1:FO_ngeoms)
    WRITE(38,*)  'FO_Mu1            = ',FO_Mu1(1:FO_ngeoms)
    WRITE(38,*)  'FO_NCrit          = ',FO_NCrit(1:FO_ngeoms)
    WRITE(38,*)
    do lay = 1, FO_nlayers
       do fine = 1, FO_nfineinput
          WRITE(38,*)  ' LAY = ',LAY,' FINE = ',FINE,&
          ' FO_xfine(LAY,FINE,1:FO_NGEOMS) = ',&
            FO_xfine(LAY,FINE,1:FO_NGEOMS)
       enddo
    enddo
    WRITE(38,*)
    do lay = 1, FO_nlayers
       do fine = 1, FO_nfineinput
          WRITE(38,*)  ' LAY = ',LAY,' FINE = ',FINE,&
          ' FO_wfine(LAY,FINE,1:FO_NGEOMS) = ',&
            FO_wfine(LAY,FINE,1:FO_NGEOMS)
       enddo
    enddo
    WRITE(38,*)
    do lay = 1, FO_nlayers
       do fine = 1, FO_nfineinput
          WRITE(38,*)  ' LAY = ',LAY,' FINE = ',FINE,&
          ' FO_csqfine(LAY,FINE,1:FO_NGEOMS) = ',&
            FO_csqfine(LAY,FINE,1:FO_NGEOMS)
       enddo
    enddo
    WRITE(38,*)
    do lay = 1, FO_nlayers
       do fine = 1, FO_nfineinput
          WRITE(38,*)  ' LAY = ',LAY,' FINE = ',FINE,&
          ' FO_cotfine(LAY,FINE,1:FO_NGEOMS) = ',&
            FO_cotfine(LAY,FINE,1:FO_NGEOMS)
       enddo
    enddo
    WRITE(38,*)
    WRITE(38,*)  'FO_Raycon         = ',FO_Raycon(1:FO_ngeoms)
    WRITE(38,*)
    do lay = 0, FO_nlayers
       do geom = 1, FO_ngeoms
          WRITE(38,*)  ' LAY = ',LAY,' GEOM = ',GEOM,&
          ' FO_cota(LAY,GEOM) = ',&
            FO_cota(LAY,GEOM)
       enddo
    enddo
    WRITE(38,*)
    do lay = 0, FO_nlayers
       do lay2 = 1, FO_nlayers
          WRITE(38,*)  ' LAY = ',LAY,' LAY2 = ',LAY2,&
          ' FO_sunpaths(LAY,LAY2,1:FO_NGEOMS) = ',&
            FO_sunpaths(LAY,LAY2,1:FO_NGEOMS)
       enddo
    enddo
    WRITE(38,*)
    do lay = 0, FO_nlayers
       do geom = 1, FO_ngeoms
          WRITE(38,*)  ' LAY = ',LAY,' GEOM = ',GEOM,&
          ' FO_ntraverse(LAY,GEOM) = ',&
            FO_ntraverse(LAY,GEOM)
       enddo
    enddo
    WRITE(38,*)
    do lay = 1, FO_nlayers
       do lay2 = 1, FO_nlayers
          do fine = 1, FO_nfineinput
         WRITE(38,*)  ' LAY = ',LAY,' LAY2 = ',LAY2,' FINE = ',FINE, &
          ' FO_sunpaths_fine(LAY,LAY2,FINE,1:FO_NGEOMS) = ',&
            FO_sunpaths_fine(LAY,LAY2,FINE,1:FO_NGEOMS)
          enddo
       enddo
    enddo
    WRITE(38,*)
    do lay = 1, FO_nlayers
       do fine = 1, FO_nfineinput
         WRITE(38,*)  ' LAY = ',LAY,' FINE = ',FINE, &
          ' FO_ntraverse_fine(LAY,FINE,1:FO_NGEOMS) = ',&
            FO_ntraverse_fine(LAY,FINE,1:FO_NGEOMS)
       enddo
    enddo

    close(38)

end subroutine first_order_debug_output

end module first_order_debug_output_m

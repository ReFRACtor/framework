
! ###############################################################
! #                                                             #
! #                       VLIDORT_2p8p3                         #
! #                                                             #
! #  Vectorized LInearized Discrete Ordinate Radiative Transfer #
! #  -          --         -        -        -         -        #
! #                                                             #
! ###############################################################

! ###############################################################
! #                                                             #
! #  Authors :     Robert. J. D. Spurr (1)                      #
! #                Matt Christi                                 #
! #                                                             #
! #  Address (1) : RT Solutions, inc.                           #
! #                9 Channing Street                            #
! #                Cambridge, MA 02138, USA                     #
! #                                                             #
! #  Tel:          (617) 492 1183                               #
! #  Email :       rtsolutions@verizon.net                      #
! #                                                             #
! #  This Version :   VLIDORT_2p8p3                             #
! #  Release Date :   31 March 2021                             #
! #                                                             #
! #  Previous VLIDORT Versions under Standard GPL 3.0:          #
! #  ------------------------------------------------           #
! #                                                             #
! #      2.7   F90, released        August 2014                 #
! #      2.8   F90, released        May    2017                 #
! #      2.8.1 F90, released        August 2019                 # 
! #      2.8.2 F90, limited release May    2020                 # 
! #                                                             #
! #  Features Summary of Recent VLIDORT Versions:               #
! #  -------------------------------------------                #
! #                                                             #
! #      NEW: TOTAL COLUMN JACOBIANS         (2.4)              #
! #      NEW: BPDF Land-surface KERNELS      (2.4R)             #
! #      NEW: Thermal Emission Treatment     (2.4RT)            #
! #      Consolidated BRDF treatment         (2.4RTC)           #
! #      f77/f90 Release                     (2.5)              #
! #      External SS / New I/O Structures    (2.6)              #
! #                                                             #
! #      SURFACE-LEAVING / BRDF-SCALING      (2.7)              #
! #      TAYLOR Series / OMP THREADSAFE      (2.7)              #
! #      New Water-Leaving Treatment         (2.8)              #
! #      LBBF & BRDF-Telescoping, enabled    (2.8)              #
! #      Several Performance Enhancements    (2.8)              #
! #      Water-leaving coupled code          (2.8.1)            #
! #      Planetary problem, media properties (2.8.1)            #
! #      Doublet geometry post-processing    (2.8.2)            #
! #      Reduction zeroing, dynamic memory   (2.8.2)            #
! #                                                             #
! #  Features Summary of This VLIDORT Version                   #
! #  ----------------------------------------                   #
! #                                                             #
! #   2.8.3, released 31 March 2021.                            #
! #     ==> Green's function RT solutions (Nstokes = 1 or 3)    #
! #     ==> Sphericity Corrections using MS source terms        #
! #     ==> BRDF upgrades, including new snow reflectance       #
! #     ==> SLEAVE Upgrades, extended water-leaving treatment   #
! #                                                             #
! ###############################################################

! ###################################################################
! #                                                                 #
! # This is Version 2.8.3 of the VLIDORT_2p8 software library.      #
! # This library comes with the Standard GNU General Public License,#
! # Version 3.0, 29 June 2007. Please read this license carefully.  #
! #                                                                 #
! #      VLIDORT Copyright (c) 2003-2021.                           #
! #          Robert Spurr, RT Solutions, Inc.                       #
! #          9 Channing Street, Cambridge, MA 02138, USA.           #
! #                                                                 #
! # This file is part of VLIDORT_2p8p3 ( Version 2.8.3 )            #
! #                                                                 #
! # VLIDORT_2p8p3 is free software: you can redistribute it         #
! # and/or modify it under the terms of the Standard GNU GPL        #
! # (General Public License) as published by the Free Software      #
! # Foundation, either version 3.0 of the License, or any           #
! # later version.                                                  #
! #                                                                 #
! # VLIDORT_2p8p3 is distributed in the hope that it will be        #
! # useful, but WITHOUT ANY WARRANTY; without even the implied      #
! # warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR         #
! # PURPOSE. See the Standard GNU General Public License (GPL)      #
! # for more details.                                               #
! #                                                                 #
! # You should have received a copy of the Standard GNU General     #
! # Public License (GPL) Version 3.0, along with the VLIDORT_2p8p3  #
! # code package. If not, see <http://www.gnu.org/licenses/>.       #
! #                                                                 #
! ###################################################################

! ###############################################################
! #                                                             #
! #  This is vlidort_findpar.f. Utility routines (from VLIDORT) #
! #                                                             #
! #      xfindpar:              J. Lavagnino, SAO, 1991         #
! #      gfindpar:              J. Lavagnino, SAO, 1991         #
! #      FINDPAR_ERROR          R. Spurr, 2001                  #
! #      LCSTRING               Lavagnino                       #
! #                                                             #
! ###############################################################

!  1/31/21, Version 2.8.3. No changes here, other than removing initial vlidort_pars_m

      MODULE vbrdf_findpar_m

      private
      public :: gfindpar, findpar_error

      CONTAINS

!

!  FINDPAR input subroutines.

!----- The basic idea

!  FINDPAR is a poor man's RDPAR.  It implements some of the
! central features of RDPAR --- the ability to find a named set
! of data values in an input file, the relaxation of the need for
! data in the input file to appear in the order in which it's
! used, the provision for comments in the input file --- in
! (almost) standard FORTRAN 77.

!  Among the features of RDPAR that are not implemented in FINDPAR:
! all redirection of input (from one parameter name to another, from
! file to terminal, from file to command line); freer-format data
! specification on input; tracing options.

!----- Using FINDPAR

!  To use FINDPAR, you open your input file using ordinary FORTRAN
! procedures.  You will also use ordinary FORTRAN to read data
! values.  The only difference is that you can use the FINDPAR call
! at any point to position the file to a point after a line containing
! a parameter name.  These parameter names help to document what the
! input values are for; they also don't have to appear in the order in
! which the program uses them.  An input file might look like this:

! ! Parameters for July 16 run.

! NLAYERS
! 5

! LAYER PARAMETERS
! 5,1013,296
! .001,.001,.001,.001,.001,.001,.001,.001

! MOLECULE FLAGS
! t,t,t,t,t,t,t,t

!  By convention, ! or a space introduces a comment line.  A
! parameter name appears on a line of its own, followed by the data
! lines to be read using FORTRAN statements.  Here is a sample of
! FORTRAN code to read this parameter file:

!      OPEN (UNIT = IN_UNIT, FILE = INPUT, STATUS = 'OLD')
!      CALL FINDPAR ( IN_UNIT, 'NLAYERS' )
!      READ (IN_UNIT, *) NLAYERS
!      CALL FINDPAR ( IN_UNIT, 'LAYER PARAMETERS' )
!      READ (IN_UNIT, *) Z, PRESS, TEMP
!      READ (IN_UNIT, *) BRO, CLO, HCHO, NO2, O2, O3, OCLO, SO2
!      CALL FINDPAR ( IN_UNIT, 'MOLECULE FLAGS' )
!      READ (IN_UNIT, *) IF_BRO, IF_CLO, IF_HCHO, IF_NO2, IF_O2, IF_O3,
!     $ IF_OCLO, IF_SO2
!      CLOSE (UNIT = IN_UNIT)

!----- Details of parameter file format

! The parameter file consists of a sequence of comments and
! parameter-name/data chunks, in any order (except as noted below).

! Comments are essentially anything that doesn't happen to get matched
! as a parameter name or read as data.  However, by convention comments
! are introduced by an exclamation point, to make absolutely sure that
! they don't get mistaken for parameter names.  It is also conventional
! to ``comment out'' a parameter name that you don't want to match by
! simply indenting it one space.

! Parameter names can theoretically be almost any sequence of
! characters.  However, FINDPAR ignores trailing spaces both in the
! parameter name you give it and in lines it reads from the parameter
! file.  And, by convention, parameter names don't begin with spaces.
! Parameter names may contain embedded spaces.  Parameter name matching
! is case-insensitive but otherwise exact: if you use odd things
! like control characters in a parameter name they have to match
! exactly, too, so it's usually best to avoid control characters and
! tabs in parameter names.

! The data lines that follow a line with a parameter name are ordinary
! FORTRAN input lines.  FINDPAR has no influence on how they're handled,
! and all the ordinary rules about FORTRAN READ statements apply.

! Comments may not appear on the same lines as parameter names or data,
! and they can't appear within parameter name/data chunks.  When FINDPAR
! locates a parameter name it positions the file at the line following
! that name; your program will then use ordinary FORTRAN READ statements
! to read the data, and will probably be confused by comments that appea
! at that point.

! There's a maximum length for input lines and parameter names, set
! by the parameter MAXLINE in the code.

! By convention, a parameter name that ends in ? precedes an
! input line that contains a single logical variable (T or F).

!----- Repeated parameters

! Sometimes you want to read several parameters with the same name:
! for example, parameters for each layer of the atmosphere, where there
! might be a large number of layers.  In this case the order of
! parameters in the file does matter, since you usually just want to
! put the numbers in the file in some natural order rather than giving
! a number in each input line specifying where it appears in the input
! order.

! The normal FINDPAR approach is to organize these parameters as
! follows: begin with some parameter that appears only once in the file
! (the number of layers, for example); then follow it with the instances
! of the repeated parameter and associated data.  The only-once
! parameter can even be a dummy, with no data following; its importance
! is that reading it gets the parameter file positioned at the right
! point.

! The other problem here is knowing when a list of repeated parameters
! is over.  FINDPAR doesn't have any of the tricks for doing this that
! RDPAR has; you must either know exactly how many instances of the
! repeated parameter there are going to be (by reading some once-only
! parameter that tells you), or else you need some special numbers to
! flag the end of the list (zeros, negative numbers, etc.).  You can't
! just keep looking for the same parameter name indefinitely, because
! FINDPAR will just rewind the file for you and loop through its
! contents forever.

!----- Errors and XFINDPAR

! If FINDPAR can't find a parameter, it prints an error message
! and returns.  Your READ statements following the FINDPAR call are
! very likely to run into errors of their own at this point, since
! you'll be positioned at some random point in the file (usually at
! the beginning, in this version, but that isn't true in all cases).

! If you want to do something more intelligent about such errors,
! you can use XFINDPAR, which is a function returning a logical value:
! .true. if the parameter name was found, .false. if not; XFINDPAR
! doesn't display FINDPAR's error message when a parameter name cannot
! be found.

! Both FINDPAR and XFINDPAR can run into ordinary FORTRAN input errors
! as they read the file: no special action is taken on these---the
! system's default action, whatever that is, occurs.

!       9/10/91         John Lavagnino

!  Find a parameter name in a parameter file.

!        SUBROUTINE FINDPAR ( NUNIT, PARNAME )

!  Input arguments.

!        INTEGER         NUNIT
!        CHARACTER * (*) PARNAME

!  Local variables.

!        LOGICAL         XFINDPAR
!        EXTERNAL        XFINDPAR

!  No local variables need to be SAVEd.

!***********************************************************************

!        IF ( .NOT. XFINDPAR ( NUNIT, PARNAME ) ) THEN
!          WRITE (*, *) 'FINDPAR error'
!          WRITE (*, *) '   Unit number ', NUNIT
!          WRITE (*, *) '   Parameter name ', PARNAME
!        END IF

!        END
!
!  Find a parameter name in a paramter file: return .true. if found,
! .false. if not.

      LOGICAL FUNCTION XFINDPAR ( NUNIT, PARNAME )

        implicit none

!  Parameters.

        INTEGER         MAXLINE
        PARAMETER       ( MAXLINE = 132 )

!  Input arguments.

        INTEGER         NUNIT
        CHARACTER * (*) PARNAME

!  Local variables.

        INTEGER                 III
        INTEGER                 PARLEN
        INTEGER                 START
        INTEGER                 END
        LOGICAL                 ENDSEEN
        CHARACTER * (MAXLINE)   LINE
        CHARACTER * (MAXLINE)   NAME

!  No local variables need to be SAVEd.

!***********************************************************************

!  Determine the length of the parameter name.

        DO 10 III = LEN ( PARNAME ), 1, -1
          IF ( PARNAME ( III : III ) .NE. ' ' ) THEN
            PARLEN = III
            GO TO 20
          END IF
10      CONTINUE

!  If we get to here, then name contains nothing but blanks.
!  We just return, claiming success, in such a case.

        GO TO 500

!  If we get here, then there's a non-null parameter name; but it
!  might still be too long, in which case we always return
!  signaling failure, since we couldn't ever find such a name in the
!  file.

20      CONTINUE
        IF ( PARLEN  .GT.  MAXLINE ) GO TO 400

!  Convert the name to lower-case.

        NAME = PARNAME ( 1 : PARLEN )
        CALL LCSTRING ( NAME ( 1 : PARLEN ) )

!  Top of main loop.

        ENDSEEN = .FALSE.

100     CONTINUE

          LINE = ' '
          READ ( UNIT = NUNIT, FMT = '(A)', END = 200 ) LINE
          CALL LCSTRING ( LINE ( 1 : PARLEN ) )
          IF ( LINE ( 1 : PARLEN ) .NE. NAME ( 1 : PARLEN ) ) &
                GO TO 100

          START = PARLEN + 1
          END = MAXLINE
          IF ( START  .GT.  END ) GO TO 500
          IF ( LINE ( START : END ) .EQ. ' ' ) GO TO 500
          GO TO 100

!  End-of-file branch.

200       CONTINUE
          REWIND ( UNIT = NUNIT )
          IF ( ENDSEEN ) THEN
            GO TO 400
          ELSE
            ENDSEEN = .TRUE.
            GO TO 100
          END IF

!  End of loop: failure, no parameter name found.

400     CONTINUE
        XFINDPAR = .FALSE.
        RETURN

!  End of loop: successful location of parameter name.

500     CONTINUE
        XFINDPAR = .TRUE.
        RETURN

      END FUNCTION XFINDPAR

!

      LOGICAL FUNCTION GFINDPAR (NUNIT, PREFIX, &
                                          ERROR, PARNAME)

!  Find a parameter name in a parameter file, with optional added
! prefix on parameter name.  Error messages are displayed but
! instead of stopping the program a logical flag is returned.
!  Returns .TRUE. if some form of the parameter name was found;
! .FALSE. if not.
! ERROR and GFINDPAR are initialised for a successful search.
! ERROR is just the negation of GFINDPAR. Error messages are removed.

        implicit none

!  Input arguments.  If PREFIX is something other than a bunch of
! blanks, it will be added onto PARNAME, with a blank to separate
! it, for the FINDPAR call.  If that FINDPAR call fails, then
! PARNAME without prefix is tried instead, and only if that call fails
! is there an error message.  This means that unprefixed versions of
! the parameter name act as defaults.
!  If PREFIX is a bunch of blanks, this works much as FINDPAR does.
!  (This is assuming that PREFIX has no trailing blanks.  Probably
! it should be checking for them and trimming them if necessary.)

        INTEGER         NUNIT
        CHARACTER * (*) PREFIX
        CHARACTER * (*) PARNAME

!  Modified argument.  If all goes well, this is not changed.  If
! the parameter name can't be found, it's set to .TRUE.

        LOGICAL         ERROR

!  Local variables.

        CHARACTER * 132 LINE

!        LOGICAL         XFINDPAR
!        EXTERNAL        XFINDPAR

!***********************************************************************

!  Initialise return values.

        GFINDPAR = .TRUE.
        ERROR    = .FALSE.

!  Compose search line for first search.

        IF ( PREFIX  .EQ.  ' ' ) THEN
          LINE = PARNAME
        ELSE
          LINE = PREFIX // ' ' // PARNAME
        END IF

!  Try the parameter name with prefix.

        IF (XFINDPAR (NUNIT, LINE)) THEN
          RETURN
        END IF

!  If that doesn't work, try it without the prefix, if the prefix
! was non-null.

        IF (PREFIX  .NE.  ' ') THEN
          IF (XFINDPAR (NUNIT, PARNAME)) THEN
            RETURN
          END IF
        END IF

!  If we got here, we just couldn't find the parameter name.

        GFINDPAR = .FALSE.
        ERROR    = .TRUE.

      END FUNCTION GFINDPAR

!

      SUBROUTINE FINDPAR_ERROR &
           (  ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  include file

      USE vlidort_pars_m, only : MAX_MESSAGES, VLIDORT_SERIOUS

      implicit none

!  subroutine input arguments

      LOGICAL       ERROR
      CHARACTER*(*) PAR_STR

!  subroutine Output or In/Out arguments

      INTEGER       STATUS
      INTEGER       NM
      CHARACTER*(*) MESSAGES(0:MAX_MESSAGES), ACTIONS(0:MAX_MESSAGES)

      IF ( ERROR ) THEN
        NM = NM + 1
        STATUS = VLIDORT_SERIOUS
        MESSAGES(NM) = 'Cannot find string: '//adjustl(trim(PAR_STR))
        ACTIONS(NM)  = 'Check Spelling of String in input file'
      ENDIF

!  finish

      RETURN
      END SUBROUTINE FINDPAR_ERROR

!  Convert a string to lower-case.  (ASCII character set assumed.)

      SUBROUTINE LCSTRING ( STRING )

      implicit none

!  Modified argument.

      CHARACTER * (*)     STRING

!  Local variables.

      INTEGER         CCC
      INTEGER         III
      INTEGER         OFFSET

!***********************************************************************

      OFFSET = ICHAR ( 'a' ) - ICHAR ( 'A' )
      DO 10 III = 1, LEN ( STRING )
        CCC = ICHAR ( STRING ( III : III ) )
        IF ( CCC .GE. ICHAR ( 'A' ) .AND. &
                CCC .LE. ICHAR ( 'Z' ) ) THEN
          STRING ( III : III ) = CHAR ( CCC + OFFSET )
        END IF
10    CONTINUE

      END SUBROUTINE LCSTRING

!  End module

      END MODULE vbrdf_findpar_m


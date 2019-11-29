!======================================================================+
!   SPECTRAL_ANALYSIS
!   Module to derive power/amplitude spectrum from an input time series
!======================================================================+
!
!   Contact information:
!   Garry Hayman
!   Centre for Ecology and Hydrology
!   Wallingford
!   Oxfordshire, OX10 8BB
!   United Kingdom
!   
!   email: garr@ceh.ac.uk
!   http://www.ceh.ac.uk
!
!   Development history
!
!   Date      Version     Author      Description
!   11/13      1.0        GDH         Initial version developed as
!                                     FORTRAN module from Graham Weedon's
!                                     time series FORTRAN code for use
!                                     with PYTHON.                               
!
!   11/19      1.1        GDH         Release version
!
!======================================================================+
!
!   This contains subroutines and functions
!
!   Subroutines
!   -----------
!
!   TIME_SERIES_RUNFILE
!   TIME_SERIES_GRIDDED
!   TIME_SERIES_DATA_EXTRACT
!   TIME_SERIES
!   SOLVE
!   HPSORT
!
!   Functions
!   ---------
!   ROFUNC
!   RAN1
!   GASDEV
!   BRENT
!   XLS
!   BRENT2
!   XLS2
!   BRENT3
!   XLS3
!
!======================================================================+
MODULE MODULE_TIME_SERIES_DATA
!======================================================================+

!   Declaration of Parameters/Variables common to
!   all routines and modules

      IMPLICIT                                  NONE

!     Define limits for arrays:

      INTEGER,PARAMETER                         :: &
            NMAX     = 300001, &
            NMAXB    = 150000

      INTEGER                                   :: &
            NDATAT

      REAL                                      :: &
            AA,ABDEVT

      REAL,DIMENSION(NMAX)                      :: &
            XT,YT,ARR

      REAL,ALLOCATABLE                           :: &
            data_local(:),time_local(:)

      CHARACTER (len=50)                        :: &
            time_series_routine

      SAVE

END MODULE MODULE_TIME_SERIES_DATA
!======================================================================+
MODULE MODULE_TIME_SERIES_PROCEDURES
!======================================================================+

CONTAINS

!======================================================================+
   SUBROUTINE TIME_SERIES_RUNFILE(tsfile,fmt,pnyq,pray,cfix,cpower, &
              cfit,nmin,nfit0,kmk,xlowp,xhighp,xnray,xnyq,crit_freq, &
              iconf_level,miss_data)
!======================================================================+

      USE                                       MODULE_TIME_SERIES_DATA

      IMPLICIT                                  NONE

      INTEGER                                   :: &
            ios

      INTEGER,INTENT(OUT)                       :: &
            nmin,nfit0,kmk,iconf_level

      REAL,INTENT(OUT)                          :: &
            crit_freq,miss_data,xlowp,xhighp,xnray,xnyq

      CHARACTER (len=  1),INTENT(OUT)           :: &
            pnyq,pray,cfix,cpower,cfit
      CHARACTER (len=  1)                       :: &
            proc1,proc2
      CHARACTER (len= 41)                       :: &
            header
      CHARACTER (len= 60),INTENT(OUT)           :: &
            fmt
      CHARACTER (len= 80),INTENT(OUT)           :: &
            tsfile
      CHARACTER (len=100)                       :: &
            runfile

!-----------------------------------------------------------------------

      time_series_routine    = 'time_series_runfile'
      WRITE(*,*) 'In routine: ',time_series_routine

!   Garry Hayman
!   November 2013
!   Removed all READ statement; variables to be passed to code by calling
!   program

      runfile='/data/grp/eow/garr/Projects/Methane/CODE/FORTRAN/'// &
            'z_TEST/time_series_runfile.dat'
      OPEN(UNIT=2,FILE=runfile,IOSTAT=ios,STATUS='OLD')

!   Garry Hayman
!   October 2010
!   Input minimum number of data points

      READ(2,'(41X,I5)'  ,IOSTAT=ios) nmin

!   Input values of lowest and highest acceptable value to remove outliers

      READ(2,'(41X,F8.2)',IOSTAT=ios) xlowp
      READ(2,'(41X,F8.2)',IOSTAT=ios) xhighp

!   Input Proceed with program?  Proceed = P, Stop = S

      READ(2,'(41A,A)'   ,IOSTAT=ios) header,proc1

!   Read in time series data from "TSFILE" as arrays T (= observation
!   time or position) and V (variable value or rock type code).

      READ(2,'(41A,(A))') header,TSFILE
      WRITE(*,'(A)') TSFILE

!   Enter the Fortran format for the location of the Observation time and variable value

      READ(2,'(41A,(A))',IOSTAT=ios) header,fmt

!   Input missing data value
 
      READ(2,'(41X,F8.2)',IOSTAT=ios) miss_data

!   Do you want to stop the program here?

      READ(2,'(41A,A)',IOSTAT=ios) header,proc2

!   Observation time/stratigraphic units

      READ(2,'(41X,I1)',IOSTAT=ios) kmk

!   Change the Rayleigh Frequency? (Y or N)

      READ(2,'(41A,A)',IOSTAT=ios) header,pray
      READ(2,*,IOSTAT=ios) xnray

!   Change the Nyquist Frequency? (Y or N)

      READ(2,'(41A,A)',IOSTAT=ios) header,pnyq
      xnyq  = 0.0
      cfix  = 'N'

!   Absolute or relative power? Abs. = A, Rel. = R

      READ(2,'(41A,A)',IOSTAT=ios) header,cpower

!   Choose the form of the spectral background

      READ(2,'(41X,I1)',IOSTAT=ios) nfit0

!   Quadratic (Q) or power law (P) fit

      READ(2,'(41A,A)') header,cfit
!
!   Input critical frequency
 
      READ(2,'(41X,F8.2)',IOSTAT=ios) crit_freq

!   Select confidence level (1=90%,2=95%,3=99%)

      READ(2,'(41X,I1)',IOSTAT=ios) iconf_level

      WRITE(*,*) crit_freq,iconf_level,miss_data, &
                 xnray,xnyq,nfit0,kmk,pray,pnyq,cpower,cfit
      WRITE(*,*) 'End of input from runfile.... '
!
!   Return to calling routine

      RETURN

   END SUBROUTINE TIME_SERIES_RUNFILE
!
!======================================================================+
   SUBROUTINE TIME_SERIES_GRIDDED(data_name,DIM_DATA,NSIZE,NTIMES, &
           nfreqs,nrays,data_all,time_all,dates,power, &
           NMIN,XLOWP,XHIGHP,xnrays,XNYQ,FMT,NFIT0,KMK,PRAY,PNYQ, &
           CFIX,CPOWER,CFIT,crit_freqs,iconf_level,miss_data, &
           FILENAME,nspect,nspect_valid,use_flag,error_flags,debug)
!======================================================================+

      USE                                       MODULE_TIME_SERIES_DATA

      IMPLICIT                                  NONE

      INTEGER,PARAMETER                         :: &
            max_dims = 10,                         &
            max_size = 200

      INTEGER,INTENT(IN)                        :: &
            nsize,nmin,nfit0,kmk,ntimes,nfreqs,iconf_level,nrays,use_flag

      INTEGER,INTENT(IN)                        :: &
            DIM_DATA(max_dims+1),DATES(NTIMES,3)

      INTEGER,INTENT(OUT)                       :: &
            nspect,nspect_valid,error_flags(nsize,nfreqs)

      INTEGER                                   :: &
            IPOWER,JPOWER,IPEAK,NDATA,ISIZE,ifreq,itime,MN,error_flag, &
            IDIM,IDIM1,IDIM2,IDIM3,IDIM4,NDIMS,NDIM1,NDIM2,NDIM3,NDIM4, &
            IFLAG_MON,IFLAG_SEA

      REAL,INTENT(IN)                           :: &
            XLOWP,XHIGHP,MISS_DATA

      REAL,INTENT(INOUT)                        :: &
            XNYQ

      REAL,INTENT(IN)                           :: &
            DATA_ALL(NSIZE,NTIMES),TIME_ALL(NTIMES), &
            xnrays(nrays),crit_freqs(nfreqs)

      REAL,INTENT(OUT)                          :: &
            POWER(NSIZE,nfreqs)

      REAL                                      :: &
            xnray,bandwidth,variance,freq,freq_min,freq_max,peak_power, &
            OUTPUT(NMAXB,7),DATA_MOD(NMAX),TIME_MOD(NMAX), &
            output_red(nfreqs,7)
!
      CHARACTER (len= 1),INTENT(INOUT)          :: &
            PNYQ,PRAY,CFIX,CPOWER,CFIT,DEBUG
      CHARACTER (len=60),INTENT(IN)             :: &
            FMT
      CHARACTER (len=20),INTENT(IN)             :: &
            data_name
      CHARACTER (len=150),INTENT(IN)            :: &
            FILENAME
      CHARACTER (len=1)                         :: &
            CFLAG
      CHARACTER (len=30)                        :: &
            GRID_REF(NSIZE)

!-----------------------------------------------------------------------

      time_series_routine    = 'time_series_gridded'
      WRITE(*,*) 'In routine: ',time_series_routine

      WRITE(*,*) DEBUG

      IF (debug == 'Y') THEN
          WRITE(*,'(A11,I2,I6)')     'USE_FLAG = ',USE_FLAG,NTIMES

          DO itime=1,NTIMES
             WRITE(*,'(3(I4))') DATES(itime,1),DATES(itime,2), &
                                DATES(itime,3)
          END DO
      END IF

!   Allocate arrays: time_local,data_local

      ALLOCATE(time_local(ntimes))
      ALLOCATE(data_local(ntimes))

      WRITE(*,*) FILENAME
      WRITE(CFLAG,'(I1)') USE_FLAG
      OPEN(UNIT=7,FILE=TRIM(FILENAME)//'_F'//CFLAG//'.spc')
      OPEN(UNIT=8,FILE=TRIM(FILENAME)//'_F'//CFLAG//'_red.spc')
      IF (NSIZE <= max_size) OPEN(UNIT=9,FILE=TRIM(FILENAME)//'_F'//CFLAG//'.dat')

!   DATA, NUMBER and TIME contain the data, number of data points
!   and the time points

      NDIMS         = DIM_DATA(1)
      IF (NDIMS .GT. max_dims) THEN
         WRITE(*,*) '*** Exceeded maximum number of dimensions ***',NDIMS,max_dims
         STOP
      ENDIF

!   Assign grid square indices for ease

      IF (NDIMS .EQ. 2) THEN

         NDIM1     = DIM_DATA(3)

         DO IDIM1  = 1,NDIM1
            IDIM   = IDIM1
!           WRITE(*,'(I6)') IDIM1
            WRITE(GRID_REF(IDIM),'(I6)') IDIM1
         ENDDO           ! Loop over IDIM1

      ELSE IF (NDIMS .EQ. 3) THEN

         NDIM1     = DIM_DATA(3)
         NDIM2     = DIM_DATA(4)

         DO IDIM1  = 1,NDIM1
            DO IDIM2  = 1,NDIM2
               IDIM   = (IDIM1-1)*NDIM2+IDIM2
!              WRITE(*,'(I6,I6,I6)') IDIM,IDIM1,IDIM2
               WRITE(GRID_REF(IDIM),'(I6,I6)') IDIM1,IDIM2
            ENDDO        ! Loop over IDIM2
         ENDDO           ! Loop over IDIM1

      ELSE IF (NDIMS .EQ. 4) THEN

         NDIM1     = DIM_DATA(3)
         NDIM2     = DIM_DATA(4)
         NDIM3     = DIM_DATA(5)

         DO IDIM1  = 1,NDIM1
            DO IDIM2  = 1,NDIM2
               DO IDIM3  = 1,NDIM3
                  IDIM   = ((IDIM1-1)*NDIM2+(IDIM2-1))*NDIM3+IDIM3
!                 WRITE(*,'(I6,I6,I6,I6)') IDIM,IDIM1,IDIM2,IDIM3
                  WRITE(GRID_REF(IDIM),'(I6,I6,I6)') IDIM1,IDIM2,IDIM3
               ENDDO     ! Loop over IDIM3
            ENDDO        ! Loop over IDIM2
         ENDDO           ! Loop over IDIM1

      ENDIF

!   Loop over length of vector - defined by ndata

      WRITE(*,*) SHAPE(GRID_REF),NSIZE,NTIMES

      nspect         = 0
      nspect_valid   = 0
      time_local     = time_all(:)
      xnray          = xnrays(1)

      DO 100 ISIZE=1,NSIZE

         WRITE(7,1000) TRIM(GRID_REF(ISIZE))
         IF (debug == 'Y') WRITE(*,1000) TRIM(GRID_REF(ISIZE))

!   Extract data

         xnyq           = 0.0
         variance       = 0.0
         bandwidth      = 0.0
         output_red     = 0.0
         data_local     = data_all(ISIZE,:)

         CALL TIME_SERIES_DATA_EXTRACT(NTIMES,ndata, &
              DATA_MOD,TIME_MOD,DATES,IFLAG_MON,IFLAG_SEA, &
              MISS_DATA,DEBUG)
!        WRITE(*,*) NMIN,ndata

!   Select option depending on ndata
!   (1) If no data, set to missing data

         IF (ndata .EQ. 0) THEN
!
            ERROR_FLAG           = -999
            WRITE(7,1010) ERROR_FLAG
            POWER(isize,:)       = MISS_DATA
            ERROR_FLAGS(isize,:) = ERROR_FLAG 

         ELSE IF (ndata .LT. NMIN) THEN

            ERROR_FLAG           = -2
            WRITE(7,1010) ERROR_FLAG
            POWER(isize,:)       = -2.0
            ERROR_FLAGS(isize,:) = ERROR_FLAG

         ELSE IF (ndata .GE. NMIN) THEN

!   Derive time-series if ndata > NMIN
!   Assign Rayleigh frequency

            IF (nrays > 1) xnray = xnrays(isize)

            CALL TIME_SERIES(ndata,NMIN,DATA_MOD,TIME_MOD,XLOWP,XHIGHP,XNRAY,XNYQ, &
                 FMT,NFIT0,KMK,PRAY,PNYQ,CFIX,CPOWER,CFIT, &
                 MN,bandwidth,variance,OUTPUT,ERROR_FLAG)

            WRITE(7,1010) ERROR_FLAG

            nspect               = nspect+1
            ERROR_FLAGS(isize,:) = ERROR_FLAG

            IF (error_flag == 0) THEN

               nspect_valid  = nspect_valid+1

               WRITE(7,1015) ndata,bandwidth,variance

               DO IPOWER=1,MN
                     WRITE(7,1020) (OUTPUT(IPOWER,JPOWER),JPOWER=1,7)
               ENDDO   ! Loop over MN

               DO ifreq=1,nfreqs

                  freq_min      = crit_freqs(ifreq)-bandwidth/2.0
                  freq_max      = crit_freqs(ifreq)+bandwidth/2.0
                  peak_power    = 0.0

                  DO IPOWER=1,MN

                     freq          = OUTPUT(IPOWER,1)
                     IF (freq >= freq_min .and. freq <= freq_max .and. &
                         OUTPUT(IPOWER,3) >= peak_power) THEN 
                           peak_power    = OUTPUT(IPOWER,3)
                           ipeak         = ipower
                           output_red(ifreq,:) = OUTPUT(ipeak,:)                 
                     END IF

                  ENDDO   ! Loop over MN

!         Assign power

                  IF ((peak_power >= OUTPUT(ipeak,4+iconf_level)) .AND. &
                      ((use_flag .EQ. 0)                          .OR.  &
                       (use_flag .EQ. 1 .AND. IFLAG_MON .EQ. 1)   .OR.  &
                       (use_flag .EQ. 2 .AND. IFLAG_SEA .EQ. 1))) THEN 
                     IF (cpower == 'A' .or. cpower == 'a') THEN
                         POWER(isize,ifreq)  = peak_power
                     ELSE
                         POWER(isize,ifreq)  = SQRT(0.5*peak_power*variance/0.147)
                     END IF
                  ELSE IF (use_flag .EQ. 1 .AND. IFLAG_MON .EQ. 0) THEN
                     POWER(isize,ifreq)       = -2
                     ERROR_FLAGS(isize,ifreq) = -3
                  ELSE IF (use_flag .EQ. 2 .AND. IFLAG_SEA .EQ. 0) THEN
                     POWER(isize,ifreq)       = -2
                     ERROR_FLAGS(isize,ifreq) = -4
                  ELSE
                     POWER(isize,ifreq)       = -1
                     ERROR_FLAGS(isize,ifreq) = -1
                  END IF

               ENDDO   ! Loop over NFREQS

            ELSE

               DO ifreq=1,nfreqs
                  POWER(isize,ifreq)       = -1
                  ERROR_FLAGS(isize,ifreq) = -1
               ENDDO   ! Loop over NFREQS
                  
            ENDIF

         ENDIF

         DO ifreq=1,nfreqs
            IF (debug == 'Y') WRITE(*,1010) ERROR_FLAGS(isize,ifreq)
            WRITE(8,1025) GRID_REF(ISIZE),ndata,error_flags(isize,ifreq), &
                          variance,crit_freqs(ifreq),bandwidth, &
                          (output_red(ifreq,JPOWER),JPOWER=1,7)
         ENDDO   ! Loop over NFREQS

!   Output input data for debug purposes

         IF (NSIZE <= max_size) THEN
            WRITE(9,1000) TRIM(GRID_REF(ISIZE))
            WRITE(9,1030) ndata,xnray,xnyq
            DO itime=1,ndata
               WRITE(9,FMT) time_mod(itime),data_mod(itime)
            END DO
         END IF

  100 ENDDO       ! Loop over NSIZE

      CLOSE(7)

!   Deallocate arrays: time_local,data_local

      DEALLOCATE(time_local)
      DEALLOCATE(data_local)

!   FORMAT statements

 1000 FORMAT('Grid reference        = ',A30)
 1010 FORMAT('Error Flag            = ',I6)
 1015 FORMAT('Number of data points = ',I6, &
            /'Bandwidth             = ',F14.8, &
            /'Variance              = ',F14.8)
 1020 FORMAT(F11.5,1X,F11.4,5(1X,F11.6))
 1025 FORMAT(A20,2I6,F12.4,9(1X,F11.5))
 1030 FORMAT('Number of data points = ',I6, &
            /'Rayleigh frequency    = ',F14.8, &
            /'Nyquist  frequency    = ',F14.8)

!   Return to calling routine
!
      RETURN

   END SUBROUTINE TIME_SERIES_GRIDDED
!======================================================================+
   SUBROUTINE TIME_SERIES_DATA_EXTRACT(NTIMES_IN,NTIMES_OUT, &
              DATA_MOD,TIME_MOD,DATES,IFLAG_MON,IFLAG_SEA, &
              MISS_DATA,DEBUG)
!======================================================================+

      USE                                       MODULE_TIME_SERIES_DATA

      IMPLICIT                                  NONE

      INTEGER,INTENT(IN)                        :: &
            ntimes_in,dates(ntimes_in,3)

      INTEGER,INTENT(OUT)                       :: &
            ntimes_out,iflag_mon,iflag_sea

      INTEGER                                   :: &
            itime,iyear,imonth,iday,iseason,months(12),seasons(4)

      INTEGER,DIMENSION(12)                     :: &
            mon2season = (/ 1,1,2,2,2,3,3,3,4,4,4,1 /)

      REAL,INTENT(IN)                           :: &
            MISS_DATA

      REAL,INTENT(OUT)                          :: &
            DATA_MOD(NMAX),TIME_MOD(NMAX)

      CHARACTER (len= 1),INTENT(IN)             :: &
            DEBUG

      CHARACTER (len=10)                        :: &
            TEMP_DATE

!-----------------------------------------------------------------------

      time_series_routine    = 'time_series_data_extract'
      IF (DEBUG == 'Y') WRITE(*,*) 'In routine: ',time_series_routine

!   DATA and TIME contain the data and the time points at a particular location

      IF (DEBUG == 'Y') WRITE(*,*) ntimes_in

      ntimes_out = 0
      DATA_MOD   = MISS_DATA
      TIME_MOD   = MISS_DATA
      MONTHS     = 0
      SEASONS    = 0
      iflag_mon  = 1
      iflag_sea  = 1

!   Loop over data, removing any missing data

      DO 100 itime=ntimes_in,1,-1

         iyear      = DATES(itime,1)
         imonth     = DATES(itime,2)
         iday       = DATES(itime,3)
         iseason    = mon2season(imonth)
         IF (debug == 'Y') WRITE(*,*) itime,iyear,imonth,iday,iseason

         IF (data_local(itime) /= MISS_DATA) THEN
             ntimes_out = ntimes_out+1
             MONTHS(imonth)       = MONTHS(imonth)+1
             SEASONS(iseason)     = SEASONS(iseason)+1
             DATA_MOD(ntimes_out) = data_local(itime)
             TIME_MOD(ntimes_out) = time_local(itime)
             IF (DEBUG == 'Y') WRITE(*,*) itime,ntimes_out, &
                        TIME_MOD(ntimes_out),DATA_MOD(ntimes_out), &
                        time_local(itime),data_local(itime)
         END IF

  100 END DO

      IF (DEBUG == 'Y') THEN
         WRITE(*,*) 'Months'
         WRITE(*,*) MONTHS
         WRITE(*,*) 'Seasons'
         WRITE(*,*) SEASONS
      END IF

!   See whether data in all months/seasons

      DO 110 imonth=1,12
         IF (MONTHS(imonth) == 0)   iflag_mon = 0
  110 END DO

      DO 120 iseason=1,4
         IF (SEASONS(iseason) == 0) iflag_sea = 0
  120 END DO

!   Return to calling routine
!
      RETURN

   END SUBROUTINE TIME_SERIES_DATA_EXTRACT
!======================================================================+
   SUBROUTINE TIME_SERIES(NDATA,NMIN,XINPUT,TINPUT,XLOWP,XHIGHP, &
           XNRAY,XNYQ,FMT,NFIT0,KMK,PRAY,PNYQ,CFIX,CPOWER,CFIT, &
           MN,BW,VARX1,OUTPUT,ERROR_FLAG)
!======================================================================+

      USE                                       MODULE_TIME_SERIES_DATA
      IMPLICIT                                  NONE

!   References to FUNCTION names not needed as accessed through module

      REAL,EXTERNAL                             :: &
            ROFUNC,GASDEV

      DOUBLE PRECISION,EXTERNAL                 :: &
            BRENT,BRENT2,BRENT3,XLS,XLS2,XLS3

!     Declare local integer variables and arrays:

      INTEGER,INTENT(IN)                        :: &
            NMIN,NFIT0

      INTEGER,INTENT(OUT)                       :: &
            MN,ERROR_FLAG

      REAL,INTENT(IN)                           :: &
            XLOWP,XHIGHP,XINPUT(NMAX),TINPUT(NMAX)

      REAL,INTENT(INOUT)                        :: &
            XNRAY,XNYQ

      REAL,INTENT(OUT)                          :: &
            BW,VARX1

      REAL,DIMENSION(NMAXB,7),INTENT(OUT)       :: &
            OUTPUT

      INTEGER                                   :: &
            I,J,K,L,N,IOS,IOS2,NDATA,NSTAT,NSTOP,NSTART,NINT, &
            NREV,NIRR,NREVINT,NVARY,NORIG,NOUT,NOUT1,NOUTN,NEW,NDOF, &
            NARSIM,NMU,NCC,NSS,NFIT,NSPECB,NWMAX,NWMIN,NWIND,NWIND2, &
            NWIND2B,IR,MID,NROWS,KMK,MUDEST,MEDPOS,MEDPOS2,MEDPOS2I, &
            NMED,IORD1,NP1,NSIM,NARSIMC,IDUM,NBIASC,NST,NTOT,NTOT2,IA, &
            ID,NMIS,CTLIM,CTLIM1,CTLIM2,CTLIM3,PSTOP,NFLIM,NLIM1,NLIM2, &
            NPOS,NNEG,PSTOP1,OUTI

      REAL                                      :: &
            Y(NMAXB),XMEDIAN(NMAXB),XLFREQ(NMAXB),XLSPBIAS(NMAXB), &
            XLP(NMAXB),ARP2(NMAXB),XLARSP(NMAXB),BIASF(NMAXB),XLOGA(NMAXB), &
            XOUT(NMAXB),PSMTH(NMAXB),ARSPAVE(NMAXB),AC(NMAXB),BC(NMAXB), &
            CIH(10),CIW(10),BWH(10),BWV(10),B(20,20),C(20),XP(20), &
            ROUT(NMAX),P(NMAXB),FREQ(NMAXB),X(NMAX),X2(NMAX),SROH, &
            SYINTA,SYINTB,AVE,CCOS,CC,CWTAU,PNOW,SSIN,SS,SUMC,SSDIFF,EP, &
            SYINTD,SUMCY,SUMS,SUMSH,SUMSY,SWTAU,VARX,WTAU,STAVE,STDIF, &
            XMAX,XMIN,YY,ANSTART,XNNYQ,XNEW,AMN,AN,DEL,A,BB,BIT, &
            CHISQ,ZINTMIN,ANWIND,ANWIND2,ANWIND2B,XXP,YYP,XNSIM, &
            XLSIN,XLSMIN,ZINT1,STVAL,ZINT,ZINTMAX,ZINTVA,ZINTVB,XVMIN, &
            XVMAX,STRATD,XMAXT,XMINT,ANORIG,XFINT,XFINT1,XF4,XF2, &
            XNYN,XNDOF,SX,SY,SXY,SXX,SIGB,B1,B2,F,F1,F2,XINTERC,SLOPE, &
            ABDEV,XFIT,AL,TAU,TAUMIN,SCALT,STDIFF,STDEV,XLSOF, &
            SDIFF,AI,TOTP,TOTP1,ZTEST,SUM,SOLMEAN,SO,SOF,TEMP,SELCT, &
            TMP,DIFF,ANTOT,ANTOT2,TSLEN,XPERIOD,XMINPER,XMAXPER,PA,PB, &
            PC,TEMPER,TLOGP,TLFRE,TEMPBK,TEMP90,TEMP95,TEMP99,ANMED, &
            XMAXLP,XMINLP,DIFFMAX,SOF2,SOFNEW,PROP,XLIM

      DOUBLE PRECISION                          :: &
            ARG,WTEMP,WI(NMAX),WPI(NMAX),WPR(NMAX),AMIN, &
            WR(NMAX),DUMMA,DUMMB,XSCAL(NMAX),TSCAL(NMAX), &
            DUM1,DUM2,DUM3,DUM4,DUM5, &
            AR11,AR12,AR13,XLOUT(NMAXB),XFR(NMAXB),SOFA,XPI,XNYQ1, &
            SOFT,GUESS,SOFIT,SOLMIN,SOLMAX,SOLINT,SOLMID, &
            SOLMI,SOLMD,SOLMA,XNYQNEW,XLARP2(NMAXB),YINTA(NMAXB), &
            YINTB(NMAXB),YINTD(NMAXB),RANGE,DUMMIN,XPROP,WEIGHT,XLNYQ, &
            XLFR(NMAXB),ROH,ROHT,ROHB,ROH2,ROHF,ROHI,ROHN,ROHF2,R90, &
            R95,R99,XLSUM,XARP2,XLP2,XLMIN,XLREV,BRENTSUM,XINT,STRAT1, &
            STRAT(NMAX),XINTFIX,FACTOR

      DOUBLE PRECISION,PARAMETER                :: &
            AR1    = 0.367879441E+00, TOL  = 3.0E-08, TOL2 = 1.0E-06, &
            TWOPID = 6.2831853071795865

      CHARACTER (len= 1),INTENT(INOUT)          :: &
            PNYQ,PRAY,CFIX,CPOWER,CFIT
      CHARACTER (len=41)                        :: &
            HEADER
      CHARACTER (len=60),INTENT(IN)             :: &
            FMT
      CHARACTER (len=80)                        :: &
            TSFILE

!-----------------------------------------------------------------------

!     RUNFILE='/data/grp/eow/garr/Projects/Methane/CODE/FORTRAN/'// &
!              'z_POWER/spectral_runfile.dat'

!     OPEN(UNIT=2,FILE=RUNFILE,IOSTAT=IOS,STATUS='OLD')

      time_series_routine    = 'time_series'
      WRITE(*,*) 'In routine: ',time_series_routine
      WRITE(*,*) fmt
      WRITE(*,*) xlowp,xhighp,nfit0,cfit,kmk,xnray,xnyq,pray,pnyq

!   Garry Hayman
!   November 2013
!   Commented out all READ statement; variables to be passed to code
!   by calling program
!
!   January 2014
!   PYTHON code written to comment out INPUT/OUTPUT-related statements
!   (WRITE, FORMAT, OPEN and CLOSE)

!   Replaced uncommented STOP with RETURN and error code in ERROR_FLAG

      ERROR_FLAG = 0
      TSFILE = 'From Python-FORTRAN'
      WRITE(*,2)
      WRITE(*,*)
      WRITE(*,*)'***            SPECTRA.F                  ***'
      WRITE(*,*)'***    Lomb-Scargle power spectrum for    ***'
      WRITE(*,*)'*** regularly or irregularly spaced data. ***'
      WRITE(*,*)'***      Graham P. Weedon, 2009           ***'                                               
      WRITE(*,*)
      WRITE(*,2)
      WRITE(*,*)
      WRITE(*,*)'Data file requirements:'
      WRITE(*,*)'a: Each line = Observation time, Variable value.'
      WRITE(*,*)'b: The time between observations can be fixed or'
      WRITE(*,*)'   variable.'
      WRITE(*,*)'c: Data must be ascii aligned columns (fixed format),'
      WRITE(*,*)'   but with NO TEXT headers.'

!   Garry Hayman
!   October 2010
!   Input minimum number of data points

!     READ(2,'(41X,I5)',IOSTAT=IOS)N
      WRITE(*,5) NMIN

    2 FORMAT('  *************************************************** ')
    5 FORMAT(' d: The time series must have > ',I4, &
            ' and < 300,000 points')
!
      WRITE(*,*)'e: The latest observation should be the first line.'
      WRITE(*,*)
      WRITE(*,*)'Program input needed:'
      WRITE(*,*)'a: The time series filename.'
      WRITE(*,*)'b: The FORTRAN format of the data.'
      WRITE(*,*)'c: The sample interval and units'
      WRITE(*,*)'   (e.g. day, year, kyr, m).'
      WRITE(*,*)
      WRITE(*,*)'Program output:'
      WRITE(*,*)'To file: spectral data are written to "spectrum"'
      WRITE(*,*)'   and to "logspec".'

!   Garry Hayman
!   October 2010
!   Input values of lowest and highest acceptable value to remove outliers

      WRITE(*,*)
      WRITE(*,*) 'Minimum and maximum acceptable values '
!     READ(2,'(41X,F8.2)',IOSTAT=IOS)XLOWP
      WRITE(*,'(F8.2)')XLOWP
!     READ(2,'(41X,F8.2)',IOSTAT=IOS)XHIGHP
      WRITE(*,'(F8.2)')XHIGHP

!   Garry Hayman
!   January 2014
!   Commented out DO loop 10 as assumed will proceed

!     CTLIM=0
!     DO 10
!        CTLIM=CTLIM+1
!        IF(CTLIM > 5)EXIT                   ! EXIT DO 10
!        WRITE(*,*)
!        WRITE(*,*)' Proceed with program?  Proceed = P, Stop = S:'
!        READ(2,'(41A,A)',IOSTAT=IOS)HEADER,PROC
!        WRITE(*,'(A)') PROC
!        IF(PROC == 'P'.OR.PROC == 'p')EXIT   ! EXIT DO 10
!        IF(PROC == 'S'.OR.PROC == 's')STOP
!        IF(IOS /= 0)THEN
!           WRITE(*,*)
!           WRITE(*,*)'  *** Enter P or S ***'
!           CYCLE                              ! CYCLE DO 10
!        ENDIF
! 10  ENDDO

!   Time series data passed to SUBROUTINE as arrays TINPUT (= observation
!   time or position) and XINPUT (variable value or rock type code).

!   Garry Hayman
!   January 2014
!   Commented out DO loop 25 as no longer reading from datafile

!     DO 25
! 15     WRITE(*,*)
!        WRITE(*,*)' Enter time series filename:'
!        WRITE(*,*)' (Nb file must be in current directory).'
!        WRITE(*,*)' (Nb filename must be no more than 80 characters).'
!        WRITE(*,*)'  Example file = bmbase'
!        WRITE(*,*)
!        READ(2,'(41A,(A))',ERR=15)HEADER,TSFILE
!        WRITE(*,'(A)') TSFILE
!        OPEN(UNIT=1,FILE=TSFILE,IOSTAT=IOS,STATUS='OLD')
!        IF(IOS /= 0)THEN
!           WRITE(*,*)
!           WRITE(*,20)TSFILE
! 20        FORMAT('  *** Error opening ',A,'***')
!           WRITE(*,*)' *** (Filename must not have more than ***'
!           WRITE(*,*)' *** 3 characters after point.         ***'
!           WRITE(*,*)
!           WRITE(*,*)' *** Try entering filename again or    ***'
!           WRITE(*,*)' *** use <Ctrl> & <C> to end program.  ***'
!           CYCLE                          ! CYCLE DO 25
!        ELSE
!           EXIT                           ! EXIT DO 25
!        ENDIF
! 25  ENDDO

!   Garry Hayman
!   January 2014
!   Commented out DO loop 30 as no longer reading from datafile

!     CTLIM=0
!     DO 30
!        CTLIM=CTLIM+1
!        IF(CTLIM > 5)THEN
!           WRITE(*,*)
!           WRITE(*,*)' *** Format not set ***'
!           STOP                           ! STOP
!        ENDIF
!        WRITE(*,*)
!        WRITE(*,*)' Enter the Fortran format for the location of the'
!        WRITE(*,*)' Observation time and variable value.'
!        WRITE(*,*)
!        WRITE(*,*)' eg.   bmbase uniform spacing:  (18X,F5.2,3X,F5.2)'
!        WRITE(*,*)' eg. bmbase irregular spacing:  (9X,F7.4,10X,F5.2)'
!        WRITE(*,*)
!        READ(2,'(41A,(A))',IOSTAT=IOS)HEADER,FMT
!        WRITE(*,'(A)')FMT
!        IF(IOS /= 0)THEN
!           WRITE(*,*)
!           WRITE(*,*)' *** Format entry problem (type brackets). ***'
!           CYCLE                           ! CYCLE DO 30
!        ELSE
!           EXIT                            ! EXIT DO 30
!        ENDIF
! 30  ENDDO

      NROWS=0
      DO 50 I=1,NDATA

!        NROWS=NROWS+1
!        READ(1,FMT,END=55,IOSTAT=IOS)STRAT(NROWS),X(NROWS)
!        IF(IOS /= 0)THEN
!           WRITE(*,*)
!           WRITE(*,*)'      *** Problem with format statement ***'
!           WRITE(*,*)'      *** Memo: Remove text headings    ***'
!           STOP
!        ENDIF

!   Garry Hayman
!   October 2010
!   Remove line if outlier

         IF(XINPUT(I) < XLOWP .OR. XINPUT(I) > XHIGHP) THEN
            WRITE(*,51) I
   51       FORMAT('Outlier removed at row',I6)
            WRITE(*,FMT) TINPUT(I),XINPUT(I)
         ELSE
            NROWS=NROWS+1
            X(NROWS)     = XINPUT(I)
            STRAT(NROWS) = TINPUT(I)
         ENDIF

   50 ENDDO

      N       = NROWS
      NSTART  = N
      ANSTART = N
      NFIT    = NFIT0

!   Check for irregularly spaced data. Compare observation times
!   for first two points with all others. If variation exceeds +/-5.0%
!   then data are considered to be irregularly spaced.

   60 FORMAT('  Obs. time      Line')
   61 FORMAT(2X,F10.5,2X,I6)          
!
      NINT=0
      NREV=0
      NIRR=0
      NREVINT=0
      ZINT1=STRAT(1)-STRAT(2)
!   ZINT1 GT 0.0
      IF(ZINT1 > 0.0)THEN
        STVAL=1.0
        ZINTMIN=ZINT1
        ZINTMAX=ZINT1
        ZINTVA=ZINT1*1.05
        ZINTVB=ZINT1*0.95
        DO 64 I=2,N-1
          ZINT=STRAT(I)-STRAT(I+1)      
          IF(ZINT <= 0.0)NREV=NREV+1
          IF(NREV == 1)THEN
            NREVINT=NREVINT+1
            IF(NREVINT == 1)THEN
             OPEN(UNIT=11,FILE='tsrev')
             WRITE(11,60)
            ENDIF
          ENDIF
          IF(ZINT <= 0.0)WRITE(11,61)STRAT(I+1),I
          IF(ZINT > ZINTVA)THEN 
            NINT=NINT+1
            IF(ZINTMAX < ZINT)ZINTMAX=ZINT
          ENDIF
          IF(ZINT < ZINTVB)THEN
            NINT=NINT+1
            IF(ZINTMIN > ZINT)ZINTMIN=ZINT
          ENDIF
   64    ENDDO
!   ZINT1 LE 0.0
      ELSE
        STVAL=-1.0
        ZINT1=STRAT(2)-STRAT(1)
        ZINTMIN=ZINT1
        ZINTMAX=ZINT1
        ZINTVA=ZINT1*1.05
        ZINTVB=ZINT1*0.95
        DO 65 I=2,N-1
          ZINT=STRAT(I+1)-STRAT(I)      
          IF(ZINT <= 0.0)NREV=NREV+1
          IF(NREV == 1)THEN
            NREVINT=NREVINT+1
            IF(NREVINT == 1)THEN
             OPEN(UNIT=11,FILE='tsrev')
             WRITE(11,60)
            ENDIF
          ENDIF
          IF(ZINT <= 0.0)WRITE(11,61)STRAT(I+1),I
          IF(ZINT > ZINTVA)THEN 
            NINT=NINT+1
            IF(ZINTMAX < ZINT)ZINTMAX=ZINT
          ENDIF
          IF(ZINT < ZINTVB)THEN
            NINT=NINT+1
            IF(ZINTMIN > ZINT)ZINTMIN=ZINT
          ENDIF
   65    ENDDO
      ENDIF

      IF(NINT >= 1)THEN
        NIRR=1
        IF(NREV > 0)THEN
          CLOSE(UNIT=11)
          WRITE(*,*)
          WRITE(*,66)NREV
   66      FORMAT('  PROBLEM: There are ',I7)
          WRITE(*,*)' equal or reversed observation times !'
          WRITE(*,*)
          WRITE(*,*)' Consult output file "tsrev" for list'
          WRITE(*,*)' of equal or reversed observation times'
          WRITE(*,*)
          error_flag = -1
          RETURN
        ENDIF
        WRITE(*,*) 
        WRITE(*,*)'  *************************************************'
        WRITE(*,*)'  *** It appears that the data have irregularly ***'
        WRITE(*,67)NINT
   67    FORMAT('   *** spaced values/gaps at ',I7,'      levels.***')
        WRITE(*,68)ZINTMIN
   68    FORMAT('   *** Minimum interval detected = ',F12.6,'  ***')
        WRITE(*,69)ZINTMAX
   69    FORMAT('   *** Maximum interval detected = ',F12.6,'  ***')
        WRITE(*,*)'  ***                                     ***'
        WRITE(*,*)'  ***  If the data should be uniformly spaced,  ***'
        WRITE(*,*)'  ***  you should interpolate the data, check   ***'
        WRITE(*,*)'  ***  the format, or correct any observation   ***'
        WRITE(*,*)'  ***  time errors before spectral analysis.    ***'
        WRITE(*,*)'  *************************************************'
      ENDIF

      WRITE(*,*)
      WRITE(*,*)' The first & last lines of data have been read as:'
      WRITE(*,*)
   70  FORMAT('  Position(1) = ',F13.6,' Variable(1) = ',F11.5)
   71  FORMAT('  Position(N) = ',F13.6,' Variable(N) = ',F11.5)
      WRITE(*,70)STRAT(1),X(1)
      WRITE(*,71)STRAT(N),X(N)
   72  FORMAT('  Number of observations = ',I7)
   73  FORMAT('  *** Number of observation times >',I7,' ***')
      WRITE(*,*)
      WRITE(*,72)N
      IF(N > NMAX-1)THEN
        WRITE(*,*)
        WRITE(*,73)NMAX-1
        WRITE(*,*)' *** Shorten dataset &/or check the format ***'
        error_flag = -2
        RETURN
      ENDIF
      IF(N < NMIN)THEN
        WRITE(*,*)
        WRITE(*,*)' *** The number of observation times is < NMIN ***'
        WRITE(*,*)' *** Lengthen data set &/or check the format ***'
        error_flag = -3
        RETURN
      ENDIF

!  Check there is variability in the data - otherwise least deviation

!  fitting will result in an infinite loop.

      NVARY=0
      DO 75 I=2,N
        IF(X(1) /= X(I))NVARY=1
   75  ENDDO
      IF(NVARY == 0)THEN
        WRITE(*,*)
        WRITE(*,*)'    **************************************'
        WRITE(*,*)'    ***                            ***'
        WRITE(*,*)'    ***  No variability in the data !  ***'
        WRITE(*,*)'    ***                            ***'
        WRITE(*,*)'    **************************************'
        WRITE(*,*)
        error_flag = -4
        RETURN
      ENDIF

!  Option to leave program for data file editing if incorrect No. points.

      XVMIN=X(1)
      XVMAX=X(1)
      DO 80 I=1,N
        IF(X(I) < XVMIN)XVMIN=X(I)
        IF(X(I) > XVMAX)XVMAX=X(I)
   80  END DO
   81  FORMAT('  Maximum variable value = ',F10.5)
   82  FORMAT('  Minimum variable value = ',F10.5)

!   Garry Hayman
!   January 2014
!   Commented out DO loop 85 as assumed will proceed

!     CTLIM=0
!     DO 85
!        CTLIM=CTLIM+1
!        IF(CTLIM > 5)EXIT                     ! EXIT DO 85
!        WRITE(*,*)
!        WRITE(*,81)XVMAX
!        WRITE(*,82)XVMIN
!        WRITE(*,*)
!        WRITE(*,*)' Do you want to stop the program here?'
!        WRITE(*,*)' (e.g. if the no. of points is incorrect)'
!        IF(NIRR == 0)THEN
!           WRITE(*,*)' (or the data are meant to be uniformly spaced)'
!        ENDIF
!        WRITE(*,*)' Proceed          = P'
!        WRITE(*,*)' Stop the program = S'
!        READ(2,'(41A,A)',IOSTAT=IOS)HEADER,LMN
!        WRITE(*,'(A)')LMN
!        IF(LMN == 'P'.OR.LMN == 'p')EXIT      ! EXIT DO 85
!        IF(LMN == 'S'.OR.LMN == 's')STOP
!        IF(IOS /= 0)THEN
!           WRITE(*,*)
!           WRITE(*,*)' *** Enter P or S ***'
!           CYCLE                              ! CYCLE DO 85
!        ENDIF
! 85  ENDDO

!   Garry Hayman
!   January 2014
!   Commented out DO loops 86-89 as using irregularly-spaced data

! Option to treat slightly irregularly-spaced data as regularly-spaced.

!     IF(NIRR == 1.AND.N > 5000.AND.(ZINTMIN/ZINTMAX > 0.75))THEN
!        CTLIM=0
!        DO 86
!           CTLIM=CTLIM+1
!           IF(CTLIM > 5)EXIT                     ! EXIT DO 86
!           WRITE(*,*)
!           WRITE(*,*)' Treat the data as fixed spacing? (Y or N):'
!           READ(*,'(A)',IOSTAT=IOS)FIX
!           IF(IOS /= 0)THEN
!              WRITE(*,*)
!              WRITE(*,*)' *** Enter Y or N ***'
!              CYCLE                              ! CYCLE DO 87
!           ENDIF
!           IF(FIX == 'N'.OR.FIX == 'n')EXIT      ! EXIT 86
!           IF(FIX == 'Y'.OR.FIX == 'y')THEN
!              NIRR=0
!              EXIT                               ! EXIT DO 86
!           ENDIF
! 86     ENDDO
!        IF(FIX == 'Y'.OR.FIX == 'y')THEN
!           CTLIM=0
!           DO 87
!              CTLIM=CTLIM+1
!              IF(CTLIM > 5)EXIT                    ! EXIT DO 87
!              WRITE(*,*)
!              WRITE(*,*)' Enter the fixed sample interval:'
!              READ(*,*,IOSTAT=IOS)XINTFIX
!              IF(IOS == 0.AND.XINTFIX > 0.0d0)THEN
!                 EXIT                              ! EXIT DO 87
!              ELSE
!                 WRITE(*,*)
!                 WRITE(*,*)' *** Enter a positive number ***'
!                 CYCLE                             ! CYCLE DO 87
!              ENDIF
!  87       ENDDO
!           XINT=XINTFIX
!           ZINTMIN=XINTFIX
!           ZINTMAX=XINTFIX
!           IF(STRAT(1) < STRAT(2))THEN
!              STRAT1=STRAT(1)
!              DO 88 I=2,NROWS
!                 STRAT(I)=STRAT1+XINT
!  88          ENDDO
!           ELSE
!              DO 89 I=2,NROWS
!                 STRAT(I)=STRAT1-XINT
!  89          ENDDO
!           ENDIF
!        ENDIF
!     ENDIF

!   Observation time/stratigraphic units.

      WRITE(*,*)
      WRITE(*,*)' Indicate time/stratigraphic units:'
      WRITE(*,*)'    cm=1, kyr=3,     yr(Age)=5, Year(AD)=7'
      WRITE(*,*)'     m=2, Myr=4, Day of year=6,  Day=>Yr=8'
!     READ(2,'(41X,I1)',IOSTAT=IOS)KMK
      WRITE(*,'(I1)')KMK
      IF(IOS == 0.AND.(KMK == 1.OR.KMK == 2.OR.KMK == 3.OR. &
         KMK == 4.OR.KMK == 5.OR.KMK == 6.OR.KMK == 7.OR. &
         KMK == 8))THEN
      ENDIF
      IF(KMK == 4)THEN
         WRITE(*,*)
         WRITE(*,*)' NOTE:'
         WRITE(*,*)' Frequency will be reported in: Cycles per kyr'
         WRITE(*,*)'             Period will be reported in: kyr'
      ENDIF
      IF(KMK == 8)THEN
         WRITE(*,*)
         WRITE(*,*)' NOTE:'
         WRITE(*,*)' Frequency will be reported in: Cycles per year'
         WRITE(*,*)'             Period will be reported in: Year'
      ENDIF

!   Check time/age scale increases upwards.

      STRATD=STRAT(1)-STRAT(2)
      IF(KMK == 6.OR.KMK == 7.OR.KMK == 8.AND.STRATD < 0.0)THEN
        WRITE(*,*)
        WRITE(*,*)' ***           WARNING              ***'
        WRITE(*,*)' *** Implications for cross spectrum.  ***'
        WRITE(*,*)' *** Time values increasing down file. ***'
        WRITE(*,*)' *** Reverse order of time series.     ***'
      ENDIF
      IF(STRATD < 0.0)STRATD=-1.0*STRATD

!   For centimetres scale convert sample interval to metres.
!   For Myr scale convert sample interval to kyr.

      IF(KMK == 7)KMK=5
      IF(KMK == 1)THEN
        ZINTMIN=ZINTMIN/100.0
        STRATD=STRATD/100.0
        DO 97 I=1,N
          STRAT(I)=STRAT(I)/100.0
   97    ENDDO
      ENDIF
      IF(KMK == 4)THEN
        ZINTMIN=ZINTMIN*1000.0
        STRATD=STRATD*1000.0
        DO 98 I=1,N
          STRAT(I)=STRAT(I)*1000.0
   98    ENDDO
      ENDIF
      IF(KMK == 8)THEN
        ZINTMIN=ZINTMIN/365.25
        STRATD=STRATD/365.25
        DO 99 I=1,N
          STRAT(I)=STRAT(I)/365.25
   99    ENDDO
      ENDIF 
      XINT=STRATD
      XMAXT=STRAT(1)
      XMINT=STRAT(N)

!   Find total time/stratigraphic range of the data.

      XMAX=STRAT(1)
      XMIN=STRAT(1)
      DO 100 J=1,N
        IF(STRAT(J) > XMAX)XMAX=STRAT(J)
        IF(STRAT(J) < XMIN)XMIN=STRAT(J)

  100  ENDDO
      STDIF=XMAX-XMIN
      STAVE=0.5*(XMAX+XMIN)

      WRITE(*,*)
      WRITE(*,*)'   *************** Spectrum **********************'
      WRITE(*,*)'   *** Calculated using linear detrending,     ***'
      WRITE(*,*)'   *** tapering, then the Lomb-Scargle method  ***'
      WRITE(*,*)'   *** for regularly and irregularly spaced    ***'
      WRITE(*,*)'   *** data. Smoothing uses three applications ***'
      WRITE(*,*)'   *** of a discrete Hanning spectral window.  ***'
      WRITE(*,*)'   ***                                    ***'
      WRITE(*,*)'   *** The spectrum will be normalised for     ***'
      WRITE(*,*)'   *** comparison of different variables (i.e. ***'
      WRITE(*,*)'   *** the spectral estimates will be divided  ***'
      WRITE(*,*)'   *** by the sum of the power).             ***'
      WRITE(*,*)'   ***********************************************'
      WRITE(*,*)
      AN=N
      ANORIG=AN
      NORIG=N

!   Determine average sample interval (XINT) for irregularly spaced data.

      IF(NIRR >= 1.OR.CFIX == 'Y'.OR.CFIX == 'y')THEN
        IF(XMAXT-XMINT < 0.0)THEN
          XINT=(XMINT-XMAXT)/(AN-1)
        ELSE
          XINT=(XMAXT-XMINT)/(AN-1)
        ENDIF
        WRITE(*,*)
  101    FORMAT('      Minimum sample interval = ',F13.6)
        WRITE(*,101)ZINTMIN
      ENDIF
  102 FORMAT('         Mean sample interval = ',F13.6)
      WRITE(*,102)XINT

!   For SI = constant: define Nyquist frequency and frequency
!           interval using SI (XINT) and find Rayleigh frequency (XFINT).
!   For SI = variable: define Nyquist frequency using minimum
!           sample interval (ZINTMIN), but use Rayleigh frequency
!           determined earlier (XFINT1).

      IF(NIRR == 0)THEN
        XNYQ=1.0/(2.0*XINT)
      ELSE
        XNYQ=1.0/(2.0*ZINTMIN)
      ENDIF
      XNYQ1=XNYQ
      XFINT=XNYQ/(ANSTART/2.0)
      XFINT1=XFINT
      NOUT=int((XNYQ/XFINT1)+0.00001)
      NOUT1=NOUT
      WRITE(*,*) NOUT,XNYQ,XFINT1,XNYQ/XFINT1,ANSTART

!   Change Rayleigh frequency.

  103 FORMAT(' Indicated Rayleigh Frequency = ',F13.6)
      WRITE(*,*)
      WRITE(*,103)XFINT1
      WRITE(*,*)
  104 FORMAT('  Indicated Nyquist Frequency = ',F13.6)
      WRITE(*,104)XNYQ
      WRITE(*,*)
      WRITE(*,*)' Change the Rayleigh Frequency? (Y or N):'
!     READ(2,'(41A,A)',IOSTAT=IOS)HEADER,PRAY
      WRITE(*,'(A)')PRAY
      IF(PRAY == 'N'.OR.PRAY == 'n')THEN
        XNRAY=XFINT1
      ENDIF
      IF(PRAY == 'Y'.OR.PRAY == 'y')THEN
        XF4=XFINT1*4.0
        XF2=XFINT1*0.5
        WRITE(*,*)      
        WRITE(*,*)' Enter choice of Rayleigh Frequency:'
  106   FORMAT('  (between ',F13.6,'  &',F13.6,')')
        WRITE(*,106)XF2,XF4
!       READ(2,*,IOSTAT=IOS2)XNRAY
        XFINT=XNRAY
        XFINT1=XNRAY
        NOUT=NOUT1
        XNYQ=NOUT*XNRAY
        XNNYQ=XNYQ
        XNYQNEW=XNYQ
        PNYQ='N'
        WRITE(*,*) XNRAY,XNYQ,NOUT,NOUT1
      ENDIF

!   Set Nyquist Frequency.
      
      IF(PRAY == 'y'.OR.PRAY == 'Y')THEN
        WRITE(*,*)
        WRITE(*,104)XNYQ
      ENDIF
      WRITE(*,*)
      XNYQNEW=XNYQ
      WRITE(*,*)'  Change the Nyquist Frequency? (Y or N):'
!     READ(2,'(41A,A)',IOSTAT=IOS)HEADER,PNYQ
      WRITE(*,'(A)',IOSTAT=IOS)PNYQ
      IF(PNYQ == 'N'.OR.PNYQ == 'n')THEN
         XNNYQ=XNYQ
      ENDIF
      IF(PNYQ == 'Y'.OR.PNYQ == 'y')THEN
         XNYN=XNYQ1/10.0
         XNNYQ=XNYQ
         IF(XNYN < XFINT*30.0) XNYN=XNYQ1/5.0
         IF(XNYN < XFINT*30.0) XNYN=XNYQ1/2.0
            XNYQNEW=XNNYQ
      ENDIF
 
!   Check new values for implied number of frequencies.

  120  NOUTN=XNNYQ/XNRAY
      NOUT=NOUTN
      IF(NOUTN > NMAXB)THEN
  121    WRITE(*,*)
        WRITE(*,*)' *** Chosen Rayleigh and Nyquist Freqs   ***'
        WRITE(*,*)' *** imply no. spectral estimates exceed ***'
  122    FORMAT('    *** the array sizes of ',I7,' ***')
        WRITE(*,122)NMAXB
        WRITE(*,*)' Applying automatic fix...'
!   Automatic fix for number of spectral estimates. Using chosen Rayleigh
!   frequency (XFINT=XNRAY), set NOUT to be largest possible given size
!   of program arrays.
        NOUT=NMAXB-1
        NOUTN=NOUT
        XNYQ=(NMAXB-1)*XNRAY
        XNNYQ=XNYQ
        WRITE(*,*)
  123    FORMAT('    Revised Nyquist frequency = ',F13.6)
        WRITE(*,123)XNYQ
      ENDIF

      NEW=NOUT*2
      XNEW=NEW
      NDOF=8
      XNDOF=8.0
      MN=NOUT+1
      AMN=MN

!   Specify N and AN here as needed for treatment of AR1 data later.

      N=NSTART
      AN=NSTART

!   For the Hanning window bandwidth (allowing for tapering of data and
!   using Rayleigh Frequency (XFINT), bandwidth (BW):
!          BW = 1 / ( I x 1.055 x N x sample interval).
!   This value is scaled according to choice of Rayleigh frequency (XNRAY).

      BW=1.0/(0.225586*1.055*ANSTART*XINT)
      BW=(XNRAY/XFINT)*BW
      WRITE(*,*)
  124  FORMAT('          Rayleigh Frequency = ',F13.6)
      WRITE(*,124)XFINT1
  125  FORMAT('        Resolution Bandwidth = ',F13.6)
      WRITE(*,125)BW

      NARSIM=0

!   Label for jump from AR1 series generation (Nb. NARSIM=1)

  200  CONTINUE

!   Linear detrending of data using robust least deviation line fitting 
!   routine MEDFIT of Press et al. (1992, p699).
!   XT = STRAT position, YT = Variable.

      NDATA=N
      SX=0.0
      SY=0.0
      SXY=0.0
      SXX=0.0
      DO 220 J=1,N
        XT(J)=STRAT(J)
        YT(J)=X(J)
        SX=SX+STRAT(J)
        SY=SY+X(J)
        SXY=SXY+STRAT(J)*X(J)
        SXX=SXX+STRAT(J)**2
  220  ENDDO
      NDATAT=NDATA
      DEL=NDATA*SXX-SX**2
      AA=(SXX*SY-SX*SXY)/DEL
      BB=(NDATA*SXY-SX*SY)/DEL
      CHISQ=0.0
      DO 221 J=1,NDATA
        CHISQ=CHISQ+(X(J)-(AA+BB*STRAT(J)))**2
  221  ENDDO
      SIGB=SQRT(CHISQ/DEL)
      B1=BB
      F1=ROFUNC(B1)
      B2=BB+SIGN(3.0*SIGB,F1)
      F2=ROFUNC(B2)
      CTLIM=0
      DO 222
        CTLIM=CTLIM+1
        IF(CTLIM > N)WRITE(*,*)' *** EXIT DO 222 CTLIM>N ***'
        IF(CTLIM > N)EXIT                 ! EXIT DO 222
        IF(F1*F2 > 0.0)THEN
          BB=2.0*B2-B1
          B1=B2
          F1=F2
          B2=BB
          F2=ROFUNC(B2)
          CYCLE                         ! CYCLE DO 222
        ENDIF
        EXIT                            ! EXIT DO 222
  222  ENDDO
      SIGB=0.01*SIGB
      CTLIM=0
      DO 223
        CTLIM=CTLIM+1 
        IF(CTLIM > N)WRITE(*,*)' *** EXIT DO 223 CTLIM>N ***'
        IF(CTLIM > N)EXIT                 ! EXIT DO 223
        IF(ABS(B2-B1) > SIGB)THEN
          BB=0.5*(B1+B2)
          IF(BB == B1.OR.BB == B2)EXIT      ! EXIT DO 223
          F=ROFUNC(BB)
          IF(F*F1 >= 0.0)THEN
            F1=F
            B1=BB
          ELSE
            F2=F
            B2=BB
          ENDIF
          CYCLE                         ! CYCLE DO 223
        ENDIF
        EXIT                            ! EXIT DO 223
  223  ENDDO
  224  XINTERC=AA
      SLOPE=BB
      ABDEV=ABDEVT/NDATA

!   End of line fitting.

!   Find the fitted line and subtract this from the observed variable
!   values in order to execute linear detrend.

      DO 226 I=1,N
        XFIT=(SLOPE*STRAT(I))+XINTERC
        X(I)=X(I)-XFIT
  226  ENDDO

!   If first pass through:
!   Calculate the average, variance (p611, Press et al., 1992)
!   and lag-1 autocorrelation of the time series .

      IF(NARSIM == 0)THEN
        AVE=0.0
        DO 227 I=1,N
          AVE=AVE+X(I)
  227    ENDDO
        AVE=AVE/ANSTART
        VARX=0.0
        EP=0.0
        SSDIFF=0.0
        DO 228 I=1,N      
          SSDIFF=X(I)-AVE
          EP=EP+SSDIFF
          VARX=VARX+SSDIFF*SSDIFF
  228    ENDDO
        VARX=(VARX-EP**2/ANSTART)/(ANSTART-1.0)
        VARX1=VARX
        STDEV=SQRT(VARX)
!   Calculate raw autocorrelation lag(1) following p238-239, Davis, 1973.   
        I=2      
        SX=0.0
        SY=0.0
        SXY=0.0
        L=N-I+1
        DO 229 J=1,L
          K=J+I-1
          SX=SX+X(J)
          SY=SY+X(K)
          SXY=SXY+X(J)*X(K)
  229    ENDDO
        AL=L
        ROH=((AL*SXY-SX*SY)/(AL*(AL-1.0)))/VARX
!
!   If ROH estimate is less than zero then search for an outlier,
!   replace with average and repeat ROH calculation.
!   This design only looks for a single outlier and retains all
!   original data in the spectrum calculations.
        WRITE(*,'(10(1PE14.6))') ROH
        IF(ROH < 0.0d0)THEN
          OUTI=0
          DO 9230 I=1,N
            XLIM=AVE-(3.0*STDEV)
            IF(X(I) < XLIM)THEN
             OUTI=I
             WRITE(*,*)
             WRITE(*,*)' *** Outlier detected at line no.:',OUTI
             EXIT                   ! EXIT DO 9230
            ENDIF
            XLIM=AVE+(3.0*STDEV)
            IF(X(I) > XLIM)THEN
             OUTI=I
             WRITE(*,*)
             WRITE(*,*)' *** Outlier detected at line no.:',OUTI
             EXIT                   ! EXIT DO 9230
            ENDIF
 9230      ENDDO
          IF(OUTI > 0)THEN
            DO 9232 I=1,N
             X2(I)=X(I)
 9232       ENDDO
            X2(OUTI)=AVE 
            I=2   
            SX=0.0
            SY=0.0
            SXY=0.0
            L=N-I+1
            DO 9234 J=1,L
             K=J+I-1
             SX=SX+X2(J)
             SY=SY+X2(K)
             SXY=SXY+X2(J)*X2(K)
 9234       ENDDO
            AL=L
            ROH=((AL*SXY-SX*SY)/(AL*(AL-1.0)))/VARX            
          ENDIF
        ENDIF
        IF(OUTI == 0.AND.ROH < 0.0d0)THEN
          WRITE(*,*)
          WRITE(*,*)' *** ROH < 0.0, but outlier search failed ***'
        ENDIF
        IF(ROH < 0.0d0)ROH=0.0d0
        ROHI=ROH
!   Estimate TAU using direct/raw estimate of RHO provided
!   it is greater than 0.0.
        SROH=ROH
        IF(ROH /= 0.0d0)THEN
          IF(NIRR /= 0)THEN
            TAUMIN=-1.0*ZINTMIN/ALOG(SROH)
          ENDIF
          TAU=-1.0*XINT/ALOG(SROH)
        ENDIF
        WRITE(*,*)
  240    FORMAT('    RHO direct (raw) estimate = ',F5.3)
        WRITE(*,240)ROH
!
!   For SI = variable estimate ROH using procedure of Mudelsee (2002).
!   Copy data to XSCAL and TSCAL, order must be reversed for STRAT
!   if not height data (i.e. reverse if Age or Depth data).
        IF(NIRR == 1)THEN
          IF(ROH == 0.0d0)THEN
            SCALT=-LOG(0.999)/XINT
          ELSE
            SCALT=-LOG(ROH)/XINT
          ENDIF  
          IF(STVAL >= 0.0)THEN
            DO 230 I=1,N
             XSCAL(I)=X(I)/SQRT(VARX)
             TSCAL(I)=-1.0*STRAT(I)*SCALT
  230       ENDDO
          ELSE
            DO 231 I=1,N
             XSCAL(I)=X(N+1-I)/SQRT(VARX)
             TSCAL(I)=(-1.0*STRAT(N+1-I))*SCALT
  231       ENDDO
          ENDIF
!   Find least square minimization three times based on three starting
!   conditions via BRENT minimization of XLS function.
          MUDEST=0
          DUM1=BRENT(-2.0D0,AR1,+2.0D0,XLS,TOL,AR11,TSCAL,XSCAL,N)
          DUM2=BRENT(AR1,0.5D0*(AR1+1.0D0),+2.0D0,XLS,TOL,AR12,TSCAL, &
         XSCAL,N)
          DUM3=BRENT(-2.0D0,0.5D0*(AR1-1.0D0),AR1,XLS,TOL,AR13,TSCAL, &
         XSCAL,N)
!   Find smallest sum of squares from the minimizations. Use corresponding
!   value of AMIN.
          IF(AR11 == -1.0D0)AR11=ROH
          IF(AR12 == -1.0D0)AR12=ROH
          IF(AR13 == -1.0D0)AR13=ROH
          DUMMIN=DUM1
          AMIN=AR11
          IF(DUMMIN > DUM2)AMIN=AR12
          IF(DUMMIN > DUM3)AMIN=AR13
          WRITE(*,'(10(1PE14.6))') DUM1,DUM2,DUM3,DUMMIN, &
                              AMIN,AR11,AR12,AR13,ROH
!   Check output of XLS and BRENT and implied Rho.
          NMU=0
          IF((ABS(AR12-AR11) > TOL2.AND.ABS(AR12-AR1) > TOL2) &
           .OR.(ABS(AR13-AR11) > TOL2.AND.ABS(AR13-AR1) > TOL2)) &
           NMU=1
          IF(NMU == 1)THEN
            WRITE(*,*)
            WRITE(*,*)' *** Least squares has multiple minima. ***'
          ENDIF
          IF(AMIN <= 0.0)THEN
            WRITE(*,*)
            WRITE(*,*)' *** Estimation problem: AMIN =< 0.0    ***'
          ENDIF
          IF(AMIN > 1.0)THEN
            WRITE(*,*)
            WRITE(*,*)' *** Estimation problem: AMIN => 1.0    ***'
          ENDIF
!         IF(NMU /= 1.OR.AMIN > 0.0.AND.AMIN <= 1.0)THEN
          IF(NMU /= 1.AND.AMIN > 0.0.AND.AMIN <= 1.0)THEN
            TAU=-1.0/SCALT*LOG(AMIN)
            ROHN=EXP(-1.0*XINT/TAU)
            IF(ROHN > 0.0.AND.ROHN <= 1.0)THEN
             MUDEST=1
             WRITE(*,*)
             WRITE(*,*)'   RHO estimate from Mudelsee (2002)'
  241         FORMAT('                  procedure = ',F5.3)
             WRITE(*,241)ROHN
            ELSE
             WRITE(*,*)
             WRITE(*,*)' *** Implied Rho outside 0.0 to 1.0      ***'
             WRITE(*,*)' *** Using Rho from direct/raw estimate. ***'
             ROHN=ROH
            ENDIF
          ELSE
            WRITE(*,*)
            WRITE(*,*)' *** Using Rho from direct/raw estimate. ***'
            ROHN=ROH
          ENDIF
          IF(ROHI > ROHN)ROH=ROHN
        ENDIF
      ENDIF

!   Data tapering using adaptation of algorithm from Bloomfield 
!   (p116, 1976). Taper first and last 5% of time series using a 
!   split cosine bell. This prevents "periodogram leakage" caused 
!   by the abrupt ends of the time series.

!   Modification of Bloomfield procedure, due to use of irregularly
!   spaced-data, is to find position of each data point 
!   relative to the full length of the data. If a data point lies
!   between proportions of 0.00 and 0.05 or between 0.95 and 1.00 
!   of the data span, the variable value is weighted (tapered).

      XPI=3.141592654
      STDIFF=STRAT(N)-STRAT(1)
      IF(STDIFF > 0.0)THEN
        BIT=STDIFF/(10.0*AN)
        DO 252 I=1,N
          XPROP=(STRAT(N)-STRAT(I)+BIT)/STDIFF
          IF(XPROP < 0.05.OR.XPROP > 0.95)THEN
            WEIGHT=0.5 - 0.5*COS(XPI*(XPROP/0.05))
            X(I)=X(I)*WEIGHT
          ENDIF
  252    ENDDO
      ELSE
        STDIFF=STRAT(1)-STRAT(N)
        BIT=STDIFF/(10.0*AN)
        DO 254 I=1,N
          XPROP=(STRAT(I)-STRAT(N)+BIT)/STDIFF
          IF(XPROP < 0.05.OR.XPROP > 0.95)THEN
            WEIGHT=0.5 - 0.5*COS(XPI*XPROP/0.05)
            X(I)=X(I)*WEIGHT
          ENDIF
  254    ENDDO
      ENDIF

!   Zero-padding: if zeroes need to be added to allow chosen
!   Rayleigh frequency to be applied out to chosen Nyquist frequency.

      IF((PNYQ == 'y'.OR.PNYQ == 'Y').OR. &
         (PRAY == 'Y'.OR.PRAY == 'y'))THEN     ! Nyq &/or Ray freqs changed
        IF(NOUT > NOUT1.OR.XNYQ > XNYQ1)THEN  ! Nyq freq &/or NOUT increased
          WRITE(*,*)' Zero padding'
          SDIFF=STRAT(1)-STRAT(2)
!   Depth or age scale:
            IF(SDIFF <= 0.0)THEN
             DO 256 I=NSTART+1,NEW
               AI=I-NSTART
               X(I)=0.0
               IF(NARSIM == 0)STRAT(I)=STRAT(NSTART)+(AI*XINT)
  256         ENDDO
!   Height or time scale:
            ELSE
             J=N
             DO 258 I=1,NEW-NSTART
               J=J+1
               AI=I-NSTART
               X(J)=0.0
               IF(NARSIM == 0)STRAT(J)=STRAT(NSTART)-(AI*XINT)
  258         ENDDO
            ENDIF                                                                                                                                       
!   Reset N following zero-padding.      
          N=NEW
          AN=NEW
        ENDIF
      ENDIF

!   Find middle of stratigraphic range after zero-padding (STAVE).

      XMAX=STRAT(1)
      XMIN=STRAT(1)
      DO 259 J=1,N
        IF(STRAT(J) > XMAX)XMAX=STRAT(J)
        IF(STRAT(J) < XMIN)XMIN=STRAT(J)
  259  ENDDO
      STAVE=0.5*(XMAX+XMIN)

!   Lomb-Scargle method ("PERIOD") of Fourier transform for irregularly 
!   spaced data (Press et al., 1992, pp572-574). This is not an FFT, so 
!   the routine is slow. As an alternative using large arrays FASPER on
!   pp575-576 of Press et al. is much faster (but large RAM is needed).
!   Number of spectral estimates, MN = N/2 +1.
!   Note here array "y" in "PERIOD" means array X and array "x" in 
!   "PERIOD" means array STRAT.

!   Calculate the average and variance of the detrended, tapered
!   ad/or zero-padded time series (p611, Press et al., 1992).

      AVE=0.0
      DO 260 I=1,N
        AVE=AVE+X(I)
  260  ENDDO
      AVE=AVE/ANSTART
      VARX=0.0
      EP=0.0
      SSDIFF=0.0
      DO 262 I=1,N      
        SSDIFF=X(I)-AVE
        EP=EP+SSDIFF
        VARX=VARX+SSDIFF*SSDIFF
  262  ENDDO
      VARX=(VARX-EP**2/ANSTART)/(ANSTART-1.0)

!   Start transform and calculation of frequencies. PNOW defines frequency.
!   Note XFINT defined differently to Press et al. algorithm (i.e. not
!   PNOW=1.0/STDIF).

      PNOW=XFINT
      AN=N
      DUMMA=-2.0
      DUMMB=0.5
      DO 263 J=1,N
        ARG=TWOPID*((STRAT(J)-STAVE)*PNOW)
        WPR(J)=DUMMA*SIN(DUMMB*ARG)**2
        WPI(J)=SIN(ARG)
        WR(J)=COS(ARG)
        WI(J)=WPI(J)
  263  ENDDO
      DO 266 I=2,MN
        FREQ(I)=PNOW
        SUMSH=0.0
        SUMC=0.0
        DO 264 J=1,N
          CCOS=WR(J)
          SSIN=WI(J)
          SUMSH=SUMSH+SSIN*CCOS
          SUMC=SUMC+(CCOS-SSIN)*(CCOS+SSIN)
  264    ENDDO
        WTAU=0.5*ATAN2(2.0*SUMSH,SUMC)
        SWTAU=SIN(WTAU)
        CWTAU=COS(WTAU)
        SUMS=0.0
        SUMC=0.0
        SUMSY=0.0
        SUMCY=0.0
        NCC=0
        NSS=0
        DO 265 J=1,N
          SSIN=WI(J)
          CCOS=WR(J)
          SS=SSIN*CWTAU-CCOS*SWTAU
          CC=CCOS*CWTAU+SSIN*SWTAU
          SUMS=SUMS+SS**2
          SUMC=SUMC+CC**2
          YY=X(J)-AVE
          SUMSY=SUMSY+YY*SS
          SUMCY=SUMCY+YY*CC
          WTEMP=WR(J)
          WR(J)=(WR(J)*WPR(J)-WI(J)*WPI(J))+WR(J)
          WI(J)=(WI(J)*WPR(J)+WTEMP*WPI(J))+WI(J)
  265    ENDDO
!       IF(CPOWER == 'R'.OR.CPOWER == 'r') THEN
           P(I)=0.5*(SUMCY**2/SUMC+SUMSY**2/SUMS)/VARX
!       ELSE
!          P(I)=0.5*(SUMCY**2/SUMC+SUMSY**2/SUMS)
!       END IF
        PNOW=PNOW+XFINT
  266  ENDDO

      IF(NARSIM > 0)GO TO 410

!   Write raw periodogram to file.
  268  FORMAT('  Frequency       Power')
  269  FORMAT(2X,F10.6,1X,F15.8)
      OPEN(UNIT=10,FILE='specpow')
      WRITE(10,268)
      DO 270 I=2,MN
        WRITE(10,269)FREQ(I),P(I)
  270  END DO
      CLOSE(UNIT=10)

!   Smooth spectrum using a 3-point Hanning window (0.25, 0.5, 0.25, 
!   end points: 0.5, 0.5). 
!   Degrees of freedom = 2/I, where I = sum of squares of final window
!   weights (ie. 0.015625, 0.09375, 0.234375, 0.3125, 0.234375, 0.09375,
!   0.015625. I = 0.225586).
!   Since data were tapered by a total of 10% this needs correction factor
!   (Bloomfield, 1976 pp193-194). Hence DOF = 2.0 / 1.055 x 0.225586 = 8.4
!   Effective degrees of freedom = 8 (used for confidence levels).
!   Zero-padding effect means DOF = 2.0 x N / 1.055 x 0.225586 x NEW
!   where N is number of data points and NEW is number of time series
!   points after zero-padding.

      IF(NEW > NSTART)THEN
        XNDOF=(2.0*ANSTART)/(1.055*0.225586*XNEW)
        NDOF=XNDOF

      ENDIF

      DO 273 K=1,3
        PSMTH(1)=P(1)
        PSMTH(2)=(0.5*P(2))+(0.5*P(3))
        PSMTH(MN)=(0.5*P(MN-1))+(0.5*P(MN))
        DO 271 I=3,MN-1
          PSMTH(I)=(0.25*P(I-1))+(0.5*P(I))+(0.25*P(I+1))
  271    ENDDO
        DO 272 I=1,MN
          P(I)=PSMTH(I)
  272    ENDDO
  273  ENDDO

!  Calculate normalized spectral estimates to allow comparison of
!  spectra for time series with different average amplitudes and
!  location of spectral background.

      WRITE(*,*)
      WRITE(*,*)' Absolute or relative power? Abs. = A, Rel. = R:'
!     READ(2,'(41A,A)',IOSTAT=IOS)HEADER,CPOWER
      WRITE(*,'(A)') CPOWER

      TOTP=0.0
      DO 274 I=2,MN
         TOTP=TOTP+(PSMTH(I))
  274  ENDDO

      IF(NARSIM == 0)THEN
        TOTP1=TOTP
      ENDIF

      DO 275 I=2,MN
        PSMTH(I)=(PSMTH(I))/TOTP
  275  ENDDO
!
      WRITE(*,*)
      WRITE(*,*)' Choose the form of the spectral background:' 
      WRITE(*,*)'                    White noise model = 0'
      WRITE(*,*)'                         AR(1) model = 1' 
      WRITE(*,*)'                        Quadratic fit = 2'
      IF(NIRR == 0)THEN
         WRITE(*,*)'                        Power law fit = 3'
!        READ(2,'(41X,I1)',IOSTAT=IOS)NFIT
         WRITE(*,'(I1)')NFIT
         IF(IOS /= 0.OR.(NFIT /= 0.AND.NFIT /= 1.AND.NFIT /= 2.AND. &
           NFIT /=3))THEN
           WRITE(*,*)
           WRITE(*,*)' *** Enter 0, 1, 2 or 3 ***'
         ENDIF
      ELSE
!        READ(2,'(41X,I1)',IOSTAT=IOS)NFIT
         WRITE(*,'(I1)')NFIT
         IF(IOS /= 0.OR.(NFIT /= 0.AND.NFIT /= 1.AND.NFIT /= 2))THEN
            WRITE(*,*)
            WRITE(*,*)' *** Enter 0, 1, or 2 ***'
         ENDIF
      ENDIF

!   Power law fitting for confidence levels.

!   Put Log power values into XLOGA, put frequency values into XLFREQ

  280  IF(NFIT == 3)THEN
        WRITE(*,*)
        WRITE(*,*)' *****************************************'
        WRITE(*,*)' *** Spectral background location via: ***'
        WRITE(*,*)' *** Power-law fitting               ***'
        WRITE(*,*)' *****************************************'
        DO 282 I=2,MN
          XLOGA(I)=ALOG10(PSMTH(I))
          Y(I)=XLOGA(I)
          XLFREQ(I)=ALOG10(FREQ(I))
  282    ENDDO
!   Start robust least deviation line fitting routine of Press et al.
!   (1992, p699) for power law fit to log-log spectrum.
        NDATA=MN-1
        SX=0.0
        SY=0.0
        SXY=0.0
        SXX=0.0
        DO 284 J=2,MN
          XT(J)=XLFREQ(J)
          YT(J)=Y(J)
          SX=SX+XLFREQ(J)
          SY=SY+Y(J)
          SXY=SXY+XLFREQ(J)*Y(J)
          SXX=SXX+XLFREQ(J)**2
  284    ENDDO
        NDATAT=NDATA
        DEL=NDATA*SXX-SX**2
        AA=(SXX*SY-SX*SXY)/DEL
        BB=(NDATA*SXY-SX*SY)/DEL
        CHISQ=0.0
        DO 286 J=1,NDATA
          CHISQ=CHISQ+(Y(J)-(AA+BB*XLFREQ(J)))**2
  286    ENDDO
        SIGB=SQRT(CHISQ/DEL)
        B1=BB
        F1=ROFUNC(B1)
        B2=BB+SIGN(3.0*SIGB,F1)
        F2=ROFUNC(B2)
        CTLIM=0
        DO 288
          CTLIM=CTLIM+1
          IF(CTLIM > N)THEN
            WRITE(*,*)' *** EXIT DO 288 CTLIM>N ***'
          ENDIF
          IF(CTLIM > N)EXIT                  ! EXIT DO 288
          IF(F1*F2 > 0.0)THEN
            BB=2.0*B2-B1
            B1=B2
            F1=F2
            B2=BB
            F2=ROFUNC(B2)
            CYCLE                           ! CYCLE DO 288
          ENDIF
          EXIT                              ! EXIT DO 288
  288    ENDDO
        SIGB=0.01*SIGB
        CTLIM=0
        DO 290
          CTLIM=CTLIM+1
          IF(CTLIM > N)THEN
            WRITE(*,*)' *** EXIT DO 290 CTLIM>N ***'
          ENDIF
          IF(CTLIM > N)EXIT                  ! EXIT DO 223
          IF(ABS(B2-B1) > SIGB)THEN
            BB=0.5*(B1+B2)
            IF(BB == B1.OR.BB == B2)EXIT       ! EXIT DO 290
            F=ROFUNC(BB)
            IF(F*F1 >= 0.0)THEN
             F1=F
             B1=BB
            ELSE
             F2=F
             B2=BB
            ENDIF
          ELSE
            CYCLE                           ! CYCLE DO 290
          ENDIF
          EXIT                              ! EXIT DO 290
  290    ENDDO
        XINTERC=AA
        SLOPE=BB
        ABDEV=ABDEVT/NDATA
!   End of line fitting.
  292    FORMAT(' Intercept =',F12.8,'  Slope =',F12.8)
        WRITE(*,292)XINTERC,SLOPE
        IF(SLOPE > 0.0)THEN
          WRITE(*,*)
          WRITE(*,*)' *** Power law fitting failed    ***'
          WRITE(*,*)' *** (positive slope).          ***'
          WRITE(*,*)' *** Switching to Quadratic fit. ***'
          NFIT=2
          GO TO 324
        ENDIF
!   Find the background spectrum (i.e. the power-law/fitted line).
        DO 294 I=2,MN
          XLOUT(I)=(SLOPE*XLFREQ(I))+XINTERC
  294    ENDDO
        DO 296 I=2,MN
          XOUT(I)=10.0**XLOUT(I)
  296    ENDDO
      ENDIF

!   Calculation of confidence levels based on fit of log spectrum.
!   However, if one or more spectral estimates equal zero in spectral
!   background fitting interval, then estimates are reset to very
!   low levels artificially.

      IF(NFIT == 0.OR.NFIT == 1.OR.NFIT == 2)THEN
        ZTEST=0.0
        DO 300 I=2,MN
          IF(PSMTH(I) <= 0.0000000000000001)THEN
            PSMTH(I)=0.0000000001
            ZTEST=ZTEST+1.0
          ENDIF
  300    ENDDO
      ENDIF

!   Find log of power values and put into array Y excluding value for
!   Freq = 0.0 since this is the spectrum for detrended data.
!   Find mean log10 power value (SOLMEAN).
      NSPECB=0
!   Re-set spectra limits if changes to Nyquist frequency needed.
  304  IF(NFIT == 0.OR.NFIT == 1.OR.NFIT == 2)THEN
        IF(NIRR == 0.OR.NSPECB == 1.OR.NSPECB == 2)THEN
          IF(XNYQ /= XNYQNEW)THEN
            XNYQ=XNYQNEW
            NOUT=NOUTN
            MN=NOUT+1
            AMN=MN
          ENDIF        
        ENDIF
        J=1
        SUM=0.0
        DO 306 I=2,MN
          XLOGA(I)=ALOG10(PSMTH(I))
          SUM=XLOGA(I)+SUM
          J=J+1
  306    ENDDO
        SOLMEAN=SUM/(AMN-1.0)
      ENDIF

!   Median smoothing for AR1 or white noise fitting.

      IF(NFIT == 0.OR.NFIT == 1)THEN
        ROH2=ROH*ROH
        SO=VARX/(1.0-ROH2)
!   Function SELECT from Press et al.
!   1992, p334, to find median.
!   Note: allowing for proportion of maximum Nyquist frequency used in
!        calculation of maximum width of median smoothing window.
!   For SI = variable, determine the number of spectral estimates
!   out to Nyquist from SIMEAN (XINT). Temporarily set MN at this limit.
        NWMAX=NOUT1/4
        NWMIN=2.0*(BW/(1/(AN*ZINTMIN)))
        NWIND=NWMAX/3
        CTLIM=0
        DO 308
          CTLIM=CTLIM+1
          IF(CTLIM > NWMIN+1)THEN
            WRITE(*,*)' *** EXIT DO 308 CTLIM>NWMIN+1 ***'
          ENDIF
          IF(CTLIM > NWMIN+1)EXIT            ! EXIT DO 308
          IF(NWIND <= NWMIN)THEN
            NWIND=NWIND+1
            CYCLE                          ! CYCLE DO 308
          ENDIF
          EXIT                             ! EXIT DO 308
  308    ENDDO
        CTLIM=0
        DO 310
          CTLIM=CTLIM+1
          IF(CTLIM > N)THEN
            WRITE(*,*)' *** EXIT DO 310 CTLIM>N ***'
          ENDIF
          IF(CTLIM > N)EXIT                  ! EXIT DO 310
          ANWIND=NWIND
          NWIND2=NWIND/2
          ANWIND2=NWIND2+0.6
          NWIND2=ANWIND2
          ANWIND2B=(ANWIND/2.0)+0.6
          NWIND2B=ANWIND2B
          IF(NWIND2B == NWIND2)THEN
            NWIND=NWIND+1
            CYCLE                          ! CYCLE DO 310
          ENDIF
          EXIT                             ! EXIT DO 310
  310    ENDDO
        NSTAT=2
        NSTOP=1+NWIND
!   Position of median value in each XLP window = MEDPOS.
!   Position of median value in XMEDIAN list = MEDPOS2.
        MEDPOS=(NWIND+1)/2
        MEDPOS2=NSTAT+MEDPOS-1 
        MEDPOS2I=MEDPOS2
!   Put XLOGA data into XLP array according to position of median
!   smoothing window (defined by NSTAT AND NSTOP).
        CTLIM=0
        DO 319
          CTLIM=CTLIM+1
          IF(CTLIM > N)THEN
            WRITE(*,*)' *** EXIT DO 319 CTLIM>N ***'
          ENDIF
          IF(CTLIM > N)EXIT                  ! EXIT DO 319
          K=1
          DO 312 I=NSTAT,NSTOP
            XLP(K)=XLOGA(I)
            K=K+1
  312      ENDDO
!      Function SELECT (p334 of Press et al 1992)
!      using N=NWIND, ARR(I)=XLP(I), K=MEDPOS
          L=1
          IR=NWIND
          CTLIM1=0
          DO 318
            CTLIM1=CTLIM1+1
            IF(CTLIM1 > N)THEN
             WRITE(*,*)' *** EXIT DO 318 CTLIM1>N ***'
            ENDIF
            IF(CTLIM1 > N)EXIT              ! EXIT DO 318      
            IF(IR-L <= 1)THEN
             IF(IR-L == 1)THEN
               IF(XLP(IR) < XLP(L))THEN
                 TEMP=XLP(L)
                 XLP(L)=XLP(IR)
                 XLP(IR)=TEMP
               ENDIF
             ENDIF
             SELCT=XLP(MEDPOS)
             EXIT                         ! EXIT DO 318
            ELSE
             MID=(L+IR)/2
             TEMP=XLP(MID)
             XLP(MID)=XLP(L+1)
             XLP(L+1)=TEMP
             IF(XLP(L+1) > XLP(IR))THEN
               TEMP=XLP(L+1)
               XLP(L+1)=XLP(IR)
               XLP(IR)=TEMP
             ENDIF
             IF(XLP(L) > XLP(IR))THEN
               TEMP=XLP(L)
               XLP(L)=XLP(IR)
               XLP(IR)=TEMP
             ENDIF
             IF(XLP(L+1) > XLP(L))THEN
               TEMP=XLP(L+1)
               XLP(L+1)=XLP(L)
               XLP(L)=TEMP
             ENDIF
             I=L+1
             J=IR
             A=XLP(L)
             CTLIM2=0
             DO 317
               CTLIM2=CTLIM2+1
               IF(CTLIM2 > N)THEN
                 WRITE(*,*)' *** EXIT DO 317 CTLIM2>N ***'
               ENDIF
                IF(CTLIM2 > N)EXIT          ! EXIT DO 317
               CTLIM3=0
               DO 314
                 CTLIM3=CTLIM3+1
                 IF(CTLIM3 > N)THEN
                  WRITE(*,*)' *** EXIT DO 314 CTLIM3>N ***'
                 ENDIF
                 IF(CTLIM3 > N)EXIT         ! EXIT DO 314
                 I=I+1
                 IF(XLP(I) < A)CYCLE        ! CYCLE DO 314
                 EXIT                     ! EXIT DO 314
  314           ENDDO
               CTLIM3=0
               DO 316
                 CTLIM3=CTLIM3+1
                 IF(CTLIM3 > N)THEN
                  WRITE(*,*)' *** EXIT DO 316 CTLIM3>N ***'
                 ENDIF
                 IF(CTLIM3 > N)EXIT        ! EXIT DO 316
                 J=J-1
                 IF(XLP(J) > A)CYCLE       ! CYCLE DO 316
                 EXIT                    ! EXIT DO 316
  316           ENDDO
               IF(J >= I)THEN
                 TEMP=XLP(I)
                 XLP(I)=XLP(J)
                 XLP(J)=TEMP
                 CYCLE                   ! CYCLE DO 317 
               ENDIF
               EXIT                      ! EXIT DO 317
  317         ENDDO
             XLP(L)=XLP(J)
             XLP(J)=A
             IF(J >= MEDPOS)IR=J-1
             IF(J <= MEDPOS)L=I
            ENDIF
            CYCLE                        ! CYCLE DO 318
  318      ENDDO
!   End of Median sorting.
          XMEDIAN(MEDPOS2)=XLP(MEDPOS)
          NSTAT=NSTAT+1
          NSTOP=NSTOP+1
          MEDPOS2=MEDPOS2+1
          IF(NSTOP < MN)CYCLE              ! CYCLE DO 319
          EXIT                           ! EXIT DO 319
  319    ENDDO
!   For AR(1) model put MEDIAN values in array Y 
!   and associated FREQ values into array XFR.
        NMED=MN-1-(2*(MEDPOS-1))-1
        K=MEDPOS2I
        DO 321 I=1,NMED
          XFR(I)=FREQ(K)
          XLFR(I)=ALOG10(FREQ(K))
          Y(I)=XMEDIAN(K)
!          IF(I > NLIM1.AND.I <= NLIM2)Y(I)=XMEDIAN(NLIM1)
          K=K+1
  321    ENDDO
!  AR(1) fitting: use logs of median smoothed data (Y(I)) and associated
!  frequencies (XFR(I)).
        DO 322 I=1,NMED
          XOUT(I)=10.0**Y(I)
          XLOUT(I)=Y(I)
  322    ENDDO
        OPEN(UNIT=30,FILE='logmedsmth')
  323    FORMAT(1X,F12.7,1X,F12.8)
        DO 1323 I=1,NMED
          WRITE(30,323)XFR(I),XLOUT(I)
 1323    ENDDO
        CLOSE(UNIT=30)      
      ENDIF

!   For Quadratic fit (NFIT=2) put XLOGA data in array Y.
  324  IF(NFIT == 2)THEN
        WRITE(*,*)
        WRITE(*,*)
        WRITE(*,*)' *****************************************'
        WRITE(*,*)' *** Spectral background location via: ***'
        WRITE(*,*)' *** Quadratic fitting               ***'
        WRITE(*,*)' *****************************************'
        J=1
        DO 326 I=2,MN 
          XFR(J)=FREQ(I)
          Y(J)=XLOGA(I)
          J=J+1
  326    ENDDO
!   Start polynomial (quadratic) fit of smoothed log power 
!   v frequency. Algorithm of Davis (1973) p213.
        IORD1=2+1
        DO 328 I=1,IORD1
          C(I)=0.0
          DO 327 J=1,IORD1
            B(I,J)=0.0
  327      ENDDO
  328    ENDDO
        DO 332 I=1,MN
          XP(1)=1.0
          DO 329 J=2,IORD1
            XP(J)=XP(J-1)*XFR(I+1)
  329      ENDDO
          DO 331 J=1,IORD1
            DO 330 K=1,IORD1
             B(J,K)=B(J,K)+XP(J)*XP(K)
  330       ENDDO
            C(J)=C(J)+XP(J)*Y(I)
  331      ENDDO
  332    ENDDO
!   Solve simultaneous equations.
        CALL SOLVE(B,C,IORD1,20,1.0E-06)
!   Calculate quadratic fit power (XOUT) as a function of frequency.
        XOUT(1)=0.0
        DO 334 I=1,MN-1
          XXP=1.0
          YYP=0.0
          DO 333 J=1,IORD1
            YYP=YYP+XXP*C(J)
            XXP=XXP*FREQ(I+1)
  333      ENDDO
          XLOUT(I+1)=YYP
          XOUT(I+1)=10.0**YYP
  334    ENDDO
      ENDIF

!   Test for best AR(1) spectrum fit to the log of the median
!   smoothed spectrum by varying Rho and SO (via least squares).
!   Regular-spaced data: use lowest of either Rho (initial estimate)
!   or Rho (Mudelsee method) (=ROH) as upper bound for rho during fitting
!   (i.e. vary ROH and SOFA).
!   Irregularly-spaced data: use Rho (Mudelsee method) (i.e. vary
!   SOFA only).
      IF(NFIT == 0.OR.NFIT == 1)THEN
!   If on Monte Carlo spectral run miss out screen notes.
        IF(NSPECB == 0)THEN
          WRITE(*,*)' ******************************************'
          WRITE(*,*)' *** Spectral background location via:  ***'
          IF(MUDEST == 0)THEN
            IF(NFIT == 0)THEN
             WRITE(*,*)' *** White Noise model conf. levels.    ***'
            ENDIF
            IF(NFIT == 1)THEN
             WRITE(*,*)' *** a robust-fitted AR(1) model       ***'
             WRITE(*,*)' *** for finding confidence levels      ***'
             WRITE(*,*)' *** (Mann & Lees, 1996, Clim. Change). ***'
            ENDIF

          ELSE
            WRITE(*,*)' *** fitting of AR(1) model based on    ***'
            WRITE(*,*)' *** best estimate of the lag-1        ***'
            WRITE(*,*)' *** correlation coefficient & finding  ***'
            WRITE(*,*)' *** best-fit mean log spectrum level   ***'
            WRITE(*,*)' *** in to find confidence levels.      ***'
          ENDIF
          WRITE(*,*)' ******************************************'
        ENDIF
        IF((XLOUT(NMED)+XLOUT(2)) == 0.0D0)THEN
          WRITE(*,*)
          WRITE(*,*)' *** Problem initialising mid spectrum ***'
          error_flag = -5
          RETURN
        ENDIF
!   Find maximum and minimum log median-smoothed power value (in XLOUT)
        XMAXLP=XLOUT(2)
        XMINLP=XLOUT(NMED) 
        DO 335 I=2,NMED
          IF(XLOUT(I) > XMAXLP)XMAXLP=XLOUT(I)
          IF(XLOUT(I) < XMINLP)XMINLP=XLOUT(I)
  335    ENDDO
        SOLMID=(XMAXLP+XMINLP)/2.0D0
        RANGE=(XMAXLP-XMINLP)/2.0D0
        IF(NFIT == 0)RANGE=XMAXLP-SOLMID
        SOLMAX=SOLMID+RANGE
        SOLMIN=SOLMID-RANGE
        SOF=0.0D0
        SOLINT=0.01D0
        IF(ROH <= 0.0D0)NFIT=0
!   Using AR1 fit and  uniformly/regularly-spaced data vary ROH and SOFA.
        IF(NFIT == 1.AND.NIRR == 0)THEN
          SOFIT=SOLMAX+SOLINT
          ROHB=ROH
          XLSIN=1000000000000.0
          XLSMIN=XLSIN
  338      FORMAT('          Testing LogSO = ',F7.4,' to ',F7.4)
  339      FORMAT('                  ROH =  0.000  to  ',F5.3)
          WRITE(*,*)
          WRITE(*,338)SOLMAX,SOLMIN
          WRITE(*,339)ROH
          WRITE(*,*)
!   Try repeatedly reducing SOFIT to find a better fit:
          CTLIM=0
          DO 340
            CTLIM=CTLIM+1
            IF(CTLIM > N)THEN
             WRITE(*,*)' *** EXIT DO 340 CTLIM>N ***'
            ENDIF
            IF(CTLIM > N)EXIT               ! EXIT DO 340
            SOFIT=SOFIT-SOLINT    
            IF(SOFIT < SOLMIN)EXIT           ! EXIT DO 340
            SOFA=10.0**SOFIT
!   Call minimization procedure (using functions BRENT2 and XLS2) to find
!   best-fit value for rho given current SOFA.
            GUESS=ROHB-0.1D0
            IF(GUESS < 0.0D0)GUESS=ROHB/2.0D0
            XLNYQ=LOG10(XNYQ)
            DUM4=BRENT2(0.0D0,GUESS,ROHB,XLS2,TOL,ROHT,XLOUT,XLFR, &
           NMED,SOFA,XLNYQ)

!   Test for minimum least squares given current value of SOFA
            IF(ROHT /= 2.0D0)THEN
             IF(ROHT < 0.0D0)THEN
               WRITE(*,*)' ROHT <= 0.0'
               EXIT                       ! EXIT DO 340
             ENDIF
             IF(ROHT > ROHB)THEN
               WRITE(*,*)' ROHT > ROHB'
               EXIT                       ! EXIT DO 340
             ENDIF
             IF(DUM4 < XLSMIN)THEN
               XLSMIN=DUM4
               ROHF=ROHT
               SOF=SOFA
               CYCLE                      ! CYCLE DO 340
             ENDIF
             IF(DUM4 /= XLSIN)CYCLE         ! CYCLE DO 340
            ENDIF
            EXIT                          ! EXIT DO 340
  340      ENDDO
          IF(XLSMIN == XLSIN)THEN
            WRITE(*,*)' *** Warning: Least Sqs min = Init. value! ***'
            error_flag = -6
            RETURN
          ENDIF
          IF(ROHT < 0.0D0.OR.ROHT > 1.0D0)THEN
            WRITE(*,*)' *** Problem with Rho fitting ***'
            error_flag = -7
            RETURN
          ENDIF
        ENDIF
!   Find best fit value of SOF given ROH (for regularly-spaced data
!   and white noise fit or irregularly spaced data).
        IF(NFIT == 0.OR.(NIRR == 1.AND.NFIT == 1))THEN
          SOLMI=10.0**SOLMIN
          SOLMA=10.0**SOLMAX
          SOLMD=10.0**SOLMID
  350      FORMAT('          Testing LogSO= ',F7.4,' to ',F7.4)
  351      FORMAT('             using ROH=  ',F5.3)
          IF(ROH <= 0.0D0)ROH=0.0D0
          ROHB=ROH
          IF(NFIT == 0)ROHB=0.0D0
          WRITE(*,350)SOLMAX,SOLMIN
          WRITE(*,351)ROHB
          WRITE(*,*)
!   Call minimization procedure (using functions BRENT3 and XLS3) to find
!   best-fit value for SOFA given ROH.
          XLNYQ=LOG10(XNYQ)
          DUM5=BRENT3(SOLMI,SOLMD,SOLMA,XLS3,TOL,SOFT,XLOUT, &
         XLFR,NMED,ROHB,XLNYQ)
!   Output is least squares best fit of SOF.
          ROHF=ROHB
          IF(SOFT == 20.0D0)THEN
            SOF=SOLMID
          ELSE
            SOF=SOFT
          ENDIF
        ENDIF
!   Check BRENT result for false XLS minimum.
        NPOS=0
        NNEG=0
        XLSUM=0.0d0
        DO 362 I=1,NMED
          XARP2=(SOF*(1.0-(ROHF*ROHF)))/(1.0-2*ROHF*COS((XFR(I)/XNYQ)* &
           XPI)+(ROHF*ROHF))
          XLP2=LOG10(XARP2)
          XLSUM=XLSUM+((XLP2-XLOUT(I))*(XLP2-XLOUT(I)))
          IF(XLP2-XLOUT(I) > 0.0)NPOS=NPOS+1
          IF(XLP2-XLOUT(I) < 0.0)NNEG=NNEG+1
  362    ENDDO
        BRENTSUM=XLSUM
!   If BRENT fitted background spectrum (almost) entirely above median smoothed
!   spectrum then drop fitted background to find true XLS minimum.
        XLMIN=BRENTSUM 
        IF(SOF <= 0.0)THEN
          WRITE(*,*)' *** SOF negative or zero! ***'
          error_flag = -8
          RETURN
        ENDIF      
        XLSOF=LOG10(SOF)
        SOF2=0.0
        PROP=NPOS/NMED
        IF(PROP > 0.7)THEN
          DO 364
            XLSOF=XLSOF-0.1
            IF(XLSOF < SOLMIN)EXIT         ! EXIT DO 364
            SOF2=10.0**XLSOF
            XLSUM=0.0d0
            DO 363 I=1,NMED
             XARP2=(SOF2*(1.0-(ROHF*ROHF)))/(1.0-2*ROHF*COS((XFR(I)/ &
              XNYQ)*XPI)+(ROHF*ROHF))
             XLP2=LOG10(XARP2)
             XLSUM=XLSUM+((XLP2-XLOUT(I))*(XLP2-XLOUT(I)))            
  363       ENDDO
            IF(XLSUM < XLMIN)THEN
             XLMIN=XLSUM
             SOFNEW=SOF2
            ENDIF
            CYCLE                        ! CYCLE DO 364
  364      ENDDO  
        ENDIF
!   If BRENT fitted background spectrum (almost) entirely below median smoothed
!   spectrum then raise fitted background to find true XLS minimum.
        PROP=NNEG/NMED
        IF(PROP > 0.7)THEN
          DO 366
            XLSOF=XLSOF+0.1
            IF(XLSOF > SOLMAX)EXIT         ! EXIT DO 366
            SOF2=10.0**XLSOF
            XLSUM=0.0d0
            DO 365 I=1,NMED
             XARP2=(SOF2*(1.0-(ROHF*ROHF)))/(1.0-2*ROHF*COS((XFR(I)/ &
              XNYQ)*XPI)+(ROHF*ROHF))
             XLP2=LOG10(XARP2)
             XLSUM=XLSUM+((XLP2-XLOUT(I))*(XLP2-XLOUT(I)))            
  365       ENDDO
            IF(XLSUM < XLMIN)THEN
             XLMIN=XLSUM
             SOFNEW=SOF2
            ENDIF
            CYCLE                        ! CYCLE DO 366
  366      ENDDO  
        ENDIF
        IF(XLMIN < BRENTSUM)THEN
          SOF=SOFNEW
          WRITE(*,*)' Spectral background from fixing BRENT fit'
        ENDIF
!   If final fitted ROH (= ROHF) is not raw or Mudelesee ROH then 
!   check simple alternative to fit from BRENT output by using fixed
!   initial ROH (=ROHI).
        XLREV=XLMIN
        IF(ROHF /= ROHI)THEN
          XLSOF=SOLMAX+0.1
          DO 368 
            XLSOF=XLSOF-0.1
            IF(XLSOF < SOLMIN)EXIT         ! EXIT DO 368
            SOF2=10.0**XLSOF
            XLSUM=0.0d0
            DO 367 I=1,NMED
             XARP2=(SOF2*(1.0-(ROHI*ROHI)))/(1.0-2*ROHI*COS((XFR(I)/ &
              XNYQ)*XPI)+(ROHI*ROHI))
             XLP2=LOG10(XARP2)
             XLSUM=XLSUM+((XLP2-XLOUT(I))*(XLP2-XLOUT(I)))            
  367       ENDDO
            IF(XLSUM < XLMIN)THEN
             XLMIN=XLSUM
             SOFNEW=SOF2
            ENDIF
            CYCLE                       ! CYCLE DO 368            
  368      ENDDO
        ENDIF
        IF(XLMIN < XLREV)THEN
          SOF=SOFNEW
          ROHF=ROHI
          WRITE(*,*)'Spectral background from ROH raw gives best fit'
        ENDIF
        IF(SOF > 0.0D0)XLSOF=LOG10(SOF)
  370    FORMAT('  Fitted spectral background used: Log SOF = ',F7.4)
  371    FORMAT('                                  ROH =  ',F5.3)
        WRITE(*,*)
        WRITE(*,370)XLSOF
        WRITE(*,371)ROHF
!   Having found best fit now calculate best fit AR(1) spectrum [for
!   irregularly spaced data this means out to Nyquist based on minimum
!   sample interval (instead of mean SI Nyquist used in fitting)]. 
        ROHF2=ROHF*ROHF
        DO 375 I=2,MN
          ARP2(I)=SOF*(1.0-ROHF2)/(1.0-2*ROHF*COS((FREQ(I)/XNYQ)*XPI) &
           +ROHF2)
          XLARP2(I)=ALOG10(ARP2(I))
  375    ENDDO
      ENDIF

!   Using a power law or quadratic fit.   

  380  IF(NFIT == 2.OR.NFIT == 3)THEN
        DO 382 I=2,MN
          XLARP2(I)=XLOUT(I)
          ARP2(I)=XOUT(I)
  382    ENDDO

      ENDIF

!   Start bias correction for spectral estimates based on Lomb-Scargle
!   algorithm (larger correction for higher frequencies) following
!   Schulz & Mudelsee (2002) procedure.

!   Reference AR1 spectrum is the confidence level background found
!   earlier i.e. ARP2(I) AND XLARP2(I).

!   Calculate 200 AR(1) time series using estimated value of lag-1
!   autocorrelation coefficient (ROH).

!   Calculate NSTART Gaussian random numbers.

      IF(NIRR == 1.AND.NSPECB == 0.AND.(NFIT == 0.OR.NFIT == 1))THEN
        NP1=NSTART+1      
        NARSIMC=0
        IDUM=-14
      ENDIF
  400  CONTINUE
      IF(NIRR == 1.AND.NSPECB == 0.AND.(NFIT == 0.OR.NFIT == 1))THEN
        DO 402 I=1,NP1
          ROUT(I)=GASDEV(IDUM)
  402    ENDDO
        SUM=0.0
        DO 404 I=1,NP1-1
          SUM=ROUT(I)+ROHF*SUM
          X(I)=SUM
  404    ENDDO
!   Calculate spectrum for each of the 200 simulated AR1 time series.
        NARSIM=1
        NARSIMC=NARSIMC+1
        IF(NARSIMC == 1)THEN
          WRITE(*,*)
          WRITE(*,*)' *****************************************'
          WRITE(*,*)' *** Generation of 200 AR1 time series ***'
          WRITE(*,*)' *** & corresponding power spectra for ***'
          WRITE(*,*)' *** bias correction of Lomb-Scargle   ***'
          WRITE(*,*)' *** spectral estimates [Schulz and    ***'
          WRITE(*,*)' *** Mudelsee (2002) Computers &      ***'
          WRITE(*,*)' *** Geosciences pp421-426.]          ***'
          WRITE(*,*)' *****************************************'
          WRITE(*,*)
        ENDIF
      ENDIF
!   Go back to generate periodogram
      IF(NARSIMC > 0)GO TO 200
  410  CONTINUE
      IF(NIRR == 1.AND.NSPECB == 0.AND.(NFIT == 0.OR.NFIT == 1))THEN
!   Add latest simulated AR1 spectrum to sum.
        DO 411 I=2,MN
          ARSPAVE(I)=ARSPAVE(I)+P(I)
  411    ENDDO
        NSIM=200
  412    FORMAT(' No. of series & spectra generated =    5/',I3)
  413    FORMAT(' No. of series & spectra generated =   50/',I3)
  414    FORMAT(' No. of series & spectra generated =  100/',I3)
  415    FORMAT(' No. of series & spectra generated =  150/',I3)
  416    FORMAT(' No. of series & spectra generated =  200/',I3) 
        IF(NARSIMC == 5)THEN
          WRITE(*,412)NSIM
        ENDIF
        IF(NARSIMC == 50)THEN
          WRITE(*,413)NSIM
        ENDIF
        IF(NARSIMC == 100)THEN
          WRITE(*,414)NSIM
        ENDIF
        IF(NARSIMC == 150)THEN

          WRITE(*,415)NSIM
        ENDIF
        IF(NARSIMC == 200)THEN
          WRITE(*,416)NSIM
        ENDIF
!   Loop back to get next simulated AR1 time series.
        IF(NARSIMC <= NSIM-1)GO TO 400
!   Find average spectrum for NSIM simulations.
        XNSIM=NSIM
        DO 417 I=2,MN
          ARSPAVE(I)=ARSPAVE(I)/XNSIM
  417    ENDDO
!   Smooth average spectrum using P as a holding array.  
        DO 420 K=1,3
          P(1)=ARSPAVE(1)
          P(2)=(0.5*ARSPAVE(2))+(0.5*ARSPAVE(3))
          P(MN)=(0.5*ARSPAVE(MN-1))+(0.5*ARSPAVE(MN))
          DO 418 I=3,MN-1
            P(I)=(0.25*ARSPAVE(I-1))+(0.5*ARSPAVE(I))+ &
            (0.25*ARSPAVE(I+1))
  418      ENDDO
          DO 419 I=1,MN
            ARSPAVE(I)=P(I)
  419      ENDDO
  420    ENDDO
!   Normalise average spectrum.
        TOTP=0.0
        DO 422 I=2,MN
           TOTP=TOTP+ARSPAVE(I)
  422    ENDDO
!
        NBIASC=0
        DO 424 I=2,MN
          IF(ARSPAVE(I) <= 0.0)NBIASC=1
  424    ENDDO
        IF(NBIASC /= 0)THEN
          WRITE(*,*)
          WRITE(*,*)' *** Bias correction step failed ! ***'
          WRITE(*,*)' *** (negative power in reference  ***'
          WRITE(*,*)' *** spectrum).                  ***'
          DO 426 I=2,MN
            ARSPAVE(I)=0.0
            XLARSP(I)=0.0
  426      ENDDO
        ELSE
          DO 428 I=2,MN
            TEMP= ARSPAVE(I)
            ARSPAVE(I)=ARSPAVE(I)/TOTP
            XLARSP(I)=ALOG10(ARSPAVE(I))
            WRITE(6,*)I,TEMP,ARSPAVE(I),XLARSP(I)
  428      ENDDO
        ENDIF
!   Store original biased estimate of spectrum.
        DO 430 I=2,MN
          XLSPBIAS(I)=XLOGA(I)
  430    ENDDO      
!   Calculate and apply bias correction function to original spectrum.
        IF(NBIASC == 0)THEN
          DO 432 I=2,MN
            BIASF(I)=XLARSP(I)-XLARP2(I)
            XLOGA(I)=(ALOG10(PSMTH(I)))-BIASF(I)
            AC(I)=AC(I)/SQRT(10.0**BIASF(I))
            BC(I)=BC(I)/SQRT(10.0**BIASF(I))
            PSMTH(I)=10.0**XLOGA(I)
  432      ENDDO
!   Re-normalise AR-average spectrum.
          TOTP=0.0
          DO 434 I=2,MN
            TOTP=TOTP+ARSPAVE(I)
  434      ENDDO
!
          DO 436 I=2,MN
            ARSPAVE(I)=ARSPAVE(I)/TOTP
            XLARSP(I)=ALOG10(ARSPAVE(I))
  436      ENDDO
!   Re-normalise bias-corrected spectrum.
          TOTP=0.0
          DO 438 I=2,MN
            TOTP=TOTP+(PSMTH(I))
  438      ENDDO
          DO 440 I=2,MN
            PSMTH(I)=(PSMTH(I))/TOTP
            XLOGA(I)=ALOG10(PSMTH(I))
  440      ENDDO
        ENDIF
!   Re-check position of background is correctly fitted.
        NSPECB=1
        WRITE(*,*)
        WRITE(*,*)'   Re-checking background fit'
        WRITE(*,*)
        GO TO 304
      ENDIF

!   If AR1 fit is inappropriate then switch to Quadratic or power-law fit.

      IF(NFIT == 1)THEN
        NMIS=0
        AMN=MN
!   Check high frequency fit.
        TMP=AMN/1.25
        NST=TMP
        DIFFMAX=0.0
        DO 452 I=NST,MN-1
          DIFF=XLARP2(I)-XLOGA(I)
          IF(DIFF > DIFFMAX)DIFFMAX=DIFF
          IF(DIFF > 1.5.OR.DIFF < -1.5)NMIS=NMIS+1
  452    ENDDO
!   Check low frequency fit.
        DO 454 I=3,MN/5
          DIFF=XLARP2(I)-XLOGA(I)
          IF(DIFF > DIFFMAX)DIFFMAX=DIFF
          IF(DIFF > 1.0.OR.DIFF < -1.0)NMIS=NMIS+1
  454    ENDDO
        NTOT=(MN-(AMN/1.25))+(MN/5)-3
        ANTOT=NTOT
        ANTOT2=0.4*ANTOT
        NTOT2=ANTOT2
        IF(NMIS >= NTOT2 .AND. NMIS>0)THEN
          NSPECB=2
          CTLIM=0
          DO 455
            CTLIM=CTLIM+1
            IF(CTLIM > 5)EXIT              !  EXIT DO 455
            WRITE(*,*)
            WRITE(*,*)' NMIS = ',NMIS
            WRITE(*,*)' DIFFMAX = ',DIFFMAX
            WRITE(*,*)
            WRITE(*,*)'   *** AR1 fit inappropriate ***'
            WRITE(*,*)' Choose switch to Quadratic fit (Q)'
            WRITE(*,*)'             or Power law fit (P):'
            WRITE(*,*)  
!           READ(2,'(41A,A)',ERR=455)HEADER,CFIT
            WRITE(*,'(A)') CFIT 
            IF(CFIT == 'Q'.OR.CFIT == 'q')THEN
             NFIT=2
             GO TO 324
            ENDIF
            IF(CFIT == 'P'.OR.CFIT == 'p')THEN
             NFIT=3
             GO TO 280
            ENDIF
            WRITE(*,*)
            WRITE(*,*)'   *** Enter Q or P ***'
            CYCLE                        ! CYCLE DO 455
  455      ENDDO
        ENDIF
      ENDIF     

!   Write biased data spectrum, average spectrum and reference
!   spectrum to file.

      OPEN(UNIT=14,FILE='specx')
  460  FORMAT('   Variance = ',F15.8)
      WRITE(14,460)VARX1
  461  FORMAT('  Frequency  Corrected spec.         AC            BC' &
              'Ref. spec.  Time series (past)  AveAR1 spec.    Bias spec.')
  462  FORMAT(2X,F10.6,1X,F15.8,1X,F15.8,1X,F15.8,1X,F15.8,1X, &
      F15.8,1X,F15.8,1X,F15.8)
      IF(NSPECB == 1.OR.NSPECB == 2)THEN
        WRITE(14,461)
        DO 463 I=2,MN
          WRITE(14,462)FREQ(I),PSMTH(I),AC(I),BC(I),ARP2(I), &
         ROUT(I),ARSPAVE(I),XLSPBIAS(I)
  463    ENDDO
      ENDIF
  464  FORMAT('   Frequency     Spectrum            AC            BC' &
              'Background  Time series (part)')
  465  FORMAT(2X,F10.6,1X,F15.8,1X,F15.8,1X,F15.8,1X,F15.8,1X,F15.8)
      IF(NSPECB == 0)THEN
        WRITE(14,464)
        DO 466 I=2,MN
          WRITE(14,465)FREQ(I),PSMTH(I),AC(I),BC(I),ARP2(I),ROUT(I)
  466    ENDDO
      ENDIF
      CLOSE(UNIT=14)
!
!   CPOWER to select:
!    = R: relative power
!    = A: absolute amplitude
!
!   Calculate approx. upper confidence levels based on chi-squared
!   distribution, given effective degrees of freedom (Priestley, 1981)

      DO 467 I=2,MN
        YINTA(I)=ARP2(I)*2.292
        YINTB(I)=ARP2(I)*2.927
        YINTD(I)=ARP2(I)*4.860
  467  ENDDO
!
!   Code added to convert relative power to absolute amplitude
!   Garry Hayman
!   Centre for Ecology and HYdrology
!   March 2011
!
      FACTOR  = 0.5*VARX1/0.147
      WRITE(*,*) VARX1,FACTOR
!
      IF(CPOWER == 'A'.OR.CPOWER == 'a') THEN
         PSMTH = SQRT(PSMTH*FACTOR)
         ARP2  = SQRT(ARP2*FACTOR)
         YINTA = SQRT(YINTA*FACTOR)
         YINTB = SQRT(YINTB*FACTOR)
         YINTD = SQRT(YINTD*FACTOR)
      ENDIF
!
!   Code to set-up OUTPUT array for python
!
      DO 468 I=2,MN
        OUTPUT(I-1,1) = FREQ(I)
        OUTPUT(I-1,2) = 1.0/FREQ(I)
        OUTPUT(I-1,3) = PSMTH(I)
        OUTPUT(I-1,4) = ARP2(I)
        OUTPUT(I-1,5) = YINTA(I)
        OUTPUT(I-1,6) = YINTB(I)
        OUTPUT(I-1,7) = YINTD(I)
  468  ENDDO
!
!   Analysis for significant spectral peaks and writing spectrum to file.

      R90=YINTA(2)/ARP2(2)
      R95=YINTB(2)/ARP2(2)
      R99=YINTD(2)/ARP2(2)

  471  FORMAT('Intercept = ',F10.5,'  Slope = ',F10.5)
  472  FORMAT('Spectral background located using power law fit')
  475  FORMAT('Spectral background located using White Noise model')
  476  FORMAT('Spectral background located using robust AR(1) model fit,' &
              'Rho-fit = ',F5.3,'   Log10(SOF) =',F6.3)
  478  FORMAT('Spectral background located using a quadratic fit')
  479  FORMAT('    Frequency   Wavelength  Maximum   Minimum    Rel./Abs. ' &
              'Confidence')
  480  FORMAT('    Frequency     Period    Maximum   Minimum    Rel./Abs. ' &
              'Confidence')
  481  FORMAT('  Cycles per yr    Year      Period   Wlength      Power ' &
              'Level')
  482  FORMAT('  Cycles per day    day      Period   Wlength      Power ' &
              'Level')
  483  FORMAT('  Cycles per m      m      Wlength   Wlength      Power ' &
              'Level')
  484  FORMAT('  Cycles per kyr    kyr      Period    Period      Power ' &
              'Level')
  485  FORMAT(1X,F10.5,1X,F11.4,1X,F10.6,1X,F11.6,1X,F11.6,1X,F11.6,1X, &
      F11.6)
  486  FORMAT(1X,F10.5,1X,F11.4,1X,F10.6,1X,F11.6,1X,F11.6,1X,F11.6,1X, &
      F11.6,1X,'>90%')
  487  FORMAT(1X,F10.5,1X,F11.4,1X,F10.6,1X,F11.6,1X,F11.6,1X,F11.6,1X, &
      F11.6,1X,'>95%')
  488  FORMAT(1X,F10.5,1X,F11.4,1X,F10.6,1X,F11.6,1X,F11.6,1X,F11.6,1X, &
      F11.6,1X,'>99%')
  489  FORMAT('Spectral method = Lomb-Scargle Method, Variance = ',F18.7)
  490  FORMAT('Data filename = ',A50)
 1490  FORMAT('Number of data points =',I7 &
            /'BW =',F13.6 &
            /'DoF =',I2 &
            /'Rho-raw = ',F5.3)
 1491  FORMAT('Number of data points =',I7 &
            /'BW =',F13.6 &
            /'DoF =',I2 &
            /'Rho-raw = ',F5.3 &
            /'Rho-Mud = ',F5.3)
  491  FORMAT('90%CL=',F6.3,'xBakgd 95%CL=',F6.3,'xBakgd 99%CL=',F6.3, &
              'xBakgd')
  492  FORMAT('90%CL=',F6.3,'xBakgd 95%CL=',F6.3,'xBakgd 99%CL=',F6.3, &
              'xBakgd  Tau =',F10.4,' Year')
  493  FORMAT('90%CL=',F6.3,'xBakgd 95%CL=',F6.3,'xBakgd 99%CL=',F6.3, &
              'xBakgd  Tau =',F10.4,' Day')
  494  FORMAT('90%CL=',F6.3,'xBakgd 95%CL=',F6.3,'xBakgd 99%CL=',F6.3, &
              'xBakgd  Tau =',F10.4,' kyr')
  495  FORMAT('90%CL=',F6.3,'xBakgd 95%CL=',F6.3,'xBakgd 99%CL=',F6.3, &
              'xBakgd  Tau =',F10.4,' m')
  496  FORMAT('Frequency    Period   Smoothed   Backgnd  90%Conf.  95%' &
              'Conf.  99%Conf. Peak')
  497  FORMAT('Frequency  Wavelnth   Smoothed   Backgnd  90%Conf.  95%' &
              'Conf.  99%Conf. Peak')
  498  FORMAT('Cycles per yr  Year      Power     Power     Level     ' &
              'Level     Level  conf.')
  499  FORMAT('Cycles per day  Day      Power     Power     Level     ' &
              'Level     Level  conf.')
  500  FORMAT('Cycles per m     m      Power     Power     Level      ' &
              'Level     Level  conf.')
  501  FORMAT('Cycles per kyr  kyr      Power     Power     Level     ' &
              'Level     Level  conf.')
  502  FORMAT(1X,F12.7,1X,F11.5,1X,F12.8,1X,F12.8,1X,F12.8,1X,F12.8,1X, &
      F12.8)
  503  FORMAT(1X,F12.7,1X,F11.5,1X,F12.8,1X,F12.8,1X,F12.8,1X,F12.8,1X, &
      F12.8,1X,'>90%')
  504  FORMAT(1X,F12.7,1X,F11.5,1X,F12.8,1X,F12.8,1X,F12.8,1X,F12.8,1X, &
      F12.8,1X,'>95%')
  505  FORMAT(1X,F12.7,1X,F11.5,1X,F12.8,1X,F12.8,1X,F12.8,1X,F12.8,1X, &
      F12.8,1X,'>99%')
  506  FORMAT('    Frequency      Period    Smoothed      Backgnd      ', &
      '90%Conf.     95%Conf.     99%Conf. Peak')
  507  FORMAT('    Frequency  Wavelength    Smoothed      Backgnd      ', &
      '90%Conf.     95%Conf.     99%Conf. Peak')
  508  FORMAT('Cycles per yr     Year       Power       Power        ', &
      'Level       Level       Level   conf.')
  509  FORMAT('Cycles per day     Day       Power       Power        ', &
      'Level       Level       Level   conf.')
  510  FORMAT('Cycles per m       m        Power       Power         ', &
      'Level       Level       Level   conf.')
  511  FORMAT('Cycles per kyr     kyr       Power       Power        ', &
      'Level       Level       Level   conf.')

!   IA > 0 Means significant peaks are present.

      IA=0
      ID=0
      DO 515 I=2,MN
        IF(PSMTH(I) >= YINTA(I))IA=IA+1
        IF(PSMTH(I) >= YINTD(I))ID=ID+1
  515  ENDDO
      IF(IA > 0)THEN
        WRITE(*,*)
        WRITE(*,*)'       *** Significant Spectral Peaks ***'
        WRITE(*,*)
        IF(KMK == 1.OR.KMK == 2)THEN
          WRITE(*,479)
        ELSE
          WRITE(*,480)
        ENDIF
        IF(KMK == 1.OR.KMK == 2)WRITE(*,483)  
        IF(KMK == 3.OR.KMK == 4)WRITE(*,484)
        IF(KMK == 5.OR.KMK == 7.OR.KMK == 8)WRITE(*,481)
        IF(KMK == 6)WRITE(*,482)
      ENDIF
      IF(IA == 0)THEN
        WRITE(*,*)' *** No spectral peaks exceed confidence levels ***'
      ENDIF
      OPEN(UNIT=3,FILE='spectrum')
      WRITE(3,490)TSFILE
      IF(NIRR == 0)THEN
         WRITE(3,1490)NROWS,BW,NDOF,ROHI
      ELSE
         WRITE(3,1491)NROWS,BW,NDOF,ROHI,ROHN
      ENDIF
      WRITE(3,489)VARX1
      IF(NFIT == 0)WRITE(3,475)
      IF(NFIT == 3)THEN
        WRITE(3,472)
        WRITE(3,471)XINTERC,SLOPE
      ENDIF
      IF(NFIT == 2)WRITE(3,478)
      IF(NFIT == 1)WRITE(3,476)ROHF,XLSOF
      IF(NFIT /= 1)THEN
        WRITE(3,491)R90,R95,R99
      ENDIF

      IF(NFIT == 1)THEN
        IF(KMK == 1.OR.KMK == 2)WRITE(3,495)R90,R95,R99,TAU
        IF(KMK == 3.OR.KMK == 4)WRITE(3,494)R90,R95,R99,TAU
        IF(KMK == 5.OR.KMK == 7.OR.KMK == 8)THEN
          WRITE(3,492)R90,R95,R99,TAU
        ENDIF
        IF(KMK == 6)WRITE(3,493)R90,R95,R99,TAU
      ENDIF
      IF(KMK == 1.OR.KMK == 2)THEN
        IF(N < 10000)THEN       
          WRITE(3,497)
          WRITE(3,500)
        ELSE
          WRITE(3,507)
          WRITE(3,510)
        ENDIF
      ELSE
        IF(N < 10000)THEN
          WRITE(3,496)
        IF(KMK == 3.OR.KMK == 4)WRITE(3,501)
        IF(KMK == 5.OR.KMK == 7.OR.KMK == 8)WRITE(3,498)
        IF(KMK == 6)WRITE(3,499)
        ELSE
          WRITE(3,506)
        IF(KMK == 3.OR.KMK == 4)WRITE(3,511)
        IF(KMK == 5.OR.KMK == 7.OR.KMK == 8)WRITE(3,508)
        IF(KMK == 6)WRITE(3,509)
        ENDIF
      ENDIF

!  Significant peak search logic:

!   Search for a spectral estimate which exceeds one or more confidence 
!   levels and has a larger clearance of the confidence level than its
!   immediate neighbours.

!   ID > 0 Means one or more spectral estimates exceed the 99% level.
!   IC > 0 Means one or more spectral estimates exceed the 98% level.

      TSLEN=AN*XINT

!   Check if first frequency value (FREQ(2)) is a significant peak.
      XPERIOD=1.0/FREQ(2)
      XMINPER=1.0/(FREQ(2)+(BW/2.0))
      XMAXPER=1.0/(FREQ(2)-(BW/2.0))
      PSTOP1=0
! First frequency >99%?
      IF(PSMTH(2) >= YINTD(2))THEN
        IF(PSMTH(2) > PSMTH(3))THEN
          PSTOP1=1
          IF(XMAXPER < 0.0)THEN
  520       FORMAT(2X,F10.4,5X,F8.3,1X,'>',F8.3,2X,F8.3,4X,F8.5, &
           4X,'>99%')
            WRITE(*,520)FREQ(2),XPERIOD,TSLEN,XMINPER,PSMTH(2)
          ELSE
  521       FORMAT(2X,F10.4,5X,F8.3,2X,F8.3,2X,F8.3,4X,F8.5,4X, &
           '>99%')
            WRITE(*,521)FREQ(2),XPERIOD,XMAXPER,XMINPER,PSMTH(2)
          ENDIF
          IF(N < 10000)THEN
            WRITE(3,488)FREQ(2),XPERIOD,PSMTH(2),ARP2(2),YINTA(2) &
            ,YINTB(2),YINTD(2)
          ELSE
            WRITE(3,505)FREQ(2),XPERIOD,PSMTH(2),ARP2(2),YINTA(2) &
            ,YINTB(2),YINTD(2)
          ENDIF
        ENDIF
      ENDIF
! First frequency >95%?
      IF(PSMTH(2) >= YINTB(2).AND.PSTOP1 == 0)THEN
        IF(PSMTH(2) > PSMTH(3))THEN
          PSTOP1=1
          IF(XMAXPER < 0.0)THEN
  522       FORMAT(2X,F10.4,5X,F8.3,1X,'>',F8.3,2X,F8.3,4X,F8.5, &
           4X,'>95%')
            WRITE(*,522)FREQ(2),XPERIOD,TSLEN,XMINPER,PSMTH(2)
          ELSE
  523       FORMAT(2X,F10.4,5X,F8.3,2X,F8.3,2X,F8.3,4X,F8.5,4X, &
           '>95%')
            WRITE(*,523)FREQ(2),XPERIOD,XMAXPER,XMINPER,PSMTH(2)
          ENDIF
          IF(N < 10000)THEN
            WRITE(3,487)FREQ(2),XPERIOD,PSMTH(2),ARP2(2),YINTA(2) &
            ,YINTB(2),YINTD(2)
          ELSE
            WRITE(3,504)FREQ(2),XPERIOD,PSMTH(2),ARP2(2),YINTA(2) &
            ,YINTB(2),YINTD(2)            
          ENDIF
        ENDIF
      ENDIF
! First frequency >90%?
      IF(PSMTH(2) >= YINTA(2).AND.PSTOP1 == 0)THEN
        IF(PSMTH(2) > PSMTH(3))THEN
          PSTOP1=1
!   Commented out = screen writing peaks >90% CL.
!            IF(XMAXPER < 0.0)THEN
!  524         FORMAT(2X,F10.4,5X,F8.3,1X,'>',F8.3,2X,F8.3,4X,F8.5,
!     *        4X,'>90%')
!             WRITE(*,524)FREQ(2),XPERIOD,TSLEN,XMINPER,PSMTH(2)
!            ELSE
!  525         FORMAT(2X,F10.4,5X,F8.3,2X,F8.3,2X,F8.3,4X,F8.5,4X,
!     *        '>90%')
!             WRITE(*,525)FREQ(2),XPERIOD,XMAXPER,XMINPER,PSMTH(2)
!            ENDIF
          IF(N < 10000)THEN
            WRITE(3,486)FREQ(2),XPERIOD,PSMTH(2),ARP2(2),YINTA(2) &
            ,YINTB(2),YINTD(2)
          ELSE
            WRITE(3,503)FREQ(2),XPERIOD,PSMTH(2),ARP2(2),YINTA(2) &
            ,YINTB(2),YINTD(2)
          ENDIF
        ENDIF
      ENDIF
! First frequency not significant
      IF(PSTOP1 == 0)THEN
        IF(N < 10000)THEN
          WRITE(3,485)FREQ(2),XPERIOD,PSMTH(2),ARP2(2),YINTA(2) &
           ,YINTB(2),YINTD(2)
        ELSE
          WRITE(3,502)FREQ(2),XPERIOD,PSMTH(2),ARP2(2),YINTA(2) &
           ,YINTB(2),YINTD(2)
        ENDIF
      ENDIF

!   Check all other frequency values for significant spectral peaks.

      DO 530 I=3,MN
        PSTOP=0
        XPERIOD=1.0/(FREQ(I))
        XMINPER=1.0/(FREQ(I)+(BW/2.0))
        XMAXPER=1.0/(FREQ(I)-(BW/2.0))
!   99% Confidence level.
        IF(PSMTH(I) >= YINTD(I).AND.PSTOP == 0)THEN    
          PA=PSMTH(I)-YINTD(I)
          PB=PSMTH(I-1)-YINTD(I-1)
          PC=PSMTH(I+1)-YINTD(I+1)
          IF(PA > PB.AND.PA > PC)THEN
            PSTOP=1
            IF(XMAXPER >= 0.0)THEN
             WRITE(*,521)FREQ(I),XPERIOD,XMAXPER,XMINPER,PSMTH(I)
            ELSE
             WRITE(*,520)FREQ(I),XPERIOD,TSLEN,XMINPER,PSMTH(I)
            ENDIF
            IF(N < 10000)THEN
             WRITE(3,488)FREQ(I),XPERIOD,PSMTH(I),ARP2(I),YINTA(I) &
              ,YINTB(I),YINTD(I)
            ELSE
             WRITE(3,505)FREQ(I),XPERIOD,PSMTH(I),ARP2(I),YINTA(I) &
              ,YINTB(I),YINTD(I)
            ENDIF
          ENDIF
          IF(PSTOP == 0)THEN
            IF(PA <= PB.OR.PA <= PC.OR.PSMTH(I) <= &
              PSMTH(I-1).OR.PSMTH(I) <= PSMTH(I+1))THEN
             PSTOP=1
             IF(N < 10000)THEN
               WRITE(3,485)FREQ(I),XPERIOD,PSMTH(I),ARP2(I),YINTA(I) &
                ,YINTB(I),YINTD(I)
             ELSE
               WRITE(3,502)FREQ(I),XPERIOD,PSMTH(I),ARP2(I),YINTA(I) &
                ,YINTB(I),YINTD(I)
             ENDIF
            ENDIF
          ENDIF
        ENDIF
!   95% Confidence level.
        IF(PSMTH(I) >= YINTB(I).AND.PSTOP == 0)THEN      
          PA=PSMTH(I)-YINTB(I)
          PB=PSMTH(I-1)-YINTB(I-1)
          PC=PSMTH(I+1)-YINTB(I+1)
          IF(PA > PB.AND.PA > PC)THEN
            PSTOP=1
            IF(XMAXPER >= 0.0)THEN
             WRITE(*,523)FREQ(I),XPERIOD,XMAXPER,XMINPER,PSMTH(I)
            ELSE
             WRITE(*,522)FREQ(I),XPERIOD,TSLEN,XMINPER,PSMTH(I)
            ENDIF
            IF(N < 10000)THEN
             WRITE(3,487)FREQ(I),XPERIOD,PSMTH(I),ARP2(I),YINTA(I) &
              ,YINTB(I),YINTD(I)
            ELSE
             WRITE(3,504)FREQ(I),XPERIOD,PSMTH(I),ARP2(I),YINTA(I) &
              ,YINTB(I),YINTD(I)            
            ENDIF
          ENDIF
          IF(PSTOP == 0)THEN
            IF(PA <= PB.OR.PA <= PC.OR.PSMTH(I) <= PSMTH(I-1) &
            .OR.PSMTH(I) <= PSMTH(I+1))THEN
             PSTOP=1
             IF(N < 10000)THEN
               WRITE(3,485)FREQ(I),XPERIOD,PSMTH(I),ARP2(I),YINTA(I) &
                ,YINTB(I),YINTD(I)
             ELSE
               WRITE(3,502)FREQ(I),XPERIOD,PSMTH(I),ARP2(I),YINTA(I) &
                ,YINTB(I),YINTD(I)
             ENDIF
            ENDIF
          ENDIF
        ENDIF
!   90% Confidence level.
        IF(PSMTH(I) >= YINTA(I).AND.PSTOP == 0)THEN      
          PA=PSMTH(I)-YINTA(I)
          PB=PSMTH(I-1)-YINTA(I-1)
          PC=PSMTH(I+1)-YINTA(I+1)
          IF(PA > PB.AND.PA > PC)THEN
            PSTOP=1
            IF(XMAXPER >= 0.0)THEN
             IF(N < 10000)THEN
               WRITE(3,486)FREQ(I),XPERIOD,PSMTH(I),ARP2(I),YINTA(I) &
                ,YINTB(I),YINTD(I)
             ELSE
               WRITE(3,503)FREQ(I),XPERIOD,PSMTH(I),ARP2(I),YINTA(I) &
                ,YINTB(I),YINTD(I)
             ENDIF
            ENDIF
          ENDIF
          IF(PSTOP == 0)THEN
            IF(PA <= PB.OR.PA <= PC.OR.PSMTH(I) <= PSMTH(I-1) &
            .OR.PSMTH(I) <= PSMTH(I+1))THEN
             PSTOP=1
             IF(N < 10000)THEN
               WRITE(3,485)FREQ(I),XPERIOD,PSMTH(I),ARP2(I),YINTA(I) &
                ,YINTB(I),YINTD(I)
             ELSE
               WRITE(3,502)FREQ(I),XPERIOD,PSMTH(I),ARP2(I),YINTA(I) &
                ,YINTB(I),YINTD(I)
             ENDIF
            ENDIF
          ENDIF
        ENDIF
!   No confidence levels exceeded.
        IF(PSTOP == 0)THEN
          IF(PSMTH(I) < YINTA(I))THEN
            IF(N < 10000)THEN
             WRITE(3,485)FREQ(I),XPERIOD,PSMTH(I),ARP2(I),YINTA(I) &
              ,YINTB(I),YINTD(I)
            ELSE
             WRITE(3,502)FREQ(I),XPERIOD,PSMTH(I),ARP2(I),YINTA(I) &
              ,YINTB(I),YINTD(I)
            ENDIF   
          ENDIF
        ENDIF
  530  ENDDO
      CLOSE(UNIT=3)
      IF(ZTEST > 0.0)THEN 
        WRITE(*,*)' NOTE: Power values close to or equalling zero'
        WRITE(*,*)'       have been reset to 0.000000001'
      ENDIF

!   Write out Log Power v Freq and Log power v log freq spectrum.

      OPEN(UNIT=12,FILE='logspec')
      WRITE(12,490)TSFILE
      IF(NIRR == 0)THEN
         WRITE(12,1490)NROWS,BW,NDOF,ROHI
      ELSE
         WRITE(12,1491)NROWS,BW,NDOF,ROHI,ROHN
      ENDIF
      WRITE(12,489)VARX1
      IF(NFIT == 0)WRITE(12,475)    
      IF(NFIT == 1)WRITE(12,476)ROHF,XLSOF
      IF(NFIT == 2)WRITE(12,478)
      IF(NFIT == 3)WRITE(12,472)
  540  FORMAT('90%CL=',F6.3,'+Bakgd 95%CL=',F6.3,'+Bakgd 99%CL=',F6.3, &
      '+Bakgd')
      WRITE(12,540)R90,R95,R99
  541  FORMAT('Frequency    Period/Wlen    Log Freq.   Log Per/Wl', &
      '   Log Power    Log Bkgd   Log 90%CL   Log 95%CL', &
      '   Log 99%CL')
  542  FORMAT('Frequency    Period/Wlen    Log Freq.   Log Per/Wl', &
      '   Log Power    Log Bkgd   Log 90%CL   Log 95%CL', &
      '   Log 99%CL')
      IF(N < 10000)THEN
        WRITE(12,541)
      ELSE
        WRITE(12,542)
      ENDIF
  544  FORMAT(1X,F10.6,1X,F12.4,1X,F12.8,1X,F12.8,1X,F11.8,1X,F11.8,1X, &
      F11.8,1X,F11.8,1X,F11.8)
  545  FORMAT(1X,F12.7,1X,F13.5,1X,F12.8,1X,F12.8,1X,F11.8,1X,F11.8,1X, &
      F11.8,1X,F11.8,1X,F11.8)
      DO 550 I=2,MN
        TEMPER=1.0/FREQ(I)
        TLOGP=ALOG10(TEMPER)
        TLFRE=ALOG10(FREQ(I))
        TEMPBK=XLARP2(I)
        SYINTA=YINTA(I)
        SYINTB=YINTB(I)
        SYINTD=YINTD(I)
        TEMP90=ALOG10(SYINTA)
        TEMP95=ALOG10(SYINTB)
        TEMP99=ALOG10(SYINTD)
        IF(N < 10000)THEN
          WRITE(12,544)FREQ(I),TEMPER,TLFRE,TLOGP,XLOGA(I),TEMPBK, &
         TEMP90,TEMP95,TEMP99
        ELSE
          WRITE(12,545)FREQ(I),TEMPER,TLFRE,TLOGP,XLOGA(I),TEMPBK, &
         TEMP90,TEMP95,TEMP99
        ENDIF
  550  ENDDO
      CLOSE(UNIT=12)
      WRITE(*,*)
      WRITE(*,*)'     ***************************************'
      WRITE(*,*)'     *** The spectral data are listed in ***'
      WRITE(*,*)'     *** "spectrum" and in "logspec".    ***'
      WRITE(*,*)'     ***************************************'
  600 WRITE(*,*)

      MN = MN-1

      RETURN
      END SUBROUTINE TIME_SERIES

END MODULE MODULE_TIME_SERIES_PROCEDURES
!======================================================================+
   SUBROUTINE SOLVE(A,B,N,N1,ZERO)
!======================================================================+
!
!   Subroutine from Davis, 1973, p144
!   Used with quadratic fit to log spectrum.
!
!----------------------------------------------------------------------+

      IMPLICIT                                  NONE

      INTEGER,INTENT(IN)                        :: &
            N,N1

      INTEGER                                   :: &
            I,J,K

      REAL,DIMENSION(N1,N1),INTENT(INOUT)       :: &
            A

      REAL,DIMENSION(N1),INTENT(INOUT)          :: &
            B

      REAL,INTENT(IN)                           :: &
            ZERO

      REAL                                      :: &
            DIV,RATIO

      DO 60 I=1,N
        DIV=A(I,I)
        IF(ABS(DIV)-ZERO)70,70,10
   10    DO 20 J=1,N
          A(I,J)=A(I,J)/DIV
   20    ENDDO
        B(I)=B(I)/DIV
        DO 50 J=1,N
          IF(I-J)30,50,30
   30      RATIO=A(J,I)
          DO 40 K=1,N
            A(J,K)=A(J,K)-RATIO*A(I,K)
   40      ENDDO
          B(J)=B(J)-RATIO*B(I)
   50    ENDDO
   60 ENDDO

      RETURN

   70 WRITE(*,*)
      WRITE(*,*)' Matrix problem during polynomial fit to spectrum'

      END SUBROUTINE SOLVE

!======================================================================+
   FUNCTION ROFUNC(B)
!======================================================================+

!   Function from Press et al 1992, p700.
!   NB Uses HPSORT in place of SELECT.

!----------------------------------------------------------------------+

      USE                                       MODULE_TIME_SERIES_DATA

      IMPLICIT                                  NONE

      REAL                                      :: &
            ROFUNC

      REAL,INTENT(IN)                           :: &
            B

      INTEGER                                   :: &
            J,NDATA,N1,NML,NMH

      REAL                                      :: &
            D,ABDEV,SUM

      NDATA = NDATAT
      N1    = NDATA+1
      NML   = N1/2
      NMH   = N1-NML

      DO 10 J=1,NDATA
        ARR(J)=YT(J)-B*XT(J)
   10 ENDDO

      CALL HPSORT(NDATA,ARR)

      AA=0.5*(ARR(NML)+ARR(NMH))
      SUM=0.0
      ABDEV=0.0
      DO 20 J=1,NDATA
        D=YT(J)-(B*XT(J)+AA)
        ABDEV=ABDEV+ABS(D)
        SUM=SUM+XT(J)*SIGN(1.0,D)
   20  ENDDO

      ROFUNC=SUM

      RETURN
      END FUNCTION ROFUNC

!======================================================================+
   SUBROUTINE HPSORT(NDATA,ARR)
!======================================================================+
!
!    Subroutine HPSORT from Press et al 1992, p329
!
!----------------------------------------------------------------------+

      IMPLICIT                                  NONE

      INTEGER,INTENT(IN)                        :: &
            NDATA

      REAL,DIMENSION(300),INTENT(INOUT)         :: &
            ARR

      INTEGER                                   :: &
            I,J,L,IR

      REAL                                      :: &
            RRA

      L=NDATA/2+1
      IR=NDATA

      DO 20
        IF(L > 1)THEN
          L=L-1
          RRA=ARR(L)
        ELSE
          RRA=ARR(IR)
          ARR(IR)=ARR(1)
          IR=IR-1
          IF(IR == 1)THEN
            ARR(1)=RRA
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
        DO 10
          IF(J <= IR)THEN
            IF(J < IR)THEN
             IF(ARR(J) < ARR(J+1))J=J+1
            ENDIF
            IF(RRA < ARR(J))THEN
             ARR(I)=ARR(J)
             I=J
             J=J+J
            ELSE
             J=IR+1
            ENDIF
            CYCLE                        ! CYCLE DO 10
          ENDIF
          EXIT                           ! EXIT DO 10
   10    ENDDO
        ARR(I)=RRA
        CYCLE                            ! CYCLE DO 20
   20  ENDDO
      END SUBROUTINE HPSORT
!
!======================================================================+
   FUNCTION RAN1(IDUM)
!======================================================================+
!
!   Function from Press et al 1992, p271
!
!----------------------------------------------------------------------+

      IMPLICIT                                  NONE

      REAL                                      :: &
            RAN1

      INTEGER, INTENT(INOUT)                    :: &
            IDUM

      INTEGER                                   :: &
            J,K,IY

      INTEGER,PARAMETER                         :: &
            IA   = 16807, IM=2147483647, IQ=127773, IR=2836, &
            NTAB =    32, NDIV=1+(IM-1)/NTAB

      INTEGER,DIMENSION(NTAB)                   :: &
            IV

      REAL, PARAMETER                           :: &
            AM=1./IM, EPS=1.2E-7, RNMX=1.-EPS

      SAVE IV,IY

      DATA IV /NTAB*0/, IY /0/

      IF(IDUM <= 0.OR.IY == 0)THEN
         IDUM=MAX(-IDUM,1)
         DO 10 J=NTAB+8,1,-1
            K=IDUM/IQ
            IDUM=IA*(IDUM-K*IQ)-IR*K
            IF(IDUM < 0)IDUM=IDUM+IM
            IF(J <= NTAB)IV(J)=IDUM
   10    ENDDO
         IY=IV(1)
      ENDIF

      K=IDUM/IQ
      IDUM=IA*(IDUM-K*IQ)-IR*K
      IF(IDUM < 0)IDUM=IDUM+IM
      J=1+IY/NDIV
      IY=IV(J)
      IV(J)=IDUM
      RAN1=MIN(AM*IY,RNMX)

      RETURN
      END FUNCTION RAN1

!======================================================================+
   FUNCTION GASDEV(IDUM)
!======================================================================+
!
!   Function from Press et al 1992, p280
!
!----------------------------------------------------------------------+

      IMPLICIT                                  NONE

      REAL                                      :: &
            GASDEV

!   References to FUNCTION names not needed as accessed through module

      REAL,EXTERNAL                             :: &
            RAN1

      INTEGER, INTENT(INOUT)                    :: &
            IDUM

      INTEGER                                   :: &
            ISET

      REAL                                      :: &
            FAC,GSET,RSQ,V1,V2

      SAVE ISET,GSET

      DATA ISET/0/

      IF(ISET == 0)THEN
        DO 10
          V1=2.0*RAN1(IDUM)-1
          V2=2.0*RAN1(IDUM)-1
          RSQ=V1**2+V2**2
          IF(RSQ >= 1.0.OR.RSQ == 0.0)CYCLE  ! CYCLE DO 10
          EXIT                           ! EXIT DO 10
   10    ENDDO

        FAC=SQRT(-2.0*LOG(RSQ)/RSQ)
        GSET=V1*FAC
        GASDEV=V2*FAC
        ISET=1
      ELSE
        GASDEV=GSET
        ISET=0
      ENDIF

      RETURN
      END FUNCTION GASDEV

!======================================================================+
   FUNCTION BRENT(AX,BX,CX,F,TOL,XMIN,XFUNC,YFUNC,NFUNC)
!======================================================================+

      IMPLICIT                                  NONE

      DOUBLE PRECISION                          :: &
            BRENT

      DOUBLE PRECISION,EXTERNAL                 :: &
            F

      INTEGER,INTENT(IN)                        :: &
            NFUNC

      DOUBLE PRECISION,INTENT(IN)               :: &
            AX,BX,CX,TOL

      DOUBLE PRECISION,INTENT(IN),DIMENSION(NFUNC) :: &
            XFUNC,YFUNC

      DOUBLE PRECISION,INTENT(INOUT)            :: &
            XMIN

      INTEGER                                   :: &
            ITER,IERR

      INTEGER,PARAMETER                         :: &
            ITMAX = 100

      DOUBLE PRECISION,PARAMETER                :: &
            CGOLD = 0.3819660E+00,ZEPS = 1.0E-18

      DOUBLE PRECISION                          :: &
            A,B,D,E,ETEMP,FU,FV,FW,FX,P,Q,R,TOL1,TOL2,U,V,W,X,XM

      A=MIN(AX,CX)
      B=MAX(AX,CX)
      V=BX
      W=V
      X=V
      E=0.D0
      FX=F(X,XFUNC,YFUNC,NFUNC)
      FV=FX
      FW=FX
      DO 11 ITER=1,ITMAX
        XM=0.5D0*(A+B)
        TOL1=TOL*ABS(X)+ZEPS
        TOL2=2.D0*TOL1
        IF(ABS(X-XM) <= (TOL2-.5D0*(B-A))) GOTO 3
        IF(ABS(E) > TOL1) THEN
          R=(X-W)*(FX-FV)
          Q=(X-V)*(FX-FW)
          P=(X-V)*Q-(X-W)*R
          Q=2.D0*(Q-R)
          IF(Q > 0.D0) P=-P
          Q=ABS(Q)
          ETEMP=E
          E=D
          IF(ABS(P) >= ABS(.5D0*Q*ETEMP).OR.P <= Q*(A-X).OR.P >= Q*(B- &
      X))GOTO 1
          D=P/Q
          U=X+D
          IF(U-A < TOL2 .OR. B-U < TOL2) D=SIGN(TOL1,XM-X)
          GOTO 2
        ENDIF
    1    IF(X >= XM) THEN
          E=A-X
        ELSE
          E=B-X
        ENDIF
        D=CGOLD*E
    2    IF(ABS(D) >= TOL1) THEN
          U=X+D
        ELSE
          U=X+SIGN(TOL1,D)
        ENDIF
        FU=F(U,XFUNC,YFUNC,NFUNC)
        IF(FU <= FX) THEN
          IF(U >= X) THEN
            A=X
          ELSE
            B=X
          ENDIF
          V=W
          FV=FW
          W=X
          FW=FX
          X=U
          FX=FU
        ELSE
          IF(U < X) THEN
            A=U
          ELSE
            B=U
          ENDIF
          IF(FU <= FW .OR. W == X) THEN
            V=W
            FV=FW
            W=U
            FW=FU
          ELSE IF(FU <= FV .OR. V == X .OR. V == W) THEN
            V=U
            FV=FU
          ENDIF
        ENDIF
   11  ENDDO
      IERR = 1
      WRITE(*,*) ' FUNCTION BRENT: Exceeding maximum iterations'
      XMIN=-1.0D0
      GOTO 4
    3  XMIN=X
    4  BRENT=FX
!
      RETURN
      END FUNCTION BRENT

!======================================================================+
   FUNCTION XLS(A,T,X,N)
!======================================================================+

      IMPLICIT                                  NONE

      DOUBLE PRECISION                          :: &
            XLS

      INTEGER,INTENT(IN)                        :: &
            N

      DOUBLE PRECISION,INTENT(IN)               :: &
            A

      DOUBLE PRECISION,INTENT(IN),DIMENSION(N)  :: &
            T,X

      INTEGER                                   :: &
            I

      XLS=0.0D0
      DO 10 I=2,N
        XLS=XLS+(X(I)-X(I-1)*DSIGN(1.0D0,A)*DABS(A)**(T(I)-T(I-1)))** &
      2.0D0
   10  ENDDO

      RETURN
      END FUNCTION XLS

!======================================================================+
   FUNCTION BRENT2(AX,BX,CX,XLS2,TOL,XMIN,XLOUT,XFR,NMED,SOFA,XNYQ)
!======================================================================+

      IMPLICIT                                  NONE

      DOUBLE PRECISION                          :: &
            BRENT2

      DOUBLE PRECISION,EXTERNAL                 :: &
            XLS2

      INTEGER,INTENT(IN)                        :: &
            NMED

      DOUBLE PRECISION,INTENT(IN)               :: &
            AX,BX,CX,TOL,SOFA,XNYQ

      DOUBLE PRECISION,INTENT(INOUT)            :: &
            XMIN

      DOUBLE PRECISION,INTENT(IN),DIMENSION(NMED) :: &
            XLOUT,XFR

      INTEGER                                   :: &
            ITER,IERR

      INTEGER,PARAMETER                         :: &
            ITMAX = 100

      DOUBLE PRECISION                          :: &
            CGOLD = 0.3819660E+00,zeps = 1.0E-18

      DOUBLE PRECISION                          :: &
            A,B,D,E,ETEMP,FU,FV,FW,FX,P,Q,R,TOL1,TOL2, &
            U,V,W,X,XM

      A     = MIN(AX,CX)
      B     = MAX(AX,CX)
      V     = BX
      W     = V
      X     = V
      E     = 0.D0
      FX    = XLS2(X,XLOUT,XFR,NMED,SOFA,XNYQ)
      FV    = FX
      FW    = FX

      DO 11 ITER=1,ITMAX
        XM=0.5D0*(A+B)
        TOL1=TOL*ABS(X)+ZEPS
        TOL2=2.D0*TOL1
        IF(ABS(X-XM) <= (TOL2-.5D0*(B-A))) GOTO 3
        IF(ABS(E) > TOL1) THEN
          R=(X-W)*(FX-FV)
          Q=(X-V)*(FX-FW)
          P=(X-V)*Q-(X-W)*R
          Q=2.D0*(Q-R)
          IF(Q > 0.D0) P=-P
          Q=ABS(Q)
          ETEMP=E
          E=D
          IF(ABS(P) >= ABS(.5D0*Q*ETEMP).OR.P <= Q*(A-X).OR.P >= Q*(B- &
      X))GOTO 1
          D=P/Q
          U=X+D
          IF(U-A < TOL2 .OR. B-U < TOL2) D=SIGN(TOL1,XM-X)
          GOTO 2
        ENDIF
    1    IF(X >= XM) THEN
          E=A-X
        ELSE
          E=B-X
        ENDIF
        D=CGOLD*E
    2    IF(ABS(D) >= TOL1) THEN
          U=X+D
        ELSE
          U=X+SIGN(TOL1,D)
        ENDIF
        FU=XLS2(U,XLOUT,XFR,NMED,SOFA,XNYQ)
        IF(FU <= FX) THEN
          IF(U >= X) THEN
            A=X
          ELSE
            B=X
          ENDIF
          V=W
          FV=FW
          W=X
          FW=FX
          X=U
          FX=FU
        ELSE
          IF(U < X) THEN
            A=U
          ELSE
            B=U
          ENDIF
          IF(FU <= FW .OR. W == X) THEN
            V=W
            FV=FW
            W=U
            FW=FU
          ELSE IF(FU <= FV .OR. V == X .OR. V == W) THEN
            V=U
            FV=FU
          ENDIF
        ENDIF
   11  ENDDO
      IERR = 1
      WRITE(*,*) ' FUNCTION BRENT2: EXCEEDING MAXIMUM ITERATIONS'
      XMIN=2.0D0
      GOTO 4
    3  XMIN=X
    4  BRENT2=FX

      RETURN
      END FUNCTION BRENT2

!======================================================================+
   FUNCTION XLS2(ROHT,XLOUT,XFR,NMED,SOFA,XNYQ)
!======================================================================+

      IMPLICIT                                  NONE

      DOUBLE PRECISION                          :: &
            XLS2

      INTEGER,INTENT(IN)                        :: &
            NMED

      DOUBLE PRECISION,INTENT(IN)               :: &
            ROHT,SOFA,XNYQ

      DOUBLE PRECISION,INTENT(IN),DIMENSION(NMED) :: &
            XLOUT,XFR

      INTEGER                                   :: &
            I

      REAL                                      :: &
            SXARP2
 
      DOUBLE PRECISION                          :: &
            ROHT2,XARP2,XLP2,XPI

      XPI=3.141592654
      XLS2=0.0D0
      ROHT2=ROHT*ROHT
      DO 10 I=1,NMED
        XARP2=(SOFA*(1.0-ROHT2))/(1.0-2*ROHT*COS((XFR(I)/XNYQ)*XPI)+ &
         ROHT2)
        SXARP2=XARP2
        XLP2=ALOG10(SXARP2)
        XLS2=XLS2+((XLP2-XLOUT(I))*(XLP2-XLOUT(I)))
   10  ENDDO
!
      RETURN
      END FUNCTION XLS2
!
!======================================================================+
   FUNCTION BRENT3(AX,BX,CX,XLS3,TOL,XMIN,XLOUT,XFR,NMED,ROHB,XNYQ)
!======================================================================+

      IMPLICIT                                  NONE

      DOUBLE PRECISION                          :: &
            BRENT3

      DOUBLE PRECISION,EXTERNAL                 :: &
            XLS3

      INTEGER,INTENT(IN)                        :: &
            NMED

      DOUBLE PRECISION,INTENT(IN)               :: &
            AX,BX,CX,TOL

      DOUBLE PRECISION,INTENT(IN),DIMENSION(NMED) :: &
            XLOUT,XFR

      DOUBLE PRECISION,INTENT(INOUT)             :: &
            XMIN,ROHB,XNYQ

      INTEGER                                    :: &
            ITER,IERR

      INTEGER,PARAMETER                          :: &
            ITMAX = 100

      DOUBLE PRECISION,PARAMETER                 :: &
            CGOLD = 0.3819660E+00,zeps = 1.0E-18

      DOUBLE PRECISION                           :: &
            A,B,D,E,ETEMP,FU,FV,FW,FX,P,Q,R,TOL1,TOL2, &
            U,V,W,X,XM

      A     = MIN(AX,CX)
      B     = MAX(AX,CX)
      V     = BX
      W     = V
      X     = V
      E     = 0.D0
      FX    = XLS3(X,XLOUT,XFR,NMED,ROHB,XNYQ)
      FV    = FX
      FW    = FX

      DO 11 ITER=1,ITMAX
        XM=0.5D0*(A+B)
        TOL1=TOL*ABS(X)+ZEPS
        TOL2=2.D0*TOL1
        IF(ABS(X-XM) <= (TOL2-.5D0*(B-A))) GOTO 3
        IF(ABS(E) > TOL1) THEN
          R=(X-W)*(FX-FV)
          Q=(X-V)*(FX-FW)
          P=(X-V)*Q-(X-W)*R
          Q=2.D0*(Q-R)
          IF(Q > 0.D0) P=-P
          Q=ABS(Q)
          ETEMP=E
          E=D
          IF(ABS(P) >= ABS(.5D0*Q*ETEMP).OR.P <= Q*(A-X).OR.P >= Q*(B- &
      X))GOTO 1
          D=P/Q
          U=X+D
          IF(U-A < TOL2 .OR. B-U < TOL2) D=SIGN(TOL1,XM-X)
          GOTO 2
        ENDIF
    1    IF(X >= XM) THEN
          E=A-X
        ELSE
          E=B-X
        ENDIF
        D=CGOLD*E
    2    IF(ABS(D) >= TOL1) THEN
          U=X+D
        ELSE
          U=X+SIGN(TOL1,D)
        ENDIF
        FU=XLS3(U,XLOUT,XFR,NMED,ROHB,XNYQ)
        IF(FU <= FX) THEN
          IF(U >= X) THEN
            A=X
          ELSE
            B=X
          ENDIF
          V=W
          FV=FW
          W=X
          FW=FX
          X=U
          FX=FU
        ELSE
          IF(U < X) THEN
            A=U
          ELSE
            B=U
          ENDIF
          IF(FU <= FW .OR. W == X) THEN
            V=W
            FV=FW
            W=U
            FW=FU
          ELSE IF(FU <= FV .OR. V == X .OR. V == W) THEN
            V=U
            FV=FU
          ENDIF
        ENDIF
   11  ENDDO
      IERR = 1
      WRITE(*,*) ' FUNCTION BRENT3: EXCEEDING MAXIMUM ITERATIONS'
      XMIN=20.0D0
      GOTO 4
    3  XMIN=X
    4  BRENT3=FX

      RETURN
      END FUNCTION BRENT3

!======================================================================+
   FUNCTION XLS3(SOFA,XLOUT,XFR,NMED,ROHB,XNYQ)
!======================================================================+

      IMPLICIT                                  NONE

      DOUBLE PRECISION                          :: &
            XLS3

      INTEGER,INTENT(IN)                        :: &
            NMED

      DOUBLE PRECISION,INTENT(IN)               :: &
            SOFA,ROHB,XNYQ

      DOUBLE PRECISION,INTENT(IN),DIMENSION(NMED) :: &
            XLOUT,XFR

      INTEGER                                   :: &
            I

      REAL                                      :: &
            SXARP2
 
      DOUBLE PRECISION                          :: &
            ROHB2,XARP2,XLP2,XPI

      XPI=3.141592654
      XLS3=0.0D0
      ROHB2=ROHB*ROHB
      DO 10 I=1,NMED
        XARP2=(SOFA*(1.0-ROHB2))/(1.0-2*ROHB*COS((XFR(I)/XNYQ)*XPI)+ &
         ROHB2)
        IF(XARP2 <= 0.0)WRITE(*,*)' *** XARP2 IN XLS3 <= 0.0 ***'
        SXARP2=XARP2
        XLP2=ALOG10(SXARP2)
        XLS3=XLS3+((XLP2-XLOUT(I))*(XLP2-XLOUT(I)))
   10  ENDDO
!
      RETURN
      END FUNCTION XLS3
!
!======================================================================+
!   References:
!======================================================================+

!   Bloomfield, P. 1976, Fourier Analysis of Time Series: an introduction.
!         Wiley, 258p.

!   Davis, J.C. 1973, Statistics and Data Analysis in Geology. Wiley 550p.

!   Mann, M.E., and Lees, J., 1996. Robust estimation of background noise
!         and signal detection in climatic time series. Climate Change
!         33 pp409-445.

!   Mudelsee, M., 2002. TAUEST: a computer program for estimating
!         persistence in unevenly spaced weather/climate time series.
!         Computers & Geosciences 28 pp69-72.

!   Press, W.H., Teukolsky, S.A., Vetterling, W.T., Flannery, B.P., 1992,
!         Numerical Recipes in FORTRAN. The Art of Scientific Computing.
!         Cambridge University Press, Second Edition, 963p.


!   Priestley, M.B. 1981, Spectral Analysis and Time Series. Academic 890p.

!   Schulz, M. and Mudelsee, M., 2002. REDFIT: estimating red-noise
!         spectra directly from unevenly spaced paleoclimatic time
!         series. Computers & Geosciences 28 pp421-426.

!======================================================================+
!END MODULE MODULE_TIME_SERIES_PROCEDURES
!======================================================================+

!======================================================================+
!   SPECTRAL_ANALYSIS
!   Module to derive power/amplitude spectrum from an input time series
!======================================================================+

!   Contact information:
!   Garry Hayman
!   Centre for Ecology and Hydrology
!   Wallingford
!   Oxfordshire, OX10 8BB
!   United Kingdom
   
!   email: garr@ceh.ac.uk
!   http://www.ceh.ac.uk

!   Development history

!   Date      Version     Author      Description
!   11/13      1.0        GDH         Initial version developed as
!                                     FORTRAN module from Graham Weedon's
!                                     time series FORTRAN code for use
!                                     with PYTHON.
!
!   11/19      1.1        GDH         Release version

!======================================================================+
MODULE MODULE_spectral_analysis
!======================================================================+

!   Declaration of Parameters/Variables common to
!   all routines and modules

      USE                                       MODULE_READ_netCDF

      IMPLICIT                                  NONE

      INTEGER                                   :: &
            ncid,ndims,nvars,nlat,nlon,ntimes,ncomp,nsize,nfreqs,nrays, &
            nspect,nspect_valid,nmin,nfit0,kmk,iconf_level,use_flag

      INTEGER                                   :: &
            dim_sizes(max_dims),dim_data(max_dims+1)

      INTEGER,ALLOCATABLE                       :: &
            dates(:,:),                            &
            error_flags(:,:)

      REAL                                      :: &
            xlowp,xhighp,xnyq,xnray,crit_freq,miss_data

      REAL,ALLOCATABLE                          :: &
            data_time_series(:,:),                 &
            time_time_series(:),                   &
            xnrays(:),                             &
            power(:,:),                            &
            crit_freqs(:)

      CHARACTER (len=  1)                       :: &
            debug,pnyq,pray,cfix,cpower,cfit

      CHARACTER (len=41)                        :: &
            header

      CHARACTER (len=60)                        :: &
            fmt

      CHARACTER (len=30)                        :: &
            var_name_data

      CHARACTER (len=80)                        :: &
            runfile, &                          ! Time Series runfile
            tsfile                              ! Time Series data file

     CHARACTER (len=150)                        :: &
            file_power                          ! Time Series output file

      SAVE

!======================================================================+
END MODULE MODULE_spectral_analysis
!======================================================================+
PROGRAM spectral_analysis_main
!======================================================================+

!   Declaration of Parameters/Variables common to
!   all routines and modules

      USE                                       MODULE_SPECTRAL_ANALYSIS
      USE                                       MODULE_TIME_SERIES_PROCEDURES

      IMPLICIT                                  NONE

      CHARACTER (len=  1)                       :: &
            option

      CHARACTER (len= 30)                       :: &
            time_series_routine

!----------------------------------------------------------------------+
!   Start of program
!-----------------------------------------------------------------------

      time_series_routine    = 'spectral_analysis_main'

!   Read runfile to set parameters

      CALL TIME_SERIES_RUNFILE(tsfile,fmt,pnyq,pray,cfix,cpower, &
           cfit,nmin,nfit0,kmk,xlowp,xhighp,xnray,xnyq, &
           crit_freq,iconf_level,miss_data)

!   Set debug

      WRITE(*,*) 'Input value for debug (Y/N): '
      READ (5,*) debug

!   Point or gridded run

      WRITE(*,*) 'Gridded or point run  (G/P): '
      READ (5,*) option

!   Get data: gridded or point
!   (a) gridded

      IF (option == 'G') THEN

        CALL spectral_analysis_gridded

!   Dimensions of original time series data

        dim_data    = 0
        dim_data(1) = 3
        dim_data(2) = ntimes
        dim_data(3) = nlat
        dim_data(4) = nlon
        WRITE(*,*) dim_data

      END IF

!   (b) point

      IF (option == 'P') THEN

        CALL spectral_analysis_point

!   Dimensions of original time series data

        dim_data    = 0
        dim_data(1) = 2
        dim_data(2) = ntimes
        dim_data(3) = nsize
        WRITE(*,*) dim_data

      END IF

!-----------------------------------------------------------------------
!   Perform time-series analysis on var_data
!-----------------------------------------------------------------------

!   Derive number of Rayleigh frequencies

      nrays       = nsize
      use_flag    = 1

      crit_freqs(1) = 1.0
      crit_freqs(2) = 2.0

!   Extract valid data from input and reshape into 2-D arrays
!   (vector x time)

      CALL TIME_SERIES_GRIDDED(var_name_data,dim_data, &
           nsize,ntimes,nfreqs,nrays, &
           data_time_series,time_time_series,dates,power, &
           nmin,xlowp,xhighp,xnrays,xnyq,fmt,nfit0,kmk,pray,pnyq, &
           cfix,cpower,cfit,crit_freqs,iconf_level,miss_data,file_power, &
           nspect,nspect_valid,use_flag,error_flags,debug)
      
      WRITE(*,*)
      WRITE(*,*) 'total and valid # of time series: ',nspect,nspect_valid
      WRITE(*,*) 'Power spectrum at ... '
      WRITE(*,'(10(1X,F10.4))') power(:,1)

!   Deallocate array

      DEALLOCATE(time_time_series)
      DEALLOCATE(data_time_series)
      DEALLOCATE(dates)
      DEALLOCATE(xnrays)
      DEALLOCATE(power)
      DEALLOCATE(error_flags)
      DEALLOCATE(crit_freqs)

!   End of program

      STOP

!======================================================================+
END PROGRAM spectral_analysis_main
!======================================================================+
SUBROUTINE spectral_analysis_gridded
!======================================================================+

!   Declaration of Parameters/Variables common to
!   all routines and modules

      USE                                       MODULE_SPECTRAL_ANALYSIS
      USE                                       MODULE_READ_netCDF

      IMPLICIT                                  NONE

      INTEGER                                   :: &
            idim,ivar,itime,ilong,ilat,igrid, &
            idx_data,idx_time,idx_ray,idx_lat,idx_lon,idx_date

      INTEGER                                   :: &
            var_ids(max_vars),nvar_dims(max_vars),var_dimids(max_vars,max_dims)

      REAL,ALLOCATABLE                          :: &
            var_data_2d(:,:),                      &
            var_data_3d(:,:,:),                    &
            var_data_4d(:,:,:,:)

      CHARACTER (len= 30)                       :: &
            dim_names(max_dims),var_names(max_vars)

      CHARACTER (len=30)                        :: &
            var_name_date,var_name_ray

      CHARACTER (len= 30)                       :: &
            time_series_routine

      CHARACTER (len=100)                       :: &
            cdfinput                            ! netCDF filename

!----------------------------------------------------------------------+
!   Start of program
!-----------------------------------------------------------------------

      time_series_routine    = 'spectral_analysis_gridded'

!   Get netCDF filename

      WRITE(*,*) 'Input netCDF Filename: '
      READ (5,*) cdfinput

      file_power = 'z_OUTPUT/time_series_power_' // &
           cdfinput(9:len_trim(cdfinput)-3)  

!   Open netCDF file
!   netCDF call to open netCDF file using netCDF FORTRAN interface

      CALL READ_netCDF_open(cdfinput,ncid)

!   Get metadata from netCDF file

      CALL READ_netCDF_info(ncid,ndims,nvars,dim_names,dim_sizes,debug)

!   Identify and assign dimensions

      DO idim=1,ndims 
        IF (dim_names(idim) == 'time')       ntimes = dim_sizes(idim)
        IF (dim_names(idim) == 'latitude')   nlat   = dim_sizes(idim)
        IF (dim_names(idim) == 'longitude')  nlon   = dim_sizes(idim)
        IF (dim_names(idim) == 'component')  ncomp  = dim_sizes(idim)
      END DO

      nsize         = nlat*nlon
      nfreqs        = 2

!   Get information on variables

      CALL READ_netCDF_varinfo(ncid,nvars,var_names,var_ids, &
           nvar_dims,var_dimids,debug)

!   List variables in netCDF file

      DO ivar = 1,nvars
        WRITE(*,'(I4,1X,A30)') ivar,var_names(ivar)
      END DO

!   Get index of netCDF data variable name

      WRITE(*,*)
      WRITE(*,*) 'Input index for variable name: '
      READ (5,*) idx_data
      var_name_data = var_names(idx_data)

!   Allocate arrays

      ALLOCATE(time_time_series(ntimes))
      ALLOCATE(data_time_series(nsize,ntimes))
      ALLOCATE(dates(ntimes,ncomp))
      ALLOCATE(xnrays(nsize))
      ALLOCATE(power(nsize,nfreqs))
      ALLOCATE(error_flags(nsize,nfreqs))
      ALLOCATE(crit_freqs(nfreqs))

      ALLOCATE(var_data_2d(nlon,nlat))
      ALLOCATE(var_data_3d(nlon,nlat,ntimes))

!   Identify and get time data  
      
      idx_time      = -1
      DO ivar = 1,nvars
         IF (var_names(ivar) == 'time') idx_time = ivar
      END DO

      IF (idx_time == -1) THEN
         WRITE(*,*) 'Variable time not present in netCDF file'
         STOP
      ELSE
         CALL READ_netCDF_get_vardata_1d_real(ncid,var_ids(idx_time), &
              ntimes,time_time_series,debug)
      END IF

!   Get data, variable id input

      CALL READ_netCDF_get_vardata_3D_real(ncid,var_ids(idx_data), &
           nlon,nlat,ntimes,var_data_3d,debug)

      IF (DEBUG == 'Y') &
         WRITE(*,*) 'Shape of var_data_3d = ',SHAPE(var_data_3d), &
                    nsize,(dim_sizes(idim),idim=1,ndims)

      data_time_series = RESHAPE(var_data_3d,(/nsize,ntimes/))

!   Get index of variable name for Rayleigh frequencies

      idx_ray       = -1
      var_name_ray  = trim(var_name_data) // '_Ray_Freq'

      DO ivar = 1,nvars
         IF (var_name_ray == var_names(ivar)) idx_ray = ivar
      END DO

      IF (idx_ray == -1) THEN
         WRITE(*,*) 'Variable of Rayleigh frequencies not present in netCDF file'
         STOP
      ELSE
         CALL READ_netCDF_get_vardata_2D_real(ncid,var_ids(idx_ray), &
              nlon,nlat,var_data_2d,debug)
         xnrays        = RESHAPE(var_data_2d,(/nsize/))
      END IF

!   Get index of variable date

      idx_date      = -1
      var_name_date = 'date'

      DO ivar = 1,nvars
         IF (var_name_date == var_names(ivar)) idx_date = ivar
      END DO

      IF (idx_date == -1) THEN
         WRITE(*,*) 'Variable date not present in netCDF file'
         STOP
      ELSE
         CALL READ_netCDF_get_vardata_2d_int(ncid,var_ids(idx_date), &
              ntimes,ncomp,dates,debug)
      END IF

      WRITE(*,*) dates(1,:),dates(ntimes,:)

!   Close netCDF file

      CALL READ_netCDF_close(ncid)
 
!   Check to see if data correct
     
  100 CONTINUE

        IF (DEBUG == 'Y') THEN 

          WRITE(*,*)
          WRITE(*,*) 'Check to see if data correct'
          WRITE(*,*) 'Input grid square co-ordinates'
          WRITE(*,*) 'To exit, enter negative value for ilong'
          WRITE(*,*) 'ilong,ilat: '
          READ (5,*) ilong,ilat

          IF (ilong < 0) THEN
            GOTO 200
          ELSE

            igrid     = (ilat-1)*dim_sizes(3)+ilong 
            WRITE(*,'(A11,3I10)')     'Elements = ',ilong,ilat,igrid
            WRITE(*,*)
            WRITE(*,'(11X,10F10.4)')   &
                 (time_time_series(itime),itime=1,ntimes)
            WRITE(*,'(11X,10F10.4)')   &
                 (var_data_3d(ilong,ilat,itime),itime=1,ntimes)
            WRITE(*,'(11X,10F10.4)')   &
                 (data_time_series(igrid,itime),itime=1,ntimes)

            GOTO 100
          END IF

        END IF

  200 CONTINUE

      DEALLOCATE(var_data_2d)
      DEALLOCATE(var_data_3d)

!   Return to calling routine

      RETURN

!======================================================================+
END SUBROUTINE spectral_analysis_gridded
!======================================================================+
SUBROUTINE spectral_analysis_point
!======================================================================+

!   Declaration of Parameters/Variables common to
!   all routines and modules

      USE                                       MODULE_SPECTRAL_ANALYSIS

      IMPLICIT                                  NONE

      INTEGER                                   :: &
            ierror,itime

      CHARACTER (len= 30)                       :: &
            time_series_routine

      CHARACTER (len=120)                       :: &
            iomsg

!----------------------------------------------------------------------+
!   Start of program
!-----------------------------------------------------------------------

      time_series_routine    = 'spectral_analysis_point'

!   Open file

      WRITE(*,*) 'Reading: ',trim(tsfile)

      file_power = 'z_OUTPUT/time_series_power_' // &
           tsfile(9:len_trim(tsfile)-4)  

      OPEN (UNIT=1, FILE=trim(tsfile), STATUS='OLD', FORM='Formatted', &
           IOSTAT=ierror, IOMSG=iomsg)

      IF (ierror > 0) WRITE(*,'(I3,A2,A150)') ierror,': ',iomsg

!   First line has number of data/time points

      READ(1,*) ntimes

      nsize         = 1
      ncomp         = 3
      nfreqs        = 2

!   Allocate arrays

      ALLOCATE(time_time_series(ntimes))
      ALLOCATE(data_time_series(nsize,ntimes))
      ALLOCATE(dates(ntimes,ncomp))
      ALLOCATE(xnrays(nsize))
      ALLOCATE(power(nsize,nfreqs))
      ALLOCATE(error_flags(nsize,nfreqs))
      ALLOCATE(crit_freqs(nfreqs))

!   Input data
!     format: 'YEAR-MM-DD          56                      0.104900'

      DO itime = 1,ntimes
        READ(1,'(I4,1X,I2,1X,I2,F12.2,F30.6)') &
          dates(itime,1),dates(itime,2),dates(itime,3), &
          time_time_series(itime),data_time_series(nsize,itime) 
      END DO
 
      CLOSE(UNIT=1)

      xnrays(1) = xnray
 
!   Check to see if data correct

      IF (DEBUG == 'Y') THEN 

        WRITE(*,*) 'Checking data:'

        DO itime = 1,ntimes
          WRITE(1,'(I4,1X,I2,1X,I2,F12.2,F30.6)') &
            dates(itime,1),dates(itime,2),dates(itime,3), &
            time_time_series(itime),data_time_series(nsize,itime) 
        END DO

      END IF

!   Return to calling routine

      RETURN

!======================================================================+
END SUBROUTINE spectral_analysis_point
!======================================================================+

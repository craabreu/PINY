/*-------------------------------------------------------------------------*/
/*       Velocity resampling routines                                      */

void sampl_vx(CLATOMS_INFO *,CLATOMS_POS *,
              SIMOPTS *,int *,int *,double *);

void sampl_vnhc(THERM_INFO *,THERM_POS *, 
                THERM_INFO *,THERM_POS *, 
                BARO *,PAR_RAHMAN *,
                ENSOPTS *,STATEPOINT *, INT_SCR *,
		int, int *, int *,double *,int,int,int,int,int );

void sampl_vx_vnhc_isok(CLATOMS_INFO *,CLATOMS_POS *,THERM_INFO *,THERM_POS *,
		INT_SCR *,SIMOPTS *,ENSOPTS *,int *,int *,double *);

void sampl_isok_restart(CLATOMS_INFO *,CLATOMS_POS *,THERM_INFO *,THERM_POS *);

/*-------------------------------------------------------------------------*/
/*       Proj_vels                                                        */

void proj_vel(CLASS *,BONDED *,GENERAL_DATA *);

void proj_vel_rollf(CLASS *,BONDED *,GENERAL_DATA *);

void proj_vel_rolli(CLASS *,BONDED *,GENERAL_DATA *);









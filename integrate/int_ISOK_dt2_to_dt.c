/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void int_dt2_to_dt_nvt_isok(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,
                            int ir_tra,int ir_tor,int ir_ter,double dt)

/*========================================================================*/
   {/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */
   int *inhc_x	   	          = class->therm_info_class.inhc_x;
   int *inhc_y		            = class->therm_info_class.inhc_y;
   int *inhc_z		            = class->therm_info_class.inhc_z;
   int i,ipart,ichain;
   int ifirst                 = 1;
   int iflag_mass             = 1;
   int tol_glob               = 1;
   int iii;
   int ix_now;
   double arg,s,sdot,a,b,rb;
   double *clatoms_vx         = class->clatoms_pos[1].vx;
   double *clatoms_vy         = class->clatoms_pos[1].vy;
   double *clatoms_vz         = class->clatoms_pos[1].vz;
   double *clatoms_fx         = class->clatoms_pos[1].fx;
   double *clatoms_fy         = class->clatoms_pos[1].fy;
   double *clatoms_fz         = class->clatoms_pos[1].fz;
   double *clatoms_mass       = class->clatoms_info.mass;
   double **therm_v1          = class->therm_class.x_nhc;
   double **therm_gkt		      = class->therm_info_class.gkt;
   int ix_respa               = general_data->timeinfo.ix_respa;
   int int_res_tra            = general_data->timeinfo.int_res_tra;
   int int_res_ter            = general_data->timeinfo.int_res_ter;
   int iconstrnt              = bonded->constrnt.iconstrnt;
   int myatm_start            = class->clatoms_info.myatm_start;
   int myatm_end              = class->clatoms_info.myatm_end;
   int mytherm_start          = class->therm_info_class.mytherm_start;
   int mytherm_end            = class->therm_info_class.mytherm_end;
   int nres_tra               = general_data->timeinfo.nres_tra;
   int nres_tor               = general_data->timeinfo.nres_tor;
   int nres_ter               = general_data->timeinfo.nres_ter;

/*==========================================================================*/
/* 0) Useful constants                                                      */

   int len_nhc            = class->therm_info_class.len_nhc;    
   double lennhc          = (double)len_nhc;      
   double LkT             = lennhc*therm_gkt[1][inhc_x[1]];
   double kinetic         = 0.0;
   double dt2             = 0.5*dt;

   /*==========================================================================*/
/* VI) Isokinetic Constraint, Force-dependent portion                      */
	  for(ipart=myatm_start;ipart<=myatm_end;ipart++){
		  a=clatoms_fx[ipart]*clatoms_vx[ipart]/LkT;
		  b=clatoms_fx[ipart]*clatoms_fx[ipart]/(clatoms_mass[ipart]*LkT);
		  rb=sqrt(b);
		  arg=dt2*rb;
                  if(arg > 0.00001) {
                     s = (1.0/rb)*sinh(arg) + (a/b)*(cosh(arg)-1.0);
                     sdot = cosh(arg) + (a/rb)*sinh(arg);
                  } else {
                     s = ((((b*a/24.0)*dt2 + b/6.0)*dt2 + 0.5*a)*dt2 + 1.0)*dt2;
                     sdot = (((b*a/6.0)*dt2 + 0.5*b)*dt2 + a)*dt2 + 1.0;
                  }
		  clatoms_vx[ipart]=(clatoms_vx[ipart] + (clatoms_fx[ipart]*s/clatoms_mass[ipart]))/sdot;
		  for(ichain=1;ichain<=len_nhc;ichain++){
			  therm_v1[ichain][inhc_x[ipart]]/=sdot;
		  }

		  a=clatoms_fy[ipart]*clatoms_vy[ipart]/LkT;
		  b=clatoms_fy[ipart]*clatoms_fy[ipart]/(clatoms_mass[ipart]*LkT);
		  rb=sqrt(b);
		  arg=dt2*rb;
                  if(arg > 0.00001) {
                     s = (1.0/rb)*sinh(arg) + (a/b)*(cosh(arg)-1.0);
                     sdot = cosh(arg) + (a/rb)*sinh(arg);
                  } else {
                     s = ((((b*a/24.0)*dt2 + b/6.0)*dt2 + 0.5*a)*dt2 + 1.0)*dt2;
                     sdot = (((b*a/6.0)*dt2 + 0.5*b)*dt2 + a)*dt2 + 1.0;
                  }
		  clatoms_vy[ipart]=(clatoms_vy[ipart] + (clatoms_fy[ipart]*s/clatoms_mass[ipart]))/sdot;
		  for(ichain=1;ichain<=len_nhc;ichain++){
			  therm_v1[ichain][inhc_y[ipart]]/=sdot;
		  }

		  a=clatoms_fz[ipart]*clatoms_vz[ipart]/LkT;
		  b=clatoms_fz[ipart]*clatoms_fz[ipart]/(clatoms_mass[ipart]*LkT);
		  rb=sqrt(b);
		  arg=dt2*rb;
                  if(arg > 0.00001) {
                     s = (1.0/rb)*sinh(arg) + (a/b)*(cosh(arg)-1.0);
                     sdot = cosh(arg) + (a/rb)*sinh(arg);
                  } else {
                     s = ((((b*a/24.0)*dt2 + b/6.0)*dt2 + 0.5*a)*dt2 + 1.0)*dt2;
                     sdot = (((b*a/6.0)*dt2 + 0.5*b)*dt2 + a)*dt2 + 1.0;
                  }
		  clatoms_vz[ipart]=(clatoms_vz[ipart] + (clatoms_fz[ipart]*s/clatoms_mass[ipart]))/sdot;
		  for(ichain=1;ichain<=len_nhc;ichain++){
			  therm_v1[ichain][inhc_z[ipart]]/=sdot;
		  }
	  }
	  /*kinetic=0.0;
	      for(ipart=myatm_start;ipart<=myatm_end;ipart++){
	    	  kinetic+=clatoms_mass[ipart]*clatoms_vx[ipart]*clatoms_vx[ipart];
	    	  kinetic+=clatoms_mass[ipart]*clatoms_vy[ipart]*clatoms_vy[ipart];
	    	  kinetic+=clatoms_mass[ipart]*clatoms_vz[ipart]*clatoms_vz[ipart];
	      }
        printf("after second v-v1 forces,dti4 %g dti2 %g dti %g ATM KE %g \n",class->therm_info_class.wdti4[1],
                class->therm_info_class.wdti2[1],class->therm_info_class.dti_nhc,
                kinetic);*/
/*==========================================================================*/

/*==========================================================================*/

/* VII) Isokinetic Constraint, Nose-like portion                             */

   if( (int_res_tra==0) && (int_res_ter==0) ){
       apply_NH_ISOK_par(&(class->clatoms_info),&(class->clatoms_pos[1]),
                          &(class->therm_info_class),&(class->therm_class),
                          &(class->int_scr),iflag_mass,&(class->class_comm_forc_pkg));
          /*kinetic=0.0;
        for(ipart=myatm_start;ipart<=myatm_end;ipart++){
          kinetic+=clatoms_mass[ipart]*clatoms_vx[ipart]*clatoms_vx[ipart];
          kinetic+=clatoms_mass[ipart]*clatoms_vy[ipart]*clatoms_vy[ipart];
          kinetic+=clatoms_mass[ipart]*clatoms_vz[ipart]*clatoms_vz[ipart];
        }
        printf("after nhc update,dti4 %g dti2 %g dti %g ATM KE %g \n",class->therm_info_class.wdti4[1],
                class->therm_info_class.wdti2[1],class->therm_info_class.dti_nhc,
                kinetic);   */
      }else{
       ix_now = 4;
       if((ir_tra==1)){ix_now=3;}
       if((ir_tra==1)&&(ir_tor==1)){ix_now=2;}
       if((ir_tra==1)&&(ir_tor==1)&&(ir_ter==1)){ix_now=1;}
       if(ix_respa>=ix_now){
           apply_NH_ISOK_par(&(class->clatoms_info),&(class->clatoms_pos[1]),
                              &(class->therm_info_class),&(class->therm_class),
                              &(class->int_scr),iflag_mass,&(class->class_comm_forc_pkg));
             /* kinetic=0.0;
        for(ipart=myatm_start;ipart<=myatm_end;ipart++){
          kinetic+=clatoms_mass[ipart]*clatoms_vx[ipart]*clatoms_vx[ipart];
          kinetic+=clatoms_mass[ipart]*clatoms_vy[ipart]*clatoms_vy[ipart];
          kinetic+=clatoms_mass[ipart]*clatoms_vz[ipart]*clatoms_vz[ipart];
        }
        printf("after nhc update,dti4 %g dti2 %g dti %g ATM KE %g \n",class->therm_info_class.wdti4[1],
                class->therm_info_class.wdti2[1],class->therm_info_class.dti_nhc,
                kinetic);   */  
      }
   }

/*--------------------------------------------------------------------------*/
  }/*end routine*/
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void int_0_to_dt2_nvt_isok(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,
                           int ir_tra,int ir_tor,int ir_ter,double dt)

/*========================================================================*/
    {/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

  int i,ipart,ichain,ifirst=1;
  int iflag_mass = 1;
  double dt2,tol_glob=0.0;
  int iii;
  int ix_now;
  double *random;
  double cst1, cst2, cst3;
  double a,b,rb,s,sdot,arg;   /* Isokinetic Constraint Variables */
  int *inhc_x                = class->therm_info_class.inhc_x;
  int *inhc_y		     = class->therm_info_class.inhc_y;
  int *inhc_z		     = class->therm_info_class.inhc_z;
  double gamma		     = 1.0/(20.0*dt);   /* Stochastic Term */
/*   double gamma		     = 1.0/(20.0*60.0*dt);  Stochastic Term */
  int *iseed		     = &class->vel_samp_class.iseed;
  int *iseed2		     = &class->vel_samp_class.iseed2;
  double *qseed		     = &class->vel_samp_class.qseed;
  double *clatoms_xold       = class->clatoms_info.xold;
  double *clatoms_yold       = class->clatoms_info.yold;
  double *clatoms_zold       = class->clatoms_info.zold;
  double *clatoms_x          = class->clatoms_pos[1].x;
  double *clatoms_y          = class->clatoms_pos[1].y;
  double *clatoms_z          = class->clatoms_pos[1].z;
  double *clatoms_vx         = class->clatoms_pos[1].vx;
  double *clatoms_vy         = class->clatoms_pos[1].vy;
  double *clatoms_vz         = class->clatoms_pos[1].vz;
  double *clatoms_fx         = class->clatoms_pos[1].fx;
  double *clatoms_fy         = class->clatoms_pos[1].fy;
  double *clatoms_fz         = class->clatoms_pos[1].fz;
  double *clatoms_mass       = class->clatoms_info.mass;
  double **therm_mass_nhc    = class->therm_info_class.mass_nhc;
  double **therm_v1          = class->therm_class.x_nhc;
  double **therm_v2          = class->therm_class.v_nhc;
  double **therm_gkt	     = class->therm_info_class.gkt;
  double *int_scr_atm_kin    = class->int_scr.atm_kin;			  

  int myatm_start            = class->clatoms_info.myatm_start;
  int myatm_end              = class->clatoms_info.myatm_end;
  int mytherm_start          = class->therm_info_class.mytherm_start;
  int mytherm_end            = class->therm_info_class.mytherm_end;
  int iconstrnt              = bonded->constrnt.iconstrnt;
  int myid_forc              = class->communicate.myid_forc;
  int int_res_tra            = general_data->timeinfo.int_res_tra;
  int int_res_ter            = general_data->timeinfo.int_res_ter;
  int ix_respa               = general_data->timeinfo.ix_respa;
  int len_nhc 		     = class->therm_info_class.len_nhc;	  
  double lennhc		     = (double)len_nhc;			
  double LkT                 = lennhc*therm_gkt[1][inhc_x[1]];
  double L_L1		     = lennhc/(lennhc+1.0);
  double stoch_e             = exp(-gamma*dt);
  double mass_factor         = 1.0;
  double sigma 		     = sqrt(therm_gkt[1][inhc_x[1]]*(1.0-exp(-2.0*gamma*dt))
		  	  	  	  	    /(mass_factor*therm_mass_nhc[1][inhc_x[1]]));
  double kinetic             = 0.0;
  double vpotnhc	     = 0.0;
  int natm_tot  	     = class->clatoms_info.natm_tot;

/* print isokinetic constraint value*/
if (general_data->timeinfo.itime%50==0){
	for(ipart=1;ipart<=1;ipart++){
	  kinetic += clatoms_mass[ipart]*clatoms_vx[ipart]*clatoms_vx[ipart];
	  kinetic += clatoms_mass[ipart]*clatoms_vy[ipart]*clatoms_vy[ipart];
	  kinetic += clatoms_mass[ipart]*clatoms_vz[ipart]*clatoms_vz[ipart];
  		for(ichain=1;ichain<=len_nhc;ichain++){
	  		vpotnhc += (therm_mass_nhc[ichain][inhc_x[ipart]]
                   		*therm_v1[ichain][inhc_x[ipart]]*therm_v1[ichain][inhc_x[ipart]]);
	  		vpotnhc += (therm_mass_nhc[ichain][inhc_y[ipart]]
                   		*therm_v1[ichain][inhc_y[ipart]]*therm_v1[ichain][inhc_y[ipart]]);
	  		vpotnhc += (therm_mass_nhc[ichain][inhc_z[ipart]]
                   		*therm_v1[ichain][inhc_z[ipart]]*therm_v1[ichain][inhc_z[ipart]]);
}
}
printf("Average IsoK Constraint: %0.16g",(kinetic + (L_L1*vpotnhc))/(3.0*therm_gkt[1][1])  );

}

/*==========================================================================*/
/* 0) Useful constants                                                      */


  dt2 = dt/2.0;

/*==========================================================================*/
/* I)Save positions                                                         */

  for(ipart=myatm_start;ipart<=(myatm_end);ipart++){
    clatoms_xold[ipart] = clatoms_x[ipart];
    clatoms_yold[ipart] = clatoms_y[ipart];
    clatoms_zold[ipart] = clatoms_z[ipart];
  }/*endfor*/

/*==========================================================================*/
/* II) Isokinetic Constraint, Nose-like portion                             */
  if( (int_res_tra==0) && (int_res_ter==0) ){
      apply_NH_ISOK_par(&(class->clatoms_info),&(class->clatoms_pos[1]),
                         &(class->therm_info_class),&(class->therm_class),
                         &(class->int_scr),iflag_mass,&(class->class_comm_forc_pkg));
           /* kinetic=0.0 ;
           for(ipart=myatm_start;ipart<=myatm_end;ipart++){
          kinetic+=clatoms_mass[ipart]*clatoms_vx[ipart]*clatoms_vx[ipart];
          kinetic+=clatoms_mass[ipart]*clatoms_vy[ipart]*clatoms_vy[ipart];
          kinetic+=clatoms_mass[ipart]*clatoms_vz[ipart]*clatoms_vz[ipart];
        }
        printf("after nhc update,dti4 %g dti2 %g dti %g ATM KE %g \n",class->therm_info_class.wdti4[1],
                class->therm_info_class.wdti2[1],class->therm_info_class.dti_nhc,
                kinetic);*/
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
                kinetic);*/
      }
  }



/*==========================================================================*/
/* III) Random Force Generation and Application to {v_2j}                   */

  for(ipart=myatm_start;ipart<=myatm_end;ipart++){
	  gaussran(len_nhc,iseed,iseed2,qseed,int_scr_atm_kin);
	  for(ichain=1;ichain<=len_nhc;ichain++){
		  therm_v2[ichain][inhc_x[ipart]]*=stoch_e;
		  therm_v2[ichain][inhc_x[ipart]]+=(sigma*int_scr_atm_kin[ichain]);
	  }
	  gaussran(len_nhc,iseed,iseed2,qseed,int_scr_atm_kin);
	  for(ichain=1;ichain<=len_nhc;ichain++){
		  therm_v2[ichain][inhc_y[ipart]]*=stoch_e;
		  therm_v2[ichain][inhc_y[ipart]]+=(sigma*int_scr_atm_kin[ichain]);

	  }
	  gaussran(len_nhc,iseed,iseed2,qseed,int_scr_atm_kin);
	  for(ichain=1;ichain<=len_nhc;ichain++){
		  therm_v2[ichain][inhc_z[ipart]]*=stoch_e;
		  therm_v2[ichain][inhc_z[ipart]]+=(sigma*int_scr_atm_kin[ichain]);

	  }
  }


/*==========================================================================*/
/* IV) Isokinetic Constraint, Force-dependent portion                       */

    	  for(ipart=myatm_start;ipart<=myatm_end;ipart++){
    		  a=clatoms_fx[ipart]*clatoms_vx[ipart]/(LkT);
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

    		  a=clatoms_fy[ipart]*clatoms_vy[ipart]/(LkT);
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

    		  a=clatoms_fz[ipart]*clatoms_vz[ipart]/(LkT);
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
   	  /* kinetic=0.0;
	      for(ipart=myatm_start;ipart<=myatm_end;ipart++){
	    	  kinetic+=clatoms_mass[ipart]*clatoms_vx[ipart]*clatoms_vx[ipart];
	    	  kinetic+=clatoms_mass[ipart]*clatoms_vy[ipart]*clatoms_vy[ipart];
	    	  kinetic+=clatoms_mass[ipart]*clatoms_vz[ipart]*clatoms_vz[ipart];
	      }
        printf("after first v-v1 forces,dti4 %g dti2 %g dti %g ATM KE %g \n",class->therm_info_class.wdti4[1],
                class->therm_info_class.wdti2[1],class->therm_info_class.dti_nhc,
                kinetic);*/
/*==========================================================================*/
/* V) Evolve positions                                                     */


  	  for(ipart=myatm_start;ipart<=(myatm_end);ipart++){
  		  clatoms_x[ipart] += clatoms_vx[ipart]*dt;
  		  clatoms_y[ipart] += clatoms_vy[ipart]*dt;
  		  clatoms_z[ipart] += clatoms_vz[ipart]*dt;
  	  }/*endfor*/


/*==========================================================================*/


/*==========================================================================*/
/* i) Recalculate positions of ghost atoms                                  */

  get_ghost_pos(&(class->clatoms_info),&(class->clatoms_pos[1]),
                &(class->ghost_atoms));


/*--------------------------------------------------------------------------*/
  }/*end routine*/
/*==========================================================================*/

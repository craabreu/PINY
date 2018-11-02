/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: int_NVT_ISOK                                 */
/*                                                                          */
/* This subprogram integrates the system using Vel Verlet                   */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

/*#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_gen.h"
#include "../proto_defs/proto_integrate_md_entry.h"
#include "../proto_defs/proto_integrate_md_local.h"
#include "../proto_defs/proto_intra_con_entry.h"
#include "../proto_defs/proto_energy_ctrl_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h" */
#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_gen.h"
#include "../proto_defs/proto_integrate_md_entry.h"
#include "../proto_defs/proto_integrate_md_local.h"
#include "../proto_defs/proto_integrate_pimd_local.h"
#include "../proto_defs/proto_intra_con_entry.h"
#include "../proto_defs/proto_energy_ctrl_entry.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_vel_sampl_class_local.h"
#include "../proto_defs/proto_communicate_wrappers.h"


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void int_NVT_ISOK(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data)
/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */ 
#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_gen.h"
	int iii;
    int i,ipart,iflag,ifirst=1;
    int exit_flag=0;
    double dt,dti,dt2,tol_glob;
    int natm_tot,ichain,inhc;
    int iflag_mass = 1;
    double ann_rate = general_data->simopts.ann_rate;
    int     num_nhc = class->therm_info_class.num_nhc;
    int     len_nhc = class->therm_info_class.len_nhc;
    int ir_tra = 1;
    int ir_tor = 1;
    int ir_ter = 1;
    int anneal_opt = general_data->simopts.anneal_opt;
    double anneal_target_temp = general_data->simopts.ann_target_temp;
    double **therm_gkt      = class->therm_info_class.gkt;
    double **therm_mass_nhc = class->therm_info_class.mass_nhc;
    double *class_clatoms_xold = class->clatoms_info.xold;
    double *class_clatoms_yold = class->clatoms_info.yold;
    double *class_clatoms_zold = class->clatoms_info.zold;
    double *class_clatoms_x    = class->clatoms_pos[1].x;
    double *class_clatoms_y    = class->clatoms_pos[1].y;
    double *class_clatoms_z    = class->clatoms_pos[1].z;
    double *class_clatoms_vx   = class->clatoms_pos[1].vx;
    double *class_clatoms_vy   = class->clatoms_pos[1].vy;
    double *class_clatoms_vz   = class->clatoms_pos[1].vz;
    double *class_clatoms_fx   = class->clatoms_pos[1].fx;
    double *class_clatoms_fy   = class->clatoms_pos[1].fy;
    double *class_clatoms_fz   = class->clatoms_pos[1].fz;
    double *class_clatoms_mass = class->clatoms_info.mass;

    double cst1,cst2,cst3;
    double **therm_v2_nhc     = class->therm_class.v_nhc;
    double **therm_v1_nhc     = class->therm_class.x_nhc;
    int *therm_inhc_x         = class->therm_info_class.inhc_x;
    int *therm_inhc_y         = class->therm_info_class.inhc_y;
    int *therm_inhc_z         = class->therm_info_class.inhc_z;
    double lennhc	      = (double)len_nhc;	     /* For when I need LkT    */

/*==========================================================================*/
/* 0) Useful constants                                                      */

    natm_tot = class->clatoms_info.natm_tot;
    (general_data->timeinfo.int_res_tra)    = 0;
    (general_data->timeinfo.int_res_ter)    = 0;
    dt  = (general_data->timeinfo.dt);
    dti = dt;
    dt2 = dt/2.0;
    class->therm_info_class.dt_nhc  = dt;
    class->therm_info_class.dti_nhc = dt/( (double)(
                             class->therm_info_class.nres_nhc) );
    set_yosh(class->therm_info_class.nyosh_nhc,
             class->therm_info_class.dti_nhc,class->therm_info_class.wdti,
             class->therm_info_class.wdti2,class->therm_info_class.wdti4,
             class->therm_info_class.wdti8,
             class->therm_info_class.wdti16);
    zero_constrt_iters(&(general_data->stat_avg));
    general_data->timeinfo.exit_flag = 0;


/*==========================================================================*/
/* 1) Evolve system from t=0 to dt/2                                        */

    int_0_to_dt2_nvt_isok(class,bonded,general_data,ir_tra,ir_tor,ir_ter,dt);

/*==========================================================================*/
/* 2) Get the new energy/force                                              */

    (class->energy_ctrl.iget_full_inter)= 1;
    (class->energy_ctrl.iget_res_inter) = 0;
    (class->energy_ctrl.iget_full_intra)= 1;
    (class->energy_ctrl.iget_res_intra) = 0;

    energy_control(class,bonded,general_data);


/*==========================================================================*/
/* 3) Evolve system from dt/2 to dt                                         */

    int_dt2_to_dt_nvt_isok(class,bonded,general_data,ir_tra,ir_tor,ir_ter,dt);

/*#define CHECK_CONSTRNT*/
#ifdef CHECK_CONSTRNT
    for(ipart=1;ipart<=natm_tot;ipart++){
    	cst1=0.0;	cst2=0.0;
    	for(ichain=1;ichain<=len_nhc;ichain++){
    		cst1+=therm_mass_nhc[ichain][therm_inhc_x[ipart]]
    		    		   *therm_v1_nhc[ichain][therm_inhc_x[ipart]]*therm_v1_nhc[ichain][therm_inhc_x[ipart]];
    	}
    	cst2+=(class_clatoms_mass[ipart]*class_clatoms_vx[ipart]*class_clatoms_vx[ipart]);
    	cst1*=lennhc/(lennhc+1.0);
    	cst3=(cst1+cst2)/therm_gkt[1][therm_inhc_x[ipart]];
    printf("x %d \t %.14g \n",ipart,cst3);
    	cst1=0.0;	cst2=0.0;
    	for(ichain=1;ichain<=len_nhc;ichain++){
    		cst1+=therm_mass_nhc[ichain][therm_inhc_y[ipart]]
    		    		   *therm_v1_nhc[ichain][therm_inhc_y[ipart]]*therm_v1_nhc[ichain][therm_inhc_y[ipart]];
    	}
    	cst2+=(class_clatoms_mass[ipart]*class_clatoms_vy[ipart]*class_clatoms_vy[ipart]);
    	cst1*=lennhc/(lennhc+1.0);
    	cst3=(cst1+cst2)/therm_gkt[1][therm_inhc_y[ipart]];
    printf("y %d \t %.14g \n",ipart,cst3);
    	cst1=0.0;	cst2=0.0;
    	for(ichain=1;ichain<=len_nhc;ichain++){
    		cst1+=therm_mass_nhc[ichain][therm_inhc_z[ipart]]
    		    		   *therm_v1_nhc[ichain][therm_inhc_z[ipart]]*therm_v1_nhc[ichain][therm_inhc_z[ipart]];
    	}
    	cst2+=(class_clatoms_mass[ipart]*class_clatoms_vz[ipart]*class_clatoms_vz[ipart]);
    	cst1*=lennhc/(lennhc+1.0);
    	cst3=(cst1+cst2)/therm_gkt[1][therm_inhc_z[ipart]];
    printf("z %d \t %.14g \n",ipart,cst3);
    }getchar();
#endif


/*==========================================================================*/
/* 4) Scale by annealing factor                                          */

  iflag=0;
  if(anneal_opt == 1){
    anneal_class(class,ann_rate,iflag,iflag_mass,anneal_target_temp,&exit_flag);
    general_data->timeinfo.exit_flag = exit_flag;
  }/*endif*/

/*==========================================================================*/
/* 5) Finalize                                                              */

   int_final_class(class,bonded,general_data,iflag);

/*==========================================================================*/

#ifdef DEBUG_GLENN
    for(iproc=0;iproc<np_forc;iproc++){
      Barrier(comm_forc);
      if(myid_forc==iproc){
       printf("1 %.12g %.12g %.12g %.12g %.12g %.12g %.12g %.12g %.12g\n",
                  x[1],y[1],z[1],vx[1],vy[1],vz[1],fx[1],fy[1],fz[1]);
       printf("n %.12g %.12g %.12g %.12g %.12g %.12g %.12g %.12g %.12g\n",
                  x[n],y[n],z[n],vx[n],vy[n],vz[n],fx[n],fy[n],fz[n]);
      }
    }
    if(myid_forc==0){scanf("%d",&iii);}
    Barrier(comm_forc);
#endif
/*--------------------------------------------------------------------------*/
/*end routine*/}
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void apply_NH_ISOK_par(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos, 
                       THERM_INFO *therm_info_class,THERM_POS *therm_class, 
                       INT_SCR *int_scr,int iflag_mass,
                       CLASS_COMM_FORC_PKG *class_comm_forc_pkg)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

    int ipart,inhc,ichain;              /* Num: for loop counters */
    int iresn,iyosh;                    /* Num: for loop counters */
    double arg,arg2,aa;               	/* Num: scalar temps      */
    int iii;
    int natm_tot,num_nhc;
    double cst1, cst2, cst3;

/* Define local pointers                                          */
      double *int_scr_atm_kin   = int_scr->atm_kin;
      double *int_scr_sc        = int_scr->sc;
      double *int_scr_sc_temp   = int_scr->sc_temp;
      int *therm_inhc_x         = therm_info_class->inhc_x;
      int *therm_inhc_y         = therm_info_class->inhc_y;
      int *therm_inhc_z         = therm_info_class->inhc_z;
      double *clatoms_vx        = clatoms_pos->vx;
      double *clatoms_vy        = clatoms_pos->vy;
      double *clatoms_vz        = clatoms_pos->vz;
      double **therm_f_nhc      = therm_class->f_nhc;
      double **therm_v2_nhc     = therm_class->v_nhc;
      double **therm_v1_nhc     = therm_class->x_nhc;
      double **therm_gkt      	= therm_info_class->gkt;
      double **therm_mass_nhc 	= therm_info_class->mass_nhc;
      double *therm_wdti2    	= therm_info_class->wdti2;
      double *therm_wdti4    	= therm_info_class->wdti4;
      double *therm_wdti8     	= therm_info_class->wdti8;
      int len_nhc				= therm_info_class->len_nhc; /* Num: length of chains  */
      double lennhc	        	= (double)len_nhc;	     /* For when I need LkT    */
      double *clatoms_mass		= clatoms_pos->mass;
      int mytherm_start         = therm_info_class->mytherm_start;
      int mytherm_end           = therm_info_class->mytherm_end;
      int myatm_start           = clatoms_info->myatm_start;
      int myatm_end             = clatoms_info->myatm_end;
      int *map_share            = therm_info_class->map_share;
      int num_nhc_share         = therm_info_class->num_nhc_share;
      int np_forc               = class_comm_forc_pkg->num_proc;
      MPI_Comm comm_forc        = class_comm_forc_pkg->comm;
      MPI_Comm world            = class_comm_forc_pkg->world;

      if(iflag_mass==1){
       clatoms_mass    = clatoms_info->mass;
     }else{
       clatoms_mass    = clatoms_pos->mass;
     }/*endif*/

      len_nhc=therm_info_class->len_nhc;
      double len_fac=((double)(len_nhc))/(((double)(len_nhc))+1.0);

/* I) Apply the nhc evolution operator using RESPA                         */


   for(iresn=1;iresn<=therm_info_class->nres_nhc;iresn++){
    for(iyosh=1;iyosh<=therm_info_class->nyosh_nhc;iyosh++){
/*--------------------------------------------------------------------------*/
/*  1) Calculate forces {G(v_1j)} and apply to {v_2j}                      */
      for(ipart=myatm_start;ipart<=myatm_end;ipart++){
    	  for(ichain=1;ichain<=len_nhc;ichain++){
              therm_f_nhc[ichain][therm_inhc_x[ipart]] =
            		  therm_mass_nhc[ichain][therm_inhc_x[ipart]]*therm_v1_nhc[ichain][therm_inhc_x[ipart]]*
            		  therm_v1_nhc[ichain][therm_inhc_x[ipart]]-therm_gkt[ichain][therm_inhc_x[ipart]];
              therm_v2_nhc[ichain][therm_inhc_x[ipart]]+=
            		  therm_wdti4[iyosh]*therm_f_nhc[ichain][therm_inhc_x[ipart]]/therm_mass_nhc[ichain][therm_inhc_x[ipart]];

              therm_f_nhc[ichain][therm_inhc_y[ipart]] =
            		  therm_mass_nhc[ichain][therm_inhc_y[ipart]]*therm_v1_nhc[ichain][therm_inhc_y[ipart]]*
            		  therm_v1_nhc[ichain][therm_inhc_y[ipart]]-therm_gkt[ichain][therm_inhc_y[ipart]];
              therm_v2_nhc[ichain][therm_inhc_y[ipart]]+=
            		  therm_wdti4[iyosh]*therm_f_nhc[ichain][therm_inhc_y[ipart]]/therm_mass_nhc[ichain][therm_inhc_y[ipart]];

              therm_f_nhc[ichain][therm_inhc_z[ipart]] =
            		  therm_mass_nhc[ichain][therm_inhc_z[ipart]]*therm_v1_nhc[ichain][therm_inhc_z[ipart]]*
            		  therm_v1_nhc[ichain][therm_inhc_z[ipart]]-therm_gkt[ichain][therm_inhc_z[ipart]];
              therm_v2_nhc[ichain][therm_inhc_z[ipart]]+=
            		  therm_wdti4[iyosh]*therm_f_nhc[ichain][therm_inhc_z[ipart]]/therm_mass_nhc[ichain][therm_inhc_z[ipart]];

    	  /*endfor chains and endfor atoms*/}   }
 /*     for(ipart=myatm_start;ipart<=myatm_end;ipart++){
        printf("%f \t %f \n",therm_v1_nhc[1][therm_inhc_x[ipart]],therm_v2_nhc[1][therm_inhc_x[ipart]]);
        }getchar();*/

/*--------------------------------------------------------------------------*/
/*  2) Calculate and apply Nose-like portion of isokinetic constraint       */


      for(ipart=myatm_start;ipart<=myatm_end;ipart++){
    	  aa = 0.0;
    	  arg=(clatoms_vx[ipart]*clatoms_vx[ipart]*clatoms_mass[ipart]);
       	  for(ichain=1;ichain<=len_nhc;ichain++){
       	  aa+=(len_fac*therm_mass_nhc[ichain][therm_inhc_x[ipart]]*therm_v1_nhc[ichain][therm_inhc_x[ipart]]
    				 *therm_v1_nhc[ichain][therm_inhc_x[ipart]]
    				 *exp(-2.0*therm_wdti2[iyosh]*therm_v2_nhc[ichain][therm_inhc_x[ipart]]));
       	  }
       	  arg +=aa;
       	  arg2=sqrt(lennhc*therm_gkt[1][therm_inhc_x[ipart]]/arg);
    	  clatoms_vx[ipart]*=arg2;
    	  for(ichain=1;ichain<=len_nhc;ichain++){
    	  therm_v1_nhc[ichain][therm_inhc_x[ipart]]*=arg2
    			  *exp(-therm_wdti2[iyosh]*therm_v2_nhc[ichain][therm_inhc_x[ipart]]);
    	  }

    	  aa = 0.0;
          arg=(clatoms_vy[ipart]*clatoms_vy[ipart]*clatoms_mass[ipart]);
    	  for(ichain=1;ichain<=len_nhc;ichain++){
          aa+=(len_fac*therm_mass_nhc[ichain][therm_inhc_y[ipart]]*therm_v1_nhc[ichain][therm_inhc_y[ipart]]
    		    	  *therm_v1_nhc[ichain][therm_inhc_y[ipart]]
    		   		  *exp(-2.0*therm_wdti2[iyosh]*therm_v2_nhc[ichain][therm_inhc_y[ipart]]));
    	  }

    	  arg+=aa;
          arg2=sqrt(lennhc*therm_gkt[1][therm_inhc_y[ipart]]/arg);
    	  clatoms_vy[ipart]*=arg2;
    	  for(ichain=1;ichain<=len_nhc;ichain++){
    	  therm_v1_nhc[ichain][therm_inhc_y[ipart]]*=arg2
    			  *exp(-therm_wdti2[iyosh]*therm_v2_nhc[ichain][therm_inhc_y[ipart]]);
    	  }

    	  aa = 0.0;
          arg=(clatoms_vz[ipart]*clatoms_vz[ipart]*clatoms_mass[ipart]);
    	  for(ichain=1;ichain<=len_nhc;ichain++){
          aa+=(len_fac*therm_mass_nhc[ichain][therm_inhc_z[ipart]]*therm_v1_nhc[ichain][therm_inhc_z[ipart]]
    		    	  *therm_v1_nhc[ichain][therm_inhc_z[ipart]]
    		   		  *exp(-2.0*therm_wdti2[iyosh]*therm_v2_nhc[ichain][therm_inhc_z[ipart]]));
    	  }

    	  arg+=aa;
          arg2=sqrt(lennhc*therm_gkt[1][therm_inhc_z[ipart]]/arg);
    	  clatoms_vz[ipart]*=arg2;
    	  for(ichain=1;ichain<=len_nhc;ichain++){
    	  therm_v1_nhc[ichain][therm_inhc_z[ipart]]*=arg2
    			  *exp(-therm_wdti2[iyosh]*therm_v2_nhc[ichain][therm_inhc_z[ipart]]);
    	  }
      }/*endfor ipart */
     /* for(ipart=myatm_start;ipart<=myatm_end;ipart++){
        printf("%f \t %f \n",therm_v1_nhc[1][therm_inhc_x[ipart]],therm_v2_nhc[1][therm_inhc_x[ipart]]);
        }getchar();*/
/*--------------------------------------------------------------------------*/
/*  3) Calculate forces {G(v_1j)} and apply to {v_2j}                      */


      for(ipart=myatm_start;ipart<=myatm_end;ipart++){
    	  for(ichain=1;ichain<=len_nhc;ichain++){
              therm_f_nhc[ichain][therm_inhc_x[ipart]] =
            		  therm_mass_nhc[ichain][therm_inhc_x[ipart]]*therm_v1_nhc[ichain][therm_inhc_x[ipart]]*
            		  therm_v1_nhc[ichain][therm_inhc_x[ipart]]-therm_gkt[ichain][therm_inhc_x[ipart]];
              therm_v2_nhc[ichain][therm_inhc_x[ipart]]+=
			therm_wdti4[iyosh]*therm_f_nhc[ichain][therm_inhc_x[ipart]]/therm_mass_nhc[ichain][therm_inhc_x[ipart]];

              therm_f_nhc[ichain][therm_inhc_y[ipart]] =
            		  therm_mass_nhc[ichain][therm_inhc_y[ipart]]*therm_v1_nhc[ichain][therm_inhc_y[ipart]]*
            		  therm_v1_nhc[ichain][therm_inhc_y[ipart]]-therm_gkt[ichain][therm_inhc_y[ipart]];
              therm_v2_nhc[ichain][therm_inhc_y[ipart]]+=
			therm_wdti4[iyosh]*therm_f_nhc[ichain][therm_inhc_y[ipart]]/therm_mass_nhc[ichain][therm_inhc_y[ipart]];

              therm_f_nhc[ichain][therm_inhc_z[ipart]] =
            		  therm_mass_nhc[ichain][therm_inhc_z[ipart]]*therm_v1_nhc[ichain][therm_inhc_z[ipart]]*
            		  therm_v1_nhc[ichain][therm_inhc_z[ipart]]-therm_gkt[ichain][therm_inhc_z[ipart]];
              therm_v2_nhc[ichain][therm_inhc_z[ipart]]+=
			therm_wdti4[iyosh]*therm_f_nhc[ichain][therm_inhc_z[ipart]]/therm_mass_nhc[ichain][therm_inhc_z[ipart]];

    	  /*endfor chains and endfor atoms*/}}

    }} /* endfor iyosh, iresn */
  /* for(ipart=myatm_start;ipart<=myatm_end;ipart++){
     printf("%f \t %f \n",therm_v1_nhc[1][therm_inhc_x[ipart]],therm_v2_nhc[1][therm_inhc_x[ipart]]);
     }getchar();*/


     printf("%.12g %.12 %.12g\n",
             therm_v1_nhc[1][therm_inhc_z[1]]
             therm_v1_nhc[2][therm_inhc_z[1]]
             therm_v1_nhc[3][therm_inhc_z[1]]); getchar();

   }/*end routine*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void init_NH_ISOK_par(GENERAL_DATA *general_data,CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
                      THERM_INFO *therm_info_class,THERM_POS *therm_class, 
                      INT_SCR *int_scr, int iflag_mass,
                      CLASS_COMM_FORC_PKG *class_comm_forc_pkg, VEL_SAMP_CLASS *vel_samp_class)
/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"
    int ipart,inhc,ichain;   /* Num: for loop counters */
    int natm_tot			   = clatoms_info->natm_tot;;
    int num_nhc,iii;
    double cst1,cst2,cst3,arg;
    double kin_tot,nhckin;
    int *therm_inhc_x		   = therm_info_class->inhc_x;
    int *therm_inhc_y		   = therm_info_class->inhc_y;
    int *therm_inhc_z		   = therm_info_class->inhc_z;
    int *iseed			       = &vel_samp_class->iseed;
    int *iseed2				   = &vel_samp_class->iseed2;
    double *qseed			   = &vel_samp_class->qseed;
    int len_nhc 			   = therm_info_class->len_nhc;	  /* Length of Nose-like chain */
    double lennhc	      	   = (double)len_nhc;		  /* For when I need LkT    */
    double L_L1				   = lennhc/(lennhc+1.0);
    double *class_clatoms_mass = clatoms_info->mass;
    double gamma     		   = 2.0;  /*Stochastic Term*/
    double therm_tau		   = 5.0;  /* Timescale for Thermostat Masses */
    double *int_scr_atm_kin    = int_scr->atm_kin;
    double **therm_v1_nhc      = therm_class->x_nhc;
    double **therm_v2_nhc      = therm_class->v_nhc;
    double **therm_mass_nhc    = therm_info_class->mass_nhc;
    double **therm_gkt         = therm_info_class->gkt;
    double *class_clatoms_vx   = clatoms_pos->vx;
    double *class_clatoms_vy   = clatoms_pos->vy;
    double *class_clatoms_vz   = clatoms_pos->vz;
    int myatm_start 		   = clatoms_info->myatm_start;
    int myatm_end 			   = clatoms_info->myatm_end;
    int mytherm_start 		   = therm_info_class->mytherm_start;
    int mytherm_end 		   = therm_info_class->mytherm_end;
    int num_nhc_share 		   = therm_info_class->num_nhc_share;
    int *map_share   		   = therm_info_class->map_share;
    int np_forc      		   = class_comm_forc_pkg->num_proc;
    MPI_Comm comm_forc  	   = class_comm_forc_pkg->comm;


/*--------------------------------------------------------------------------*/
/* I) Assign Thermostat Masses  (Mass[v_{1j}]=Mass[v_{2j}])*/

 /*   for(ipart=myatm_start;ipart<=myatm_end;ipart++){
    	for(ichain=1;ichain<=len_nhc;ichain++){
    	  therm_mass_nhc[ichain][therm_inhc_x[ipart]]=therm_gkt[ichain][therm_inhc_x[ipart]]
    	                                              *therm_tau*therm_tau;
    	  therm_mass_nhc[ichain][therm_inhc_y[ipart]]=therm_gkt[ichain][therm_inhc_y[ipart]]
    	                                              *therm_tau*therm_tau;
    	  therm_mass_nhc[ichain][therm_inhc_z[ipart]]=therm_gkt[ichain][therm_inhc_z[ipart]]
    	                                              *therm_tau*therm_tau;
    	 }
    }*/


/*--------------------------------------------------------------------------*/
/* II) Randomly Assign {v_{2j}} and {v_{2j}}  */

 /*   for(ipart=myatm_start;ipart<=myatm_end;ipart++){

    	  gaussran(len_nhc,iseed,iseed2,qseed,int_scr_atm_kin);
    	  for(ichain=1;ichain<=len_nhc;ichain++){
    		  therm_v2_nhc[ichain][therm_inhc_x[ipart]]
    		    =sqrt(therm_gkt[ichain][therm_inhc_x[ipart]]/therm_mass_nhc[ichain][therm_inhc_x[ipart]])
    		    *int_scr_atm_kin[ichain];
    	  }


    	  gaussran(len_nhc,iseed,iseed2,qseed,int_scr_atm_kin);
    	  for(ichain=1;ichain<=len_nhc;ichain++){
    		  therm_v2_nhc[ichain][therm_inhc_y[ipart]]
    		    =sqrt(therm_gkt[ichain][therm_inhc_y[ipart]]/therm_mass_nhc[ichain][therm_inhc_y[ipart]])
    		    *int_scr_atm_kin[ichain];
    	  }


    	  gaussran(len_nhc,iseed,iseed2,qseed,int_scr_atm_kin);
    	  for(ichain=1;ichain<=len_nhc;ichain++){
    		  therm_v2_nhc[ichain][therm_inhc_z[ipart]]
    		    =sqrt(therm_gkt[ichain][therm_inhc_z[ipart]]/therm_mass_nhc[ichain][therm_inhc_z[ipart]])
    		    *int_scr_atm_kin[ichain];
    	  }
    }


/*--------------------------------------------------------------------------*/
/* III) Assign {v and v_{1j}} to Satisfy Isokinetic Constraint  */
/*    for(ipart=myatm_start;ipart<=myatm_end;ipart++){
  	  gaussran(len_nhc,iseed,iseed2,qseed,int_scr_atm_kin);
  	  class_clatoms_vx[ipart]=int_scr_atm_kin[1];
  	  gaussran(len_nhc,iseed,iseed2,qseed,int_scr_atm_kin);
    	for(ichain=1;ichain<=len_nhc;ichain++){
    		therm_v1_nhc[ichain][therm_inhc_x[ipart]]=int_scr_atm_kin[ichain];
    	}
    	cst1=class_clatoms_vx[ipart]*class_clatoms_vx[ipart];
    	for(ichain=1;ichain<=len_nhc;ichain++){
    		cst1+=therm_v1_nhc[ichain][therm_inhc_x[ipart]]*therm_v1_nhc[ichain][therm_inhc_x[ipart]];
    	}
    	cst1=sqrt(cst1);
    	class_clatoms_vx[ipart]/=cst1;
    	for(ichain=1;ichain<=len_nhc;ichain++){
    		therm_v1_nhc[ichain][therm_inhc_x[ipart]]/=cst1;
    	}

    	class_clatoms_vx[ipart]/=sqrt(class_clatoms_mass[ipart]/(lennhc*therm_gkt[1][therm_inhc_x[ipart]]));
    	for(ichain=1;ichain<=len_nhc;ichain++){
    		therm_v1_nhc[ichain][therm_inhc_x[ipart]]/=sqrt(therm_mass_nhc[ichain][therm_inhc_x[ipart]]/((lennhc+1.0)*therm_gkt[1][therm_inhc_x[ipart]]));
    	}
    	cst1=class_clatoms_mass[ipart]*class_clatoms_vx[ipart]*class_clatoms_vx[ipart];
    	for(ichain=1;ichain<=len_nhc;ichain++){
    		cst1+=L_L1*therm_mass_nhc[ichain][therm_inhc_x[ipart]]*therm_v1_nhc[ichain][therm_inhc_x[ipart]]*therm_v1_nhc[ichain][therm_inhc_x[ipart]];
    	}



    	  gaussran(len_nhc,iseed,iseed2,qseed,int_scr_atm_kin);
    	  class_clatoms_vy[ipart]=int_scr_atm_kin[1];
    	  gaussran(len_nhc,iseed,iseed2,qseed,int_scr_atm_kin);
      	for(ichain=1;ichain<=len_nhc;ichain++){
      		therm_v1_nhc[ichain][therm_inhc_y[ipart]]=int_scr_atm_kin[ichain];
      	}
      	cst1=class_clatoms_vy[ipart]*class_clatoms_vy[ipart];
      	for(ichain=1;ichain<=len_nhc;ichain++){
      		cst1+=therm_v1_nhc[ichain][therm_inhc_y[ipart]]*therm_v1_nhc[ichain][therm_inhc_y[ipart]];
      	}
      	cst1=sqrt(cst1);
      	class_clatoms_vy[ipart]/=cst1;
      	for(ichain=1;ichain<=len_nhc;ichain++){
      		therm_v1_nhc[ichain][therm_inhc_y[ipart]]/=cst1;
      	}
      	class_clatoms_vy[ipart]/=sqrt(class_clatoms_mass[ipart]/(lennhc*therm_gkt[1][therm_inhc_y[ipart]]));
      	for(ichain=1;ichain<=len_nhc;ichain++){
      		therm_v1_nhc[ichain][therm_inhc_y[ipart]]/=sqrt(therm_mass_nhc[ichain][therm_inhc_y[ipart]]/((lennhc+1.0)*therm_gkt[1][therm_inhc_y[ipart]]));
      	}
      	cst1=class_clatoms_mass[ipart]*class_clatoms_vy[ipart]*class_clatoms_vy[ipart];
      	for(ichain=1;ichain<=len_nhc;ichain++){
      		cst1+=therm_mass_nhc[ichain][therm_inhc_y[ipart]]*therm_v1_nhc[ichain][therm_inhc_y[ipart]]*therm_v1_nhc[ichain][therm_inhc_y[ipart]];
      	}


    	  gaussran(len_nhc,iseed,iseed2,qseed,int_scr_atm_kin);
    	  class_clatoms_vz[ipart]=int_scr_atm_kin[1];
    	  gaussran(len_nhc,iseed,iseed2,qseed,int_scr_atm_kin);
      	for(ichain=1;ichain<=len_nhc;ichain++){
      		therm_v1_nhc[ichain][therm_inhc_z[ipart]]=int_scr_atm_kin[ichain];
      	}
      	cst1=class_clatoms_vz[ipart]*class_clatoms_vz[ipart];
      	for(ichain=1;ichain<=len_nhc;ichain++){
      		cst1+=therm_v1_nhc[ichain][therm_inhc_z[ipart]]*therm_v1_nhc[ichain][therm_inhc_z[ipart]];
      	}
      	cst1=sqrt(cst1);
      	class_clatoms_vz[ipart]/=cst1;
      	for(ichain=1;ichain<=len_nhc;ichain++){
      		therm_v1_nhc[ichain][therm_inhc_z[ipart]]/=cst1;
      	}
      	class_clatoms_vz[ipart]/=sqrt(class_clatoms_mass[ipart]/(lennhc*therm_gkt[1][therm_inhc_z[ipart]]));
      	for(ichain=1;ichain<=len_nhc;ichain++){
      		therm_v1_nhc[ichain][therm_inhc_z[ipart]]/=sqrt(therm_mass_nhc[ichain][therm_inhc_z[ipart]]/((lennhc+1.0)*therm_gkt[1][therm_inhc_z[ipart]]));
      	}
      	cst1=class_clatoms_mass[ipart]*class_clatoms_vx[ipart]*class_clatoms_vx[ipart];
      	for(ichain=1;ichain<=len_nhc;ichain++){
      		cst1+=therm_mass_nhc[ichain][therm_inhc_z[ipart]]*therm_v1_nhc[ichain][therm_inhc_z[ipart]]*therm_v1_nhc[ichain][therm_inhc_z[ipart]];
      	}


    }
*/

    for(ipart=1;ipart<=natm_tot;ipart++){
      kin_tot += class_clatoms_mass[ipart]*class_clatoms_vx[ipart]*class_clatoms_vx[ipart];
      kin_tot += class_clatoms_mass[ipart]*class_clatoms_vy[ipart]*class_clatoms_vy[ipart];
      kin_tot += class_clatoms_mass[ipart]*class_clatoms_vz[ipart]*class_clatoms_vz[ipart];
      for(ichain=1;ichain<=len_nhc;ichain++){
          nhckin += therm_mass_nhc[ichain][therm_inhc_x[ipart]]
                          *therm_v1_nhc[ichain][therm_inhc_x[ipart]]*therm_v1_nhc[ichain][therm_inhc_x[ipart]];
          nhckin += therm_mass_nhc[ichain][therm_inhc_y[ipart]]
                          *therm_v1_nhc[ichain][therm_inhc_y[ipart]]*therm_v1_nhc[ichain][therm_inhc_y[ipart]];
          nhckin += therm_mass_nhc[ichain][therm_inhc_z[ipart]]
                          *therm_v1_nhc[ichain][therm_inhc_z[ipart]]*therm_v1_nhc[ichain][therm_inhc_z[ipart]];
      }

    }


#define CHECK_CONSTRNT_off
#ifdef CHECK_CONSTRNT
    for(ipart=myatm_start;ipart<=myatm_end;ipart++){
    	cst1=0.0;	cst2=0.0;
    	for(ichain=1;ichain<=len_nhc;ichain++){
    		cst1+=therm_mass_nhc[ichain][therm_inhc_x[ipart]]
    		    		   *therm_v1_nhc[ichain][therm_inhc_x[ipart]]*therm_v1_nhc[ichain][therm_inhc_x[ipart]];
    	}
    	cst2+=(class_clatoms_mass[ipart]*class_clatoms_vx[ipart]*class_clatoms_vx[ipart]);
    	cst1*=lennhc/(lennhc+1.0);
    	cst3=(cst1+cst2)/therm_gkt[1][therm_inhc_x[ipart]];
    printf("%d \t %f \n",ipart,cst3);
    }getchar();
#endif

/*--------------------------------------------------------------------------*/
/*end routine*/}
/*==========================================================================*/










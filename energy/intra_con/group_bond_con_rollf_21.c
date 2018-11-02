/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                   Module: grp_bond_con.c                                 */
/*                                                                          */
/* This routine controls the atom based group constraint routines           */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../proto_defs/proto_intra_con_local.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#include "../proto_defs/proto_energy_ctrl_entry.h"

#define NCON_21 1
#define NAT_21 2


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void shake_21_rollf(GRP_BOND_CON *grp_bond_con,
                    CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
                    PTENS *ptens,double dt,double *aiter,
                    PAR_RAHMAN *par_rahman,int ifirst,CELL *cell,
                    CLASS_COMM_FORC_PKG *class_comm_forc_pkg)

/*==========================================================================*/
/*        Begin Routine                                                     */
   {/*Begin Routine*/
/*=======================================================================*/
/*         Local Variable declarations                                   */
  
#include "../typ_defs/typ_mask.h"

 int i,iter,igrp,*ind1,*ind2,jtyp;
 int iii,ktemp;

 double r12sq,r12osq,r12_r12o;
 double rmu12;
 double dts;

 double a,b,c,dd; /* coefficients for quadratic equation */
 double dxn12,dyn12,dzn12;
 double pnorm;

 double *xlam1;
 double *rm1,*rm2;
 double *dx12,*dy12,*dz12;
 double *dxo12,*dyo12,*dzo12;
 double *dxo120,*dyo120,*dzo120;
 double *x1,*x2,*y1,*y2,*z1,*z2,*xo1,*xo2,*yo1,*yo2,*zo1,*zo2;
 double *dij1;
 double *p11,*p12,*p13,*p23,*p33,*p22;

/* Local pointers */

  double *clatoms_mass         = clatoms_info->mass;
  double *clatoms_xold         = clatoms_info->xold;
  double *clatoms_yold         = clatoms_info->yold;
  double *clatoms_zold         = clatoms_info->zold;
  double *clatoms_x            = clatoms_pos->x;
  double *clatoms_y            = clatoms_pos->y;
  double *clatoms_z            = clatoms_pos->z;
  double *clatoms_vx           = clatoms_pos->vx;
  double *clatoms_vy           = clatoms_pos->vy;
  double *clatoms_vz           = clatoms_pos->vz;

  int ngrp                     = grp_bond_con->num_21;
  int *grp_bond_con_j1_21      = grp_bond_con->j1_21;
  int *grp_bond_con_j2_21      = grp_bond_con->j2_21;
  int *grp_bond_con_jtyp_21    = grp_bond_con->jtyp_21;
  double **grp_bond_con_eq_21  = grp_bond_con->eq_21;
  double **grp_bond_con_al_21  = grp_bond_con->al_21;

  double *pvten_inc            = ptens->pvten_inc;
  double *pvten_tmp            = ptens->pvten_tmp;
  double *pvten_tmp2           = ptens->pvten_tmp_res;

  double *roll_mtv             = par_rahman->roll_mtv;
  double *fgmat_p              = par_rahman->fgmat_p;
  double *vgmat                = par_rahman->vgmat;
  double roll_scg              = par_rahman->roll_scg;
  double mass_hm               = par_rahman->mass_hm;

  int iperd                    = cell->iperd;
  int hmat_cons_typ            = cell->hmat_cons_typ;
  int hmat_int_typ             = cell->hmat_int_typ;

  int np_forc                  = class_comm_forc_pkg->num_proc;
  int myid_forc                = class_comm_forc_pkg->myid;
  MPI_Comm comm_forc           = class_comm_forc_pkg->comm;

/*=======================================================================*/

  if(ngrp > 0){
     xlam1= dvector(1,ngrp);

     rm1 = dvector(1,ngrp);
     rm2 = dvector(1,ngrp);
 
     dx12 = dvector(1,ngrp);
     dy12 = dvector(1,ngrp);
     dz12 = dvector(1,ngrp);

     dxo12 = dvector(1,ngrp);
     dyo12 = dvector(1,ngrp);
     dzo12 = dvector(1,ngrp);

     dxo120 = dvector(1,ngrp);
     dyo120 = dvector(1,ngrp);
     dzo120 = dvector(1,ngrp);

      dij1= dvector(1,ngrp);

      x1= dvector(1,ngrp);
      y1= dvector(1,ngrp);
      z1= dvector(1,ngrp);
      x2= dvector(1,ngrp);
      y2= dvector(1,ngrp);
      z2= dvector(1,ngrp);

      xo1= dvector(1,ngrp);
      yo1= dvector(1,ngrp);
      zo1= dvector(1,ngrp);
      xo2= dvector(1,ngrp);
      yo2= dvector(1,ngrp);
      zo2= dvector(1,ngrp);
       
       p11= dvector(1,ngrp);
       p12= dvector(1,ngrp);
       p13= dvector(1,ngrp);
       p22= dvector(1,ngrp);
       p23= dvector(1,ngrp);
       p33= dvector(1,ngrp);

      ind1 = (int *) calloc((ngrp+1),sizeof(int));
      ind2 = (int *) calloc((ngrp+1),sizeof(int));
  }/*endif*/
/*=======================================================================*/
/*=======================================================================*/

 dts = dt*dt;
 pnorm = 2.0/dts;
 *aiter = 0.0;
  pvten_tmp[1] = 0.0;
  pvten_tmp[2] = 0.0;
  pvten_tmp[3] = 0.0;
  pvten_tmp[4] = 0.0;
  pvten_tmp[5] = 0.0;
  pvten_tmp[6] = 0.0;
  pvten_tmp[7] = 0.0;
  pvten_tmp[8] = 0.0;
  pvten_tmp[9] = 0.0;

/* assign positions and masses */

if(ifirst == 2){
   for(igrp=1;igrp <= ngrp ; igrp++) {
     grp_bond_con_al_21[1][igrp] = 0.0;
   }
 }


 for(igrp=1;igrp <= ngrp ; igrp++) {
   ind1[igrp] = grp_bond_con_j1_21[igrp];
   ind2[igrp] = grp_bond_con_j2_21[igrp];
 }

 for(igrp=1;igrp <= ngrp ; igrp++) {
   ktemp = ind1[igrp];
   x1[igrp] = clatoms_x[ktemp];
   y1[igrp] = clatoms_y[ktemp];
   z1[igrp] = clatoms_z[ktemp];
 }

 for(igrp=1;igrp <= ngrp ; igrp++) {
   ktemp = ind2[igrp];
   x2[igrp] = clatoms_x[ktemp];
   y2[igrp] = clatoms_y[ktemp];
   z2[igrp] = clatoms_z[ktemp];
 }


 for(igrp=1;igrp <= ngrp ; igrp++) {
   ktemp = ind1[igrp];
   xo1[igrp] = clatoms_xold[ktemp];
   yo1[igrp] = clatoms_yold[ktemp];
   zo1[igrp] = clatoms_zold[ktemp];
   rm1[igrp] = 1.0/clatoms_mass[ktemp];
 }

 for(igrp=1;igrp <= ngrp ; igrp++) {
   ktemp = ind2[igrp];
   xo2[igrp] = clatoms_xold[ktemp];
   yo2[igrp] = clatoms_yold[ktemp];
   zo2[igrp] = clatoms_zold[ktemp];
   rm2[igrp] = 1.0/clatoms_mass[ktemp];
 }


 for(igrp=1;igrp <= ngrp ; igrp++) {
   jtyp = grp_bond_con_jtyp_21[igrp];
   dij1[igrp] = grp_bond_con_eq_21[1][jtyp];
 }

 for(igrp=1;igrp <= ngrp ; igrp++) {
    dx12[igrp] = x1[igrp]-x2[igrp];
    dy12[igrp] = y1[igrp]-y2[igrp];
    dz12[igrp] = z1[igrp]-z2[igrp];
  }

 for(igrp=1;igrp <= ngrp ; igrp++) {
    dxo120[igrp] = xo1[igrp]-xo2[igrp];
    dyo120[igrp] = yo1[igrp]-yo2[igrp];
    dzo120[igrp] = zo1[igrp]-zo2[igrp];
  }

 for(igrp=1;igrp <= ngrp ; igrp++) {
  dxo12[igrp] = (dxo120[igrp])*roll_mtv[1]
              + (dyo120[igrp])*roll_mtv[2]
              + (dzo120[igrp])*roll_mtv[3];
  dyo12[igrp] = (dxo120[igrp])*roll_mtv[4]
              + (dyo120[igrp])*roll_mtv[5]
              + (dzo120[igrp])*roll_mtv[6];
  dzo12[igrp] = (dxo120[igrp])*roll_mtv[7]
              + (dyo120[igrp])*roll_mtv[8]
              + (dzo120[igrp])*roll_mtv[9];
 }

/*======================================================================== */
/* determine value of multiplier using quadratic equation */
/* note that only the positive root has physical significance */

if(ifirst==2||ifirst==0){
 for(igrp=1;igrp <= ngrp ; igrp++) {

    rmu12= rm1[igrp] + rm2[igrp];
    r12osq = dxo12[igrp]*dxo12[igrp] + dyo12[igrp]*dyo12[igrp] 
           + dzo12[igrp]*dzo12[igrp];
    r12sq = dx12[igrp]*dx12[igrp] + dy12[igrp]*dy12[igrp] 
          + dz12[igrp]*dz12[igrp];
    r12_r12o = dx12[igrp]*dxo12[igrp]+dy12[igrp]*dyo12[igrp]
             + dz12[igrp]*dzo12[igrp];

    a= r12osq*rmu12*rmu12;
    b= 2.0*r12_r12o*rmu12;
    c= r12sq - dij1[igrp]*dij1[igrp]; 
    dd= b*b - 4.0*a*c;

   xlam1[igrp] = b - sqrt(dd);
   xlam1[igrp] /= (2.0*a);

 }
}else{
  for(igrp=1;igrp <= ngrp ; igrp++) {
    xlam1[igrp] = grp_bond_con_al_21[1][igrp];
    grp_bond_con_al_21[1][igrp] = 0.0;
  }
}

/* position update */
/* scaled distance differences */
#ifndef NO_PRAGMA
#pragma IVDEP
#endif
 for(igrp=1;igrp <= ngrp ; igrp++) {
  ktemp = ind1[igrp];
  clatoms_x[ktemp] -= (xlam1[igrp]*dxo12[igrp])*rm1[igrp];
  clatoms_y[ktemp] -= (xlam1[igrp]*dyo12[igrp])*rm1[igrp];
  clatoms_z[ktemp] -= (xlam1[igrp]*dzo12[igrp])*rm1[igrp];
 }

#ifndef NO_PRAGMA
#pragma IVDEP
#endif
 for(igrp=1;igrp <= ngrp ; igrp++) {
  ktemp = ind2[igrp];
  clatoms_x[ktemp] += (xlam1[igrp]*dxo12[igrp]*rm2[igrp]);
  clatoms_y[ktemp] += (xlam1[igrp]*dyo12[igrp]*rm2[igrp]);
  clatoms_z[ktemp] += (xlam1[igrp]*dzo12[igrp]*rm2[igrp]);
 }

/* Velocity update */
/* unscaled differences in positions */

#ifndef NO_PRAGMA
#pragma IVDEP
#endif
 for(igrp=1;igrp <= ngrp ; igrp++) {
  ktemp = ind1[igrp];
  clatoms_vx[ktemp]-=(xlam1[igrp]*dxo120[igrp])*rm1[igrp]/dt;
  clatoms_vy[ktemp]-=(xlam1[igrp]*dyo120[igrp])*rm1[igrp]/dt;
  clatoms_vz[ktemp]-=(xlam1[igrp]*dzo120[igrp])*rm1[igrp]/dt;
 }

#ifndef NO_PRAGMA
#pragma IVDEP
#endif
 for(igrp=1;igrp <= ngrp ; igrp++) {
  ktemp = ind2[igrp];
  clatoms_vx[ktemp] += (xlam1[igrp]*dxo120[igrp]*rm2[igrp])/dt;
  clatoms_vy[ktemp] += (xlam1[igrp]*dyo120[igrp]*rm2[igrp])/dt;
  clatoms_vz[ktemp] += (xlam1[igrp]*dzo120[igrp]*rm2[igrp])/dt;
 }


/* Pressure tensor update */
/* Compute difference vectors  unscaled old distances */
 for(igrp=1;igrp <= ngrp ; igrp++) {
   dxn12 = dxo120[igrp];
   dyn12 = dyo120[igrp];
   dzn12 = dzo120[igrp];

   p11[igrp] = xlam1[igrp]*dxn12*dxn12;
   p22[igrp] = xlam1[igrp]*dyn12*dyn12;
   p33[igrp] = xlam1[igrp]*dzn12*dzn12;
   p12[igrp] = xlam1[igrp]*dxn12*dyn12;
   p13[igrp] = xlam1[igrp]*dxn12*dzn12;
   p23[igrp] = xlam1[igrp]*dyn12*dzn12;
 }/*end for*/

#ifndef NO_PRAGMA
#pragma IVDEP
#endif
 for(igrp=1;igrp <= ngrp ; igrp++) {
    pvten_tmp[1] -= (p11[igrp]*pnorm); 
    pvten_tmp[2] -= (p12[igrp]*pnorm);
    pvten_tmp[3] -= (p13[igrp]*pnorm);
    pvten_tmp[4] -= (p12[igrp]*pnorm);
    pvten_tmp[5] -= (p22[igrp]*pnorm);
    pvten_tmp[6] -= (p23[igrp]*pnorm);
    pvten_tmp[7] -= (p13[igrp]*pnorm);
    pvten_tmp[8] -= (p23[igrp]*pnorm);
    pvten_tmp[9] -= (p33[igrp]*pnorm);
 }/*end for*/


/* Save multiplier */

 for(igrp=1;igrp <= ngrp; igrp++) {
   grp_bond_con_al_21[1][igrp] += xlam1[igrp];
 } /* end for igrp */
    

/*=======================================================================*/
/*  IV)Allreduce pvten_tmp     */

  if(np_forc > 1 ){
   for(i=1;i<=9;i++){
    pvten_tmp2[i] = pvten_tmp[i];
   }/*endfor*/
   Allreduce(&(pvten_tmp2[1]), &(pvten_tmp[1]),9,MPI_DOUBLE,
                   MPI_SUM,0,comm_forc);
  }/*endif*/

  pvten_inc[1] += pvten_tmp[1];
  pvten_inc[2] += pvten_tmp[2];
  pvten_inc[3] += pvten_tmp[3];
  pvten_inc[4] += pvten_tmp[4];
  pvten_inc[5] += pvten_tmp[5];
  pvten_inc[6] += pvten_tmp[6];
  pvten_inc[7] += pvten_tmp[7];
  pvten_inc[8] += pvten_tmp[8];
  pvten_inc[9] += pvten_tmp[9];

 if(ifirst == 0){
   constr_cell_mat(iperd,hmat_cons_typ,hmat_int_typ,pvten_tmp);
   for(i=1;i<=9;i++){
     fgmat_p[i] += pvten_tmp[i];
     vgmat[i]   += (pvten_tmp[i]*roll_scg*0.5*dt/mass_hm);
   }/*endfor*/
 }

/* free locally assigned memory */
 if(ngrp > 0){
     free_dvector(xlam1,1,ngrp);

     free_dvector(rm1,1,ngrp);
     free_dvector(rm2,1,ngrp);

     free_dvector(dx12,1,ngrp);
     free_dvector(dy12,1,ngrp);
     free_dvector(dz12,1,ngrp);

     free_dvector(dxo12,1,ngrp);
     free_dvector(dyo12,1,ngrp);
     free_dvector(dzo12,1,ngrp);

     free_dvector(dxo120,1,ngrp);
     free_dvector(dyo120,1,ngrp);
     free_dvector(dzo120,1,ngrp);

    free_dvector(dij1,1,ngrp);

    free_dvector(x1,1,ngrp);
    free_dvector(y1,1,ngrp);
    free_dvector(z1,1,ngrp);

    free_dvector(x2,1,ngrp);
    free_dvector(y2,1,ngrp);
    free_dvector(z2,1,ngrp);

    free_dvector(xo1,1,ngrp);
    free_dvector(yo1,1,ngrp);
    free_dvector(zo1,1,ngrp);

    free_dvector(xo2,1,ngrp);
    free_dvector(yo2,1,ngrp);
    free_dvector(zo2,1,ngrp);

    free_dvector(p11,1,ngrp);
    free_dvector(p12,1,ngrp);
    free_dvector(p13,1,ngrp);
    free_dvector(p22,1,ngrp);
    free_dvector(p23,1,ngrp);
    free_dvector(p33,1,ngrp);

    free(ind1);
    free(ind2);
 }/*endif*/

/*=======================================================================*/
   } /* end routine */
/*=======================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void rattle_21_rollf(GRP_BOND_CON *grp_bond_con,
                     CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
                     PTENS *ptens,double dt,
                     PAR_RAHMAN *par_rahman,int ifirst,CELL *cell,
                     CLASS_COMM_FORC_PKG *class_comm_forc_pkg)

/*==========================================================================*/
/*        Begin Routine                                                     */
    {/* Begin routine */
/*=======================================================================*/
/*         Local Variable declarations                                   */

#include "../typ_defs/typ_mask.h"

  int igrp,*ind1,*ind2,jtyp;
  int i,n,ktemp;

  double pnorm;
  double roll_mtvvi[10],rolli_by_vg[10];
  double det_roll_mtvv;

  double roll_sci;

  double r12sq,r12_v12,rmu12;
  double *x1,*x2,*y1,*y2,*z1,*z2;
  double *vx1,*vx2,*vy1,*vy2,*vz1,*vz2;
  double *p11,*p12,*p13,*p22,*p23,*p33;
  double *rm1,*rm2;
  double *dx12,*dy12,*dz12;
  double *dvx12,*dvy12,*dvz12;
  double *xlam1;

/* Local pointers */

  double *clatoms_mass         = clatoms_info->mass;
  double *clatoms_roll_sc      = clatoms_info->roll_sc;
  double *clatoms_x            = clatoms_pos->x;
  double *clatoms_y            = clatoms_pos->y;
  double *clatoms_z            = clatoms_pos->z;
  double *clatoms_vx           = clatoms_pos->vx;
  double *clatoms_vy           = clatoms_pos->vy;
  double *clatoms_vz           = clatoms_pos->vz;

  int ngrp                     = grp_bond_con->num_21;
  int *grp_bond_con_j1_21      = grp_bond_con->j1_21;
  int *grp_bond_con_j2_21      = grp_bond_con->j2_21;
  int *grp_bond_con_jtyp_13    = grp_bond_con->jtyp_21;
  double **grp_bond_con_al_21  = grp_bond_con->al_21;

  double *pvten_inc            = ptens->pvten_inc;
  double *pvten_tmp            = ptens->pvten_tmp;
  double *pvten_tmp2           = ptens->pvten_tmp_res;

  double *roll_mtvv            = par_rahman->roll_mtvv;
  double *vgmat_g              = par_rahman->vgmat_g;
  double *fgmat_p              = par_rahman->fgmat_p;
  double roll_scg              = par_rahman->roll_scg;
  double mass_hm               = par_rahman->mass_hm;

  int iperd                    = cell->iperd;
  int hmat_cons_typ            = cell->hmat_cons_typ;
  int hmat_int_typ             = cell->hmat_int_typ;

  int np_forc                  = class_comm_forc_pkg->num_proc;
  int myid_forc                = class_comm_forc_pkg->myid;
  MPI_Comm comm_forc           = class_comm_forc_pkg->comm;

/*=======================================================================*/
/* assign local memory */

 if(ngrp > 0){
    x1= dvector(1,ngrp);
    y1= dvector(1,ngrp);
    z1= dvector(1,ngrp);
    x2= dvector(1,ngrp);
    y2= dvector(1,ngrp);
    z2= dvector(1,ngrp);

    vx1= dvector(1,ngrp);
    vy1= dvector(1,ngrp);
    vz1= dvector(1,ngrp);
    vx2= dvector(1,ngrp);
    vy2= dvector(1,ngrp);
    vz2= dvector(1,ngrp);

  p11= dvector(1,ngrp);
  p12= dvector(1,ngrp);
  p13= dvector(1,ngrp);
  p22= dvector(1,ngrp);
  p23= dvector(1,ngrp);
  p33= dvector(1,ngrp);

  rm1= dvector(1,ngrp);
  rm2= dvector(1,ngrp);

  dx12= dvector(1,ngrp);
  dy12= dvector(1,ngrp);
  dz12= dvector(1,ngrp);

  dvx12= dvector(1,ngrp);
  dvy12= dvector(1,ngrp);
  dvz12= dvector(1,ngrp);

 xlam1= dvector(1,ngrp);
 ind1 = (int *) calloc((ngrp+1),sizeof(int));
 ind2 = (int *) calloc((ngrp+1),sizeof(int));
 }/*endif*/

/*=======================================================================*/
/*=======================================================================*/

 pnorm = 2.0/dt;
 n=3;
 gethinv(roll_mtvv,roll_mtvvi,&det_roll_mtvv,n);
 matmul_tt(roll_mtvvi,vgmat_g,rolli_by_vg,n);

 pvten_tmp[1] = 0.0;
 pvten_tmp[2] = 0.0;
 pvten_tmp[3] = 0.0;
 pvten_tmp[4] = 0.0;
 pvten_tmp[5] = 0.0;
 pvten_tmp[6] = 0.0;
 pvten_tmp[7] = 0.0;
 pvten_tmp[8] = 0.0;
 pvten_tmp[9] = 0.0;


 for(igrp=1;igrp <= ngrp ; igrp++) {
 
  ind1[igrp] = grp_bond_con_j1_21[igrp];
  ind2[igrp] = grp_bond_con_j2_21[igrp];
}

 for(igrp=1;igrp <= ngrp ; igrp++) {
  rm1[igrp] = 1.0/clatoms_mass[ind1[igrp]];
  rm2[igrp] = 1.0/clatoms_mass[ind2[igrp]];
}

 for(igrp=1;igrp <= ngrp ; igrp++) {
  ktemp = ind1[igrp];
  x1[igrp] = clatoms_x[ktemp];
  y1[igrp] = clatoms_y[ktemp];
  z1[igrp] = clatoms_z[ktemp];
} /*end for*/

 for(igrp=1;igrp <= ngrp ; igrp++) {
  ktemp = ind2[igrp];
  x2[igrp] = clatoms_x[ktemp];
  y2[igrp] = clatoms_y[ktemp];
  z2[igrp] = clatoms_z[ktemp];
} /*end for*/


 for(igrp=1;igrp <= ngrp ; igrp++) {
  ktemp = ind1[igrp];
  roll_sci=1.0/clatoms_roll_sc[ktemp];/*all roll scales the same in same cons*/
  vx1[igrp] = clatoms_vx[ktemp]
        + (x1[igrp]*rolli_by_vg[1] +y1[igrp]*rolli_by_vg[2]
         + z1[igrp]*rolli_by_vg[3]) *roll_sci;
  vy1[igrp] = clatoms_vy[ktemp]
        + (x1[igrp]*rolli_by_vg[4] +y1[igrp]*rolli_by_vg[5]
         + z1[igrp]*rolli_by_vg[6]) *roll_sci;
  vz1[igrp] = clatoms_vz[ktemp]
        + (x1[igrp]*rolli_by_vg[7] +y1[igrp]*rolli_by_vg[8]
         + z1[igrp]*rolli_by_vg[9]) *roll_sci;
} 

 for(igrp=1;igrp <= ngrp ; igrp++) {
  ktemp = ind2[igrp];
  roll_sci=1.0/clatoms_roll_sc[ktemp];/*all roll scales the same in same cons*/
  vx2[igrp] = clatoms_vx[ktemp]
        + (x2[igrp]*rolli_by_vg[1] +y2[igrp]*rolli_by_vg[2]
         + z2[igrp]*rolli_by_vg[3]) *roll_sci;
  vy2[igrp] = clatoms_vy[ktemp]
        + (x2[igrp]*rolli_by_vg[4] +y2[igrp]*rolli_by_vg[5]
         + z2[igrp]*rolli_by_vg[6]) *roll_sci;
  vz2[igrp] = clatoms_vz[ktemp]
        + (x2[igrp]*rolli_by_vg[7] +y2[igrp]*rolli_by_vg[8]
         + z2[igrp]*rolli_by_vg[9]) *roll_sci;

}


/* Define useful constants */
 for(igrp=1;igrp <= ngrp ; igrp++) {
  
    dx12[igrp] = x1[igrp]-x2[igrp];
    dy12[igrp] = y1[igrp]-y2[igrp];
    dz12[igrp] = z1[igrp]-z2[igrp];

    dvx12[igrp] = vx1[igrp]-vx2[igrp];
    dvy12[igrp] = vy1[igrp]-vy2[igrp];
    dvz12[igrp] = vz1[igrp]-vz2[igrp];
 }


/* compute lagrange multiplier */
 for(igrp=1;igrp <= ngrp ; igrp++) {
    r12sq = dx12[igrp]*dx12[igrp] + dy12[igrp]*dy12[igrp] + dz12[igrp]*dz12[igrp];
    r12_v12 = dx12[igrp]*dvx12[igrp] + dy12[igrp]*dvy12[igrp] + dz12[igrp]*dvz12[igrp];
    rmu12= rm1[igrp] + rm2[igrp];
   
    xlam1[igrp]= r12_v12/(r12sq*rmu12);
    
 }

#ifndef NO_PRAGMA
#pragma IVDEP
#endif
 for(igrp=1;igrp <= ngrp ; igrp++) {
   double xlam_1;
   double dx_12;
   double dy_12;
   double dz_12;
   int ktemp1,ktemp2;

    xlam_1= xlam1[igrp];

    dx_12= dx12[igrp]; 
    dy_12= dy12[igrp];
    dz_12= dz12[igrp];

    ktemp1=ind1[igrp]; 
    ktemp2=ind2[igrp]; 

   clatoms_vx[ktemp1] -= (xlam_1*dx_12)*rm1[igrp];
   clatoms_vy[ktemp1] -= (xlam_1*dy_12)*rm1[igrp];
   clatoms_vz[ktemp1] -= (xlam_1*dz_12)*rm1[igrp];

   clatoms_vx[ktemp2] += xlam_1*dx_12*rm2[igrp];
   clatoms_vy[ktemp2] += xlam_1*dy_12*rm2[igrp];
   clatoms_vz[ktemp2] += xlam_1*dz_12*rm2[igrp];


/* Pressure Tensor update */

    p11[igrp] = xlam_1*dx_12*dx_12;
    p22[igrp] = xlam_1*dy_12*dy_12;
    p33[igrp] = xlam_1*dz_12*dz_12;
    p12[igrp] = xlam_1*dx_12*dy_12;
    p13[igrp] = xlam_1*dx_12*dz_12;
    p23[igrp] = xlam_1*dy_12*dz_12;
}/*end for*/

#ifndef NO_PRAGMA
#pragma IVDEP
#endif
 for(igrp=1;igrp <= ngrp ; igrp++) {
    pvten_tmp[1] -= (p11[igrp]*pnorm); 
    pvten_tmp[2] -= (p12[igrp]*pnorm);
    pvten_tmp[3] -= (p13[igrp]*pnorm);
    pvten_tmp[4] -= (p12[igrp]*pnorm);
    pvten_tmp[5] -= (p22[igrp]*pnorm);
    pvten_tmp[6] -= (p23[igrp]*pnorm);
    pvten_tmp[7] -= (p13[igrp]*pnorm);
    pvten_tmp[8] -= (p23[igrp]*pnorm);
    pvten_tmp[9] -= (p33[igrp]*pnorm);

 }/* end for igrp */

/* Save multiplier */
 for(igrp=1;igrp <= ngrp; igrp++) {
  grp_bond_con_al_21[1][igrp] = xlam1[igrp];
 }

/*=======================================================================*/
/*  IV)Allreduce pvten_tmp     */

  if(np_forc > 1 ){
   for(i=1;i<=9;i++){
    pvten_tmp2[i] = pvten_tmp[i];
   }/*endfor*/
   Allreduce(&(pvten_tmp2[1]), &(pvten_tmp[1]),9,MPI_DOUBLE,
                   MPI_SUM,0,comm_forc);
  }/*endif*/

  pvten_inc[1] += pvten_tmp[1];
  pvten_inc[2] += pvten_tmp[2];
  pvten_inc[3] += pvten_tmp[3];
  pvten_inc[4] += pvten_tmp[4];
  pvten_inc[5] += pvten_tmp[5];
  pvten_inc[6] += pvten_tmp[6];
  pvten_inc[7] += pvten_tmp[7];
  pvten_inc[8] += pvten_tmp[8];
  pvten_inc[9] += pvten_tmp[9];

 if(ifirst == 0){
   constr_cell_mat(iperd,hmat_cons_typ,hmat_int_typ,pvten_tmp);
   for(i=1;i<=9;i++){
     fgmat_p[i] += pvten_tmp[i];
     vgmat_g[i] += (pvten_tmp[i]*roll_scg*0.5*dt/mass_hm);
  }/*endfor*/
 }/*endif*/

/* free locally assigned memory */  
 if(ngrp > 0){
   free_dvector(x1,1,ngrp);
   free_dvector(y1,1,ngrp);
   free_dvector(z1,1,ngrp);

   free_dvector(x2,1,ngrp);
   free_dvector(y2,1,ngrp);
   free_dvector(z2,1,ngrp);

   free_dvector(vx1,1,ngrp);
   free_dvector(vy1,1,ngrp);
   free_dvector(vz1,1,ngrp);

   free_dvector(vx2,1,ngrp);
   free_dvector(vy2,1,ngrp);
   free_dvector(vz2,1,ngrp);

   free_dvector(p11,1,ngrp);
   free_dvector(p12,1,ngrp);
   free_dvector(p13,1,ngrp);
   free_dvector(p22,1,ngrp);
   free_dvector(p23,1,ngrp);
   free_dvector(p33,1,ngrp);

   free_dvector(rm1,1,ngrp);
   free_dvector(rm2,1,ngrp);

   free_dvector(dx12,1,ngrp); 
   free_dvector(dy12,1,ngrp); 
   free_dvector(dz12,1,ngrp); 

   free_dvector(dvx12,1,ngrp); 
   free_dvector(dvy12,1,ngrp); 
   free_dvector(dvz12,1,ngrp); 

   free_dvector(xlam1,1,ngrp);

   free(ind1);
   free(ind2);
 }/*endif*/

/*=======================================================================*/
   } /* end routine */
/*=======================================================================*/





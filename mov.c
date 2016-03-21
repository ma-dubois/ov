/********************************************************************************
*      Copyright (C) 1995 Greg Landrum
*
*  This file is part of yaehmop.
*
*   This is free software.
* 
*  Permission is granted to modify, or otherwise fold, spindle, and mutilate this
*    code provided all copyright notices are left intact.
*
*  This code may be distributed to your heart's content, in whatever form,
*    provided no fee is charged for the distribution, all copyright notices are
*    left intact, and the source is distributed (without fee) along with any
*    binaries to anyone who requests it.
*
*  There are, of course, no warranties at all on this program.
*
********************************************************************************/



/********************************************************************************
*
*     this file contains stuff for evaluating overlaps
*
*  created:  greg landrum  August 1993
*
********************************************************************************/

/***
  Edit History:

  March '98: WG
  - f orbital overlaps added

***/
/* #include "bind.h" */

#define real double
#define abfns abfns_
#define lovlap lovlap_

/********************************************************************************
*
*                   Procedure mov
*
* Arguments: sigma,pi,delta,phi: pointers to reals
*                 which1,which2: int
*                          dist: real
*        q_num1,q_num2,l1,l2,nn: int
*                         atoms: atom_type
*
*
* Returns: none
*
* Action:  does whatever MOV did in the original program
*
*   comments will follow the clue when I get one
*
********************************************************************************/
void mov(sigma,pi,delta,phi,dist,q_num1,q_num2,
	 l1,l2,exp11,exp12,coeff11,coeff12,exp21,exp22,coeff21,coeff22)
  real *sigma,*pi,*delta,*phi,dist;
  real exp11,exp12,coeff11,coeff12,exp21,exp22,coeff21,coeff22;
  int q_num1,q_num2,l1,l2;
{
  int i,j,num_zeta1,num_zeta2;
  real coeff_1,coeff_2,sk1,sk2,r;
  
  real A_fn_values[30], B_fn_values[30];
  real ang_ind_overlap[4];
  
  int loopvar,max, m=0, nn;  /* definitions for abfns.f and lovlap.f */
  max = q_num1 + q_num2;

  if(l1>l2){
    nn=l2;
  }
  else{
    nn=l1;
  }

  /* initialize the components of the overlap to zero */
  *sigma=*pi=*delta=*phi=0.0;
  ang_ind_overlap[0]=ang_ind_overlap[1]=ang_ind_overlap[2]=ang_ind_overlap[3]=0.0;
  
  /* figure out whether or not we are using double zeta f'ns */

  if( coeff12 >0 ) num_zeta1 = 2;
  else num_zeta1 = 1;
  if( coeff22 >0 ) num_zeta2 = 2;
  else num_zeta2 = 1;

  /* get the exponents */

    sk1 = exp11;
    sk2 = exp21;

  /* deal with zeta1 - zeta1 overlap */

    coeff_1 = coeff11;
    coeff_2 = coeff21;

  /* call the routine to evaluate the A & B functions */

  abfns(A_fn_values,B_fn_values,&sk1,&sk2,&dist,&l1,&l2,&m,&q_num1,&q_num2,&max);

  /* test print of A and B functions */
  /* 
    printf("zeta1[%f]-zeta1[%f] call ...\n",sk1,sk2);
    loopvar=0;
    while(loopvar <= max)
      {
	printf("A(%d)= %.16f \t B(%d)= %f\n",loopvar,A_fn_values[loopvar],loopvar,B_fn_values[loopvar]);
	loopvar++;
      }
  */
  
  /* call the routine to evaluate sigma,pi,... overlaps in the local ref frame */
  
  for(i=0;i<=nn;i++)
    {
      m=i;
      lovlap(&(ang_ind_overlap[i]),A_fn_values,B_fn_values,&sk1,&sk2,&dist,&l1,&l2,&m,&q_num1,&q_num2,&max);
    }
  
  /* add in the contributions we have found thus far */
  
  *sigma += coeff_1*coeff_2*ang_ind_overlap[0];
  *pi += coeff_1*coeff_2*ang_ind_overlap[1];
  *delta += coeff_1*coeff_2*ang_ind_overlap[2];
  *phi += coeff_1*coeff_2*ang_ind_overlap[3];
  
  /* now do zeta1 - zeta2 overlap if applicable */
  
  if( num_zeta2 == 2){
      sk2 = exp22;
      coeff_2 = coeff22;

    abfns(A_fn_values,B_fn_values,&sk1,&sk2,&dist,&l1,&l2,&m,&q_num1,&q_num2,&max);

    for(i=0;i<=nn;i++){
      m=i;
      lovlap(&(ang_ind_overlap[i]),A_fn_values,B_fn_values,&sk1,&sk2,&dist,&l1,&l2,&m,&q_num1,&q_num2,&max);
    }
    
    *sigma += coeff_1*coeff_2*ang_ind_overlap[0];
    *pi += coeff_1*coeff_2*ang_ind_overlap[1];
    *delta += coeff_1*coeff_2*ang_ind_overlap[2];
    *phi += coeff_1*coeff_2*ang_ind_overlap[3];
  }
  
  /* now do zeta2 - zeta2 */
  
  if( num_zeta1 == 2 && num_zeta2 == 2 ){
      sk1 = exp12;
      coeff_1 = coeff12;
    
    abfns(A_fn_values,B_fn_values,&sk1,&sk2,&dist,&l1,&l2,&m,&q_num1,&q_num2,&max);
    
    for(i=0;i<=nn;i++){
      m=i;
      lovlap(&(ang_ind_overlap[i]),A_fn_values,B_fn_values,&sk1,&sk2,&dist,&l1,&l2,&m,&q_num1,&q_num2,&max);
    }
    
    *sigma += coeff_1*coeff_2*ang_ind_overlap[0];
    *pi += coeff_1*coeff_2*ang_ind_overlap[1];
    *delta += coeff_1*coeff_2*ang_ind_overlap[2];
    *phi += coeff_1*coeff_2*ang_ind_overlap[3];
  }
    
  /* finally do zeta2 - zeta1 */
    
  if( num_zeta1 == 2 ){
	   sk2 = exp21;
	   coeff_2 = coeff21;
	   sk1 = exp12;
	   coeff_1 = coeff12;
      
     abfns(A_fn_values,B_fn_values,&sk1,&sk2,&dist,&l1,&l2,&m,&q_num1,&q_num2,&max);
      
     for(i=0;i<=nn;i++){
	   m=i;
	   lovlap(&(ang_ind_overlap[i]),A_fn_values,B_fn_values,&sk1,&sk2,&dist,&l1,&l2,&m,&q_num1,&q_num2,&max);
     }
     
     *sigma += coeff_1*coeff_2*ang_ind_overlap[0];
     *pi += coeff_1*coeff_2*ang_ind_overlap[1];
     *delta += coeff_1*coeff_2*ang_ind_overlap[2];
     *phi += coeff_1*coeff_2*ang_ind_overlap[3];
  }
}

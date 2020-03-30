#include "consistency.h"
/*
double density_correction_MLS(Particle** particles, Setup* setup, int n_p){
  // double rho = particles[j]->rho
  // double m = particles[i]->m
  // kernel = setup->kernel

  return(0);
}*/

xy* correct_grad(xy *current_grad, Particle *p, double kh, Kernel kernel){
	double m11 = 0; double m12 = 0; double m21 = 0; double m22 = 0;
	ListNode *current = p->neighborhood->head;
	while(current != NULL){
		Particle *j = current->v;

		double r = sqrt(pow(p->pos->x - j->pos->x, 2) + pow(p->pos->y - j->pos->y, 2));
		double derivativeKernel = derivative_kernel(r, kh, kernel);

		m11 -= (j->m/j->rho)*derivativeKernel*(1/r)*pow(p->pos->x - j->pos->x,2);
		m12 -= (j->m/j->rho)*derivativeKernel*(1/r)*(p->pos->x - j->pos->x)*(p->pos->y - j->pos->y);
		m21 -= (j->m/j->rho)*derivativeKernel*(1/r)*(p->pos->x - j->pos->x)*(p->pos->y - j->pos->y);
		m22 -= (j->m/j->rho)*derivativeKernel*(1/r)*pow(p->pos->y - j->pos->y,2);
		current = current->next;
	}
	double det = m11*m22 - m12*m21;

	return xy_new((m22*current_grad->x - m21*current_grad->y)/det, (-m12*current_grad->x + m22*current_grad->y)/det);
}

void Corrective_Smoothed_Particle_Method(Particle *p,Particle_derivatives *dp, Setup *setup){
    double num1=0.0;
    double denom1=0.0;
    double denom2=0.0;
    double denom3=0.0;
        //Correction of the density of our particle.
    ListNode *current = p->neighborhood->head;
    while(current != NULL){
        Particle *j = current->v;
        num1 += j->rho*eval_kernel(p->pos,j->pos,setup->kh,setup->kernel)*j->m/j->rho;
        denom1 += eval_kernel(p->pos,j->pos,setup->kh,setup->kernel)*j->m/j->rho;
        current = current->next;
    }

    //p->rho->x = num1/denom1;
        //Correction of the gradient of pressure of our particle
    
    current = p->neighborhood->head;
    while(current != NULL){
        Particle *j = current->v;
        denom2+=eval_kernel(p->pos,j->pos,setup->kh,XXX)*(j->pos->x-p->pos->x)*j->m/j->rho;
        denom3+=eval_kernel(p->pos,j->pos,setup->kh,XXX)*(j->pos->y-p->pos->y)*j->m/j->rho;
        current = current->next;
        
    }
    dp->grad_P->x = dp->grad_P->x/denom2;
    dp->grad_P->y = dp->grad_P->y/denom3;
}
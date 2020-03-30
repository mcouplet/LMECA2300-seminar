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
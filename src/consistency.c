#include "consistency.h"




double get_M0_local(Particle* pi, double kh, Kernel kernel){
	double M = 0;
	ListNode* current = pi->neighborhood->head;
	while(current != NULL){
		Particle* pj = current->v;
		double W = eval_kernel(pi->pos, pj->pos, kh, kernel);
		double V = pj->m/pj->rho;
		M += W*V;
		current = current-> next;
	}
	return M;
}

double get_M1_local(Particle* pi, double kh, Kernel kernel){
	double M = 0;
	ListNode* current = pi->neighborhood->head;
	while(current != NULL){
		Particle* pj = current->v;
		double W = eval_kernel(pi->pos, pj->pos, kh, kernel);
		double V = pj->m/pj->rho;

		double d = sqrt(squared(pi->pos->x-pj->pos->x) + squared(pi->pos->y-pj->pos->y));
		double h = kh/2;
		double Xi_Xj = d/h; //Pas sur
		M += W*V*Xi_Xj;
		current = current-> next;
	}
	return M;
}
void get_M0(Particle** p, int n_p, double kh, Kernel kernel){
	double max = 0;
	for(int i = 0; i < n_p ; i++){
		Particle* pi = p[i];
		double M0 = get_M0_local(pi,kh,kernel);
		if (max <= abs(M0)){
			max = abs(M0);
		}
	}
	printf("Largest M0 error : %f\n", max);
}



xy* correct_grad(xy *current_grad, Particle *p, double kh, Kernel kernel){
	double m11 = 0; double m12 = 0; double m21 = 0; double m22 = 0;
	ListNode *current = p->neighborhood->head;
	while(current != NULL){
		Particle *j = current->v;
		double r = sqrt(pow(p->pos->x - j->pos->x, 2) + pow(p->pos->y - j->pos->y, 2));
		if(r != 0){
			double derivativeKernel = derivative_kernel(r, kh, kernel);
			m11 -= (j->m/j->rho)*derivativeKernel*(1/r)*pow(p->pos->x - j->pos->x,2);
			m12 -= (j->m/j->rho)*derivativeKernel*(1/r)*(p->pos->x - j->pos->x)*(p->pos->y - j->pos->y);
			m21 -= (j->m/j->rho)*derivativeKernel*(1/r)*(p->pos->x - j->pos->x)*(p->pos->y - j->pos->y);
			m22 -= (j->m/j->rho)*derivativeKernel*(1/r)*pow(p->pos->y - j->pos->y,2);
		}
		current = current->next;
	}
	double det = m11*m22 - m12*m21;
	return xy_new((m22*current_grad->x - m21*current_grad->y)/det, (-m12*current_grad->x + m22*current_grad->y)/det);
}
void density_correction_MLS(Particle** p, int n_p, double kh, Kernel kernel){
	double* density_corr = (double*)malloc(n_p * sizeof(double));
	for(int i = 0; i < n_p ; i++){
		Particle* pi = p[i];
		density_corr[i] = density_correction_MLS_local(pi,kh,kernel);
	}
	for(int i = 0; i < n_p ; i++){
		Particle* pi = p[i];
		pi->rho = density_corr[i];
	}
	free(density_corr);
}
double density_correction_MLS_local(Particle* pi, double kh, Kernel kernel){
// We only compute beta once since it does not depend on other particles.
  double beta[3] = {0};
	double A[3][3];
	for(int i = 0; i < 3 ; i++){
		for(int j = 0; j < 3 ; j++){
			A[i][j] = 0;
		}
	}
	get_A(pi,kh,kernel,A);
	get_beta(A,beta);

  double num = 0;
  double den = 0;

  ListNode* current = pi->neighborhood->head;
  while(current != NULL){
    Particle* pj = current->v;
    double W_MLS = get_W_MLS(pi,pj,kh,kernel,beta);
    num += W_MLS*pj->m;
    den += W_MLS*pj->m/pj->rho;
    current = current-> next;
  }
  // Correction
  return (num/den);
}

double get_W_MLS(Particle* pi, Particle* pj, double kh, Kernel kernel, double* beta){
// Return the kernel function corrected by mean least squares
  double W = eval_kernel(pi->pos, pj->pos, kh, kernel);
  double xi = pi->pos->x;   double xj = pj->pos->x;
  double yi = pi->pos->y;   double yj = pj->pos->y;

  double W_MLS = (beta[0] + beta[1]*(xi-xj) + beta[2]*(yi-yj))*W;
  return W_MLS;
}
//
void get_A(Particle* pi, double kh, Kernel kernel, double A[3][3]){
// Return the matrix A which has to be inversed and multiply to obtain beta
  ListNode* current = pi->neighborhood->head;
  while(current != NULL){
    Particle* pj = current->v;
    double W = eval_kernel(pi->pos, pj->pos, kh, kernel);
    double k = W*pj->m/pj->rho;
    double xixj = pi->pos->x - pj->pos->x;
    double yiyj = pi->pos->y - pj->pos->y;
    // Matrix filling
    A[0][0] += k;           A[0][1] += k*xixj;     			A[0][2] += k*yiyj;
    A[1][0] += k*A[0][1];   A[1][1] += k*pow(xixj,2);   A[1][2] += k*xixj*yiyj;
    A[2][0] += k*A[0][2];   A[2][1] += k*A[1][2];  			A[2][2] += k*pow(yiyj,2);
    // go to the next neighbor
    current= current->next;
  }
}
//
void get_beta(double A[3][3], double* beta){
  // B = [1,0,0]^T
  double C11 =   A[0][0]*(A[1][1]*A[2][2] - A[0][2]*A[2][1]);
  double C21 = - A[0][2]*(A[0][1]*A[2][2] - A[1][2]*A[2][0]);
  double C31 =   A[0][2]*(A[1][0]*A[2][1] - A[1][1]*A[2][0]);

  double det_A = C11 + C21 + C31;

	beta[0] = C11/det_A;
	beta[1] = C21/det_A;
	beta[2] = C31/det_A;
}

// void Corrective_Smoothed_Particle_Method(Particle *p,Particle_derivatives *dp, double kh, Kernel kernel){
//     double num1=0.0;
//     double denom1=0.0;
//     double denom2=0.0;
//     double denom3=0.0;
//         //Correction of the density of our particle.
//     ListNode *current = p->neighborhood->head;
//     while(current != NULL){
//         Particle *j = current->v;
//         num1 += j->rho*eval_kernel(p->pos,j->pos,kh,kernel)*j->m/j->rho;
//         denom1 += eval_kernel(p->pos,j->pos,kh,kernel)*j->m/j->rho;
//         current = current->next;
//     }
//
//     p->rho = num1/denom1;
//         //Correction of the gradient of pressure of our particle
//
//     current = p->neighborhood->head;
//     while(current != NULL){
//         Particle *j = current->v;
//         denom2+=eval_kernel(p->pos,j->pos,kh,XXX)*(j->pos->x-p->pos->x)*j->m/j->rho;
//         denom3+=eval_kernel(p->pos,j->pos,kh,XXX)*(j->pos->y-p->pos->y)*j->m/j->rho;
//         current = current->next;
//
//     }
//     dp->grad_P->x = dp->grad_P->x/denom2;
//     dp->grad_P->y = dp->grad_P->y/denom3;
// }

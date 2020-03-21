#include "derivatives.h"

double compute_div(Particle * particle, xy_getter get, Kernel kernel, double kh) {
    double div = 0;
    Particle *pi = particle;
    xy *fi = get(pi);
    ListNode *node = pi->neighborhood->head;
    while(node != NULL) {
        Particle *pj = node->v;
        xy *grad_W = grad_kernel(pi->pos, pj->pos, kh, kernel);
        xy *fj = get(pj);
        div += (pj->m / pi->rho) *
            ((fj->x - fi->x) * grad_W->x + (fj->y - fi->y) * grad_W->y);
        free(grad_W);
        node = node->next;
    }
    return div;
}

xy * compute_grad(Particle * particle, scalar_getter get, Kernel kernel, double kh) {
    xy *grad = xy_new(0,0);
    Particle *pi = particle;
    double fi = get(pi);
    ListNode *node = pi->neighborhood->head;
    while(node != NULL) {
        Particle *pj = node->v;
        xy *grad_W = grad_kernel(pi->pos, pj->pos, kh, kernel);
        double fj = get(pj);
        grad->x += pi->rho * pj->m * (fi/squared(pi->rho) + fj/squared(pj->rho)) * grad_W->x; // sign is not the same as in the def...
        grad->y += pi->rho * pj->m * (fi/squared(pi->rho) + fj/squared(pj->rho)) * grad_W->y;
        free(grad_W);
        node = node->next;
    }
    return grad;
}

double compute_lapl(Particle *particle, scalar_getter get, Kernel kernel, double kh) {
    double lapl = 0;
    Particle *pi = particle;
    double fi = get(pi);
    ListNode *node = pi->neighborhood->head;
    while(node != NULL) {
        Particle *pj = node->v;
        xy *grad_W = grad_kernel(pi->pos, pj->pos, kh, kernel);
        double fj = get(pj);
        double d2 = squared(pi->pos->x - pj->pos->x) + squared(pi->pos->y - pj->pos->y); // squared distance between particles
        xy *DXij = xy_new((pi->pos->x - pj->pos->x) / d2, (pi->pos->y - pj->pos->y) / d2); // Delta X_{ij}
        lapl += 2 * pj->m / pj->rho * (fi - fj) * (DXij->x * grad_W->x + DXij->y * grad_W->y);
        free(grad_W);
        free(DXij);
        node = node->next;
    }
    return lapl;
}

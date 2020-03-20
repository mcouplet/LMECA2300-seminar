#include "derivatives.h"

double compute_div(Particle * particle, getter get, Kernel kernel, double kh) {
    double div = 0;
    Particle *pi = particle;
    ListNode *node = pi->neighborhood->head;
    while(node != NULL) {
        Particle *pj = node->v;
        xy *grad_W = grad_kernel(pi->pos, pj->pos, kh, kernel);
        div += (pj->m / pi->density) *
            ((pj->v->x - pi->v->x) * grad_W->x + (pj->v->y - pi->v->y) * grad_W->y);
        node = node->next;
    }
    return div;
}

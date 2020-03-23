#ifndef UTILS_H
#define UTILS_H
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define M_PI 3.14159265358979323846

typedef struct List List; // Generic list
typedef struct ListNode ListNode;
typedef struct xy xy;
typedef struct Residual Residual;
typedef void(*destructor)(void *object);

struct List {
	int n; // list size
	ListNode *head; // pointer to first node for traversal
	ListNode *tail; // pointer to last node for append
};

struct ListNode {
	void *v; // node content
	ListNode *next; // pointer to next node
};

struct xy {
	double x;
	double y;
};

struct Residual {
	double mass_eq;
	double momentum_x_eq;
	double momentum_y_eq;
};

List* List_new();
void List_append(List *l, void *v);
void List_free(List*, destructor);

xy* xy_new(double x, double y);
void xy_reset(xy *p);
double rand_interval(double a, double b);

Residual* residual_new();

double squared(double x);
double norm(xy *v);

xy* map_to_circle(xy* pos_square);
xy* generate_circle(int k, int n, int nb, double radius);

#endif

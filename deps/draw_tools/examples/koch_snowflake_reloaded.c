#include "draw_tools.h"
#include <math.h>

#define TRANSITION_TIME 2.0

#define SQRT3_6 0.288675135 // sqrt(3)/6

// create p1, p2 and p3 from the line segment p0->p4, with interpolation factor 's' between 0 and 1
static void divide(GLfloat* p0, GLfloat* p1, GLfloat* p2, GLfloat* p3, GLfloat* p4, float s)
{
	float d[2] = {p4[0] - p0[0], p4[1] - p0[1]};

	for (int i=0; i<2; i++) {
		p1[i] = p0[i] + s*1./3.*d[i] + (1-s)*0.25*d[i];
		p2[i] = p0[i] + 0.5*d[i];
		p3[i] = p0[i] + s*2./3.*d[i] + (1-s)*0.75*d[i];
	}

	p2[0] -= s*SQRT3_6*d[1];
	p2[1] += s*SQRT3_6*d[0];
}

// update the properties (width and outline width) of the pointset (given width)
static void points_properties_update(points_t* pointset, float segWidth)
{
	points_set_width(pointset, segWidth);
	points_set_outline_width(pointset, 0.66*segWidth);
}

// 2x more indices 0 x 2x 3x instead of 0 2x 4x 8x ...
static void indices_update(GLuint* indices, size_t nSegment, size_t maxSegment)
{
	size_t subStep = maxSegment/nSegment;

	indices[0] = 0;
	for (size_t i=0; i<nSegment+1; i++) {
		indices[i+1] = i*subStep;
	}
	indices[nSegment+2] = maxSegment;
}

// update the coordinates, with interpolation factor 's' between 0 and 1
static void coords_update(GLfloat* coords, GLuint* indices, size_t nSegment, float s)
{
	for (size_t i=0; i<nSegment; i++) {
		divide(&coords[2*indices[4*i+1]],
		       &coords[2*indices[4*i+2]],
		       &coords[2*indices[4*i+3]],
		       &coords[2*indices[4*i+4]],
		       &coords[2*indices[4*i+5]],
		       s);
	}
}

// smooth second derative (start and stop progressively)
static inline float smoothstep(float x) {
	return x * x * (3 - 2 * x);
}


int main(int argc, char *argv[])
{
	window_t* window = window_new(0,0, argv[0]);
	window_set_color(window, (float[4]) {0.3, 0.3, 0.3, 1.0});

	const int maxIters = 5;
	const size_t maxSegment = 1ULL<<(2*maxIters); // 4^maxIter with bits operations
	const size_t maxPoints = maxSegment + 1;

	// allocations
	GLfloat* coords = malloc(2*sizeof(GLfloat)*maxPoints);
	GLuint* indices = malloc(sizeof(GLuint)*(maxPoints+2));
	points_t* pointset = points_new(NULL, maxPoints, GL_STATIC_DRAW);
	order_t* order = order_new(NULL, maxPoints+2, GL_DYNAMIC_DRAW);

	// initialization
	size_t nSegment = 1;
	float segWidth = 0.5;

	coords[0] = -1.0;
	coords[1] = 0.0;
	coords[2*maxSegment] = 1.0;
	coords[2*maxSegment+1] = 0.0;
	points_set_color(pointset, (float[4]) {0.7, 0.5, 0.0, 1.0});

	// iterations
	for (int i=0; i<maxIters-1; i++) {
		indices_update(indices, 4*nSegment, maxSegment);
		order_update(order, indices, 4*nSegment + 3);

		// animation (here we recompute the points each times)
		double tbegin = window_get_time(window);
		double tnow = tbegin;
		while(tnow - tbegin < TRANSITION_TIME) {
			// interpolation factor
			// float s = (tnow - tbegin)/TRANSITION_TIME;
			float s = smoothstep((tnow - tbegin)/TRANSITION_TIME);
			
			coords_update(coords, indices, nSegment, s);
			points_update(pointset, coords, maxPoints);
			points_properties_update(pointset, (1-s)*segWidth + s*segWidth/3);

			curve_draw_with_order(window, pointset, order);
			window_update(window);

			tnow = window_get_time(window);

			if(window_is_closed(window))
				goto end_of_program;
		}

		nSegment*=4;
		segWidth*=1./3.;
	}

	text_t* end = text_new("Max. iteration level reached", GL_STATIC_DRAW);
	text_set_pos(end, -0.9, -0.9);
	text_set_space_type(end, PIXEL_WIDTHS_WITHOUT_TRANSLATIONS);
	text_set_scale(end, 48, 48);

	points_properties_update(pointset, segWidth);
	points_update(pointset, coords, maxPoints);

	while(!window_is_closed(window)) {
		curve_draw_with_order(window, pointset, order);
		text_draw(window, end);
		window_update_and_wait_events(window);
	}

	text_delete(end);

end_of_program:

	order_delete(order);
	free(indices);
	points_delete(pointset);
	free(coords);
	window_delete(window);

	return EXIT_SUCCESS;
}

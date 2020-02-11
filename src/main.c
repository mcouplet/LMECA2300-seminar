#include "BOV.h"
#include <math.h>

// see stringification process
#define xstr(s) str(s)
#define str(s) #s
#define NPOINTS 10000

/* for v, a floating point value between 0 and 1, this function fills color with
 * the improved jet colormap color corresponding to v */
static void colormap(float v, float color[3])
{
	float v1 = 3.5*(v-0.7);
	float v2 = 1.25*v;
	float v3 = fminf(0.5,v)*2.0;
    
    color[0] = 0.8f*(-v1*v1+1.0f);
    color[1] = 0.8f*(6.0f*v2*v2*(1.0f-v2));
    color[2] = 0.8f*(5.5f*v3*(1.0f-v3)*(1.0f-v3));
}


static void fillData(GLfloat (*data)[8])
{
	float rmax = 100.0*sqrtf(2.0f);
	for(int i=0; i<NPOINTS; i++) {
		data[i][0] = rand() * 200.0 / RAND_MAX - 100.0; // x (rand between -100 and 100)
		data[i][1] = rand() * 200.0 / RAND_MAX - 100.0; // y (rand between -100 and 100)
		data[i][2] = 0; // speed x (not used by default visualization)
		data[i][3] = 0; // speed y (not used by default visualization)
		float r = sqrt(data[i][0]*data[i][0] + data[i][1]*data[i][1]);
		colormap(r/rmax, &data[i][4]); // fill color
		data[i][7] = 0.8f; // transparency
	}
}


int main()
{
	bov_window_t* window = bov_window_new(1024, 780, "ANM Project: SPH");
	bov_window_set_color(window, (GLfloat[]){0.9f, 0.85f, 0.8f, 0.0f});

	GLfloat (*data)[8] = malloc(sizeof(data[0])*NPOINTS);
	fillData(data);

	/* send data to GPU, and receive reference to those data in a points object */
	bov_points_t *particles = bov_particles_new(data, NPOINTS, GL_STATIC_DRAW);

	/* setting particles appearance */
	bov_points_set_width(particles, 0.02);
	bov_points_set_outline_width(particles, 0.0025);

	/* fit particles in a [-0.8 0.8]x[-0.8 0.8] square */
	bov_points_scale(particles, (GLfloat[2]) {0.008, 0.008});
	bov_points_set_pos(particles, (GLfloat[2]) {0.0, -0.1});

	/* we got 0.2 at the top to write something */
	bov_text_t* msg =  bov_text_new((unsigned char[]){
			"Rendering " xstr(NPOINTS) " particles"
		},
		GL_STATIC_DRAW);
	bov_text_set_pos(msg, (GLfloat[2]){-0.95, 0.82});
	bov_text_set_fontsize(msg, 0.1);

	while(!bov_window_should_close(window)){
		bov_particles_draw(window, particles);
		// bov_points_draw(window, particles, 0, BOV_TILL_END);

		bov_text_draw(window, msg);

		bov_window_update(window);
	}

	bov_text_delete(msg);
	bov_points_delete(particles);
	free(data);
	bov_window_delete(window);

	return EXIT_SUCCESS;
}
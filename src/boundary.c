#include "boundary.h"

Boundary* Boundary_new(double xleft, double xright, double ybottom, double ytop){
	Boundary* boundary = (Boundary*) malloc(sizeof(Boundary));
	boundary->xleft = xleft;
	boundary->xright = xright;
	boundary->ybottom = ybottom;
	boundary->ytop = ytop;
	return boundary;
}
void Boundary_free(Boundary* boundary){
	free(boundary);
}

void center_reflection_right(Particle* pi, double CR, double Rp, double d){
		pi->pos->x -= (1+CR)*(Rp - d);
}
void center_reflection_left(Particle* pi, double CR, double Rp, double d){
		pi->pos->x += (1+CR)*(Rp - d);
}
void center_reflection_top(Particle* pi, double CR, double Rp, double d){
		pi->pos->y -= (1+CR)*(Rp - d);
}
void center_reflection_bottom(Particle* pi, double CR, double Rp, double d){
		pi->pos->y += (1+CR)*(Rp - d);
}
void velocity_reflection_vertical(Particle* pi, double CR, double CF){
	double vpN = pi->v->x;
	double vpT = pi->v->y;
	pi->v->x = -vpN*CR;
	pi->v->y = (1-CF)* vpT;
}
void velocity_reflection_horizontal(Particle* pi, double CR, double CF){
	double vpN = pi->v->y;
	double vpT = pi->v->x;
	pi->v->x = vpN*CR;
	pi->v->y = -(1-CF)* vpT;
}
void reflective_boundary(Particle** p, int n_p, double CR, double CF, Boundary* boundary,double Rp){
	// We have just computed the time integration. We correct the positions of the particles
	for(int i = 0; i < n_p ; i++){
		// For each particle we check if its position is close to a wall
		Particle* pi = p[i];
		double xright = boundary->xright;	double xleft = boundary->xleft;
		double ytop = boundary->ytop;			double ybottom = boundary->ybottom;
		double dxright,dxleft,dytop,dybottom;

		// Collision test
		int collision, collision_right, collision_left, collision_top, collision_bottom;
		collision = 0; collision_right = 0; collision_left = 0; collision_top = 0; collision_bottom = 0;

		if(pi->pos->x > boundary->xright - Rp){
			collision_right = 1;
			collision = 1;
			dxright = fabs(pi->pos->x - xright);
			center_reflection_right(pi,CR,Rp,dxright);
			velocity_reflection_vertical(pi,CR,CF);
		}
		if(pi->pos->x < boundary->xleft + Rp){
			collision_left = 1;
			collision = 1;
			dxleft = fabs(pi->pos->x - xleft);
			center_reflection_left(pi,CR,Rp,dxleft);
			velocity_reflection_vertical(pi,CR,CF);
		}
		if(pi->pos->y > boundary->ytop - Rp){
			collision_top = 1;
			collision = 1;
			dytop = fabs(pi->pos->y - ytop);
			center_reflection_top(pi,CR,Rp,dytop);
			velocity_reflection_horizontal(pi,CR,CF);
		}
		if(pi->pos->y < boundary->ybottom + Rp){
			collision_bottom = 1;
			collision = 1;
			dybottom = fabs(pi->pos->y - ytop);
			center_reflection_bottom(pi,CR,Rp,dytop);
			velocity_reflection_horizontal(pi,CR,CF);
		}
	}
}
// void reflective_boundary(GLfloat(*data)[8], GLfloat(*coord)[2], double timestep, double xmin, double xmax, double ymin, double ymax, double iter)
// {
// 	for (int i = 0; i < NPTS; i++) {
//
// 		double maxspeed = 1000;
//
// 		if (iter == 0)
// 		{
// 			data[i][2] = rand() *maxspeed / RAND_MAX - 1.0;
// 			data[i][3] = rand() *maxspeed / RAND_MAX - 1.0;
// 		}
//
// 		//UPDATE COORD WITH EULEUR'S METHOD
// 		data[i][0] += data[i][2] * timestep;
// 		data[i][1] += data[i][3] * timestep;
// 		coord[i][0] = data[i][0];
// 		coord[i][1] = data[i][1];
//
// 		double d_x;
//
// 		//RIGTH BOUNDARY
// 		if(data[i][0] > xmax-particule_radius)
// 		{
// 			data[i][2] = -CR * data[i][2];
// 			data[i][3] = (1 - CF) * data[i][3];
// 			while (data[i][0] > xmax-particule_radius)
// 			{
// 				d_x = abs(data[i][0] - xmax);
// 				data[i][0] = data[i][0] - (1 + CR)*(particule_radius - d_x);
// 			//	data[i][0] = data[i][0] - (d_x + particule_radius)*(1 + CR);
// 			}
// 		}
//
// 		//LEFT BOUNDARY
// 		if (data[i][0] < xmin+particule_radius)
// 		{
// 			data[i][2] = -CR * data[i][2];
// 			data[i][3] = (1 - CF) * data[i][3];
// 			while (data[i][0] < xmin+particule_radius)
// 			{
// 				d_x = abs(data[i][0] - xmin);
// 				data[i][0] = data[i][0] + (1 + CR)*(particule_radius - d_x);
// 			//	data[i][0] = data[i][0] + (d_x + particule_radius)*(1 + CR);
// 			}
// 		}
//
// 		double d_y;
//
// 		//UP BOUNDARY
// 		if (data[i][1] > ymax-particule_radius)
// 		{
// 			data[i][3] = -CR * data[i][3];
// 			data[i][2] = (1 - CF) * data[i][2];
// 			while (data[i][1] > ymax-particule_radius)
// 			{
// 				d_y = abs(data[i][1] - ymax);
// 				data[i][1] = data[i][1] - (1 + CR)*(particule_radius - d_y);
// 				//data[i][1] = data[i][1] - (d_y + particule_radius)*(1 + CR);
// 			}
// 		}
//
// 		//LEFT BOUNDARY
// 		if (data[i][1] < ymin+particule_radius)
// 		{
// 			data[i][3] = -CR * data[i][3];
// 			data[i][2] = (1 - CF) * data[i][2];
// 			while (data[i][1] < ymin+particule_radius)
// 			{
// 				d_y = abs(data[i][1] - ymin);
// 			    data[i][1] = data[i][1] + (1 + CR)*(particule_radius - d_y);
// 				//data[i][1] = data[i][1] + (d_y + particule_radius)*(1 + CR);
// 			}
// 		}
// 	}
// }
//
// void fill_fictitious_particles(GLfloat(*data_fict)[8], GLfloat(*coord_fict)[2], double minTemp, double maxTemp, float(*sol)[5])
// {
// 	int j = 0;
// 	for (int i = 0; i < NPTS_FICT / 4; i++) //LEFT BOUNDARY
// 	{
// 		    data_fict[i][0] = -100;
// 			data_fict[i][1] = -100+j*(2*100/(NPTS_FICT/4));
// 			coord_fict[i][0] = data_fict[i][0];
// 			coord_fict[i][1] = data_fict[i][1];
// 			j += 1;
// 	}
// 	j = 0;
// 	for (int i = NPTS_FICT/4; i < 2*NPTS_FICT / 4; i++) //UP BOUNDARY
// 	{
// 		data_fict[i][0] = -100+j* (2 * 100 / (NPTS_FICT / 4));
// 		data_fict[i][1] = 100;
// 		coord_fict[i][0] = data_fict[i][0];
// 		coord_fict[i][1] = data_fict[i][1];
// 		j += 1;
// 	}
// 	j = 0;
// 	for (int i = 2*NPTS_FICT / 4; i < 3 * NPTS_FICT / 4; i++) //RIGHT BOUNDARY
// 	{
// 		data_fict[i][0] = 100;
// 		data_fict[i][1] = 100-j* (2 * 100 / (NPTS_FICT / 4));
// 		coord_fict[i][0] = data_fict[i][0];
// 		coord_fict[i][1] = data_fict[i][1];
// 		j += 1;
// 	}
// 	j = 0;
// 	for (int i = 3 * NPTS_FICT / 4; i < 4 * NPTS_FICT / 4; i++) //DOWN BOUNDARY
// 	{
// 		data_fict[i][0] = 100-j* (2 * 100 / (NPTS_FICT / 4));
// 		data_fict[i][1] = -100;
// 		coord_fict[i][0] = data_fict[i][0];
// 		coord_fict[i][1] = data_fict[i][1];
// 		j += 1;
// 	}
//
// 	//SOURCETEMPS
// 	for (int i = 0; i < NPTS_FICT ; i++)
// 	{
// 		tempToColor(maxTemp, &data_fict[i][4], minTemp, maxTemp); // fill color
// 		sol[i+NPTS][0] = maxTemp;
// 		data_fict[i][7] = 1; // transparency
// 	}
// }
//
// void repulsive_boundary(GLfloat(*data)[8], GLfloat(*data_fict)[8], GLfloat(*coord)[2], GLfloat(*coord_fict)[2] ,double timestep, double xmin, double xmax, double ymin, double ymax, double time)
// {
//
// 	double maxspeed = 100;
// 	double psi = 0.00002;
// 	double r_0 = 4; //cut off distance
// 	double d_x;
// 	double d_xi;
// 	double d_y;
// 	double d_yi;
// 	double F_x;
// 	double F_y;
// 	double acc_x;
// 	double acc_y;
//
// 	for (int i = 0; i < NPTS; i++){
//
// 		if (time == 0)
// 		{
// 			data[i][2] = rand() *maxspeed / RAND_MAX - 1.0;
// 			data[i][3] = rand() *maxspeed / RAND_MAX - 1.0;
// 		}
//
// 		//UPDATE COORD WITH EULEUR'S METHOD
// 		data[i][0] += data[i][2] * timestep;
// 		data[i][1] += data[i][3] * timestep;
// 		coord[i][0] = data[i][0];
// 		coord[i][1] = data[i][1];
//
// 		d_x = min(abs(data[i][0] - xmax), abs(data[i][0] - xmin));
// 		d_y = min(abs(data[i][1] - ymax), abs(data[i][1] - ymin));
// 		F_x = 0;
// 		F_y = 0;
//
// 		if (r_0 / d_x > 1.0 || r_0 / d_y > 1.0)
// 		{
// 			for (int j = 0; j < NPTS_FICT; j++)
// 			{
// 				d_xi = abs(data[i][0] - data_fict[j][0]);
// 				d_yi = abs(data[i][1] - data_fict[j][1]);
// 				if (r_0 / d_xi > 1.0 && d_xi !=0)
// 				{
// 					F_x -= psi * (pow(r_0 / d_xi, 12) - pow(r_0 / d_xi, 4))*((data_fict[j][0]-data[i][0])/ pow(d_xi,2));
// 					//printf("force x = %f \n", psi * (pow(r_0 / d_xi, 12) - pow(r_0 / d_xi, 4))*((data_fict[j][0] - data[i][0]) / pow(d_xi, 2)));
// 				}
// 				if (r_0 / d_yi > 1.0 && d_yi !=0)
// 				{
// 					F_y -= psi * (pow(r_0 / d_yi, 12) - pow(r_0 / d_yi, 4))*((data_fict[j][1]-data[i][1]) / pow(d_yi,2));
// 				//	printf("force y = %f \n", psi * (pow(r_0 / d_yi, 12) - pow(r_0 / d_yi, 4))*((data_fict[j][1] - data[i][1]) / pow(d_yi, 2)));
// 				}
// 			}
//
// 				acc_x = F_x / MASS;
// 				acc_y = F_y / MASS;
//
// 				//printf("acc x : %f \n",acc_x);
// 				//printf("acc y : %f \n", acc_y);
//
// 				data[i][2] += acc_x * timestep;
// 				data[i][0] += data[i][2] * timestep;
// 				//printf("speed before : %f \n",data[i][3]);
// 				data[i][3] += acc_y * timestep;
// 				//printf("speed after : %f \n", data[i][3]);
// 				data[i][1] += data[i][3] * timestep;
// 		}
// 	}

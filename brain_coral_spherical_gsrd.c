//	Spherical Gray-Scott Reaction-Diffusion Brain Coral Simulation

//Compile:
//	cc brain_coral_spherical_gsrd.c -o ~/brain_coral_spherical_gsrd -lm -lX11 -pthread 

//Run with the defaults:
//	~/brain_coral_spherical_gsrd

//Run with specific parameters:
//	~/brain_coral_spherical_gsrd 20000 0.0256 0.0544 0.2

//  tab stop @ 8
/*****************************************************************************/
#define MESH_DETAIL	99	//0 for undivided icosahedron
#define SUBDIV	(1+MESH_DETAIL)	//1 for undivided icosahedron
#define T_NUM_	(SUBDIV*SUBDIV)	//triangulation number
#define POINTS	(T_NUM_*10 + 2)	
#define FACETS	(T_NUM_*20)	

#define N_THRDS 	8
#include <pthread.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <errno.h>		//	ESRCH
#include <stdio.h>		//	printf
#include <stdlib.h>		//	atoi atof
#include <math.h>		//	sqrt
const	double	PI = M_PI;

/*****************************************************************************/
typedef struct
{
	double	x, y, z;	//vertex coordinates
	double	area;		//vertex area enclosed by centers of inscribed circles within neighboring triangles
	double	norm;		//area | mass normalization (the smaller the vertex area, the faster the change [in temparature | concentration])
	double	cotan;		//negative of twice the diagonal of cotangent matrix (-2*Cii)
	double	u, Lu;		//Laplacian subject to area | mass normalization
	double	v, Lv;		//two-component reaction–diffusion system
	int	NFcnt, NF[6];	//neighboring facets
	int	NPcnt, NP[6];	//neighboring points
}
Point;

/*****************************************************************************/
typedef struct
{
	int	A, B, C;	//vertex points of triangle ABC
	double	a2;		//square of side length a for segment BC =(xB-xC)^2+(yB-yC)^2+(zB-zC)^2
	double	b2;		//square of side length b for segment CA =(xC-xA)^2+(yC-yA)^2+(zC-zA)^2
	double	c2;		//square of side length c for segment AB =(xA-xB)^2+(yA-yB)^2+(zA-zB)^2
	double	T;		//area of triangle ABC =SQRT((a2+b2+c2)^2-2*(a2^2+b2^2+c2^2))/4
	double	cot_alpha;	//cotangent of angle CAB =(-a2+b2+c2)/4/T
	double	cot_beta_;	//cotangent of angle ABC =(+a2-b2+c2)/4/T
	double	cot_gamma;	//cotangent of angle BCA =(+a2+b2-c2)/4/T

	double	uC_uB_cotA;	//(uC - uB) * cotangent alpha
	double	uA_uC_cotB;	//(uA - uC) * cotangent beta_
	double	uB_uA_cotC;	//(uB - uA) * cotangent gamma

	double	vC_vB_cotA;	//(vC - vB) * cotangent alpha
	double	vA_vC_cotB;	//(vA - vC) * cotangent beta_
	double	vB_vA_cotC;	//(vB - vA) * cotangent gamma
}
Facet;

/*****************************************************************************/
typedef struct
{
	double	r;		//radial (from center)
	double	t;		//theta (inclination from north pole/90-latitude|elevation)
	double	p;		//phi (azimuth/longitude)
}
Spherical_Coords;		//physics convention
Spherical_Coords
	icosahedron_sphord[ 12 ];
Point	icosahedron_points[ 12 ];
Facet	icosahedron_facets[ 20 ]=
{	{10	,0	,2	},	{0	,1	,2	},	{3	,2	,1	},	{3	,1	,11	}
,	{10	,2	,4	},	{2	,3	,4	},	{5	,4	,3	},	{5	,3	,11	}
,	{10	,4	,6	},	{4	,5	,6	},	{7	,6	,5	},	{7	,5	,11	}
,	{10	,6	,8	},	{6	,7	,8	},	{9	,8	,7	},	{9	,7	,11	}
,	{10	,8	,0	},	{8	,9	,0	},	{1	,0	,9	},	{1	,9	,11	}
};//	NP:north circumpolar		NQ:north equatorial		SQ:south equatorial		SP:south circumpolar		

/*****************************************************************************/
Point spherial_to_cartesian( Spherical_Coords coord )
{
	Point	P;
	double	r = coord.r;	//radial (from center)
	double	t = coord.t;	//theta (from north pole/90-latitude)
	double	p = coord.p;	//phi (longitude/prime meridian)
	double	r_sin_t = r*sin(t);
		P.z 	= r*cos(t);		//z = r cos_theta
		P.y 	= r_sin_t*sin(p);	//y = r sin_theta sin_phi
		P.x 	= r_sin_t*cos(p);	//x = r sin_theta cos_phi
	return	P;
}

/*****************************************************************************/
void init_icosahedron()
{
	double	PI5 = PI/5;			//36 degrees
	double	t = atan(2);			//theta
	double	r = 1;
	Point*P 		= icosahedron_points;
	Spherical_Coords*S 	= icosahedron_sphord;
	S[ 0].r=r	;S[ 0].t=t	;S[ 0].p=0	;P[ 0]=spherial_to_cartesian(S[ 0]);
	S[ 1].r=r	;S[ 1].t=PI-t	;S[ 1].p=-PI5	;P[ 1]=spherial_to_cartesian(S[ 1]);
	S[ 2].r=r	;S[ 2].t=t	;S[ 2].p=-PI5*2	;P[ 2]=spherial_to_cartesian(S[ 2]);
	S[ 3].r=r	;S[ 3].t=PI-t	;S[ 3].p=-PI5*3	;P[ 3]=spherial_to_cartesian(S[ 3]);
	S[ 4].r=r	;S[ 4].t=t	;S[ 4].p=-PI5*4	;P[ 4]=spherial_to_cartesian(S[ 4]);
	S[ 5].r=r	;S[ 5].t=PI-t	;S[ 5].p=+PI	;P[ 5]=spherial_to_cartesian(S[ 5]);
	S[ 6].r=r	;S[ 6].t=t	;S[ 6].p=+PI5*4	;P[ 6]=spherial_to_cartesian(S[ 6]);
	S[ 7].r=r	;S[ 7].t=PI-t	;S[ 7].p=+PI5*3	;P[ 7]=spherial_to_cartesian(S[ 7]);
	S[ 8].r=r	;S[ 8].t=t	;S[ 8].p=+PI5*2	;P[ 8]=spherial_to_cartesian(S[ 8]);
	S[ 9].r=r	;S[ 9].t=PI-t	;S[ 9].p=+PI5	;P[ 9]=spherial_to_cartesian(S[ 9]);
	S[10].r=r	;S[10].t=0	;S[10].p=0	;P[10]=spherial_to_cartesian(S[10]);	//North Pole
	S[11].r=r	;S[11].t=PI	;S[11].p=0	;P[11]=spherial_to_cartesian(S[11]);	//South Pole
}

/*****************************************************************************/
static	Point	point[ POINTS ];
static	Facet	facet[ FACETS ];

/*****************************************************************************/
double	dot( Point A, Point B )
{
	return	( A.x * B.x
		+ A.y * B.y
		+ A.z * B.z );
}

/*****************************************************************************/
Point	SLERP( Point p0, Point p1, double t, double cos_omega, double sin_omega, double omega )	//spherical linear interpolation
{
	if( t == 0.0 )	return	p0;
	if( t == 1.0 )	return	p1;
	double	prop1 = sin(    t *omega ) / sin_omega;		/* t	in linear interpolation (lerp) */
	double	prop0 = sin( (1-t)*omega ) / sin_omega;		/* 1-t	in linear interpolation (lerp) */
	Point	P;
	P.x = p0.x * prop0 + p1.x * prop1;
	P.y = p0.y * prop0 + p1.y * prop1;
	P.z = p0.z * prop0 + p1.z * prop1;
	return	P;
}

/*****************************************************************************/
Point	slerp( Point p0, Point p1, double t )	//spherical linear interpolation
{
	if( t == 0.0 )	return	p0;
	if( t == 1.0 )	return	p1;
	double	cos_omega = dot( p0, p1 );
	double	sin_omega = sqrt( 1 - cos_omega*cos_omega );
	double	omega = acos( cos_omega );
	return	SLERP( p0, p1, t, cos_omega, sin_omega, omega );
}

/*****************************************************************************/
void geodesic_mesh()
{
	int divisions = SUBDIV;
	int T = divisions*divisions;
	for( int s = 0 ; s < 5 ; s++ )	//sector in 0,1,2,3,4
	{							//sector	0	1	2	3	4
		Point*	P0 = &icosahedron_points[ (2*s  )%10 ];	//icos		0	2	4	6	8
		Point*	P1 = &icosahedron_points[ (2*s+1)%10 ];	//icos		1	3	5	7	9
		Point*	P2 = &icosahedron_points[ (2*s+2)%10 ];	//icos		2	4	6	8	0
		Point*	P3 = &icosahedron_points[ (2*s+3)%10 ];	//icos		3	5	7	9	1
		Point*	PN = &icosahedron_points[ 10 ];		//icosahedron	N	N	N	N	N
		Point*	PS = &icosahedron_points[ 11 ];		//icosahedron	S	S	S	S	S

		/*****************************************************************************/
		for( int j = 0 ; j < divisions ; j++ )
		for( int k = 0 ; k < divisions ; k++ )	//iterate to North Pole
		{
			int i = 2*T*s + divisions*j + k;
			Point P;

			if( j == 0 && k == 0 )
				P = *P0;
			else if( j < k )
			//Circumpolar panel
			{
				Point	Pk0N = slerp( *P0, *PN, k*1.0 / divisions );
				Point	Pk02 = slerp( *P0, *P2, k*1.0 / divisions );
				P = slerp( Pk0N, Pk02, 1.0 * j/k );
			}
			else
			//Equatorial panel
			{
				Point	Pj01 = slerp( *P0, *P1, j*1.0 / divisions );
				Point	Pj02 = slerp( *P0, *P2, j*1.0 / divisions );
				P = slerp( Pj01, Pj02, 1.0 * k/j );
			}

			point[ i ] = P;
		}

		for( int j = 0 ; j < divisions ; j++ )	//iterate to South Pole
		for( int k = 0 ; k < divisions ; k++ )
		{
			int i = T + 2*T*s + divisions*j + k;
			Point P;

			if( j == 0 && k == 0 )
				P = *P1;
			else if( j < k )
			//Equatorial panel
			{
				Point	Pk12 = slerp( *P1, *P2, k*1.0 / divisions );
				Point	Pk13 = slerp( *P1, *P3, k*1.0 / divisions );
				P = slerp( Pk12, Pk13, 1.0 * j/k );
			}
			else
			//Circumpolar panel
			{
				Point	Pj1S = slerp( *P1, *PS, j*1.0 / divisions );
				Point	Pj13 = slerp( *P1, *P3, j*1.0 / divisions );
				P = slerp( Pj1S, Pj13, 1.0 * k/j );
			}

			point[ i ] = P;
		}

		/*****************************************************************************/
		for( int j = 0 ; j<2*divisions ; j++ )	//iterate to South Pole
		for( int k = 0 ; k < divisions ; k++ )	//iterate to North Pole
		{
			int E = 2*T*s + divisions*j + k;
			int N = (E + 1) % (POINTS - 2);
			int S = (E + divisions) % (POINTS - 2);
			int W = (E + divisions + 1) % (POINTS - 2);

			if( k == divisions - 1 )	//fix NW edges
			{
			    if( j < divisions )		//North Polar
			    {
				if( j == 0 )
				{
					N = POINTS-2;	//North Pole
					W = (2*T*(s+1) + (divisions-1-j)) % (POINTS - 2);
				}
				else	//j>0
				{
					N = (2*T*(s+1) + (divisions - j)) % (POINTS - 2);
					W = (2*T*(s+1) + (divisions-1-j)) % (POINTS - 2);
				}
			    }
			    else	//j>=divisions	//South Equatorial
			    {
					N = (2*T*(s+1) + (j - divisions)*divisions) % (POINTS - 2);
					W = (2*T*(s+1) + (j+1-divisions)*divisions) % (POINTS - 2);
			    }
			}

			if( j == 2*divisions - 1 )	//fix SW edge
			{
				if( k == 0 )
				{
					S = POINTS-1;	//South Pole
					W = (T + 2*T*(s+1) + (divisions-1-k)*divisions) % (POINTS - 2);
				}
				else	//k>0
				{
					S = (T + 2*T*(s+1) + (divisions - k)*divisions) % (POINTS - 2);
					W = (T + 2*T*(s+1) + (divisions-1-k)*divisions) % (POINTS - 2);
				}
			}

			facet[ 2*E ].A = E;	facet[2*E+1].A = E;
			facet[ 2*E ].B = W;	facet[2*E+1].B = S;
			facet[ 2*E ].C = N;	facet[2*E+1].C = W;
		}
	}

	point[ POINTS-2 ] = icosahedron_points[ 10 ];	//North Pole
	point[ POINTS-1 ] = icosahedron_points[ 11 ];	//South Pole
}

/*****************************************************************************/
//	Spherical Gray-Scott Reaction-Diffusion Brain Coral simulation is licensed under
//	MIT License
//
//	Copyright (c) 2025 DrT0M
//
//	Permission is hereby granted, free of charge, to any person obtaining a copy
//	of this software and associated documentation files (the "Software"), to deal
//	in the Software without restriction, including without limitation the rights
//	to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
//	copies of the Software, and to permit persons to whom the Software is
//	furnished to do so, subject to the following conditions:
//	
//	The above copyright notice and this permission notice shall be included in all
//	copies or substantial portions of the Software.
//
//	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//	IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//	FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//	AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//	LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//	OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
//	SOFTWARE.
//

/*****************************************************************************/
double	F	= 0.0256;	//feed rate
double	K	= 0.0544;	//kill rate
double	Dt	= 1.0;		//time step
double	Du	= 1.0;		//diffusion constant | rate | coefficient
double	Dv	= 0.5;		//diffusion constant | rate | coefficient
double	Scale	= 0.2;		//1:5 fine scale
double	Feature	= 1/32.;
double	SphArea = 0;
double	AvgArea = 0;
double	MinArea = 0;
double	MaxArea = 0;

/*****************************************************************************/
void init_points()
{
	for( int k = 0; k < POINTS ; k++ )
	{
		Point*p = &(point[k]);
		p->area = 0.0;
		p->cotan = 0.0;
		p->u = 1.0;
		p->v = 0.0;
		p->NPcnt = 0;	p->NFcnt = 0;
		p->NP[0] = -1;	p->NF[0] = -1;
		p->NP[1] = -1;	p->NF[1] = -1;
		p->NP[2] = -1;	p->NF[2] = -1;
		p->NP[3] = -1;	p->NF[3] = -1;
		p->NP[4] = -1;	p->NF[4] = -1;
		p->NP[5] = -1;	p->NF[5] = -1;
	}
}

/*****************************************************************************/
void init_neighbor( int q, int p, int f )
{
	int k = 0 ;
	for( k = 0 ; k < point[q].NFcnt ; k++ )
	if( point[q].NF[ k ] == f )
		break;

	if( point[q].NF[ k ] != f )
	if( point[q].NFcnt < 6 )
	{
		point[q].NF[ point[q].NFcnt ] = f;
		point[q].NFcnt++;
	}

	k = 0 ;
	for( k = 0 ; k < point[q].NPcnt ; k++ )
	if( point[q].NP[ k ] == p )
		break;

	if( point[q].NP[ k ] != p )
	if( point[q].NPcnt < 6 )
	{
		point[q].NP[ point[q].NPcnt ] = p;
		point[q].NPcnt++;
	}
}

/*****************************************************************************/
void init_neighbors( int F, int A, int B, int C )
{
	init_neighbor( A, B, F );	init_neighbor( B, A, F );
	init_neighbor( B, C, F );	init_neighbor( C, B, F );
	init_neighbor( C, A, F );	init_neighbor( A, C, F );
}

/*****************************************************************************/
void init_facets()
{
	for( int k = 0; k < FACETS ; k++ )
	{
		Facet*f = &(facet[k]);
		double	difx, dify, difz, tmp1;

		init_neighbors( k, f->A, f->B, f->C );

		//square of side length a for segment BC =(xB-xC)^2+(yB-yC)^2+(zB-zC)^2
		difx = point[f->B].x - point[f->C].x;
		dify = point[f->B].y - point[f->C].y;
		difz = point[f->B].z - point[f->C].z;
		f->a2 = difx*difx + dify*dify + difz*difz;

		//square of side length b for segment CA =(xC-xA)^2+(yC-yA)^2+(zC-zA)^2
		difx = point[f->C].x - point[f->A].x;
		dify = point[f->C].y - point[f->A].y;
		difz = point[f->C].z - point[f->A].z;
		f->b2 = difx*difx + dify*dify + difz*difz;

		//square of side length c for segment AB =(xA-xB)^2+(yA-yB)^2+(zA-zB)^2
		difx = point[f->A].x - point[f->B].x;
		dify = point[f->A].y - point[f->B].y;
		difz = point[f->A].z - point[f->B].z;
		f->c2 = difx*difx + dify*dify + difz*difz;

		//area of triangle ABC =SQRT((a2+b2+c2)^2-2*(a2^2+b2^2+c2^2))/4
		tmp1 = f->a2 + f->b2 + f->c2;
		f->T = sqrt(tmp1*tmp1 - 2*
				(f->a2*f->a2
				+f->b2*f->b2
				+f->c2*f->c2))/4;

		//cotangent of angle CAB =(-a2+b2+c2)/4/T
		f->cot_alpha = (-f->a2 +f->b2 +f->c2)/4/f->T;

		//cotangent of angle ABC =(+a2-b2+c2)/4/T
		f->cot_beta_ = (+f->a2 -f->b2 +f->c2)/4/f->T;

		//cotangent of angle BCA =(+a2+b2-c2)/4/T
		f->cot_gamma = (+f->a2 +f->b2 -f->c2)/4/f->T;

		//contribution to the vertex area of A =0.5*(T-0.25*a2*(-a2+b2+c2)/4/T)
		point[f->A].area +=  0.5*(f->T - 0.25*f->a2*(-f->a2 +f->b2 +f->c2)/4/f->T);

		//contribution to the vertex area of B =0.5*(T-0.25*b2*(+a2-b2+c2)/4/T)
		point[f->B].area +=  0.5*(f->T - 0.25*f->b2*(+f->a2 -f->b2 +f->c2)/4/f->T);

		//contribution to the vertex area of C =0.5*(T-0.25*c2*(+a2+b2-c2)/4/T)
		point[f->C].area +=  0.5*(f->T - 0.25*f->c2*(+f->a2 +f->b2 -f->c2)/4/f->T);

		//negative of twice the diagonal of cotangent matrix (-2*Cii)
		point[f->A].cotan += f->cot_gamma;	point[f->A].cotan += f->cot_beta_;
		point[f->B].cotan += f->cot_alpha;	point[f->B].cotan += f->cot_gamma;
		point[f->C].cotan += f->cot_beta_;	point[f->C].cotan += f->cot_alpha;
	}
}

/*****************************************************************************/
void init_sphere()
{
	for( int k = 0; k < POINTS ; k++ )
		SphArea += point[k].area;

	AvgArea = SphArea/POINTS;
	MinArea = SphArea;
	MaxArea = 0;

	for( int k = 0; k < POINTS ; k++ )
	{
		if( MaxArea < point[k].area )
		    MaxArea = point[k].area;
		if( MinArea > point[k].area )
		    MinArea = point[k].area;
	}

	//area | mass normalization
	for( int k = 0; k < POINTS ; k++ )
		point[k].norm = point[k].area / MinArea;
}

/*****************************************************************************/
void init_area_v()
{
#define INIT_CENTER 0		//initialization point
	if( INIT_CENTER >= 0 )
	{
		point[ INIT_CENTER ].v = 1.0;
		for( int n = 0; n < point[ INIT_CENTER ].NPcnt ; n++ )
			point[ point[ INIT_CENTER ].NP[ n ]].v = 1.0;
	}

#define LAST_CENTER -1		//second initialization point if non-negative, e.g. POINTS-1
	if( LAST_CENTER >= 0 )
	{
		point[ LAST_CENTER ].v = 1.0;
		for( int n = 0; n < point[ LAST_CENTER ].NPcnt ; n++ )
			point[ point[ LAST_CENTER ].NP[ n ]].v = 1.0;
	}
}

/*****************************************************************************/
void update(void*t)
{
	for( int k = (long)t; k < POINTS ; k+=N_THRDS )
	{
	    for( int n = 0; n < point[k].NFcnt ; n++ )
	    {
		Facet*f = &(facet[ point[k].NF[n] ]);
		if( f->A == k )
		{
			point[k].Lu += (point[f->B].u - point[f->A].u) * f->cot_gamma;
			point[k].Lu -= (point[f->A].u - point[f->C].u) * f->cot_beta_;
			point[k].Lv += (point[f->B].v - point[f->A].v) * f->cot_gamma;
			point[k].Lv -= (point[f->A].v - point[f->C].v) * f->cot_beta_;
		}
		else
		if( f->B == k )
		{
			point[k].Lu += (point[f->C].u - point[f->B].u) * f->cot_alpha;
			point[k].Lu -= (point[f->B].u - point[f->A].u) * f->cot_gamma;
			point[k].Lv += (point[f->C].v - point[f->B].v) * f->cot_alpha;
			point[k].Lv -= (point[f->B].v - point[f->A].v) * f->cot_gamma;
		}
		else
		if( f->C == k )
		{
			point[k].Lu += (point[f->A].u - point[f->C].u) * f->cot_beta_;
			point[k].Lu -= (point[f->C].u - point[f->B].u) * f->cot_alpha;
			point[k].Lv += (point[f->A].v - point[f->C].v) * f->cot_beta_;
			point[k].Lv -= (point[f->C].v - point[f->B].v) * f->cot_alpha;
		}
	    }
	}
}

/*****************************************************************************/
void update_single_threaded()
{
	for ( long t = 0 ; t < N_THRDS ; t++ )
		update( (void*)t );
}

/*****************************************************************************/
void*update_multi_threaded(void*t)
{
	update(t);
	pthread_exit(NULL);
}

/*****************************************************************************/
void iter_8(void*t)
{
	// L*u = M\C*u = diag(1/2/Ai)*2C*u
	for( int k = (long)t; k < POINTS ; k+=N_THRDS )
	{
		Point*p = &(point[k]);

		//area | mass normalization (the smaller the vertex area, the faster the change [in temparature | concentration])

		//-1 at center of stencil
		p->Lu = (p->Lu / p->cotan);
		p->Lv = (p->Lv / p->cotan);

		//area | mass normalization
		p->Lu /= p->norm;
		p->Lv /= p->norm;

	}

	for( int k = (long)t; k < POINTS ; k+=N_THRDS )
	{
		Point*p = &(point[k]);

		//two-component reaction-diffusion system
		double	uv2 = p->u * p->v * p->v;
		p->u += (Scale*Du*p->Lu - uv2 + F*(1 - p->u))*Dt;
		p->v += (Scale*Dv*p->Lv + uv2 - (K + F)*p->v)*Dt;

		//u, v bounded to [0,1]
		if( p->u < 0. )	p->u = 0.;	if( p->u > 1. )	p->u = 1.;
		if( p->v < 0. )	p->v = 0.;	if( p->v > 1. )	p->v = 1.;
	}

	for( int k = (long)t; k < POINTS ; k+=N_THRDS )
	{
		Point*p = &(point[k]);
		p->Lu = 0.0;
		p->Lv = 0.0;
	}
}

/*****************************************************************************/
void iter_8_single_threaded()
{
	for ( long t = 0 ; t < N_THRDS ; t++ )
		iter_8( (void*)t );
}

/*****************************************************************************/
void*iter_8_multi_threaded(void*t)
{
	iter_8(t);
	pthread_exit(NULL);
}

/*****************************************************************************/
void pthreaded_run(pthread_t *threads, const pthread_attr_t *attr, void *(*start_routine) (void *))
{
	for ( long t = 0 ; t < N_THRDS ; t++ )		//launch threads
	{
		int rc = pthread_create(&threads[ t ], NULL, start_routine, (void*)t);
		if( rc )
		{
			fprintf( stderr, "ERROR: thread %ld return code from pthread_create() is %d\n", t, rc );
			fflush( stderr );
			exit(-1);
		}
	}

	for ( long t = 0 ; t < N_THRDS ; t++ )		//join all threads
	{
		int rc = pthread_join( threads[ t ], NULL );
		if( rc && rc != ESRCH )
		{
			(void)	fprintf(stderr,	"0:pthread_join(threads[%ld])->%d\n", t , rc );
			fflush( stderr );
		}
	}
}

/*****************************************************************************/
typedef unsigned int Color24;

Color24 setRGB(double R1, double G1, double B1)
{
	return	(((unsigned int)(R1*255)%256)<<16)+
		(((unsigned int)(G1*255)%256)<<8 )+
		(((unsigned int)(B1*255)%256)    );
}

/*****************************************************************************/
int main(int argc, char*argv[])
{
	init_icosahedron();
	geodesic_mesh();

	/*****************************************************************************/
	if( argc > 2 )	F = atof( argv[2] );
	if( argc > 3 )	K = atof( argv[3] );
	if( argc > 4 )	Scale = atof( argv[4] );

	/*****************************************************************************/
	init_points();
	init_facets();
	init_sphere();
	init_area_v();

	/*****************************************************************************/
#define RADIUS 125
#define PIXELS 128
	//X11 graphics
	unsigned int width  = PIXELS * 4;
	unsigned int height = PIXELS * 2;
	Display*display = XOpenDisplay(NULL);
	Window  window  = XCreateSimpleWindow(display, RootWindow(display, 0), 0, 0, width, height, 1, 0, 0);
	Visual *visual  = DefaultVisual(display, 0);
	if(visual->class != TrueColor)
	{
		fprintf(stderr, "Cannot handle non-true color visual ...\n");
		exit(1);
	}
	char*image32 = (char*) calloc(width*height*4, sizeof(char));
	XImage *ximage = XCreateImage(display, visual, 24, ZPixmap, 0, image32, width, height, 32, 0);
	XStoreName(display, window, argv[0]);
	XSelectInput(display, window, ButtonPressMask | ExposureMask | StructureNotifyMask);
	XMapWindow(display, window);
	for(;;)
	{
		XEvent ev;
		XNextEvent(display, &ev);
		if(ev.type == MapNotify)	break;
	}
	//Wait for the MapNotify event
	//which means that the window has appeared on the screen.
	//initial black screen
	XPutImage(display, window, DefaultGC(display, 0), ximage, 0, 0, 0, 0, width, height);
	XFlush(display);

	/*****************************************************************************/
	//initial u, v point cloud on X11 display
	for( int k = 0; k < POINTS ; k++ )
	{
		Point*p = &(point[k]);
		int v = p->z * RADIUS;
		int h = p->y * RADIUS;
		XPutPixel(ximage,( p->x > 0 ? 3*PIXELS+h : PIXELS-h ), PIXELS-v, setRGB(p->u, p->u, p->v));
	}
	XStoreName(display, window, "Iteration 0");
	XPutImage(display, window, DefaultGC(display, 0), ximage, 0, 0, 0, 0, width, height);
	XFlush(display);

	/*****************************************************************************/
	//key parameters shown in terminal 
	printf("POINTS=%d, FACETS=%d, ", 
		POINTS   , FACETS);
	printf("Dt=%f, Du=%f, Dv=%f, F=%f, K=%f, Scale=%f\n", 
		Dt   , Du   , Dv   , F   , K   , Scale);
	printf("Iterating:\t");
	(void)fflush(stdout);

	/*****************************************************************************/
	pthread_t*threads = calloc( (size_t)N_THRDS, sizeof(pthread_t) );	//create thread array

	for( int i = 1 ; i <=( argc > 1 ? atoi( argv[1] ) : 20000 ) ; i++ )
	{
		pthreaded_run(threads, NULL, update_multi_threaded);	//	update_single_threaded();
		pthreaded_run(threads, NULL, iter_8_multi_threaded);	//	iter_8_single_threaded();

		//progress shown in terminal 
		if( !(i % 1000) )
		{
			printf("%d\t", i);
			(void)fflush(stdout);
		}

		//progress shown in X11 window
		if( !(i % 50) )
		{
			char counter[100];
			(void)	sprintf( counter, "Iteration %d", i );
			XStoreName(display, window, counter);
			for( int k = 0; k < POINTS ; k++ )
			{
				Point*p = &(point[k]);
				int v = p->z * RADIUS;
				int h = p->y * RADIUS;
				XPutPixel(ximage,( p->x > 0 ? 3*PIXELS+h : PIXELS-h ), PIXELS-v, setRGB(p->u, p->u, p->v));
			}
			XPutImage(display, window, DefaultGC(display, 0), ximage, 0, 0, 0, 0, width, height);
			XFlush(display);
		}
	}

	(void) free( threads );

	/*****************************************************************************/
	FILE*S = fopen(	"brain_coral_spherical_gsrd.scad","w");
	(void)fprintf(S,"scale(25)\npolyhedron(Spherical_GSRD_points(), Spherical_GSRD_facets());\n");
	(void)fprintf(S,"function Spherical_GSRD_points() =\n");
	for( int k = 0; k < POINTS ; k++ )
	{
		Point*p = &(point[k]);
		double	p__u = p->u*Feature + (1.0-Feature);
		(void)fprintf(S,"%c[%f\t,%f\t,%f\t]\n", (k==0 ? '[' : ','), p__u*p->x, p__u*p->y, p__u*p->z);
	}
	(void)fprintf(S,"];\n");
	(void)fprintf(S,"function Spherical_GSRD_facets() =\n");
	for( int k = 0; k < FACETS ; k++ )
	{
		Facet*f = &(facet[k]);
		(void)fprintf(S,"%c[%d\t,%d\t,%d\t]\n", (k==0 ? '[' : ','), f->A, f->B, f->C);
	}
	(void)fprintf(S,"];\n");
	(void)fclose( S );

	/*****************************************************************************/
	char keyboard_input[100];
	printf("\nPress enter to close graphics display.\n");
	fgets (keyboard_input,100,stdin); 
}

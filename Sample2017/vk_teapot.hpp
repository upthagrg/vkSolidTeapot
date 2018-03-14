#include <cmath>
#include <iostream>
#include "glm/vec2.hpp"
#include "glm/vec3.hpp"
#include "vk_teapot_data.h"

#include "vertex.hpp"

/* -- STATIC VARS: CACHES ---------------------------------------------------- */

/* General defs */
#define GLUT_SOLID_N_SUBDIV  8
#define GLUT_WIRE_N_SUBDIV   10

/* Bernstein coefficients only have to be precomputed once (number of patch subdivisions is fixed)
 * Can thus define arrays for them here, they will be filled upon first use.
 * 3rd order Bezier surfaces have 4 Bernstein coeffs.
 * Have separate caches for solid and wire as they use a different number of subdivisions
 * _0 is for Bernstein polynomials, _1 for their first derivative (which we need for normals)
 */
static float bernWire_0 [GLUT_WIRE_N_SUBDIV] [4];
static float bernWire_1 [GLUT_WIRE_N_SUBDIV] [4];
static float bernSolid_0[GLUT_SOLID_N_SUBDIV][4];
static float bernSolid_1[GLUT_SOLID_N_SUBDIV][4];

/* Teapot defs */
#define GLUT_TEAPOT_N_PATCHES       (6*4 + 4*2)                                                                     /* 6 patches are reproduced (rotated) 4 times, 4 patches (flipped) 2 times */
#define GLUT_SOLID_TEAPOT_N_VERT    GLUT_SOLID_N_SUBDIV*GLUT_SOLID_N_SUBDIV * GLUT_TEAPOT_N_PATCHES                 /* N_SUBDIV^2 vertices per patch */
#define GLUT_SOLID_TEAPOT_N_TRI     (GLUT_SOLID_N_SUBDIV-1)*(GLUT_SOLID_N_SUBDIV-1) * GLUT_TEAPOT_N_PATCHES * 2     /* if e.g. 7x7 vertices for each patch, there are 6*6 squares for each patch. Each square is decomposed into 2 triangles */

#define GLUT_WIRE_TEAPOT_N_VERT     GLUT_WIRE_N_SUBDIV*GLUT_WIRE_N_SUBDIV * GLUT_TEAPOT_N_PATCHES                   /* N_SUBDIV^2 vertices per patch */

/* Bit of caching:
 * vertex indices and normals only need to be generated once for
 * a given number of subdivisions as they don't change with scale.
 * Vertices can be cached and reused if scale didn't change.
 */
static unsigned short vertIdxsTeapotS[GLUT_SOLID_TEAPOT_N_TRI*3];
static float  normsTeapotS   [GLUT_SOLID_TEAPOT_N_VERT*3];
static float  vertsTeapotS   [GLUT_SOLID_TEAPOT_N_VERT*3];
static float  texcsTeapotS   [GLUT_SOLID_TEAPOT_N_VERT*2];
static float  lastScaleTeapotS = 0.f;
static bool initedTeapotS   = false;

static unsigned short vertIdxsTeapotW[GLUT_WIRE_TEAPOT_N_VERT*2];
static float  normsTeapotW   [GLUT_WIRE_TEAPOT_N_VERT*3];
static float  vertsTeapotW   [GLUT_WIRE_TEAPOT_N_VERT*3];
static float  lastScaleTeapotW = 0.f;
static bool initedTeapotW   = false;

int vkTeapotSize=0;

glm::vec3 cur_color = glm::vec3(1.0, 0.0, 0.0);

// *******************
// VERTEX BUFFER DATA:
// *******************

//#ifdef DEFINITION
/*
struct vertex
{
	glm::vec3	position;
	glm::vec3	normal;
	glm::vec3	color;
	glm::vec2	texCoord;
};
//#endif
*/

static void bernstein3(int i, float x, float *r0, float *r1)
{
    float invx = 1.f - x;

    /* r0: zero order coeff, r1: first deriv coeff */
    switch (i)
    {
        float temp;
    case 0:
        temp = invx*invx;
        *r0 = invx * temp;                  /* invx * invx * invx */
        *r1 = -3 * temp;                    /*   -3 * invx * invx */
        break;
    case 1:
        temp = invx*invx;
        *r0 = 3 * x * temp;                 /* 3 * x * invx * invx */
        *r1 = 3 * temp  -  6 * x * invx;    /* 3 * invx * invx  -  6 * x * invx */
        break;
    case 2:
        temp = x*x;
        *r0 = 3 * temp * invx;              /* 3 * x * x * invx */
        *r1 = 6 * x * invx  -  3 * temp;    /* 6 * x * invx  -  3 * x * x */
        break;
    case 3:
        temp = x*x;
        *r0 = x * temp;                     /* x * x * x */
        *r1 = 3 * temp;                     /* 3 * x * x */
        break;
    default:
        *r0 = *r1 = 0;
    }
}


static void pregenBernstein(int nSubDivs, float (*bern_0)[4], float (*bern_1)[4])
{
    int s,i;
    for (s=0; s<nSubDivs; s++)
    {
        float x = s/(nSubDivs-1.f);
        for (i=0; i<4; i++) /* 3rd order polynomial */
            bernstein3(i,x,bern_0[s]+i,bern_1[s]+i);
    }
}


/* based on flag either rotate patches around y axis to other 3 quadrants (flag=4) or reflect patch across x-y plane (flag=2) */
static void rotOrReflect(int flag, int nVals, int nSubDivs, float *vals)
{
    int u,i,o;

    if (flag==4)
    {
        int i1=nVals, i2=nVals*2, i3=nVals*3;
        for (o=0; o<nVals; o+=3)
        {
            /* 90° rotation */
            vals[i1+o+0] =  vals[o+2];
            vals[i1+o+1] =  vals[o+1];
            vals[i1+o+2] = -vals[o+0];
            /* 180° rotation */
            vals[i2+o+0] = -vals[o+0];
            vals[i2+o+1] =  vals[o+1];
            vals[i2+o+2] = -vals[o+2];
            /* 270° rotation */
            vals[i3+o+0] = -vals[o+2];
            vals[i3+o+1] =  vals[o+1];
            vals[i3+o+2] =  vals[o+0];
        }
    }
    else if (flag==2)
    {
        /* copy over values, reversing row order to keep winding correct, and negating z to perform the flip */
        for (u=0; u<nSubDivs; u++)  /* per row */
        {
            int off =   (nSubDivs-u-1)*nSubDivs*3;  /* read last row first from the already existing rows */
            o       = nVals + u   *nSubDivs*3;      /* write last row as first row to output */
            for (i=0; i<nSubDivs*3; i+=3, o+=3)     /* each row has nSubDivs points consisting of three values */
            {
                vals[o+0] =  vals[off+i+0];
                vals[o+1] =  vals[off+i+1];
                vals[o+2] = -vals[off+i+2];
            }
        }
    }
}


/* verts array should be initialized to 0! */
static int evalBezierWithNorm(float cp[4][4][3], int nSubDivs, float (*bern_0)[4], float (*bern_1)[4], int flag, int normalFix, float *verts, float *norms)
{
    int nVerts    = nSubDivs*nSubDivs;
    int nVertVals = nVerts*3;               /* number of values output for one patch, flag (2 or 4) indicates how many times we will write this to output */
    int u,v,i,j,o;
//int cnt = 0;
    /* generate vertices and coordinates for the patch */
    for (u=0,o=0; u<nSubDivs; u++)
    {
        for (v=0; v<nSubDivs; v++, o+=3)
        {
            /* for normals, get two tangents at the vertex using partial derivatives of 2D Bezier grid */
            float tan1[3]={0}, tan2[3]={0}, len;
            for (i=0; i<=3; i++)
            {
                float vert_0[3]={0}, vert_1[3]={0};
                for (j=0; j<=3; j++)
                {
                    vert_0[0] += bern_0[v][j] * cp[i][j][0];
                    vert_0[1] += bern_0[v][j] * cp[i][j][1];
                    vert_0[2] += bern_0[v][j] * cp[i][j][2];

                    vert_1[0] += bern_1[v][j] * cp[i][j][0];
                    vert_1[1] += bern_1[v][j] * cp[i][j][1];
                    vert_1[2] += bern_1[v][j] * cp[i][j][2];
                }

                verts[o+0] += bern_0[u][i]*vert_0[0];
				//cnt++;
				//std::cout<<"vert: " <<verts[o+0]<<std::endl;
                verts[o+1] += bern_0[u][i]*vert_0[1];
				//cnt++;
				//std::cout<<"vert: " <<verts[o+1]<<std::endl;
                verts[o+2] += bern_0[u][i]*vert_0[2];
				//cnt++;
				//std::cout<<"vert: " <<verts[o+2]<<std::endl;
				

                tan1[0] += bern_0[u][i]*vert_1[0];
                tan1[1] += bern_0[u][i]*vert_1[1];
                tan1[2] += bern_0[u][i]*vert_1[2];
                tan2[0] += bern_1[u][i]*vert_0[0];
                tan2[1] += bern_1[u][i]*vert_0[1];
                tan2[2] += bern_1[u][i]*vert_0[2];
            }
            /* get normal through cross product of the two tangents of the vertex */
            norms[o+0] = tan1[1] * tan2[2] - tan1[2] * tan2[1];
            norms[o+1] = tan1[2] * tan2[0] - tan1[0] * tan2[2];
            norms[o+2] = tan1[0] * tan2[1] - tan1[1] * tan2[0];
            len = (float)sqrt(norms[o+0] * norms[o+0] + norms[o+1] * norms[o+1] + norms[o+2] * norms[o+2]);
            norms[o+0] /= len;
            norms[o+1] /= len;
            norms[o+2] /= len;
        }
    }
//std::cout << "cnt: " <<cnt<<std::endl;
    /* Fix normal vector if needed */
    if (normalFix)
    {
        for (o=0; o<nSubDivs*3; o+=3) /* whole first row (first nSubDivs normals) is broken: replace normals for the whole row */
        {
            norms[o+0] = 0.f;
            norms[o+1] = normalFix==1? 1.f:-1.f;
            norms[o+2] = 0.f;
        }
    }

    /* now based on flag either rotate patches around y axis to other 3 quadrants (flag=4) or reflect patch across x-y plane (flag=2) */
    rotOrReflect(flag, nVertVals, nSubDivs, verts);
    rotOrReflect(flag, nVertVals, nSubDivs, norms);

    return nVertVals*flag;
}



static void fghTeaset( double scale, bool useWireMode,
                       float (*cpdata)[3], int (*patchdata)[16],
                       unsigned short *vertIdxs,
                       float *verts, float *norms, float *texcs,
                       float *lastScale, bool *inited,
                       bool needNormalFix, bool rotFlip, float zOffset,
                       int nVerts, int nInputPatches, int nPatches, int nTriangles )
{
    /* for internal use */
    int p,o;
    float cp[4][4][3];
    /* to hold pointers to static vars/arrays */
    float (*bern_0)[4], (*bern_1)[4];
    int nSubDivs;

    /* Get relevant static arrays and variables */
    bern_0      = useWireMode ? bernWire_0                : bernSolid_0;
    bern_1      = useWireMode ? bernWire_1                : bernSolid_1;
    nSubDivs    = useWireMode ? GLUT_WIRE_N_SUBDIV        : GLUT_SOLID_N_SUBDIV;

    /* check if need to generate vertices */
    if (!*inited || scale != *lastScale)
    {
        /* set vertex array to all 0 (not necessary for normals and vertex indices) */
        memset(verts,0,nVerts*3*sizeof(float));

        /* pregen Berstein polynomials and their first derivatives (for normals) */
        if (!*inited)
            pregenBernstein(nSubDivs,bern_0,bern_1);

        /* generate vertices and normals */
        for (p=0, o=0; p<nInputPatches; p++)
        {
            /* set flags for evalBezier function */
            int flag      = rotFlip?p<6?4:2:1;                  /* For teapot and teacup, first six patches get 3 copies (rotations), others get 2 copies (flips). No rotating or flipping at all for teaspoon */
            int normalFix = needNormalFix?p==3?1:p==5?2:0:0;    /* For teapot, fix normal vectors for vertices on top of lid (patch 4) and on middle of bottom (patch 6). Different flag value as different normal needed */

            /* collect control points */
            int i;
            for (i=0; i<16; i++)
            {
                /* Original code draws with a 270° rot around X axis, a scaling and a translation along the Z-axis.
                 * Incorporating these in the control points is much cheaper than transforming all the vertices.
                 * Original:
                 * glRotated( 270.0, 1.0, 0.0, 0.0 );
                 * glScaled( 0.5 * scale, 0.5 * scale, 0.5 * scale );
                 * glTranslated( 0.0, 0.0, -zOffset );  -> was 1.5 for teapot, but should be 1.575 to center it on the Z axis. Teacup and teaspoon have different offsets
                 */
                cp[i/4][i%4][0] =  cpdata[patchdata[p][i]][0]         *scale/2.f;
                cp[i/4][i%4][1] = (cpdata[patchdata[p][i]][2]-zOffset)*scale/2.f;
                cp[i/4][i%4][2] = -cpdata[patchdata[p][i]][1]         *scale/2.f;
            }

            /* eval bezier patch */
            //if (!*inited)   /* first time, generate normals as well */
                o += evalBezierWithNorm(cp,nSubDivs,bern_0,bern_1, flag, normalFix, verts+o,norms+o);
            //else            /* only need to regen vertices */
            //    o += evalBezier(cp,nSubDivs,bern_0, flag, verts+o);
        }
        *lastScale = scale;

        if (!*inited)
        {
            int r,c;
            /* generate texture coordinates if solid teapot/teacup/teaspoon */
            if (!useWireMode)
            {
                /* generate for first patch */
                for (r=0,o=0; r<nSubDivs; r++)
                {
                    float u = r/(nSubDivs-1.f);
                    for (c=0; c<nSubDivs; c++, o+=2)
                    {
                        float v = c/(nSubDivs-1.f);
                        texcs[o+0] = u;
                        texcs[o+1] = v;
                    }
                }
                /* copy it over for all the other patches */
                for (p=1; p<nPatches; p++)
                    memcpy(texcs+p*nSubDivs*nSubDivs*2,texcs,nSubDivs*nSubDivs*2*sizeof(float));
            }

            /* build vertex index array */
           // if (useWireMode)
            //{
            //}
            //else
            //{
                /* build vertex indices to draw teapot/teacup/teaspoon as triangles */
                for (p=0,o=0; p<nPatches; p++)
                {
                    int idx = nSubDivs*nSubDivs*p;
                    for (r=0; r<nSubDivs-1; r++)
                    {
                        int loc = r*nSubDivs;
                        for (c=0; c<nSubDivs-1; c++, o+=6)
                        {
                            /* ABC ACD, where B and C are one row lower */
                            int row1 = idx+loc+c;
                            int row2 = row1+nSubDivs;

                            vertIdxs[o+0] = row1+0;
                            vertIdxs[o+1] = row2+0;
                            vertIdxs[o+2] = row2+1;

                            vertIdxs[o+3] = row1+0;
                            vertIdxs[o+4] = row2+1;
                            vertIdxs[o+5] = row1+1;
                        }
                    }
                }
            //}

            *inited = true;
        }
    }

    /* draw */
	/*
    if (useWireMode)
        fghDrawGeometryWire (verts, norms,        nVerts, vertIdxs, nPatches*nSubDivs*2, nSubDivs, GL_LINE_STRIP, NULL,0,0);
    else
        fghDrawGeometrySolid(verts, norms, texcs, nVerts, vertIdxs,1,nTriangles*3);
	*/
}

void set_teapot_color(float r, float g, float b){
	cur_color = glm::vec3(r, g, b);
}

struct vertex* vkSolidTeapot( double size )
{
    //FREEGLUT_EXIT_IF_NOT_INITIALISED ( "glutSolidTeapot" );
    fghTeaset( size, false,
               cpdata_teapot, patchdata_teapot,
               vertIdxsTeapotS,
               vertsTeapotS, normsTeapotS, texcsTeapotS,
               &lastScaleTeapotS, &initedTeapotS,
               true, true, 1.575f,
               GLUT_SOLID_TEAPOT_N_VERT, GLUT_TEAPOT_N_INPUT_PATCHES, GLUT_TEAPOT_N_PATCHES, GLUT_SOLID_TEAPOT_N_TRI);
	//struct vertex
	//struct vertex* ret = new struct vertex[GLUT_SOLID_TEAPOT_N_VERT];
	//struct vertex* ret = new struct vertex[GLUT_SOLID_TEAPOT_N_TRI*3];
	struct vertex* ret = new struct vertex[9408];
	vkTeapotSize = 9408 * sizeof(struct vertex);
	//vkTeapotSize = (GLUT_SOLID_TEAPOT_N_VERT) * sizeof(struct vertex);
	//vkTeapotSize = (GLUT_SOLID_TEAPOT_N_TRI*3) * sizeof(struct vertex);
	//int track = 0;
	int track2 = 0;
	//for(int i=0; i<(GLUT_SOLID_TEAPOT_N_VERT)-2; i++ ){
	//for(int i=0; i<(GLUT_SOLID_TEAPOT_N_TRI*3)-2; i++){
	//	(ret[i]).position = glm::vec3(vertsTeapotS[track], vertsTeapotS[track+1], vertsTeapotS[track+3]);
	//	std::cout<<"point " <<i<< ": " <<vertsTeapotS[track]<<"," <<vertsTeapotS[track+1]<<","<< vertIdxsTeapotS[track+2]<<std::endl;
	//	(ret[i]).normal = glm::vec3(normsTeapotS[track], normsTeapotS[track+1], normsTeapotS[track+2]);
	//	track+=3;
		//(ret[i]).texCoord = glm::vec2(texcsTeapotS[track2], texcsTeapotS[track2+1]);
		//track2+=2;
	//	(ret[i]).color = cur_color;
		
	//}
	int track=0;
	//for(int i=0; i<2048; i++){
	//	(ret[i]).position = glm::vec3(
	//}
	//for(int i=0; i<(GLUT_SOLID_TEAPOT_N_VERT)-2; i++){
	for(int i=0; i<(GLUT_SOLID_TEAPOT_N_TRI*3); i++){
	//	for(int i=0; i<(6); i++){
			//(ret[i]).position = glm::vec3(vertsTeapotS[track], vertsTeapotS[track+1], vertsTeapotS[track+3]);
			//std::cout<<"vertIdxsTeapotS["<<i<<"]: "<<vertIdxsTeapotS[i]<<std::endl;
			(ret[i]).position = glm::vec3(vertsTeapotS[((vertIdxsTeapotS[i] * 3))], vertsTeapotS[((vertIdxsTeapotS[i] * 3)+1)], vertsTeapotS[((vertIdxsTeapotS[i] * 3)+2)]);
			(ret[i]).normal = glm::vec3(normsTeapotS[((vertIdxsTeapotS[i] * 3))], normsTeapotS[((vertIdxsTeapotS[i] * 3)+1)], normsTeapotS[((vertIdxsTeapotS[i] * 3)+2)]);
			(ret[i]).texCoord = glm::vec2(0.0,0.0);
			//std::cout << "pos: " << (ret[i]).position[0] << " " << (ret[i]).position[1] << " " << (ret[i]).position[2] << std::endl;
			(ret[i]).color = cur_color;
			//std::cout <<"i: " << i <<std::endl;
			
			
			
			//(ret[i]).position = glm::vec3(vertsTeapotS[track], vertsTeapotS[track+1], vertsTeapotS[track+2]);
			//(ret[i]).position = glm::vec3(normsTeapotS[((vertIdxsTeapotS[i] * 3))], normsTeapotS[((vertIdxsTeapotS[i] * 3)+1)], normsTeapotS[((vertIdxsTeapotS[i] * 3)+2)]);
			//(ret[i]).texCoord = glm::vec2(0.0,0.0);
			//std::cout << "pos: " << (ret[i]).position[0] << " " << (ret[i]).position[1] << " " << (ret[i]).position[2] << std::endl;
			//(ret[i]).color = cur_color;
			
	}
	//for(int i=0; i<(GLUT_SOLID_TEAPOT_N_TRI*3); i++){
	//		std::cout << "pos: "(ret[i]).position[0]<<" "<<(ret[i]).position[1]<<" "<<(ret[i]).position[2]<<std::endl;
	//		std::cout << "norm: "(ret[i]).normal[0]<<" "<<(ret[i]).normal[1]<<" "<<(ret[i]).normal[2]<<std::endl;
	//}
	//for(int i=0; i<6144; i++){
	//		std::cout<<"vert["<<i<<"]: "<<vertsTeapotS[i]<<std::endl;
	//}
	return ret;
}

int vkTeapotGetSize(){
	return vkTeapotSize;
}
				

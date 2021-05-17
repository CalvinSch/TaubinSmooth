
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <vector>
#include <set>
using namespace std;

#define ITERATIONS 100


// Parameters
float lambda = 0.3;	//Taubin parameter lambda
float mu = -0.2;	//Taubin parameter mu
char* fileName = "2CMP_noise.off";	//input OFF mesh file


typedef struct {
	float x;
	float y;
	float z;
}FLTVECT;

typedef struct {
	int a;
	int b;
	int c;
}INT3VECT;

typedef struct {
	int nv;
	int nf;
	FLTVECT *vertex;
	INT3VECT *face;
}SurFacemesh;

typedef struct Node{
    void *data;
    struct Node *next;
}llnode_t;

// Surface mesh obtained from .off file
SurFacemesh* surfmesh;

FLTVECT operator*(const float f, const FLTVECT u)
{
	FLTVECT ret;
	ret.x = f * u.x;
	ret.y = f * u.y;
	ret.z = f * u.z;
	return ret;
}

FLTVECT operator+(const FLTVECT v, const FLTVECT u)
{
	FLTVECT ret;
	ret.x = v.x + u.x;
	ret.y = v.y + u.y;
	ret.z = v.z + u.z;
	return ret;
}

FLTVECT operator-(const FLTVECT v, const FLTVECT u)
{
    FLTVECT ret;
    ret.x = v.x - u.x;
    ret.y = v.y - u.y;
    ret.z = v.z - u.z;
    return ret;

}

void readPolygon()
{
	int num, n, m;
	int a, b, c, d;
	float x, y, z;
	char line[256];
	FILE *fin;

	fin = fopen(fileName, "r");
	/* OFF format */
	while (fgets(line, 256, fin) != NULL) {
		if (line[0] == 'O' && line[1] == 'F' && line[2] == 'F')
			break;
	}
	fscanf(fin, "%d %d %d\n", &m, &n, &num);

	surfmesh = (SurFacemesh*)malloc(sizeof(SurFacemesh));
	surfmesh->nv = m;
	surfmesh->nf = n;
	surfmesh->vertex = (FLTVECT *)malloc(sizeof(FLTVECT)*surfmesh->nv);
	surfmesh->face = (INT3VECT *)malloc(sizeof(INT3VECT)*surfmesh->nf);

	for (n = 0; n < surfmesh->nv; n++) {
		fscanf(fin, "%f %f %f\n", &x, &y, &z);
		surfmesh->vertex[n].x = x;
		surfmesh->vertex[n].y = y;
		surfmesh->vertex[n].z = z;
	}

	for (n = 0; n < surfmesh->nf; n++) {
		fscanf(fin, "%d %d %d %d\n", &a, &b, &c, &d);
		surfmesh->face[n].a = b;
		surfmesh->face[n].b = c;
		surfmesh->face[n].c = d;
		if (a != 3)
			printf("Errors: reading surfmesh .... \n");
	}
	fclose(fin);

}


void writePolygon()
{
	FILE *fp = fopen("output.off", "wb");
	if (!fp)
	{
		printf("\nOpenning file was unsuccessful . . . \n");
	}
	fprintf(fp, "OFF\n");
	
	fprintf(fp, "%d %d %d\n", surfmesh->nv, surfmesh->nf, 0);

	for (int n = 0; n<surfmesh->nv; n++) 
	{
		fprintf(fp, "%f %f %f\n", surfmesh->vertex[n].x, surfmesh->vertex[n].y, surfmesh->vertex[n].z);
	}

	for (int n = 0; n<surfmesh->nf; n++) 
	{
		fprintf(fp, "3 %d %d %d\n", surfmesh->face[n].a, surfmesh->face[n].b, surfmesh->face[n].c);
	}

	fclose(fp);
}

FLTVECT calc_delta_p(int i, llnode_t *neighbor_list)
{
	//Delta P is defined as a vector average in the paper
    FLTVECT ret;
    int n = 0; //Number of vertecies compared

    //vj - vi from the paper
    FLTVECT vSub;
    vSub.x = 0;
    vSub.y = 0;
    vSub.z = 0;

    llnode_t *t = neighbor_list; //use temp linked list for traversal
    while(t != NULL)
    {
        INT3VECT currFace = *((INT3VECT *)t->data); //Recall that the linked list data was filled with faces
        vSub = (vSub + surfmesh->vertex[currFace.a]) - surfmesh->vertex[i];
        vSub = (vSub + surfmesh->vertex[currFace.b]) - surfmesh->vertex[i];
        vSub = (vSub + surfmesh->vertex[currFace.c]) - surfmesh->vertex[i];

        t = (llnode_t *)(t->next);
        n += 2; //From POV of current vertex on the triangular face, we've compared 2 other verts
    }
    
    ret = (1.0 / ((float) n)) * vSub;
    return ret;
}

void gaussian_smooth(float c, llnode_t *neighbors[])
{
    FLTVECT delta_p_list[surfmesh->nv]; //Store delta_p calculations for each vertex
    for(int i = 0; i < surfmesh->nv; i++)
    {
        //Calculate delta p given the current vertex and its neighboring vertecies
        delta_p_list[i] = calc_delta_p(i, neighbors[i]);
    }
    for(int i = 0; i < surfmesh->nv; i++)
    {
        surfmesh->vertex[i] = surfmesh->vertex[i] + (c * delta_p_list[i]);
    }
}

//Taubin Smoothing
void smooth()
{
	int iter, i;
	llnode_t *neighborhood[surfmesh->nv] = {NULL}; //Create neighborhood for each vertex (used for delta p calcs)
	for(i = 0; i < surfmesh->nf; i++)
    {
        INT3VECT *currFace = &(surfmesh->face[i]); //Keep track of current face from surfmesh
        llnode_t *t; //Temporary linked list

        t = (llnode_t *)malloc(sizeof(llnode_t)); //Malloc linkedlist
        t->data = currFace; //store face
        t->next = neighborhood[(*currFace).a]; //Set next to array at current Vertex
        neighborhood[(*currFace).a] = t; //Store linkedlist to array

        t = (llnode_t *)malloc(sizeof(llnode_t));
        t->data = currFace;
        t->next = neighborhood[(*currFace).b];
        neighborhood[(*currFace).b] = t;

        t = (llnode_t *)malloc(sizeof(llnode_t));
        t->data = currFace;
        t->next = neighborhood[(*currFace).c];
        neighborhood[(*currFace).c] = t;

        //Now, neighborhood of the current vertex has a linked list to keep track of each vertex on it's face. This is our "neighborhood" used for calculating delta P
    }

	for(iter = 0; iter < ITERATIONS; iter++)
	{
		gaussian_smooth(lambda, neighborhood); //Shrink
        gaussian_smooth(mu, neighborhood); //Inflate
	}
}


int main(int argc, char *argv[])
{
    readPolygon();
	
	smooth();

	writePolygon();

	return 0;
}

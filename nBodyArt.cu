//Optimized using shared memory and on chip memory 																																			
// nvcc nBodyArt.cu -o nBodyArt -lglut -lm -lGLU -lGL
//To stop hit "control c" in the window you launched it from.
#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>
#include <GL/glut.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <cuda.h>
using namespace std;

FILE* ffmpeg;

#define PI 3.141592654
#define BLOCK 256

// Globals to be read in from parameter file.
int NumberOfBodies;
float TotalRunTime;
float Dt;
float G;
float H;
float Epsalon;
float MassOfBody;
float DiameterOfBody;
float* Diameter;
float VelocityMax;
float Drag;
int DrawRate;

// Other Globals
int Pause;
float4 *BodyPosition, *BodyVelocity, *BodyForce;
float4 *BodyPositionGPU, *BodyVelocityGPU, *BodyForceGPU;
float4 *BodyColor;
dim3 Blocks, Grids;
int DrawTimer;
float RunTime;
int* Buffer;
int MovieOn;

// Window globals
static int Window;
int XWindowSize;
int YWindowSize; 
double Near;
double Far;
double EyeX;
double EyeY;
double EyeZ;
double CenterX;
double CenterY;
double CenterZ;
double UpX;
double UpY;
double UpZ;

// Prototyping functions
void readSimulationParameters();
void allocateMemory();
void setInitailConditions();
void drawPicture();
void nBody();
void errorCheck(const char*);
void setup();

#include "./callBackFunctions.h"

void readSimulationParameters()
{
	ifstream data;
	string name;
	
	data.open("./simulationSetup");
	
	if(data.is_open() == 1)
	{
		getline(data,name,'=');
		data >> NumberOfBodies;
		
		getline(data,name,'=');
		data >> TotalRunTime;
		
		getline(data,name,'=');
		data >> Dt;
		
		getline(data,name,'=');
		data >> G;
		
		getline(data,name,'=');
		data >> H;
		
		getline(data,name,'=');
		data >> Epsalon;
		
		getline(data,name,'=');
		data >> MassOfBody;
		
		getline(data,name,'=');
		data >> DiameterOfBody;
		
		getline(data,name,'=');
		data >> VelocityMax;
		
		getline(data,name,'=');
		data >> Drag;
		
		getline(data,name,'=');
		data >> DrawRate;
	}
	else
	{
		printf("\nTSU Error could not open simulationSetup file\n");
		exit(0);
	}
	data.close();
	
	printf("\n\n Parameter file has been read");
}

void allocateMemory()
{
	Blocks.x = BLOCK;
	Blocks.y = 1;
	Blocks.z = 1;
	
	Grids.x = (NumberOfBodies - 1)/Blocks.x + 1;
	Grids.y = 1;
	Grids.z = 1;
	
	BodyPosition = (float4*)malloc(NumberOfBodies*sizeof(float4));
	BodyVelocity = (float4*)malloc(NumberOfBodies*sizeof(float4));
	BodyForce    = (float4*)malloc(NumberOfBodies*sizeof(float4));
	BodyColor    = (float4*)malloc(NumberOfBodies*sizeof(float4));
	Diameter = (float*)malloc(NumberOfBodies*sizeof(float));
	
	cudaMalloc( (void**)&BodyPositionGPU, NumberOfBodies *sizeof(float4));
	errorCheck("cudaMalloc BodyPositionGPU");
	cudaMalloc( (void**)&BodyVelocityGPU, NumberOfBodies *sizeof(float4));
	errorCheck("cudaMalloc BodyDiameterOfBodyVelocityGPU");
	cudaMalloc( (void**)&BodyForceGPU,    NumberOfBodies *sizeof(float4));
	errorCheck("cudaMalloc BodyForceGPU");
	
	
	printf("\n\n Memory has been allocated");
}

void setInitailConditions()
{
    float dx, dy, dz, d, d2;
    int test;
	time_t t;
	
	srand((unsigned) time(&t));
	for(int i = 0; i < NumberOfBodies; i++)
	{
		Diameter[i] = ((float)rand()/(float)RAND_MAX)*0.01 - 0.001;	
	}	
	for(int i = 0; i < NumberOfBodies; i++)
	{
		
		test = 0;
		while(test == 0)
		{
			// Get random number between -1 at 1.
			BodyPosition[i].x = ((float)rand()/(float)RAND_MAX)*2.0 - 1.0;
			BodyPosition[i].y = ((float)rand()/(float)RAND_MAX)*2.0 - 1.0;
			BodyPosition[i].z = ((float)rand()/(float)RAND_MAX)*2.0 - 1.0;
			BodyPosition[i].w = ((float)rand()/(float)RAND_MAX)*2.0 - 1.0; 	//MassOfBody;
			test = 1;
			
			for(int j = 0; j < i; j++)
			{
				dx = BodyPosition[i].x-BodyPosition[j].x;
				dy = BodyPosition[i].y-BodyPosition[j].y;
				dz = BodyPosition[i].z-BodyPosition[j].z;
				d2  = dx*dx + dy*dy + dz*dz;
				d = sqrt(d2);
				
				if(d < DiameterOfBody)
				{
					test = 0;
					break;
				}
			}
			
			if(test == 1)
			{
				BodyVelocity[i].x = 0.0; //VelocityMax*((float)rand()/(float)RAND_MAX)*2.0 - 1.0;
				BodyVelocity[i].y = 0.0; //VelocityMax*((float)rand()/(float)RAND_MAX)*2.0 - 1.0;
				BodyVelocity[i].z = 0.0;  //VelocityMax*((float)rand()/(float)RAND_MAX)*2.0 - 1.0;
				BodyVelocity[i].w = 0.0;
				
				BodyColor[i].x = ((float)rand()/(float)RAND_MAX);
				BodyColor[i].y = ((float)rand()/(float)RAND_MAX);
				BodyColor[i].z = ((float)rand()/(float)RAND_MAX);
				BodyColor[i].w = 0.0;
			}
		}
	}
	printf("\n\n Initail conditions have been set.");
}

void drawPicture()
{
	//glClear(GL_COLOR_BUFFER_BIT);
	//glClear(GL_DEPTH_BUFFER_BIT);
	
	for(int i = 0; i < NumberOfBodies; i++)
	{
		glColor3d(BodyColor[i].x, BodyColor[i].y, BodyColor[i].z);
		//glColor3d(1.0, 1.0, 1.0);
		glPushMatrix();
			glTranslatef(BodyPosition[i].x, BodyPosition[i].y, BodyPosition[i].z);
			glutSolidSphere(Diameter[i]/2.0, 30, 30);
			
			//glTranslatef(0.0, 0.0, 0.0);
			//glutSolidSphere(1.0, 20, 20);
		glPopMatrix();
	}
	glutSwapBuffers();
	
	
	if(MovieOn == 1)
	{
		glReadPixels(5, 5, XWindowSize, YWindowSize, GL_RGBA, GL_UNSIGNED_BYTE, Buffer);
		fwrite(Buffer, sizeof(int)*XWindowSize*YWindowSize, 1, ffmpeg);
	}
}
                                 
__device__ float3 getBodyBodyForce(float4 p0, float4 p1, float G, float H, float Epsalon)
{
    float3 f;
    float dx = p1.x - p0.x;
    float dy = p1.y - p0.y;
    float dz = p1.z - p0.z;
    float r2 = dx*dx + dy*dy + dz*dz + Epsalon;
    float r = sqrt(r2);
    
    float force  = (G*p0.w*p1.w)/(r2) - (H*p0.w*p1.w)/(r2*r2);
    
    f.x = force*dx/r;
    f.y = force*dy/r;
    f.z = force*dz/r;
    
    return(f);
}

__global__ void getForces(float4 *pos, float4 *vel, float4 * force, float G, float H, float Epsalon, int n)
{
	int j,ii;
    float3 force_mag, forceSum;
    float4 posMe;
    __shared__ float4 shPos[BLOCK];
    int id = threadIdx.x + blockDim.x*blockIdx.x;
    
    forceSum.x = 0.0;
	forceSum.y = 0.0;
	forceSum.z = 0.0;
		
	posMe.x = pos[id].x;
	posMe.y = pos[id].y;
	posMe.z = pos[id].z;
	posMe.w = pos[id].w;
	    
    for(j=0; j < gridDim.x; j++)
    {
    	shPos[threadIdx.x] = pos[threadIdx.x + blockDim.x*j];
    	__syncthreads();
   
		#pragma unroll 32
        for(int i=0; i < blockDim.x; i++)	
        {
        	ii = i + blockDim.x*j;
		    if(ii != id && ii < n) 
		    {
		    	force_mag = getBodyBodyForce(posMe, shPos[i], G, H, Epsalon);
			    forceSum.x += force_mag.x;
			    forceSum.y += force_mag.y;
			    forceSum.z += force_mag.z;
		    }
	   	 }
	}
	if(id < n)
	{
	    force[id].x = forceSum.x;
	    force[id].y = forceSum.y;
	    force[id].z = forceSum.z;
    }
}

__global__ void moveBodies(float4 *pos, float4 *vel, float4 * force, float drag, float dt, int n)
{
    int id = threadIdx.x + blockDim.x*blockIdx.x;
    if(id < n)
    {
	    vel[id].x += ((force[id].x-drag*vel[id].x)/pos[id].w)*dt;
	    vel[id].y += ((force[id].y-drag*vel[id].y)/pos[id].w)*dt;
	    vel[id].z += ((force[id].z-drag*vel[id].z)/pos[id].w)*dt;
	
	    pos[id].x += vel[id].x*dt;
	    pos[id].y += vel[id].y*dt;
	    pos[id].z += vel[id].z*dt;
	    
	    
    }
}

void nBody()
{
	//drawPicture();
	//while(1);
	if(Pause != 1)
	{	
		getForces<<<Grids, Blocks>>>(BodyPositionGPU, BodyVelocityGPU, BodyForceGPU, G, H, Epsalon, NumberOfBodies);
		moveBodies<<<Grids, Blocks>>>(BodyPositionGPU, BodyVelocityGPU, BodyForceGPU, Drag, Dt, NumberOfBodies);
        
        DrawTimer++;
		if(DrawTimer == DrawRate) 
		{
		    cudaMemcpy( BodyPosition, BodyPositionGPU, NumberOfBodies*sizeof(float4), cudaMemcpyDeviceToHost );
			drawPicture();
			//printf("\n Time = %f", RunTime);
			DrawTimer = 0;
		}
		RunTime += Dt; 
		if(TotalRunTime < RunTime)
		{
			printf("\n\n Done\n");
			exit(0);
		}
	}
}

void errorCheck(const char *message)
{
  cudaError_t  error;
  error = cudaGetLastError();

  if(error != cudaSuccess)
  {
    printf("\n CUDA ERROR: %s = %s\n", message, cudaGetErrorString(error));
    exit(0);
  }
}

void setup()
{	
	readSimulationParameters();
	allocateMemory();
	setInitailConditions();
	cudaMemcpy( BodyPositionGPU, BodyPosition, NumberOfBodies*sizeof(float4), cudaMemcpyHostToDevice );
    cudaMemcpy( BodyVelocityGPU, BodyVelocity, NumberOfBodies*sizeof(float4), cudaMemcpyHostToDevice );
    DrawTimer = 0;
	RunTime = 0.0;
	Pause = 1;
}

int main(int argc, char** argv)
{
	setup();
	
	XWindowSize = 1000;
	YWindowSize = 1000; 
	Buffer = new int[XWindowSize*YWindowSize];

	// Clip plains
	Near = 0.2;
	Far = 30.0;

	//Direction here your eye is located location
	EyeX = 0.0;
	EyeY = 0.0;
	EyeZ = 2.0;

	//Where you are looking
	CenterX = 0.0;
	CenterY = 0.0;
	CenterZ = 0.0;

	//Up vector for viewing
	UpX = 0.0;
	UpY = 1.0;
	UpZ = 0.0;
	
	glutInit(&argc,argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_DEPTH | GLUT_RGB);
	glutInitWindowSize(XWindowSize,YWindowSize);
	glutInitWindowPosition(5,5);
	Window = glutCreateWindow("N Body");
	
	gluLookAt(EyeX, EyeY, EyeZ, CenterX, CenterY, CenterZ, UpX, UpY, UpZ);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glFrustum(-0.2, 0.2, -0.2, 0.2, Near, Far);
	glMatrixMode(GL_MODELVIEW);
	glClearColor(0.0, 0.0, 0.0, 0.0);
	
	GLfloat light_position[] = {1.0, 1.0, 1.0, 0.0};
	GLfloat light_ambient[]  = {0.0, 0.0, 0.0, 1.0};
	GLfloat light_diffuse[]  = {1.0, 1.0, 1.0, 1.0};
	GLfloat light_specular[] = {1.0, 1.0, 1.0, 1.0};
	GLfloat lmodel_ambient[] = {0.2, 0.2, 0.2, 1.0};
	GLfloat mat_specular[]   = {1.0, 1.0, 1.0, 1.0};
	GLfloat mat_shininess[]  = {10.0};
	glShadeModel(GL_SMOOTH);
	glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);
	glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
	glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lmodel_ambient);
	glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
	glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_COLOR_MATERIAL);
	glEnable(GL_DEPTH_TEST);
	
	glutDisplayFunc(Display);
	glutReshapeFunc(reshape);
	glutMouseFunc(mymouse);
	glutKeyboardFunc(KeyPressed);
	glutIdleFunc(idle);
	glutMainLoop();
	return 0;
}







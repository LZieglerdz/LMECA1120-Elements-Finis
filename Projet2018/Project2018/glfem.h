#ifndef _GLFEM_H_
#define _GLFEM_H_

#define GLFW_INCLUDE_GLU
#include <GLFW/glfw3.h>

#include "fem.h"

GLFWwindow* window2;

void    glfemDrawColorElement(float *x, float *y, double *u, int n);
void    glfemDrawElement(float *x, float *y, int n);
void    glfemDrawNodes(double* x, double* y,int n,double r);
void    glfemDrawCircle(double x, double y,double r);
void    glfemDrawDisk(double x, double y, double r);

void    glfemReshapeWindows(double r, int width, int height);

void 		glfemPlotField(femMesh *theMesh, double *u);
void 		glfemPlotMesh(femMesh *theMesh);
void 		glfemPlotEdges(femEdges *theEdges);
void 		glfemPlotBnd(femEdges *theEdges);

void    glfemMessage(char *message);
void    glfemDrawMessage(int h, int v, char *message);
void    glfemSetRasterSize(int width, int height);
GLFWwindow* glfemInit(char *windowName);

#endif

/*
 *  main.c
 *  Library for MECA1120 : Finite Elements for dummies
 *
 *  Copyright (C) 2018 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 *
 */

#include "glfem.h"

int main(void)
{
  femPoissonProblem* theProblem = femPoissonCreate("../data/meshMedium.txt");


    float omegaf;
    float muf;
    float nGrains;
    printf("enter number of grains >>");
    scanf("%f", &nGrains);
    printf("enter angular velocity between -4 and 4 >>");
    scanf("%f", &omegaf);
    printf("enter viscosity between 1e-3 and 1e8 >>");
    scanf("%f", &muf);

    int    n          = (double) nGrains;
    double radiusIn   = radIn(theProblem);
    double radiusOut  = radOut(theProblem);
    double omega      = (double) omegaf;  //En rad/s -> min -4, max 4
    double mu         = (double) muf;     // -> min 10-3, max 10-8
    double gamma      = 0.5;
    double radius     = 0.05;
    double mass       = 0.1;
    double dt         = 1e-1;
    double tEnd       = 8.0;
    double tol        = 1e-6;
    double t          = 0;
    double iterMax    = 100;

  femGrains* myGrains = femGrainsCreateSimple(n,radius,mass,radiusIn,radiusOut, gamma);

  printf("Number of elements    : %4d\n", theProblem->mesh->nElem);
  printf("Number of local nodes : %4d\n", theProblem->mesh->nLocalNode);
  printf("Number of segments    : %4d\n", theProblem->edges->nBoundary);
  printf("Number of unknowns(x) : %4d\n", theProblem->systemX->size);
  printf("Number of unknowns(y) : %4d\n", theProblem->systemY->size);

  femPoissonSolve(theProblem, omega, mu, myGrains);

  printf("Maximum value : %.4f\n", femMax(theProblem->systemX->B,theProblem->systemX->size));
  printf("Maximum value : %.4f\n", femMax(theProblem->systemY->B,theProblem->systemY->size));
  fflush(stdout);

  char theMessage[256];
  sprintf(theMessage, "Max : %.4f", femMax(theProblem->systemX->B,theProblem->systemX->size));
  sprintf(theMessage, "Max : %.4f", femMax(theProblem->systemY->B,theProblem->systemY->size));

    GLFWwindow* window = glfemInit("MECA1120 : Projet 2018 - Un Poisson dans une machine à laver");
    glfwMakeContextCurrent(window);
    int theRunningMode = 1.0;
    float theVelocityFactor = 0.25;

    do {
        int i,w,h;
        double currentTime = glfwGetTime();

        glfwGetFramebufferSize(window,&w,&h);
        glfemReshapeWindows(radiusOut,w,h);
        double u[theProblem->systemX->size];
        for (i=0 ;i < theProblem->systemX->size; i++) {
          u[i] = sqrt(pow(theProblem->systemX->B[i], 2) + pow(theProblem->systemY->B[i],2));
        }
        glfemPlotField(theProblem->mesh, u );

        for (i=0 ;i < myGrains->n; i++) {
            glColor3f(0,0,0);
            glfemDrawDisk(myGrains->x[i],myGrains->y[i],myGrains->r[i]);
          }
        glColor3f(0,0,0); glfemDrawCircle(0,0,radiusOut);
        glColor3f(0,0,0); glfemDrawCircle(0,0,radiusIn);
        char theMessage[256];
        char theVisco[256];
        char theDrag[256];
        char theVelocity[256];
        sprintf(theVelocity,"Angular Velocity = %.2f rad\n",omega);
        sprintf(theVisco,"Viscosity = %.2e Pa*s\n",mu);
        sprintf(theDrag,"Drag = %.2f",gamma);
        sprintf(theMessage,"Time = %g sec",t);
        glColor3f(0,0,0); glfemDrawMessage(20,400,theMessage);
        glColor3f(0,0,0); glfemDrawMessage(20,420,theVelocity);
        glColor3f(0,0,0); glfemDrawMessage(20,440,theVisco);
        glColor3f(0,0,0); glfemDrawMessage(20,460,theDrag);

        glfwSwapBuffers(window);
        glfwPollEvents();

        if (t < tEnd && theRunningMode == 1) {
            printf("Time = %4g : ",t);
            femGrainsUpdate(theProblem, myGrains, dt, tol, iterMax, omega, mu);
            t += dt;
          }

        while ( glfwGetTime()-currentTime < theVelocityFactor ) {
          if (glfwGetKey(window,'R') == GLFW_PRESS)
        	    theRunningMode = 1;
          if (glfwGetKey(window,'S') == GLFW_PRESS)
        	    theRunningMode = 0;
            }
    }
    while (glfwGetKey(window, GLFW_KEY_ESCAPE ) != GLFW_PRESS &&
	        (!glfwWindowShouldClose(window)));

    glfwTerminate();
    femGrainsFree(myGrains);
    exit(EXIT_SUCCESS);
}

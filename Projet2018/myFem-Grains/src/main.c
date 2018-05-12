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

  int    n = 250;
  double omega = M_PI;   //En rad/s
  double mu = 1;
  double radius    = 0.03;
  double mass      = 0.1;
  double radiusIn  = 0.5;
  double radiusOut = 2.0;
  double dt      = 1e-1;
  double tEnd    = 8.0;
  double tol     = 1e-6;
  double t       = 0;
  double iterMax = 100;
  femGrains* theGrains = femGrainsCreateSimple(n,radius,mass,radiusIn,radiusOut);

  femPoissonProblem* theProblem = femPoissonCreate("../data/meshMedium.txt");

  printf("Number of elements    : %4d\n", theProblem->mesh->nElem);
  printf("Number of local nodes : %4d\n", theProblem->mesh->nLocalNode);
  printf("Number of segments    : %4d\n", theProblem->edges->nBoundary);
  printf("Number of unknowns(x) : %4d\n", theProblem->systemX->size);
  printf("Number of unknowns(y) : %4d\n", theProblem->systemY->size);

  femPoissonSolve(theProblem, omega, mu);

  printf("Maximum value : %.4f\n", femMax(theProblem->systemX->B,theProblem->systemX->size));
  printf("Maximum value : %.4f\n", femMax(theProblem->systemY->B,theProblem->systemY->size));
  fflush(stdout);

  char theMessage[256];
  sprintf(theMessage, "Max : %.4f", femMax(theProblem->systemX->B,theProblem->systemX->size));
  sprintf(theMessage, "Max : %.4f", femMax(theProblem->systemY->B,theProblem->systemY->size));


  //  A decommenter pour obtenir l'exemple de la seance d'exercice :-)
  //  femGrains* theGrains = femGrainsCreateTiny(radiusIn,radiusOut);;

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

        for (i=0 ;i < theGrains->n; i++) {
            glColor3f(102/255,255/255,102/255);
            glfemDrawDisk(theGrains->x[i],theGrains->y[i],theGrains->r[i]); }
        glColor3f(0,0,0); glfemDrawCircle(0,0,radiusOut);
        glColor3f(0,0,0); glfemDrawCircle(0,0,radiusIn);
        char theMessage[256];
        sprintf(theMessage,"Time = %g sec",t);
        glColor3f(1,0,0); glfemDrawMessage(20,460,theMessage);
        glfwSwapBuffers(window);
        glfwPollEvents();

        if (t < tEnd && theRunningMode == 1) {
            printf("Time = %4g : ",t);
  //
  // A decommenter pour pouvoir progresser pas par pas
  //          printf("press CR to compute the next time step >>");
  //          char c= getchar();
  //
            femGrainsUpdate(theProblem, theGrains,dt,tol,iterMax);
            t += dt; }

        while ( glfwGetTime()-currentTime < theVelocityFactor ) {
          if (glfwGetKey(window,'R') == GLFW_PRESS)
        	    theRunningMode = 1;
          if (glfwGetKey(window,'S') == GLFW_PRESS)
        	    theRunningMode = 0; }

    }
    while (glfwGetKey(window, GLFW_KEY_ESCAPE ) != GLFW_PRESS &&
	        (!glfwWindowShouldClose(window)));


    glfwTerminate();
    femGrainsFree(theGrains);
    exit(EXIT_SUCCESS);
}

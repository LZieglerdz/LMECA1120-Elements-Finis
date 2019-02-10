/*
 *  main.c
 *  Library for MECA1120 : Finite Elements for dummies
 *
 *  Copyright (C) 2018 UCL-EPL : Vincent Legat
 *  All rights reserved.
 *
 */

#include "glfem.h"
#include "time.h"

int main(void)
{
	clock_t begin = clock();

    femPoissonProblem* theProblem = femPoissonCreate("../data/meshMedium.txt");

    printf("Number of elements    : %4d\n", theProblem->mesh->nElem);
    printf("Number of local nodes : %4d\n", theProblem->mesh->nLocalNode);
    printf("Number of segments    : %4d\n", theProblem->edges->nBoundary);
    printf("Number of unknowns X  : %4d\n", theProblem->systemX->size);
		printf("Number of unknowns Y  : %4d\n", theProblem->systemY->size);

    femPoissonSolve(theProblem);

    printf("Maximum value X : %.4f\n", femMax(theProblem->systemX->B,theProblem->systemX->size));
		printf("Maximum value Y : %.4f\n", femMax(theProblem->systemY->B,theProblem->systemY->size));
    fflush(stdout);

    char theMessage[256];
    sprintf(theMessage, "Max : %.4f", femMax(theProblem->systemX->B,theProblem->systemX->size));
		sprintf(theMessage, "Max : %.4f", femMax(theProblem->systemY->B,theProblem->systemY->size));

        clock_t end = clock();
       	double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	printf("execution time: %f \n", time_spent);



    GLFWwindow* window = glfemInit("MECA1120 : homework 3 ");
    glfwMakeContextCurrent(window);
    do {
        int w,h,i;
        glfwGetFramebufferSize(window,&w,&h);
        glfemReshapeWindows(theProblem->mesh,w,h);
				double u[theProblem->systemX->size];
        for (i=0 ;i < theProblem->systemX->size; i++) {
          u[i] = theProblem->systemX->B[i] + theProblem->systemY->B[i] ;
        }
        glfemPlotField(theProblem->mesh,u);
        glColor3f(1.0,0.0,0.0); glfemDrawMessage(20,460,theMessage);
        glfwSwapBuffers(window);
        glfwPollEvents();
    }
    while( glfwGetKey(window,GLFW_KEY_ESCAPE) != GLFW_PRESS &&
             glfwWindowShouldClose(window) != 1 );

    // Check if the ESC key was pressed or the window was closed

    glfwTerminate();
    femPoissonFree(theProblem);

    exit(EXIT_SUCCESS);



}

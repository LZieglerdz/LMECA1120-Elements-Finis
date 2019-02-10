/*
 *  glfem.c
 *  Library for MECA1120 : Finite Elements for dummies
 *
 *  Copyright (C) 2018 UCL-EPL : Vincent Legat
 *  All rights reserved.
 *
 *  Pour GLFW (version utilisée 3.2.1)
 *  Pour l'installation de la librairie, voir http://www.glfw.org/
 *
 */

#include "glfem.h"


void glMakeRasterFont(void);
void getColor(double value, int numberOfColors, float* R, float* G, float* B);
double glScale(double minimum, double maximum, double value);

static int gRasterH = 800;
static int gRasterV = 600;
static int numberColors = 50;

GLubyte space[] = 
    {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00};

GLubyte letters[][13] = {
    {0x00, 0x00, 0xc3, 0xc3, 0xc3, 0xc3, 0xff, 0xc3, 0xc3, 0xc3, 0x66, 0x3c, 0x18}, // A
    {0x00, 0x00, 0xfe, 0xc7, 0xc3, 0xc3, 0xc7, 0xfe, 0xc7, 0xc3, 0xc3, 0xc7, 0xfe}, // B
    {0x00, 0x00, 0x7e, 0xe7, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xe7, 0x7e}, // C
    {0x00, 0x00, 0xfc, 0xce, 0xc7, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc7, 0xce, 0xfc}, // D
    {0x00, 0x00, 0xff, 0xc0, 0xc0, 0xc0, 0xc0, 0xfc, 0xc0, 0xc0, 0xc0, 0xc0, 0xff}, // E
    {0x00, 0x00, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xfc, 0xc0, 0xc0, 0xc0, 0xff}, // F
    {0x00, 0x00, 0x7e, 0xe7, 0xc3, 0xc3, 0xcf, 0xc0, 0xc0, 0xc0, 0xc0, 0xe7, 0x7e}, // G
    {0x00, 0x00, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xff, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3}, // H
    {0x00, 0x00, 0x7e, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x7e}, // I
    {0x00, 0x00, 0x7c, 0xee, 0xc6, 0x06, 0x06, 0x06, 0x06, 0x06, 0x06, 0x06, 0x06}, // J
    {0x00, 0x00, 0xc3, 0xc6, 0xcc, 0xd8, 0xf0, 0xe0, 0xf0, 0xd8, 0xcc, 0xc6, 0xc3}, // K
    {0x00, 0x00, 0xff, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0}, // L
    {0x00, 0x00, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xdb, 0xff, 0xff, 0xe7, 0xc3}, // M
    {0x00, 0x00, 0xc7, 0xc7, 0xcf, 0xcf, 0xdf, 0xdb, 0xfb, 0xf3, 0xf3, 0xe3, 0xe3}, // N
    {0x00, 0x00, 0x7e, 0xe7, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xe7, 0x7e}, // O
    {0x00, 0x00, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xfe, 0xc7, 0xc3, 0xc3, 0xc7, 0xfe}, // P
    {0x00, 0x00, 0x3f, 0x6e, 0xdf, 0xdb, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0x66, 0x3c}, // Q
    {0x00, 0x00, 0xc3, 0xc6, 0xcc, 0xd8, 0xf0, 0xfe, 0xc7, 0xc3, 0xc3, 0xc7, 0xfe}, // R
    {0x00, 0x00, 0x7e, 0xe7, 0x03, 0x03, 0x07, 0x7e, 0xe0, 0xc0, 0xc0, 0xe7, 0x7e}, // S
    {0x00, 0x00, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0xff}, // T
    {0x00, 0x00, 0x7e, 0xe7, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3}, // U
    {0x00, 0x00, 0x18, 0x3c, 0x3c, 0x66, 0x66, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3}, // V
    {0x00, 0x00, 0xc3, 0xe7, 0xff, 0xff, 0xdb, 0xdb, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3}, // W
    {0x00, 0x00, 0xc3, 0x66, 0x66, 0x3c, 0x3c, 0x18, 0x3c, 0x3c, 0x66, 0x66, 0xc3}, // X
    {0x00, 0x00, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x3c, 0x3c, 0x66, 0x66, 0xc3}, // Y
    {0x00, 0x00, 0xff, 0xc0, 0xc0, 0x60, 0x30, 0x7e, 0x0c, 0x06, 0x03, 0x03, 0xff}  // Z
};

GLubyte lowletters[][13] = {
    {0x00, 0x00, 0x7d, 0xc3, 0xc3, 0xc3, 0x7f, 0x03, 0x7e, 0x00, 0x00, 0x00, 0x00}, // a
    {0x00, 0x00, 0xfe, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xfe, 0xc0, 0xc0, 0xc0, 0xc0}, // b
    {0x00, 0x00, 0x7f, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0x7f, 0x00, 0x00, 0x00, 0x00}, // c
    {0x00, 0x00, 0x7f, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0x7f, 0x03, 0x03, 0x03, 0x03}, // d
    {0x00, 0x00, 0x7e, 0xc0, 0xc0, 0xfe, 0xc3, 0xc3, 0x7e, 0x00, 0x00, 0x00, 0x00}, // e
    {0x00, 0x00, 0x3c, 0x18, 0x18, 0x18, 0x18, 0x18, 0x7e, 0x18, 0x18, 0x18, 0x0e}, // f
    {0x7f, 0x03, 0x7f, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0x7e, 0x00, 0x00, 0x00, 0x00}, // g
    {0x00, 0x00, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xfe, 0xc0, 0xc0, 0xc0, 0xc0}, // h
    {0x00, 0x00, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x00, 0x18, 0x18, 0x00}, // i
    {0x70, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x00, 0x18, 0x18, 0x00}, // j
    {0x00, 0x00, 0xc3, 0xc7, 0xce, 0xfc, 0xfe, 0xc7, 0xc3, 0xc0, 0xc0, 0xc0, 0xc0}, // k
    {0x00, 0x00, 0x0c, 0x1c, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18}, // l
    {0x00, 0x00, 0xdb, 0xdb, 0xdb, 0xdb, 0xdb, 0xdb, 0xfe, 0x00, 0x00, 0x00, 0x00}, // m
    {0x00, 0x00, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xfe, 0x00, 0x00, 0x00, 0x00}, // n
    {0x00, 0x00, 0x7e, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0x7e, 0x00, 0x00, 0x00, 0x00}, // o
    {0xc0, 0xc0, 0xfe, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xfe, 0x00, 0x00, 0x00, 0x00}, // p
    {0x03, 0x03, 0x7f, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0x7f, 0x00, 0x00, 0x00, 0x00}, // q
    {0x00, 0x00, 0xc0, 0xc0, 0xc0, 0xc0, 0xe0, 0xf0, 0xdf, 0x00, 0x00, 0x00, 0x00}, // r
    {0x00, 0x00, 0xfe, 0x03, 0x03, 0x7e, 0xc0, 0xc0, 0x7f, 0x00, 0x00, 0x00, 0x00}, // s
    {0x00, 0x00, 0x0e, 0x18, 0x18, 0x18, 0x18, 0x18, 0x7e, 0x18, 0x18, 0x18, 0x18}, // t
    {0x00, 0x00, 0x7f, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0x00, 0x00, 0x00, 0x00}, // u
    {0x00, 0x00, 0x18, 0x3c, 0x66, 0x66, 0xc3, 0xc3, 0xc3, 0x00, 0x00, 0x00, 0x00}, // v
    {0x00, 0x00, 0x66, 0x7e, 0xdb, 0xdb, 0xdb, 0xdb, 0xdb, 0x00, 0x00, 0x00, 0x00}, // w
    {0x00, 0x00, 0xc3, 0xe7, 0x3c, 0x18, 0x3c, 0xe7, 0xc3, 0x00, 0x00, 0x00, 0x00}, // x
    {0x7f, 0x03, 0x7f, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0x00, 0x00, 0x00, 0x00}, // y
    {0x00, 0x00, 0xff, 0xc0, 0x70, 0x1c, 0x06, 0x03, 0xff, 0x00, 0x00, 0x00, 0x00}, // z
};

GLubyte numletters[][13] = {
    {0x00, 0x00, 0x3c, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x3c}, // 0
    {0x00, 0x00, 0x3c, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x78, 0x38, 0x18}, // 1
    {0x00, 0x00, 0x7e, 0x60, 0x60, 0x60, 0x60, 0x3c, 0x06, 0x06, 0x66, 0x66, 0x3c}, // 2
    {0x00, 0x00, 0x3c, 0x66, 0x06, 0x06, 0x06, 0x1c, 0x06, 0x06, 0x06, 0x66, 0x3c}, // 3
    {0x00, 0x00, 0x06, 0x06, 0x06, 0x06, 0x06, 0x7f, 0x66, 0x36, 0x1e, 0x0e, 0x06}, // 4
    {0x00, 0x00, 0x3c, 0x66, 0x06, 0x06, 0x06, 0x7c, 0x60, 0x60, 0x60, 0x60, 0x7e}, // 5
    {0x00, 0x00, 0x3c, 0x66, 0x66, 0x66, 0x66, 0x66, 0x7c, 0x60, 0x60, 0x66, 0x3c}, // 6
    {0x00, 0x00, 0x0c, 0x0c, 0x0c, 0x0c, 0x0c, 0x1f, 0x06, 0x06, 0x06, 0x06, 0x7e}, // 7
    {0x00, 0x00, 0x3c, 0x66, 0x66, 0x66, 0x66, 0x3c, 0x66, 0x66, 0x66, 0x66, 0x3c}, // 8
    {0x00, 0x00, 0x3c, 0x66, 0x06, 0x06, 0x06, 0x3e, 0x66, 0x66, 0x66, 0x66, 0x3c}, // 9
    {0x00, 0x00, 0x00, 0x18, 0x18, 0x00, 0x00, 0x18, 0x18, 0x00, 0x00, 0x00, 0x00}, // :
    {0x00, 0x00, 0x30, 0x18, 0x18, 0x00, 0x00, 0x18, 0x18, 0x00, 0x00, 0x00, 0x00}, // ;
    {0x00, 0x00, 0x06, 0x1c, 0x30, 0x60, 0x30, 0x1c, 0x06, 0x00, 0x00, 0x00, 0x00}, // <
    {0x00, 0x00, 0x00, 0x00, 0x3c, 0x00, 0x00, 0x3c, 0x00, 0x00, 0x00, 0x00, 0x00}, // =
    {0x00, 0x00, 0x60, 0x38, 0x0c, 0x06, 0x0c, 0x38, 0x60, 0x00, 0x00, 0x00, 0x00}, // >
    {0x00, 0x00, 0x18, 0x18, 0x00, 0x18, 0x18, 0x18, 0x0c, 0x06, 0x06, 0x66, 0x3c}, // ?
};

GLubyte specialletters[][13] = {
    {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00}, // space
    {0x00, 0x00, 0x18, 0x18, 0x00, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18}, // !
    {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x24, 0x24, 0x00, 0x00}, // "
    {0x00, 0x00, 0x24, 0x24, 0x7e, 0x7e, 0x24, 0x7e, 0x7e, 0x24, 0x24, 0x00, 0x00}, // #
    {0x00, 0x00, 0x18, 0x3c, 0x5a, 0x5a, 0x1a, 0x3c, 0x58, 0x58, 0x5a, 0x3c, 0x18}, // $
    {0x00, 0x00, 0x44, 0x4a, 0x6a, 0x24, 0x30, 0x18, 0x0c, 0x24, 0x56, 0x52, 0x22}, // %
    
    {0x00, 0x00, 0x79, 0xcf, 0xc6, 0xcf, 0x79, 0x70, 0x78, 0xcc, 0xcc, 0xcc, 0x78}, // &
    {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x08, 0x08, 0x18, 0x00, 0x00}, // '
    {0x00, 0x00, 0x0c, 0x18, 0x18, 0x30, 0x30, 0x30, 0x30, 0x30, 0x18, 0x18, 0x0c}, // (
    {0x00, 0x00, 0x30, 0x18, 0x18, 0x0c, 0x0c, 0x0c, 0x0c, 0x0c, 0x18, 0x18, 0x30}, // )
    {0x00, 0x00, 0x00, 0x00, 0x10, 0x54, 0x38, 0x54, 0x10, 0x00, 0x00, 0x00, 0x00}, // *
    {0x00, 0x00, 0x00, 0x00, 0x10, 0x10, 0x7c, 0x10, 0x10, 0x00, 0x00, 0x00, 0x00}, // +
    {0x00, 0x30, 0x18, 0x18, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00}, // ,
    {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x3c, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00}, // -
    {0x00, 0x00, 0x18, 0x18, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00}, // .
    {0x00, 0x00, 0x60, 0x60, 0x30, 0x30, 0x18, 0x18, 0x18, 0x0c, 0x0c, 0x06, 0x06}, // /
};

GLuint fontOffset;

void glMakeRasterFont(void)
{
    GLuint i, j;
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);  
    glShadeModel (GL_FLAT);
    fontOffset = glGenLists (128);
    
    for (i = 0,j = 'A'; i < 26; i++,j++) {
        glNewList(fontOffset + j, GL_COMPILE);
        glBitmap(8, 13, 0.0, 2.0, 10.0, 0.0, letters[i]);
        glEndList(); }
    
    for (i = 0,j = 'a'; i < 26; i++,j++) {
        glNewList(fontOffset + j, GL_COMPILE);
        glBitmap(8, 13, 0.0, 2.0, 10.0, 0.0, lowletters[i]);
        glEndList(); }
    
    for (i = 0,j = '0'; i < 16; i++,j++) {
        glNewList(fontOffset + j, GL_COMPILE);
        glBitmap(8, 13, 0.0, 2.0, 10.0, 0.0, numletters[i]);
        glEndList(); }
    
    for (i = 0,j = ' '; i < 16; i++,j++) {
        glNewList(fontOffset + j, GL_COMPILE);
        glBitmap(8, 13, 0.0, 2.0, 10.0, 0.0, specialletters[i]);
        glEndList(); }
    glShadeModel (GL_SMOOTH);
}

void glfemSetRasterSize(int h, int v)
{
    gRasterH = h;
    gRasterV = v;
}

void glfemDrawMessage(int h, int v, char *s)
{
    // Les coordonnées négatives sont normalement admises :-)
    // On conserve le strint entier même si le début est off-screen

    int off;    
    glPushAttrib(GL_ALL_ATTRIB_BITS);   
    glDisable(GL_TEXTURE_2D);   
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    glOrtho (0.0, gRasterH, gRasterV, 0.0, -1.0, 1.0);
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    
    if (h < 0)   h = gRasterH + h - strlen(s)*10;
    if (v < 0)   v = gRasterV + v; 
    glRasterPos2i(h, v);    
    glListBase(fontOffset);
    
    if (h >= 0) {
      glRasterPos2i(h, v);
      glCallLists((GLsizei) strlen(s), GL_UNSIGNED_BYTE, (GLubyte *) s); }
    else { 
      off = (h-9)/10;
      glRasterPos2i(h - off*10, v);
      if (strlen(s)+off > 0) glCallLists((GLsizei) strlen(s)+off, GL_UNSIGNED_BYTE, (GLubyte *) &s[-off]); }
    
    glPopMatrix(); 
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPopAttrib ();
}


void getColor(double value, int numberOfColors, float* r, float* g, float* b)
{
    if (value > 1) value = 1;
    value = value * (numberOfColors);
    value = (int)(value - 0.00000001);
    value = value / (numberOfColors - 1);
    
    value = 1 - value;
    if (value < 0) value=0;
    if (value > 1) value=1; 
    *r = 3.5*(1-value)*(1-value);
    *g = (1-value)*(value)*3.5;
    *b = value*value;
}

double glScale(double minimum, double maximum, double value)
{
    if (value < minimum)        return 0;
    if (minimum == maximum)     return minimum;
    return (value - minimum) / fabs(maximum - minimum);
}

void glfemDrawNodes(double* x, double* y,int n) 
{
    int i;
    glEnable(GL_POINT_SMOOTH);
    glPointSize(10.0);
    glBegin(GL_POINTS);
    for (i = 0; i < n; i++) {      
        glVertex2f(x[i],y[i]); }
    glEnd();
}
    
void glfemDrawColorElement(double *x, double *y, double *u, int n) 
{   
    GLfloat r,g,b;
    int j;
    glBegin(GL_POLYGON);
    for (j=0; j < n; ++j) {
        getColor(u[j],numberColors,&r,&g,&b);
        glColor3f(r,g,b);
        glVertex2f(x[j],y[j]);}
    glEnd();

}

double interpolateTri(double u[3], double xsi, double eta)
{
    return u[0]*xsi + u[1]*eta + u[2]*(1.0-xsi-eta);
}


void glfemDrawColorRecursiveTriangle(femApproxProblem *theProblem, 
           double *x, double *y, double *xsi, double *eta, int iElem, int n)
{
    int i,j;    
    const int nodes[4][3] = {{0,3,5},{3,1,4},{5,4,2},{3,4,5}};
    const double xsiRef[6] = {0.0,1.0,0.0,0.5,0.5,0.0};
    const double etaRef[6] = {0.0,0.0,1.0,0.0,0.5,0.5};
    double xSub[3], xsiSub[3];
    double ySub[3], etaSub[3];
    double values[3];
        
    if (n <= 0) {
        for (j=0; j<3; j++) {
            double myValue = femApproxEval(theProblem,iElem,xsi[j],eta[j]);
            double minValue = theProblem->minValue;
            double maxValue = theProblem->maxValue;            
            values[j] = glScale(minValue,maxValue,myValue); }                 
            glfemDrawColorElement(x,y,values,3); 
     //       glColor3f(0.0,0.0,0.0); glfemDrawElement(x,y,3); 
             }
    else {
        for (i=0; i<4; i++) {
            for (j=0; j<3; j++) {
                double xsiLoc = xsiRef[nodes[i][j]];
                double etaLoc = etaRef[nodes[i][j]];
                xSub[j] = interpolateTri(x,xsiLoc,etaLoc);
                ySub[j] = interpolateTri(y,xsiLoc,etaLoc);  
                xsiSub[j] = interpolateTri(xsi,xsiLoc,etaLoc);
                etaSub[j] = interpolateTri(eta,xsiLoc,etaLoc); } 
            glfemDrawColorRecursiveTriangle(theProblem,xSub,ySub,xsiSub,etaSub,iElem,n-1); 
            }}
                
         
}


void glfemPlotFieldRecursiveTriangle(femApproxProblem *theProblem, int level)
{
    femMesh *theMesh = theProblem->mesh;   
    int i,j,*nodes;
    double x[3],y[3];
    double xsi[3]  = { 0.0, 1.0, 0.0 };
    double eta[3]  = { 0.0, 0.0, 1.0 };   
    for (i = 0; i < theMesh->nElem; ++i) {
        nodes = &(theMesh->elem[i*3]);
        for (j=0; j < 3; ++j) {
            x[j] = theMesh->X[nodes[j]];
            y[j] = theMesh->Y[nodes[j]]; }              
        glfemDrawColorRecursiveTriangle(theProblem,x,y,xsi,eta,i,level); 
        }
    glColor3f(0.0, 0.0, 0.0);
    glfemPlotMesh(theMesh);
}

void glfemDrawIntegrationNodes(femApproxProblem *theProblem)
{
    femIntegration *theRule = theProblem->rule;
    femMesh *theMesh = theProblem->mesh;   
    double x[3],y[3];

    glEnable(GL_POINT_SMOOTH);
    glPointSize(10.0);
    glBegin(GL_POINTS);
    int i,j,*nodes;
    for (i = 0; i < theMesh->nElem; ++i) {
        nodes = &(theMesh->elem[i*3]);
        for (j=0; j < 3; ++j) {
            x[j] = theMesh->X[nodes[j]];
            y[j] = theMesh->Y[nodes[j]]; }              
        for (int iInteg=0; iInteg < theRule->n; iInteg++) {    
            double xsi    = theRule->xsi[iInteg];
            double eta    = theRule->eta[iInteg];
             double xloc = x[2]*eta + x[1]*xsi + x[0]*(1.0-xsi-eta);
            double yloc = y[2]*eta + y[1]*xsi + y[0]*(1.0-xsi-eta); 
    //        printf("%g %g \n",xloc,yloc);
            glVertex2f(xloc,yloc);  }}     
    glEnd();    
}


void glfemDrawValueNodes(femApproxProblem *theProblem)
{
    femDiscrete *theSpace = theProblem->space;
    femMesh *theMesh = theProblem->mesh;   
    double x[3],y[3];


    double xsis[10], etas[10];    
    femDiscreteXsi2(theSpace,xsis,etas);

    glEnable(GL_POINT_SMOOTH);
    glPointSize(10.0);
    glBegin(GL_POINTS);
    int i,j,*nodes;
    for (i = 0; i < theMesh->nElem; ++i) {
        nodes = &(theMesh->elem[i*3]);
        for (j=0; j < 3; ++j) {
            x[j] = theMesh->X[nodes[j]];
            y[j] = theMesh->Y[nodes[j]]; }              
        for (int index = 0; index < theSpace->n; index++) {  
            double xsi = xsis[index];
            double eta = etas[index];
            double xloc = x[2]*eta + x[1]*xsi + x[0]*(1.0-xsi-eta);
            double yloc = y[2]*eta + y[1]*xsi + y[0]*(1.0-xsi-eta); 
    //        printf("%g %g \n",xloc,yloc);
            glVertex2f(xloc,yloc);  }}     
    glEnd();    
}


void glfemDrawElement(double *x, double *y, int n)  
{
    int j;
    glBegin(GL_LINE_STRIP);         
    for (j = 0; j < n; j++) {      
        glVertex2f(x[j],y[j]); }
    glVertex2f(x[0],y[0]);
    glEnd();

}

void glfemReshapeWindows(femMesh *theMesh, int w, int h)
{
    double minX  = femMin(theMesh->X,theMesh->nNode);
    double maxX  = femMax(theMesh->X,theMesh->nNode);
    double minY  = femMin(theMesh->Y,theMesh->nNode);
    double maxY  = femMax(theMesh->Y,theMesh->nNode);
    double sizeX = (maxX-minX)/1.65;
    double meanX = (maxX+minX)/2.0; 
    double sizeY = (maxY-minY)/1.65;
    double meanY = (maxY+minY)/2.0;
    
    double ratio = (GLfloat) h / (GLfloat) w;
    double size = fmax(sizeX,sizeY);
    double left = meanX - size;
    double right = meanX + size;
    double bottom = meanY- size*ratio;
    double top = meanY + size*ratio;        
        
    glViewport(0,0,w,h);    
    glClearColor( 0.9f, 0.9f, 0.8f, 0.0f );
  //  glClearColor( 1.0f, 1.0f, 1.0f, 0.0f );  // for white plot
    glClear(GL_COLOR_BUFFER_BIT);   
    glViewport(0, 0, (GLsizei) w, (GLsizei) h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(left, right, bottom, top, -5.0, 5.0);
    
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

void glfemPlotField(femMesh *theMesh, double *u)
{
    int i,j,*nodes;
    double xLoc[4];
    double yLoc[4];
    double uLoc[4];
    double uMax = femMax(u,theMesh->nNode);
    double uMin = femMin(u,theMesh->nNode);
    int nLocalNode = theMesh->nLocalNode;
    
    for (i = 0; i < theMesh->nElem; ++i) {
        nodes = &(theMesh->elem[i*nLocalNode]);
        for (j=0; j < nLocalNode; ++j) {
            xLoc[j] = theMesh->X[nodes[j]];
            yLoc[j] = theMesh->Y[nodes[j]];
            uLoc[j] = glScale(uMin,uMax,u[nodes[j]]); }              
        glfemDrawColorElement(xLoc,yLoc,uLoc,nLocalNode); }
    glColor3f(0.0, 0.0, 0.0);
    glfemPlotMesh(theMesh);
}

void glfemPlotMesh(femMesh *theMesh)
{
    int i,j,*nodes;
    double  xLoc[4];
    double  yLoc[4];
    int nLocalNode = theMesh->nLocalNode;

    for (i = 0; i < theMesh->nElem; ++i) {
        nodes = &(theMesh->elem[i*nLocalNode]);
        for (j=0; j < nLocalNode; ++j) {
            xLoc[j] = theMesh->X[nodes[j]];
            yLoc[j] = theMesh->Y[nodes[j]]; }                 
        glfemDrawElement(xLoc,yLoc,nLocalNode); }
}

void glfemMessage(char *message)
{
    glfemDrawMessage(20,460,message);   
}

GLFWwindow* glfemInit(char *theWindowName)
{
    glfwInit();
    GLFWwindow* window = glfwCreateWindow(480,480,"Simple example with graphics",NULL,NULL);    
    glfwMakeContextCurrent(window);
    glfemSetRasterSize(480,480);
    glfwSetWindowTitle(window,theWindowName);
    glShadeModel(GL_SMOOTH);
    glMakeRasterFont();
    return(window);
}

void glfemPlotEdges(femEdges *theEdges)
{
    int i;
    femMesh* theMesh = theEdges->mesh;
    glBegin(GL_LINES);
    for (i = 0; i < theEdges->nEdge; i++) {
        glVertex2d(theMesh->X[theEdges->edges[i].node[0]], theMesh->Y[theEdges->edges[i].node[0]]);
        glVertex2d(theMesh->X[theEdges->edges[i].node[1]], theMesh->Y[theEdges->edges[i].node[1]]); }
    glEnd();
}

void glfemPlotBnd(femEdges *theEdges)
{
    int i;
    femMesh* theMesh = theEdges->mesh;
    glBegin(GL_LINES);
    for (i = 0; i < theEdges->nEdge; i++) {
        if (theEdges->edges[i].elem[1] < 0) {
            glVertex2d(theMesh->X[theEdges->edges[i].node[0]], theMesh->Y[theEdges->edges[i].node[0]]);
            glVertex2d(theMesh->X[theEdges->edges[i].node[1]], theMesh->Y[theEdges->edges[i].node[1]]); }}
    glEnd();
}



#define _CRT_SECURE_NO_DEPRECATE
#include <stdio.h>
#include <stdlib.h>
#include <GL\glut.h>
#include <math.h>
#include <iostream>

#define ON 1
#define OFF 0

typedef float Matrix4x4[4][4]; 

typedef struct _point {
	float x, y, z;
} point;

typedef struct _faceStruct {
	int v1, v2, v3;
	int n1, n2, n3;
} faceStruct;


void meshReader(char *filename, int sign);
void	display(void);
void drawObject();
void drawAxis();
void draw_blue_tetraheadron();
void draw_red_rectangle();
void	resize(int x, int y);
void OnMouseMove(int x, int y);
void OnMouseDown(int button, int state, int x, int y);
void	keyboard(unsigned char key, int x, int y);
void setIdentity(Matrix4x4 m);
void matrixMultiply(Matrix4x4 m1, Matrix4x4 m2);
void matrixTranslate(float tx, float ty, float tz, Matrix4x4 m);
void matrixRotate(_point p1, _point p2, float radianAngle, Matrix4x4 m);
void matrixScale(float sx, float sy, float sz, _point fixedPt, Matrix4x4 m);
void printMatrix(Matrix4x4 m);
void arrayToMatrix(GLfloat* array, Matrix4x4 m);
void matrixToarray(Matrix4x4 m, float array[]);
void MatrixAssignment(Matrix4x4 m1, Matrix4x4 m2);
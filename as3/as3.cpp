// EECS366 Assignment 2 
//James Du, Haifeng Chen

#include "as3.h"
#define TRANSLATE_INDEX 0.5
#define ROTATION_INDEX_RADIAN 0.05
#define SCALE_INDEX 0.1


// Global variables
int window_width, window_height;    // Window dimensions
int PERSPECTIVE = ON;
int OBJECT_ON = ON;
int AXIS_ON = ON;
int degrees = 0;
float phi = 1.570796;
float theta = 0;
float rho = 10;
float centerX = 0;
float centerY = 0;
float centerZ = 0;

float eyeX = -4;
float eyeY = 0;
float eyeZ = 0;

float upX = 0;
float upY = 0;
float upZ = 1;
int oldX, oldY;
int rotateHead, zoom, moveBody;
int toggleCamera = 0;
int snapCamera = 0;
float moveX, moveY, zoomVal;
float scale = 1;
Matrix4x4 modelMatrix, localRotation, worldRotation;
point worldTranslation, distanceToOrigin; //faking as a 3d vector


int verts, faces, norms;    // Number of vertices, faces and normals in the system
point *vertList, *normList; // Vertex and Normal Lists
faceStruct *faceList;	    // Face List

							// The mesh reader itself
							// It can read *very* simple obj files
typedef struct _Vector3f {
	float x, y, z;
} Vector3f;

void meshReader(char *filename, int sign)
{
	float x, y, z, len;
	int i;
	char letter;
	point v1, v2, crossP;
	int ix, iy, iz;
	int *normCount;
	FILE *fp;

	fp = fopen(filename, "r");
	if (fp == NULL) {
		printf("Cannot open %s\n!", filename);
		return;
	}

	// Count the number of vertices and faces
	while (!feof(fp))
	{
		fscanf(fp, "%c %f %f %f\n", &letter, &x, &y, &z);
		if (letter == 'v')
		{
			verts++;
		}
		else
		{
			faces++;
		}
	}

	fclose(fp);

	printf("verts : %d\n", verts);
	printf("faces : %d\n", faces);

	// Dynamic allocation of vertex and face lists
	faceList = (faceStruct *)malloc(sizeof(faceStruct)*faces);
	vertList = (point *)malloc(sizeof(point)*verts);
	normList = (point *)malloc(sizeof(point)*verts);

	fp = fopen(filename, "r");

	// Read the veritces
	for (i = 0; i < verts; i++)
	{
		fscanf(fp, "%c %f %f %f\n", &letter, &x, &y, &z);
		vertList[i].x = x;
		vertList[i].y = y;
		vertList[i].z = z;
	}

	// Read the faces
	for (i = 0; i < faces; i++)
	{
		fscanf(fp, "%c %d %d %d\n", &letter, &ix, &iy, &iz);
		faceList[i].v1 = ix - 1;
		faceList[i].v2 = iy - 1;
		faceList[i].v3 = iz - 1;
	}
	fclose(fp);


	// The part below calculates the normals of each vertex
	normCount = (int *)malloc(sizeof(int)*verts);
	for (i = 0; i < verts; i++)
	{
		normList[i].x = normList[i].y = normList[i].z = 0.0;
		normCount[i] = 0;
	}

	for (i = 0; i < faces; i++)
	{
		v1.x = vertList[faceList[i].v2].x - vertList[faceList[i].v1].x;
		v1.y = vertList[faceList[i].v2].y - vertList[faceList[i].v1].y;
		v1.z = vertList[faceList[i].v2].z - vertList[faceList[i].v1].z;
		v2.x = vertList[faceList[i].v3].x - vertList[faceList[i].v2].x;
		v2.y = vertList[faceList[i].v3].y - vertList[faceList[i].v2].y;
		v2.z = vertList[faceList[i].v3].z - vertList[faceList[i].v2].z;

		crossP.x = v1.y*v2.z - v1.z*v2.y;
		crossP.y = v1.z*v2.x - v1.x*v2.z;
		crossP.z = v1.x*v2.y - v1.y*v2.x;

		len = sqrt(crossP.x*crossP.x + crossP.y*crossP.y + crossP.z*crossP.z);

		crossP.x = -crossP.x / len;
		crossP.y = -crossP.y / len;
		crossP.z = -crossP.z / len;

		normList[faceList[i].v1].x = normList[faceList[i].v1].x + crossP.x;
		normList[faceList[i].v1].y = normList[faceList[i].v1].y + crossP.y;
		normList[faceList[i].v1].z = normList[faceList[i].v1].z + crossP.z;
		normList[faceList[i].v2].x = normList[faceList[i].v2].x + crossP.x;
		normList[faceList[i].v2].y = normList[faceList[i].v2].y + crossP.y;
		normList[faceList[i].v2].z = normList[faceList[i].v2].z + crossP.z;
		normList[faceList[i].v3].x = normList[faceList[i].v3].x + crossP.x;
		normList[faceList[i].v3].y = normList[faceList[i].v3].y + crossP.y;
		normList[faceList[i].v3].z = normList[faceList[i].v3].z + crossP.z;
		normCount[faceList[i].v1]++;
		normCount[faceList[i].v2]++;
		normCount[faceList[i].v3]++;
	}
	for (i = 0; i < verts; i++)
	{
		normList[i].x = (float)sign*normList[i].x / (float)normCount[i];
		normList[i].y = (float)sign*normList[i].y / (float)normCount[i];
		normList[i].z = (float)sign*normList[i].z / (float)normCount[i];
	}

}



// The display function. It is called whenever the window needs
// redrawing (ie: overlapping window moves, resize, maximize)
// You should redraw your polygons here
void	display(void)
{
	// Clear the background
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	
	if (PERSPECTIVE) {
		// Perpective Projection 
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		gluPerspective(60, (GLdouble)window_width / window_height, 0.01, 10000);
		glutSetWindowTitle("Assignment 2 Template (perspective)");
		// Set the camera position, orientation and target
		/*
		eyeX = rho * cos(theta)*sin(phi);
		eyeY = rho * sin(theta)*sin(phi);
		eyeZ = rho * cos(phi);
		gluLookAt(eyeX, eyeY, eyeZ, 0, 0, 0, 0, 1, 0);
		*/
	}
	else {
		// Orthogonal Projection 
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		glOrtho(-0.25*rho, 0.25*rho, -0.25*rho, 0.25*rho, -10000, 10000);
		glutSetWindowTitle("Assignment 2 Template (orthogonal)");
	}

	//viewMatrix= viewTransform();
	if (!snapCamera) {
		centerX = 1 * cos(theta)*sin(phi) + eyeX;
		centerY = 1 * sin(theta)*sin(phi) + eyeY;
		centerZ = 1 * cos(phi) + eyeZ;
	}

	printf("center x: %f, y: %f, z: %f\n", centerX, centerY, centerZ);
	Vector3f f = { centerX - eyeX, centerY - eyeY, centerZ - eyeZ };
	float f_mag = sqrt(f.x * f.x + f.y * f.y + f.z * f.z);
	Vector3f f_n = { f.x / f_mag, f.y / f_mag, f.z / f_mag };
	Vector3f s = { f_n.y * upZ - f_n.z * upY, f_n.z * upX - f_n.x * upZ, f_n.x * upY - f_n.y * upX };
	float s_mag = sqrt(s.x * s.x + s.y * s.y + s.z * s.z);
	Vector3f s_n = { s.x / s_mag, s.y / s_mag, s.z / s_mag };
	Vector3f u = { s_n.y * f_n.z - s_n.z * f_n.y, s_n.z * f_n.x - s_n.x * f_n.z, s_n.x * f_n.y - s_n.y * f_n.x };
	GLfloat M[16] = { s_n.x, u.x, -1 * f_n.x, 0,
		s_n.y, u.y, -1 * f_n.y, 0,
		s_n.z, u.z, -1 * f_n.z, 0,
		0,   0,    0,         1 };
	eyeX = eyeX + s_n.x * moveX;
	eyeY = eyeY + s_n.y * moveX;
	eyeZ = eyeZ + s_n.z * moveX;

	eyeX = eyeX + u.x * moveY;
	eyeY = eyeY + u.y * moveY;
	eyeZ = eyeZ + u.z * moveY;

	if (!snapCamera) {
		centerX = centerX + s_n.x * moveX;
		centerY = centerY + s_n.y * moveX;
		centerZ = centerZ + s_n.z * moveX;

		centerX = centerX + u.x * moveY;
		centerY = centerY + u.y * moveY;
		centerZ = centerZ + u.z * moveY;
	}
	
	eyeX = eyeX + f_n.x * zoomVal;
	eyeY = eyeY + f_n.y * zoomVal;
	eyeZ = eyeZ + f_n.z * zoomVal;

	for (int i = 0; i < 16; i++) {
		//printf("the item is: %f\n", M[i]);
	}
	GLfloat T[16] = { 1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0, 1, 0,
		-eyeX, -eyeY, -eyeZ,  1 };
	
	glMultMatrixf(M);
	glMultMatrixf(T);
	


	GLfloat tempArray[16];
	glMatrixMode(GL_MODELVIEW);
	//printMatrix(modelMatrix);
	
	matrixMultiply(worldRotation, modelMatrix);
	
	matrixTranslate(worldTranslation, modelMatrix);

	matrixLocalRotate(localRotation, modelMatrix);

	//matrixMultiply(localRotation, modelMatrix);
	
	matrixToarray(modelMatrix, tempArray);
	

//	Matrix4x4 tempMatrix;
//	glGetFloatv(GL_MODELVIEW_MATRIX, tempArray);
//	arrayToMatrix(tempArray, tempMatrix);
	glLoadMatrixf(tempArray);
	//printMatrix(modelMatrix);

	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	drawObject();
	drawAxis();  //local axis
	
	glLoadIdentity();
	drawAxis();    //global axis

	//update and reinitiate 
	/*
	distanceToOrigin.x += worldTranslation.x;
	distanceToOrigin.y += worldTranslation.y;
	distanceToOrigin.z += worldTranslation.z;
	printf("x %f y %f z %f", distanceToOrigin.x, distanceToOrigin.y, distanceToOrigin.z);
	printMatrix(modelMatrix);
	matrixMultiply(worldRotation, distanceToOrigin);
	*/
	init();
	// (Note that the origin is lower left corner)
	// (Note also that the window spans (0,1) )
	// Finish drawing, update the frame buffer, and swap buffers
	glutSwapBuffers();
}

void drawObject() {
	if (OBJECT_ON) {
		glColor3f(1, 0, 0);
		glBegin(GL_POLYGON);
		for (int i = 0; i < faces; i++) {
			glNormal3f(normList[(faceList[i].v1)].x / scale, normList[(faceList[i].v1)].y / scale, normList[(faceList[i].v1)].z / scale); // normals
			glVertex3f(vertList[(faceList[i].v1)].x / scale, vertList[(faceList[i].v1)].y / scale, vertList[(faceList[i].v1)].z / scale);

			glNormal3f(normList[(faceList[i].v2)].x / scale, normList[(faceList[i].v2)].y / scale, normList[(faceList[i].v2)].z / scale); // normals
			glVertex3f(vertList[(faceList[i].v2)].x / scale, vertList[(faceList[i].v2)].y / scale, vertList[(faceList[i].v2)].z / scale);

			glNormal3f(normList[(faceList[i].v3)].x / scale, normList[(faceList[i].v3)].y / scale, normList[(faceList[i].v3)].z / scale); // normals
			glVertex3f(vertList[(faceList[i].v3)].x / scale, vertList[(faceList[i].v3)].y / scale, vertList[(faceList[i].v3)].z / scale);
		}
		glEnd();
	}
}



void drawAxis() {
	if (AXIS_ON) {
		// Draw a green line
		glColor3f(0, 1, 0);
		glBegin(GL_LINES);

		glVertex3f(0.0, 0.0, 0.0);
		glVertex3f(4, 0, 0);
		glEnd();

		// Draw a green line
		glColor3f(0, 0, 1);
		glBegin(GL_LINES);
		glVertex3f(0.0, 0.0, 0.0);
		glVertex3f(0, 4, 0);
		glEnd();

		// Draw a green line
		glColor3f(1, 0, 0);
		glBegin(GL_LINES);
		glVertex3f(0.0, 0.0, 0.0);
		glVertex3f(0, 0, 4);
		glEnd();
	}
}

void draw_blue_tetraheadron() {
	glColor3f(0, 0, 1);
	glBegin(GL_TRIANGLES);
	glVertex3f(0.0, 1.6, 0.0);
	glVertex3f(0.8, -0.4, 0.8);
	glVertex3f(-0.8, -0.4, 0.8);

	glVertex3f(0.0, 1.6, 0.0);
	glVertex3f(0.8, -0.4, 0.8);
	glVertex3f(0.0, -0.4, -0.8);

	glVertex3f(0.0, 1.6, 0.0);
	glVertex3f(0.0, -0.4, -0.8);
	glVertex3f(-0.8, -0.4, 0.8);

	glVertex3f(-0.8, -0.4, 0.8);
	glVertex3f(0.8, -0.4, 0.8);
	glVertex3f(0.0, -0.4, -0.8);
	glEnd();
}

void draw_red_rectangle() {
		glColor3f(1,0,0);
		glBegin(GL_POLYGON);
		glVertex3f(0.8,0.8,-0.8);
		glVertex3f(0.8,-0.8,-0.8);
		glVertex3f(-0.8,-0.8,-0.0);
		glVertex3f(-0.8,0.8,-0.0);
	   glEnd();
}

// This function is called whenever the window is resized. 
// Parameters are the new dimentions of the window
void	resize(int x, int y)
{
	glViewport(0, 0, x, y);
	window_width = x;
	window_height = y;
	if (PERSPECTIVE) {
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		gluPerspective(60, (GLdouble)window_width / window_height, 0.01, 10000);
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
	}
	printf("Resized to %d %d\n", x, y);
}


//This function is called whenever the mouse is moved with a mouse button held down.
// x and y are the location of the mouse (in window-relative coordinates)
void OnMouseMove(int x, int y) {
	if (rotateHead) {
		moveX = 0;
		moveY = 0;
		theta += (x - oldX)*0.01f;
		//prevent the camera from going upside down (limit how far up it can move)
		//if ((phi + (y - oldY) * 0.01f  > 1.570796 && y > oldY) || (phi + (y - oldY) * 0.01f < -1.570796 && y < oldY)) {
		//do nothing
		//}f
		//else {
		phi += (y - oldY)*0.01f;
		//}
		degrees = degrees + x - oldX;
		glutPostRedisplay();
	}
	else if (moveBody) {
		moveX = (x - oldX) * .01f;
		moveY = (y - oldY) * .01f;
		glutPostRedisplay();
	}
	else if (zoom) {
		moveX = 0;
		moveY = 0;
		zoomVal = (y - oldY) * 0.01f;
		glutPostRedisplay();
	}
	oldX = x;
	oldY = y;

}

//determine which mouse button is pressed
void OnMouseDown(int button, int state, int x, int y) {
	rotateHead = 0;
	moveBody = 0;
	zoomVal = 0;
	zoom = 0;
	if (button == GLUT_LEFT_BUTTON) {
		snapCamera = 0;
		oldX = x;
		oldY = y;
		rotateHead = 1;
	}
	else if (button == GLUT_RIGHT_BUTTON) {
		snapCamera = 0;
		zoom = 1;
		oldY = y;
	}
	else if (button == GLUT_MIDDLE_BUTTON) {
		snapCamera = 0;
		moveBody = 1;
		oldX = x;
		oldY = y;
	}
}
// This function is called whenever there is a keyboard input
// key is the ASCII value of the key pressed
// x and y are the location of the mouse
void keyboard(unsigned char key, int x, int y)
{
	switch (key) {
	case 'a':
		AXIS_ON = !AXIS_ON;
		break;
	case 'q':                           /* Quit */
		exit(1);
		break;
	case 'p':
	case 'P':
		// Toggle Projection Type (orthogonal, perspective)
		if (PERSPECTIVE) {
			// switch from perspective to orthogonal
			PERSPECTIVE = OFF;
		}
		else {
			// switch from orthogonal to perspective
			PERSPECTIVE = ON;
		}
		break;
	case 's':
		OBJECT_ON = !OBJECT_ON;
		break;
	case '4': //negative translate along x axis
		worldTranslation.x -= TRANSLATE_INDEX;
		break;
	case '6':
		worldTranslation.x += TRANSLATE_INDEX;
		break;
	case '8':
		worldTranslation.y += TRANSLATE_INDEX;
		break;
	case '2':
		worldTranslation.y -= TRANSLATE_INDEX;
		break;
	case '9':
		worldTranslation.z += TRANSLATE_INDEX;
		break;
	case '1':
		worldTranslation.z -= TRANSLATE_INDEX;
		break;
	case '[':
		matrixRotate(-ROTATION_INDEX_RADIAN, 'x', worldRotation);
		break;
	case ']':
		matrixRotate(ROTATION_INDEX_RADIAN, 'x', worldRotation);
		break;
	case ';':
		matrixRotate(-ROTATION_INDEX_RADIAN, 'y', worldRotation);
		break;
	case '\'':
		matrixRotate(ROTATION_INDEX_RADIAN, 'y' , worldRotation);
		break;
	case '.':
		matrixRotate(-ROTATION_INDEX_RADIAN, 'z', worldRotation);
		break;
	case '/':
		matrixRotate(ROTATION_INDEX_RADIAN, 'z', worldRotation);
		break;
	case '=':
		scale += SCALE_INDEX;
		break;
	case '-':
		scale -= SCALE_INDEX;
		break;
	case 'i':
		matrixRotate(-ROTATION_INDEX_RADIAN, 'x', localRotation);
		break;
	case 'o':
		matrixRotate(ROTATION_INDEX_RADIAN, 'x', localRotation);
		break;
	case 'k':  
		matrixRotate(-ROTATION_INDEX_RADIAN, 'y', localRotation);
		break;
	case 'l':  
		matrixRotate(ROTATION_INDEX_RADIAN, 'y', localRotation);
		break;
	case 'm': 
		matrixRotate(-ROTATION_INDEX_RADIAN, 'z', localRotation);
		break;
	case ',':
		matrixRotate(ROTATION_INDEX_RADIAN, 'z', localRotation);
		break;
	case 'c':
		toggleCamera = !toggleCamera;
		snapCamera = 1;
		if (toggleCamera) {
			centerX = 0;
			centerY = 0;
			centerZ = 0;
		}
		else {
			updateDistanceToOrigin();
			centerX = distanceToOrigin.x;
			centerY = distanceToOrigin.y;
			centerZ = distanceToOrigin.z;
		}
		glutPostRedisplay();
		break;
	
	
	default:
		printf("wtf");
		break;
	}

	// Schedule a new display event
	glutPostRedisplay();
}


/* Construct the 4 x 4 identity matrix. */
void setIdentity(Matrix4x4 m)
{
	int row, col;

	for (row = 0; row < 4; row++)
		for (col = 0; col < 4; col++)
			m[row][col] = (row == col);
}

/* Premultiply matrix m1 by matrix m2, store result in m2. */

void matrixMultiply(Matrix4x4 m1, Matrix4x4 m2)
{
	int row, col;
	Matrix4x4 matTemp;

	for (row = 0; row < 4; row++)
		for (col = 0; col < 4; col++)
			matTemp[row][col] = m1[row][0] * m2[0][col] + m1[row][1] *
			m2[1][col] + m1[row][2] * m2[2][col] +
			m1[row][3] * m2[3][col];
	for (row = 0; row < 4; row++)
		for (col = 0; col < 4; col++)
			m2[row][col] = matTemp[row][col];
}

/*
	matrix multiply to a vector, store result in p
*/
void matrixMultiply(Matrix4x4 m, point p) {
	p.x = m[0][0] * p.x + m[0][1] * p.y + m[0][2] * p.z + m[0][3];
	p.y = m[1][0] * p.x + m[1][1] * p.y + m[1][2] * p.z + m[1][3];
	p.z = m[2][0] * p.x + m[2][1] * p.y + m[2][2] * p.z + m[2][3];
}


/*  Procedure for generating 3D translation matrix.  */

void matrixTranslate(point translation ,Matrix4x4 m)
{
	Matrix4x4 matTransl3D;
	float tx = translation.x;
	float ty = translation.y;
	float tz = translation.z;

	//  Initialize translation matrix to identity.  
	setIdentity(matTransl3D);

	matTransl3D[0][3] = tx;
	matTransl3D[1][3] = ty;
	matTransl3D[2][3] = tz;

	//  Concatenate matTransl3D with composite matrix.  
	matrixMultiply(matTransl3D, m);
}

void matrixLocalRotate(Matrix4x4 localRotation, Matrix4x4 modelMatrix) {
	updateDistanceToOrigin();
	matrixTranslate(reversePoint(distanceToOrigin), modelMatrix);
//	printMatrix(modelMatrix);
	matrixMultiply(localRotation, modelMatrix);
	matrixTranslate(distanceToOrigin, modelMatrix);
	
}

void updateDistanceToOrigin() {
	distanceToOrigin.x = modelMatrix[0][3];
	distanceToOrigin.y = modelMatrix[1][3];
	distanceToOrigin.z = modelMatrix[2][3];
}

point reversePoint(point p){
	point temp;
	temp.x = -p.x;
	temp.y = -p.y;
	temp.z = -p.z;
	return temp;
}

/*  Procedure for generating a rotation matrix.  */

void matrixRotate(float radian, char axis, Matrix4x4 m) {
	float matrixX[4][4] = { 
						  {1,0,0,0},
						  {0,cos(radian),-1.0f * sin(radian),0},
						  {0,sin(radian),cos(radian),0},
						  {0,0,0,1} };
	float matrixY[4][4] = {
						  {cos(radian),0,sin(radian),0},
						  {0,1,0,0},
						  {-1.0f*sin(radian),0,cos(radian),0},
						  {0,0,0,1} };
							     
	float matrixZ[4][4] = { 
						  {cos(radian),-1.0f*sin(radian),0,0},
						  {sin(radian),cos(radian),0,0},
						  {0,0,1,0},
						  {0,0,0,1} };

	switch (axis) {
	case 'x':
		matrixMultiply(matrixX, m);
		break;
	case 'y':
		matrixMultiply(matrixY, m);
		break;
	case 'z':
		matrixMultiply(matrixZ, m);
		break;
	default:
		printf("asfdsad");
	}
}


/*
	m1 assign to m2
*/
void MatrixAssignment(Matrix4x4 m1, Matrix4x4 m2) {
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			m2[i][j] = m1[i][j];
		}
	}
}

/*
   element in matrix transfter to column major array
*/
void matrixToarray(Matrix4x4 m, float array[]) {
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			array[i + 4 * j] = m[i][j];
		}
	}
}

/*
   element in column major array transfter to matrix m
*/
void arrayToMatrix(GLfloat* array, Matrix4x4 m) {
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			m[i][j] = (float)array[i + 4 * j];
		}
	}

}

void printMatrix(Matrix4x4 m) {
	std::cout << "current matrix\n";
	std::cout << m[0][0] << " " << m[0][1] << " " << m[0][2] << " " << m[0][3] << '\n';
	std::cout << m[1][0] << " " << m[1][1] << " " << m[1][2] << " " << m[1][3] << '\n';
	std::cout << m[2][0] << " " << m[2][1] << " " << m[2][2] << " " << m[2][3] << '\n';
	std::cout << m[3][0] << " " << m[3][1] << " " << m[3][2] << " " << m[3][3] << '\n';
}

void printArray(float* a) {
	for (int i = 0; i < 16; i++) {
		printf("%f ", a[i]);
	}
	printf("\n");
}

void init() {
	setIdentity(worldRotation);
	setIdentity(localRotation);
	worldTranslation.x = 0;
	worldTranslation.y = 0;
	worldTranslation.z = 0;
	
}


void errorCheck() {
	GLenum glErr;
	glErr = glGetError();
	printf("%s\n", gluErrorString(glErr));
}

int main(int argc, char* argv[])
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
	glutCreateWindow("Assignment 2 Template (orthogonal)");
	glutDisplayFunc(display);
	glutReshapeFunc(resize);
	glutMouseFunc(OnMouseDown);
	glutMotionFunc(OnMouseMove);
	glutKeyboardFunc(keyboard);
	/*
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-2.5, 2.5, -2.5, 2.5, -10000, 10000);
	glEnable(GL_DEPTH_TEST);
	*/
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	meshReader("teapot.obj", 1);

	setIdentity(modelMatrix);
	distanceToOrigin.x = 0;
	distanceToOrigin.y = 0;
	distanceToOrigin.z = 0;
	init();
	glutMainLoop();
	return 0;
}

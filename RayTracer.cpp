/*========================================================================
* COSC 363  Computer Graphics (2018)
* Ray tracer 
* By Alex Bull
* Boilplate/start by Department of Computer Science and Software 
* Engineering, University of Canterbury,
*=========================================================================
*/
#include <iostream>
#include <cmath>
#include <vector>
#include <glm/glm.hpp>
#include "Sphere.h"
#include "SceneObject.h"
#include "Ray.h"
#include <GL/glut.h>
#include "Plane.h"
#include "TextureBMP.h"
using namespace std;

const float WIDTH = 20.0;  
const float HEIGHT = 20.0;
const float EDIST = 40.0;
const int NUMDIV = 600;
const int MAX_STEPS = 5;
const float XMIN = -WIDTH * 0.5;
const float XMAX =  WIDTH * 0.5;
const float YMIN = -HEIGHT * 0.5;
const float YMAX =  HEIGHT * 0.5;

vector<SceneObject*> sceneObjects;  //A global list containing pointers to objects in the scene

TextureBMP texture;
TextureBMP texture2;


//---The most important function in a ray tracer! ---------------------------------- 
//   Computes the colour value obtained by tracing a ray and finding its 
//     closest point of intersection with objects in the scene.
//----------------------------------------------------------------------------------
glm::vec3 trace(Ray ray, int step)
{
    float phongsConst = 20;
    float eta = 1 / 1.5;
	glm::vec3 backgroundCol(0);
	glm::vec3 light(15, 20, -20);
	glm::vec3 ambientCol(0.2); 
    glm::vec3 specCol(0);
    glm::vec3 colorSum(0);

    ray.closestPt(sceneObjects);		//Compute the closest point of intersetion of objects with the ray

    if(ray.xindex == -1) return backgroundCol;      //If there is no intersection return background colour

    glm::vec3 materialCol = sceneObjects[ray.xindex]->getColor(); //else return object's colour
    
    float a1 = -20;
    float a2 = 20;
    float b1 = -80;
    float b2 = -200;
    
    if(ray.xindex == 2){
        float texcoords = (ray.xpt.x - a1)/(a2-a1);
        float texcoordt = (ray.xpt.z - b1)/(b2-b1);
        materialCol = texture.getColorAt(texcoords, texcoordt);
    }
    
    glm::vec3 normalVector = sceneObjects[ray.xindex]->normal(ray.xpt); //normal vector at the point of intersection
    
    if(ray.xindex == 11){
        float texcoords = asin(normalVector.x)/M_PI + 0.5;
        float texcoordt = asin(normalVector.y)/M_PI + 0.5;
        materialCol = texture2.getColorAt(texcoords, texcoordt);
    }
    
    //Simple procedural texture
    if(ray.xindex == 12){
        if ((int(ray.xpt.y) + int(ray.xpt.z)) % 2 == 0){
            materialCol = glm::vec3(0, 0, 0);
        } else {
            materialCol = glm::vec3(0.8, 0, 0.8);
        }
    }
    
    //Diffuse reflections
    glm::vec3 lightVector = light - ray.xpt; // vector from the point of intersection towards the light source
    glm::vec3 lightUnitVector = glm::normalize(lightVector); // unit vector of lightVector
    float lDotn = glm::dot(lightUnitVector, normalVector); // dot product of the unit vector and the normal vector
    
    //Specular reflections
    glm::vec3 reflVector = glm::reflect(-lightUnitVector, normalVector); // reflection vector
    float rDotv = glm::dot(reflVector, -ray.dir); // dot product for spec. reflections calc.
    float rDotvf = pow(rDotv, phongsConst); // to the power of Phong's constant
    
    if(rDotv < 0){
        specCol = glm::vec3(0, 0, 0);
    } else {
        specCol = glm::vec3(1, 1, 1) * rDotvf;
    }
    
    
    //Shadows
    Ray shadow(ray.xpt, lightUnitVector);
    shadow.closestPt(sceneObjects);
    
    float lightDist = glm::distance(ray.xpt, light);
    
    if(lDotn <= 0 || (shadow.xindex > -1 && (shadow.xdist < lightDist))){
        colorSum += ambientCol * materialCol;
    } else {
        colorSum += (ambientCol * materialCol) + (lDotn * materialCol) + specCol; // sum of accumalated color values
    }
    
    //For reflective object(s)
    if(ray.xindex == 0 && step < MAX_STEPS){
        glm::vec3 reflectedDir = glm::reflect(ray.dir, normalVector);
        Ray reflectedRay(ray.xpt, reflectedDir);
        glm::vec3 reflectedCol = trace(reflectedRay, step+1);  //Recursion
        colorSum = colorSum + (1.0f*reflectedCol);
    }
    
    if(ray.xindex == 11 && step < MAX_STEPS){
        glm::vec3 reflectedDir = glm::reflect(ray.dir, normalVector);
        Ray reflectedRay(ray.xpt, reflectedDir);
        glm::vec3 reflectedCol = trace(reflectedRay, step+1);  //Recursion
        colorSum = colorSum + (0.3f*reflectedCol);
    }
    
    //for refractive object(s)
	if (ray.xindex == 9 && step < MAX_STEPS) {
        glm::vec3 m;
        glm::vec3 h;
        Ray refrRay2;
        glm::vec3 refractCol;
        
        glm::vec3 n = sceneObjects[ray.xindex]->normal(ray.xpt);
		glm::vec3 g = glm::refract(ray.dir, n, eta);
        
		Ray refrRay1(ray.xpt, g);
		refrRay1.closestPt(sceneObjects);
		
        if (refrRay1.xindex == 9){
            m = sceneObjects[refrRay1.xindex]->normal(refrRay1.xpt);
            h = glm::refract(g, -m, 1.0f/eta);
            refrRay2 = Ray(refrRay1.xpt, h);
            refractCol = trace(refrRay2, step + 1);
        } else {
            refractCol = trace(refrRay1, step + 1);
        }
        colorSum += (refractCol * 1.0f);
    }
    
    
    //for transparent object(s)
	if (ray.xindex == 10 && step < MAX_STEPS) {

        glm::vec3 refractCol;
        
		Ray refrRay(ray.xpt, ray.dir);
		refrRay.closestPt(sceneObjects);
		
        refractCol = trace(refrRay, step + 1);

        colorSum += (refractCol * 1.0f);
    }
    
    
    
    return colorSum;
    
}

//A better trace function/calls the trace function multiple times through four divisions of a cell
glm::vec3 supersampleTrace(glm::vec3 eye, glm::vec3 dir, float xp, float yp, float cellWidth, float cellHeight)
{
    float centerLeftX = xp + (cellWidth / 2);
    float centerLowerY = yp + (cellHeight / 2);
    int divisions = 4;
    float red = 0;
    float blue = 0;
    float green = 0;
    
    
    glm::vec3 divisionCols[divisions];
    
    glm::vec3 rayDir;
    
    Ray ray;
   
    rayDir = glm::vec3(centerLeftX, centerLowerY, -EDIST);
    ray = Ray(eye, rayDir);
    ray.normalize();
    divisionCols[0] = trace(ray, 1);
    rayDir = glm::vec3(centerLeftX + (cellWidth / 2), centerLowerY, -EDIST);
    ray = Ray(eye, rayDir);
    ray.normalize();
    divisionCols[1] = trace(ray, 1);
    rayDir = glm::vec3(centerLeftX + (cellWidth / 2), centerLowerY + (cellHeight / 2), -EDIST);
    ray = Ray(eye, rayDir);
    ray.normalize();
    divisionCols[2] = trace(ray, 1);
    rayDir = glm::vec3(centerLeftX, centerLowerY + (cellHeight / 2), -EDIST);
    ray = Ray(eye, rayDir);
    ray.normalize();
    divisionCols[3] = trace(ray, 1);
    
    for (int i = 0; i < divisions; i++) {
        red += divisionCols[i][0];
        green += divisionCols[i][1];
        blue += divisionCols[i][2];
    }
    
    return glm::vec3(red/divisions, green/divisions, blue/divisions);
    
}


//---The main display module -----------------------------------------------------------
// In a ray tracing application, it just displays the ray traced image by drawing
// each cell as a quad.
//---------------------------------------------------------------------------------------
void display()
{
	float xp, yp;  //grid point
	float cellX = (XMAX-XMIN)/NUMDIV;  //cell width
	float cellY = (YMAX-YMIN)/NUMDIV;  //cell height

	glm::vec3 eye(0., 0., 0.);  //The eye position (source of primary rays) is the origin

	glClear(GL_COLOR_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

	glBegin(GL_QUADS);  //Each cell is a quad.

	for(int i = 0; i < NUMDIV; i++)  	//For each grid point xp, yp
	{
		xp = XMIN + i*cellX;
		for(int j = 0; j < NUMDIV; j++)
		{
			yp = YMIN + j*cellY;

		    glm::vec3 dir(xp+0.5*cellX, yp+0.5*cellY, -EDIST);	//direction of the primary ray

		    Ray ray = Ray(eye, dir);		//Create a ray originating from the camera in the direction 'dir'
			ray.normalize();				//Normalize the direction of the ray to a unit vector
		    //glm::vec3 col = trace (ray, 1); //Trace the primary ray and get the colour value (no anti-aliasing)

            glm::vec3 col = supersampleTrace(eye, dir, xp, yp, cellX, cellY);

			glColor3f(col.r, col.g, col.b);
			glVertex2f(xp, yp);				//Draw each cell with its color value
			glVertex2f(xp+cellX, yp);
			glVertex2f(xp+cellX, yp+cellY);
			glVertex2f(xp, yp+cellY);
        }
    }

    glEnd();
    glFlush();
}


//---This function initializes the scene ------------------------------------------- 
//   Specifically, it creates scene objects (spheres, planes, cones, cylinders etc)
//     and add them to the list of scene objects.
//   It also initializes the OpenGL orthographc projection matrix for drawing the
//     the ray traced image.
//----------------------------------------------------------------------------------
void initialize()
{
    glMatrixMode(GL_PROJECTION);
    gluOrtho2D(XMIN, XMAX, YMIN, YMAX);
    glClearColor(0, 0, 0, 1);
    
    char file[] = "floor.bmp";
    texture = TextureBMP(file);
    char file2[] = "marble.bmp";
    texture2 = TextureBMP(file2);

	//-- Create a pointer to a sphere object
	Sphere *sphere1 = new Sphere(glm::vec3(-8.0, -5.0, -90.0), 10.0, glm::vec3(0.6, 0.6, 0.6));
    sceneObjects.push_back(sphere1); 
        
    Sphere *sphere2 = new Sphere(glm::vec3(5.0, 5.0, -80.0), 5.0, glm::vec3(1, 0, 0));
    sceneObjects.push_back(sphere2); 
    
    //Pointer to plane object
    Plane *plane = new Plane (
    glm::vec3(-20., -20, -80),  //Point A
    glm::vec3(20., -20, -80),   //Point B
    glm::vec3(20., -20, -200),  //Point C
    glm::vec3(-20., -20, -200), //Point D
    glm::vec3(0.5, 0.5, 0));    //Colour
    sceneObjects.push_back(plane); 
    
    // Cube
    float cx = 10;
    float cy = 5;
    float cz = -65;
    float csize = 2;
    
    Plane *bottom = new Plane (
	 glm::vec3(cx-csize, cy-csize, cz-csize),
	 glm::vec3(cx+csize, cy-csize, cz-csize), 
	 glm::vec3(cx+csize, cy-csize, cz+csize),
	 glm::vec3(cx-csize, cy-csize, cz+csize), 
	 glm::vec3(1, 0.5, 0));
     sceneObjects.push_back(bottom);
    Plane *back = new Plane (
	 glm::vec3(cx-csize, cy-csize, cz-csize),
	 glm::vec3(cx-csize, cy+csize, cz-csize), 
	 glm::vec3(cx+csize, cy+csize, cz-csize),
	 glm::vec3(cx+csize, cy-csize, cz-csize), 
	 glm::vec3(1, 0.5, 0));
     sceneObjects.push_back(back);
    Plane *right = new Plane (
	 glm::vec3(cx+csize, cy-csize, cz-csize),
	 glm::vec3(cx+csize, cy+csize, cz-csize), 
	 glm::vec3(cx+csize, cy+csize, cz+csize),
	 glm::vec3(cx+csize, cy-csize, cz+csize), 
	 glm::vec3(1, 0.5, 0));
     sceneObjects.push_back(right);
    Plane *left = new Plane (
	 glm::vec3(cx-csize, cy-csize, cz-csize),
	 glm::vec3(cx-csize, cy-csize, cz+csize), 
	 glm::vec3(cx-csize, cy+csize, cz+csize),
	 glm::vec3(cx-csize, cy+csize, cz-csize), 
	 glm::vec3(1, 0.5, 0));
     sceneObjects.push_back(left);
    Plane *front = new Plane (
	 glm::vec3(cx-csize, cy-csize, cz+csize),
	 glm::vec3(cx+csize, cy-csize, cz+csize), 
	 glm::vec3(cx+csize, cy+csize, cz+csize),
	 glm::vec3(cx-csize, cy+csize, cz+csize), 
	 glm::vec3(1, 0.5, 0));
     sceneObjects.push_back(front);
    Plane *top = new Plane (
	 glm::vec3(cx-csize, cy+csize, cz-csize),
	 glm::vec3(cx+csize, cy+csize, cz-csize), 
	 glm::vec3(cx+csize, cy+csize, cz+csize),
	 glm::vec3(cx-csize, cy+csize, cz+csize), 
	 glm::vec3(1, 0.5, 0));
     sceneObjects.push_back(top);
     
     
    Sphere *sphere3 = new Sphere(glm::vec3(0.0, -10.0, -80.0), 3.0, glm::vec3(0, 0, 0));
    sceneObjects.push_back(sphere3); 
    
    Sphere *sphere4 = new Sphere(glm::vec3(0.0, 7.0, -70.0), 3.0, glm::vec3(0.0, 0.5, 0.5));
    sceneObjects.push_back(sphere4); 
    
    Sphere *sphere5 = new Sphere(glm::vec3(10.0, -15.0, -90.0), 3.0, glm::vec3(0, 0, 0));
    sceneObjects.push_back(sphere5); 
    
    Plane *plane2 = new Plane (
        glm::vec3(-20., -20, -80),  //Point A
        glm::vec3(-20., -20, -200),   //Point B
        glm::vec3(-20., 20, -200),  //Point C
        glm::vec3(-20., 20, -90), //Point D
        glm::vec3(0, 0, 0));    //Colour
    sceneObjects.push_back(plane2); 
	


	//--Add the above to the list of scene objects.

}



int main(int argc, char *argv[]) {
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB );
    glutInitWindowSize(600, 600);
    glutInitWindowPosition(20, 20);
    glutCreateWindow("Raytracer");

    glutDisplayFunc(display);
    initialize();

    glutMainLoop();
    return 0;
}

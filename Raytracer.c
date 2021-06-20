#include <GL/glut.h>
#include <stdio.h>
#include <math.h>

typedef struct TRIANGLE{
  float triPtOne[3];
  float triPtTwo[3];
  float triPtThree[3];
  float triNormal[3];
} Triangle;

Triangle *allObjects[14];
float centerOfSphere[3] = {0.5f, 0.5f, 1.5f};

float ambient = 0.1f;

// Function Declarations
float* get_normal(Triangle *triangle);
float get_dot_product(float vectorOne[], float vectorTwo[]);
float* get_cross_product(float vectorOne[], float vectorTwo[]);
float* get_intersection(float u, float pointOne[], float pointTwo[]);
float get_u(Triangle *triangle, float pointOne[], float pointTwo[]);
float get_sphere_u(float pointOne[], float pointTwo[]);
float get_magnitude(float vector[]);
float get_dist_from_sphere(float pointOne[], float pointTwo[]);
float* get_difference(float vectorOne[], float vectorTwo[]);
float* get_sum(float vectorOne[], float vectorTwo[]);
float* get_scalar_product(float scalar, float vector[]);
float ray(Triangle *triangle, float pointOne[], float pointTwo[]);
float ray_to_light(float pointOne[], float pointTwo[]);
void drawpixel(float x,float y,float r,float g,float b);
float get_pixel_color(float pointOne[], float pointTwo[]);
float get_diffuse(float pointOne[], float pointTwo[]);
void init_objects_and_light(void);
void display(void);

// Functions
void init_objects_and_light(void){
  // Object Declarations
  // Floor
  Triangle *floorOne;
  floorOne = (Triangle *)malloc(sizeof(Triangle));
  floorOne->triPtOne[0] = 0.0f;
  floorOne->triPtOne[1] = 0.0f;
  floorOne->triPtOne[2] = 0.0f;

  floorOne->triPtTwo[0] = 0.0f;
  floorOne->triPtTwo[1] = 0.0f;
  floorOne->triPtTwo[2] = 1.0f;

  floorOne->triPtThree[0] = 1.0f;
  floorOne->triPtThree[1] = 0.0f;
  floorOne->triPtThree[2] = 1.0f;

  float *floorOneNormal = get_normal(floorOne);
  floorOne->triNormal[0] = floorOneNormal[0];
  floorOne->triNormal[1] = floorOneNormal[1];
  floorOne->triNormal[2] = floorOneNormal[2];

  Triangle *floorTwo;
  floorTwo = (Triangle *)malloc(sizeof(Triangle));
  floorTwo->triPtOne[0] = 0.0f;
  floorTwo->triPtOne[1] = 0.0f;
  floorTwo->triPtOne[2] = 0.0f;

  floorTwo->triPtTwo[0] = 1.0f;
  floorTwo->triPtTwo[1] = 0.0f;
  floorTwo->triPtTwo[2] = 0.0f;

  floorTwo->triPtThree[0] = 1.0f;
  floorTwo->triPtThree[1] = 0.0f;
  floorTwo->triPtThree[2] = 1.0f;

  float *floorTwoNormal = get_normal(floorTwo);
  floorTwo->triNormal[0] = floorTwoNormal[0];
  floorTwo->triNormal[1] = floorTwoNormal[1];
  floorTwo->triNormal[2] = floorTwoNormal[2];

  // Cube
  // Front
  Triangle *frontOne;
  frontOne = (Triangle *)malloc(sizeof(Triangle));
  frontOne->triPtOne[0] = 0.3f;
  frontOne->triPtOne[1] = 0.25f;
  frontOne->triPtOne[2] = 0.5f;

  frontOne->triPtTwo[0] = 0.1f;
  frontOne->triPtTwo[1] = 0.05f;
  frontOne->triPtTwo[2] = 0.5f;

  frontOne->triPtThree[0] = 0.1f;
  frontOne->triPtThree[1] = 0.25f;
  frontOne->triPtThree[2] = 0.5f;

  float *frontOneNormal = get_normal(frontOne);
  frontOne->triNormal[0] = frontOneNormal[0];
  frontOne->triNormal[1] = frontOneNormal[1];
  frontOne->triNormal[2] = frontOneNormal[2];

  Triangle *frontTwo;
  frontTwo = (Triangle *)malloc(sizeof(Triangle));
  frontTwo->triPtOne[0] = 0.3f;
  frontTwo->triPtOne[1] = 0.25f;
  frontTwo->triPtOne[2] = 0.5f;

  frontTwo->triPtTwo[0] = 0.3f;
  frontTwo->triPtTwo[1] = 0.05f;
  frontTwo->triPtTwo[2] = 0.5f;

  frontTwo->triPtThree[0] = 0.1f;
  frontTwo->triPtThree[1] = 0.05f;
  frontTwo->triPtThree[2] = 0.5f;

  float *frontTwoNormal = get_normal(frontTwo);
  frontTwo->triNormal[0] = frontTwoNormal[0];
  frontTwo->triNormal[1] = frontTwoNormal[1];
  frontTwo->triNormal[2] = frontTwoNormal[2];

  // Back
  Triangle *backOne;
  backOne = (Triangle *)malloc(sizeof(Triangle));
  backOne->triPtOne[0] = 0.3f;
  backOne->triPtOne[1] = 0.25f;
  backOne->triPtOne[2] = 0.7f;

  backOne->triPtTwo[0] = 0.1f;
  backOne->triPtTwo[1] = 0.25f;
  backOne->triPtTwo[2] = 0.7f;

  backOne->triPtThree[0] = 0.3f;
  backOne->triPtThree[1] = 0.05f;
  backOne->triPtThree[2] = 0.7f;

  float *backOneNormal = get_normal(backOne);
  backOne->triNormal[0] = backOneNormal[0];
  backOne->triNormal[1] = backOneNormal[1];
  backOne->triNormal[2] = backOneNormal[2];

  Triangle *backTwo;
  backTwo = (Triangle *)malloc(sizeof(Triangle));
  backTwo->triPtOne[0] = 0.1f;
  backTwo->triPtOne[1] = 0.25f;
  backTwo->triPtOne[2] = 0.7f;

  backTwo->triPtTwo[0] = 0.1f;
  backTwo->triPtTwo[1] = 0.05f;
  backTwo->triPtTwo[2] = 0.7f;

  backTwo->triPtThree[0] = 0.3f;
  backTwo->triPtThree[1] = 0.05f;
  backTwo->triPtThree[2] = 0.7f;

  float *backTwoNormal = get_normal(backTwo);
  backTwo->triNormal[0] = backOneNormal[0];
  backTwo->triNormal[1] = backTwoNormal[1];
  backTwo->triNormal[2] = backTwoNormal[2];

  // Left
  Triangle *leftOne;
  leftOne = (Triangle *)malloc(sizeof(Triangle));
  leftOne->triPtOne[0] = 0.1f;
  leftOne->triPtOne[1] = 0.05f;
  leftOne->triPtOne[2] = 0.5f;

  leftOne->triPtTwo[0] = 0.1f;
  leftOne->triPtTwo[1] = 0.05f;
  leftOne->triPtTwo[2] = 0.7f;

  leftOne->triPtThree[0] = 0.1f;
  leftOne->triPtThree[1] = 0.25f;
  leftOne->triPtThree[2] = 0.7f;

  float *leftOneNormal = get_normal(leftOne);
  leftOne->triNormal[0] = leftOneNormal[0];
  leftOne->triNormal[1] = leftOneNormal[1];
  leftOne->triNormal[2] = leftOneNormal[2];

  Triangle *leftTwo;
  leftTwo = (Triangle *)malloc(sizeof(Triangle));
  leftTwo->triPtOne[0] = 0.1f;
  leftTwo->triPtOne[1] = 0.05f;
  leftTwo->triPtOne[2] = 0.5f;

  leftTwo->triPtTwo[0] = 0.1f;
  leftTwo->triPtTwo[1] = 0.25f;
  leftTwo->triPtTwo[2] = 0.7f;

  leftTwo->triPtThree[0] = 0.1f;
  leftTwo->triPtThree[1] = 0.25f;
  leftTwo->triPtThree[2] = 0.5f;

  float *leftTwoNormal = get_normal(leftTwo);
  leftTwo->triNormal[0] = leftTwoNormal[0];
  leftTwo->triNormal[1] = leftTwoNormal[1];
  leftTwo->triNormal[2] = leftTwoNormal[2];

  // Right
  Triangle *rightOne;
  rightOne = (Triangle *)malloc(sizeof(Triangle));
  rightOne->triPtOne[0] = 0.3f;
  rightOne->triPtOne[1] = 0.25f;
  rightOne->triPtOne[2] = 0.7f;

  rightOne->triPtTwo[0] = 0.3f;
  rightOne->triPtTwo[1] = 0.05f;
  rightOne->triPtTwo[2] = 0.5f;

  rightOne->triPtThree[0] = 0.3f;
  rightOne->triPtThree[1] = 0.25f;
  rightOne->triPtThree[2] = 0.5f;

  float *rightOneNormal = get_normal(rightOne);
  rightOne->triNormal[0] = rightOneNormal[0];
  rightOne->triNormal[1] = rightOneNormal[1];
  rightOne->triNormal[2] = rightOneNormal[2];

  Triangle *rightTwo;
  rightTwo = (Triangle *)malloc(sizeof(Triangle));
  rightTwo->triPtOne[0] = 0.3f;
  rightTwo->triPtOne[1] = 0.05f;
  rightTwo->triPtOne[2] = 0.5f;

  rightTwo->triPtTwo[0] = 0.3f;
  rightTwo->triPtTwo[1] = 0.25f;
  rightTwo->triPtTwo[2] = 0.7f;

  rightTwo->triPtThree[0] = 0.3f;
  rightTwo->triPtThree[1] = 0.05f;
  rightTwo->triPtThree[2] = 0.7f;

  float *rightTwoNormal = get_normal(rightTwo);
  rightTwo->triNormal[0] = rightTwoNormal[0];
  rightTwo->triNormal[1] = rightTwoNormal[1];
  rightTwo->triNormal[2] = rightTwoNormal[2];

  // Top
  Triangle *topOne;
  topOne = (Triangle *)malloc(sizeof(Triangle));
  topOne->triPtOne[0] = 0.3f;
  topOne->triPtOne[1] = 0.25f;
  topOne->triPtOne[2] = 0.7f;

  topOne->triPtTwo[0] = 0.3f;
  topOne->triPtTwo[1] = 0.25f;
  topOne->triPtTwo[2] = 0.5f;

  topOne->triPtThree[0] = 0.1f;
  topOne->triPtThree[1] = 0.25f;
  topOne->triPtThree[2] = 0.5f;

  float *topOneNormal = get_normal(topOne);
  topOne->triNormal[0] = topOneNormal[0];
  topOne->triNormal[1] = topOneNormal[1];
  topOne->triNormal[2] = topOneNormal[2];

  Triangle *topTwo;
  topTwo = (Triangle *)malloc(sizeof(Triangle));
  topTwo->triPtOne[0] = 0.3f;
  topTwo->triPtOne[1] = 0.25f;
  topTwo->triPtOne[2] = 0.7f;

  topTwo->triPtTwo[0] = 0.1f;
  topTwo->triPtTwo[1] = 0.25f;
  topTwo->triPtTwo[2] = 0.5f;

  topTwo->triPtThree[0] = 0.1f;
  topTwo->triPtThree[1] = 0.25f;
  topTwo->triPtThree[2] = 0.7f;

  float *topTwoNormal = get_normal(topTwo);
  topTwo->triNormal[0] = topTwoNormal[0];
  topTwo->triNormal[1] = topTwoNormal[1];
  topTwo->triNormal[2] = topTwoNormal[2];

  // Bottom
  Triangle *bottomOne;
  bottomOne = (Triangle *)malloc(sizeof(Triangle));
  bottomOne->triPtOne[0] = 0.3f;
  bottomOne->triPtOne[1] = 0.05f;
  bottomOne->triPtOne[2] = 0.7f;

  bottomOne->triPtTwo[0] = 0.1f;
  bottomOne->triPtTwo[1] = 0.05f;
  bottomOne->triPtTwo[2] = 0.5f;

  bottomOne->triPtThree[0] = 0.3f;
  bottomOne->triPtThree[1] = 0.05f;
  bottomOne->triPtThree[2] = 0.5f;

  float *bottomOneNormal = get_normal(bottomOne);
  bottomOne->triNormal[0] = bottomOneNormal[0];
  bottomOne->triNormal[1] = bottomOneNormal[1];
  bottomOne->triNormal[2] = bottomOneNormal[2];

  Triangle *bottomTwo;
  bottomTwo = (Triangle *)malloc(sizeof(Triangle));
  bottomTwo->triPtOne[0] = 0.3f;
  bottomTwo->triPtOne[1] = 0.05f;
  bottomTwo->triPtOne[2] = 0.7f;

  bottomTwo->triPtTwo[0] = 0.1f;
  bottomTwo->triPtTwo[1] = 0.05f;
  bottomTwo->triPtTwo[2] = 0.7f;

  bottomTwo->triPtThree[0] = 0.1f;
  bottomTwo->triPtThree[1] = 0.05f;
  bottomTwo->triPtThree[2] = 0.5f;

  float *bottomTwoNormal = get_normal(bottomTwo);
  bottomTwo->triNormal[0] = bottomTwoNormal[0];
  bottomTwo->triNormal[1] = bottomTwoNormal[1];
  bottomTwo->triNormal[2] = bottomTwoNormal[2];

  allObjects[0] = floorOne;
  allObjects[1] = floorTwo;
  allObjects[2] = leftOne;
  allObjects[3] = leftTwo;
  allObjects[4] = rightOne;
  allObjects[5] = rightTwo;
  allObjects[6] = backOne;
  allObjects[7] = backTwo;
  allObjects[8] = frontOne;
  allObjects[9] = frontTwo;
  allObjects[10] = topOne;
  allObjects[11] = topTwo;
  allObjects[12] = bottomOne;
  allObjects[13] = bottomTwo;
}

float* get_normal(Triangle *triangle){
   float *triOne = triangle->triPtOne;
   float *triTwo = triangle->triPtTwo;
   float *triThree = triangle->triPtThree;

   float *vectorOne = get_difference(triTwo, triOne);
   float *vectorTwo = get_difference(triThree, triOne);

   float *normal = get_cross_product(vectorOne, vectorTwo);
   return normal;
}

float get_dot_product(float vectorOne[], float vectorTwo[]){
   float product = 0;
   for(int i = 0; i < 3; i++){
      product = product + (vectorOne[i] * vectorTwo[i]);
   }
   return product;
}

float* get_cross_product(float vectorOne[], float vectorTwo[]){
   float *product = calloc(3, sizeof(float));
   product[0] = (vectorOne[1]*vectorTwo[2]) - (vectorOne[2]*vectorTwo[1]);
   product[1] = (vectorOne[2]*vectorTwo[0]) - (vectorOne[0]*vectorTwo[2]);
   product[2] = (vectorOne[0]*vectorTwo[1]) - (vectorOne[1]*vectorTwo[0]);
   return product;
}

float* get_intersection(float u, float pointOne[], float pointTwo[]){
   float *difference = get_difference(pointTwo, pointOne);
   float *applyScalar = get_scalar_product(u, difference);
   float *intsersection = get_sum(pointOne, applyScalar);
   return intsersection;
}

float get_u(Triangle *triangle, float pointOne[], float pointTwo[]){
   float *normal = triangle->triNormal;
   float *trianglePoint = triangle->triPtThree;
   float numerator = get_dot_product(normal, get_difference(trianglePoint, pointOne));
   float denominator = get_dot_product(normal, get_difference(pointTwo, pointOne));
   float u = numerator / denominator;
   return u;
}

float* get_difference(float vectorOne[], float vectorTwo[]){
   float *difference = calloc(3, sizeof(float));
   difference[0] = vectorOne[0] - vectorTwo[0];
   difference[1] = vectorOne[1] - vectorTwo[1];
   difference[2] = vectorOne[2] - vectorTwo[2];
   return difference;
}

float* get_sum(float vectorOne[], float vectorTwo[]){
   float *sum = calloc(3, sizeof(float));
   sum[0] = vectorOne[0] + vectorTwo[0];
   sum[1] = vectorOne[1] + vectorTwo[1];
   sum[2] = vectorOne[2] + vectorTwo[2];
   return sum;
}

float* get_scalar_product(float scalar, float vector[]){
   float *scaled = calloc(3, sizeof(float));
   scaled[0] = vector[0] * scalar;
   scaled[1] = vector[1] * scalar;
   scaled[2] = vector[2] * scalar;
   return scaled;
}

float get_magnitude(float vector[]){
   float numUnderRoot = ((vector[0]*vector[0]) + (vector[1]*vector[1]) + (vector[2]*vector[2]));
   float magnitude = sqrt(numUnderRoot);
   return magnitude;
}

float get_dist_from_sphere(float pointOne[], float pointTwo[]){
   float xMultOne = centerOfSphere[0] - pointOne[0];
   float xMultTwo = pointTwo[0] - pointOne[0];
   float xProduct = xMultOne * xMultTwo;

   float yMultOne = centerOfSphere[1] - pointOne[1];
   float yMultTwo = pointTwo[1] - pointOne[1];
   float yProduct = yMultOne * yMultTwo;

   float zMultOne = centerOfSphere[2] - pointOne[2];
   float zMultTwo = pointTwo[2] - pointOne[2];
   float zProduct = zMultOne * zMultTwo;

   float numerator = xProduct + yProduct + zProduct;
   float* pTwoMinuspOne = get_difference(pointTwo, pointOne);
   float denominator = get_magnitude(pTwoMinuspOne) * get_magnitude(pTwoMinuspOne);
   float u = numerator/denominator;

   // Find the P value
   float* difference = get_difference(pointTwo, pointOne);
   float* applyU = get_scalar_product(u, pTwoMinuspOne);
   float* p = get_sum(pointOne, applyU);

   // Get magnitude of P-P3
   float* pMinusPthree = get_difference(p, centerOfSphere);
   float mag = get_magnitude(pMinusPthree);

   return mag;
}

float get_pixel_color(float pointOne[], float pointTwo[]){
   float diffuse = get_diffuse(pointOne, pointTwo);
   return diffuse + ambient;
}

float get_diffuse(float pointOne[], float pointTwo[]){
    float bright = 0.0f;
    float uObject = -1.0f;
    float uLight = -1.0f;
    float objectIndex = 0.0f;

    for(int i = 0; i < 14; i++){
        Triangle *current = allObjects[i];

        if(ray(allObjects[i], pointOne, pointTwo) == 1.0f){
            float u = get_u(allObjects[i], pointOne, pointTwo);
            if((uObject > u || uObject == -1.0f) && u > 0.0f){
                uObject = u;
                objectIndex = i;
            }
        }
    }

    if(ray_to_light(pointOne, pointTwo) == 1.0f){
      float xMultOne = centerOfSphere[0] - pointOne[0];
      float xMultTwo = pointTwo[0] - pointOne[0];
      float xProduct = xMultOne * xMultTwo;

      float yMultOne = centerOfSphere[1] - pointOne[1];
      float yMultTwo = pointTwo[1] - pointOne[1];
      float yProduct = yMultOne * yMultTwo;

      float zMultOne = centerOfSphere[2] - pointOne[2];
      float zMultTwo = pointTwo[2] - pointOne[2];
      float zProduct = zMultOne * zMultTwo;

      float numerator = xProduct + yProduct + zProduct;

      float* pTwoMinuspOne = get_difference(pointTwo, pointOne);
      float denominator = get_magnitude(pTwoMinuspOne) * get_magnitude(pTwoMinuspOne);
      float u = numerator/denominator;
         if((uLight > u || uLight == -1.0f) && u > 0){
             uLight = u;
         }
    }

    if(uLight > 0.0f && (uLight < uObject || uObject < 0.0f)){
       bright = 1.0f;
    }
    if(uObject > 0.0f && (uObject < uLight || uLight < 0)){
         float* intersection = get_intersection(uObject, pointOne, pointTwo);

         if(!(pointTwo[0] == centerOfSphere[0] && pointTwo[1] == centerOfSphere[1] && pointTwo[2] == centerOfSphere[2])){
             float numerator = get_diffuse(intersection, centerOfSphere);
             float* difference = get_difference(centerOfSphere, intersection);
             float denominator = get_magnitude(difference);
             bright = bright + (numerator / denominator);
          }
    }
    return bright;
}

float ray(Triangle *triangle, float pointOne[], float pointTwo[]){
   float *triOne = triangle->triPtOne;
   float *triTwo = triangle->triPtTwo;
   float *triThree = triangle->triPtThree;

   float u = get_u(triangle, pointOne, pointTwo);
   float *intsersection = get_intersection(u, pointOne, pointTwo);

   float *vOne = get_difference(triOne, intsersection);
   float *vTwo = get_difference(triTwo, intsersection);
   float *vThree = get_difference(triThree, intsersection);

   float *crossOne = get_cross_product(vOne, vTwo);
   float *crossTwo = get_cross_product(vTwo, vThree);
   float *crossThree = get_cross_product(vThree, vOne);

   float dotOne = get_dot_product(crossOne, crossTwo);
   float dotTwo = get_dot_product(crossTwo, crossThree);
   float dotThree = get_dot_product(crossThree, crossOne);

   float hit = 0.0f;

   if(dotOne > 0.0f && dotTwo > 0.0f && dotThree > 0.0f){
      hit = 1.0f;
   }
   return hit;
}

float ray_to_light(float pointOne[], float pointTwo[]){
    float hit = 0.0f;
    if(get_dist_from_sphere(pointOne, pointTwo) < 0.1f){
       hit = 1.0f;
    }
    return hit;
}

void drawpixel(float x,float y,float r,float g,float b){
#define SZ  0.004
	glBegin(GL_TRIANGLES);
		glColor3f(r,g,b);
	 	glVertex2f(-1.0+.004*(float)x,     -1.0+.004*(float)y);
	 	glVertex2f(-1.0+.004*(float)x,     -1.0+.004*(float)y+ SZ);
	 	glVertex2f(-1.0+.004*(float)x+ SZ, -1.0+.004*(float)y);

	 	glVertex2f(-1.0+.004*(float)x+ SZ, -1.0+.004*(float)y);
	 	glVertex2f(-1.0+.004*(float)x+ SZ, -1.0+.004*(float)y+ SZ);
	 	glVertex2f(-1.0+.004*(float)x,     -1.0+.004*(float)y+ SZ);
	glEnd();
}

void display(void){
   init_objects_and_light();
   glClear(GL_COLOR_BUFFER_BIT);

   float eye[3] = { 0.5f, 0.5f, -1.0f };
   float screen[3] = { 0.0f, 0.0f, 0.0f };

   for(float x = 0.0f; x < 1.0f; x = x + 0.002f){
      screen[0] = x;
      for(float y = 0.0f; y < 1.0f; y = y + 0.002f){
         screen[1] = y;

         float pixelColor = get_pixel_color(eye, screen);
         if(pixelColor > 1.0f ){
            pixelColor = 1.0f;
         }

         float scaledX = x*500;
         float scaledY = y*500;
         drawpixel(scaledX, scaledY, pixelColor, pixelColor, pixelColor);
      }
   }
	glFlush();
}

int main(int argc, char** argv){
  glutInit(&argc,argv);
	glutCreateWindow("simple");
	glutDisplayFunc(display);
	glutMainLoop();
}

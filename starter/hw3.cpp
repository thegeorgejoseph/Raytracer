/* **************************
 * CSCI 420
 * Assignment 3 Raytracer
 * Name: georgejo
 * *************************
*/
#ifdef WIN32
  #include <windows.h>
#endif

#if defined(WIN32) || defined(linux)
  #include <GL/gl.h>
  #include <GL/glut.h>
#elif defined(__APPLE__)
  #include <OpenGL/gl.h>
  #include <GLUT/glut.h>
#endif

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <limits>
#include <imageIO.h>
#include <float.h>
#include <algorithm>
#include <float.h>
#include <iostream>
#include <random>
#include "Vector.h"

#ifdef WIN32
  #define strcasecmp _stricmp
#endif

#define MAX_TRIANGLES 20000
#define MAX_SPHERES 100
#define MAX_LIGHTS 100

#define WIDTH 640
#define HEIGHT 480

char * filename = NULL;


#define MODE_DISPLAY 1
#define MODE_JPEG 2

int mode = MODE_DISPLAY;


#define fov 60.0
const bool ANTI_ALIASING =  true;
const bool MODE_SOFT = true;

using namespace std;

unsigned char buffer[HEIGHT][WIDTH][3];
double y_min = -tan((fov * 3.14159 / 180)/2);
double x_min = -((double)WIDTH/HEIGHT)*tan((fov * 3.14159 / 180)/2);
double w_screen = 2*((double)WIDTH/HEIGHT)*tan((fov * 3.14159 / 180)/2);
double h_screen = 2*tan((fov * 3.14159 / 180)/2);

struct Vertex
{
  Vector position;
  Vector color_diffuse;
  Vector color_specular;
  Vector normal;
  double shininess;
};

struct Triangle
{
  Vertex v[3];
};

struct Sphere
{
  Vector position;
  Vector color_diffuse;
  Vector color_specular;
  double shininess;
  double radius;
};

struct Light
{
  Vector position;
  Vector color;
};

Vector normalise(Vector vec){
    return vec * (1/sqrt(vec.x * vec.x + vec.y * vec.y + vec.z * vec.z));
}

struct Ray
{
    Vector origin, direction;
};

Triangle triangles[MAX_TRIANGLES];
Sphere spheres[MAX_SPHERES];
Light lights[MAX_LIGHTS];
Vector ambient_light;

double thresh,thresh2; int oneIdx,twoIndex;

int sphereNumber = 0;
int triangleNumber = 0;

random_device device;
mt19937 generator(device());
// mt9937 generator(0);
uniform_real_distribution<double> distribution(0.0f, 1);

void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel(int x,int y,unsigned char r,unsigned char g,unsigned char b);
bool intersectingShadow(const Light &light, const Vertex &point);
void checkSphereAndTriangleIntersection(const Ray &ray, double &sphereTemp, int &sphereIdx, double &triangleTemp, int &triangleIdx);
void antiAliasingOff(unsigned int x, unsigned int y);
void antiAliasingOn(unsigned int x, unsigned int y);

Vector diffuse, specular, normal, positionIntersected;
double shininess;
Vector cameraVec = Vector(0,0,0);
bool doesIntersectSphere, doesIntersectTriangle;

int lightNumber = 0;

Vector getLightingModel(const Ray &ray, int times)
{
  

  //check if spheres and triangles intersect
  checkSphereAndTriangleIntersection(ray,thresh,oneIdx,thresh2,twoIndex);
  
  //check if there is no intersection of spheres or triangles then return bg
  if (doesIntersectSphere == false && doesIntersectTriangle == false)
    return Vector(0.823,0.74,0.74);
  
  // if sphere intersection
  Vertex intersectionPosition;
  if ((doesIntersectSphere && !doesIntersectTriangle) || (doesIntersectSphere && doesIntersectTriangle && thresh<thresh2))
  {
    Vector position = ray.origin + ray.direction * thresh;
    diffuse = spheres[oneIdx].color_diffuse;
    specular = spheres[oneIdx].color_specular;
    normal = ray.origin + ray.direction * thresh - spheres[oneIdx].position;
    shininess = spheres[oneIdx].shininess;
    intersectionPosition = {
            position,
            diffuse,
            specular,
            normalise(normal),
            shininess
    };
  }
  // if triangle intersection
  else
  {
    double ratio[3];
    Vector position = ray.origin + ray.direction * thresh2;
    int &index = twoIndex;
    double (&r)[3] = ratio;
    Vector xy = triangles[index].v[1].position - triangles[index].v[0].position;
    Vector xz = triangles[index].v[2].position - triangles[index].v[0].position;
    Vector px = triangles[index].v[0].position - position;
    Vector py = triangles[index].v[1].position - position;
    Vector pz = triangles[index].v[2].position - position;
    
    double xyC_Area = sqrt(Vector(xy.y * xz.z - xy.z * xz.y, xy.z * xz.x - xy.x * xz.z, xy.x * xz.y - xy.y * xz.x).x * Vector(xy.y * xz.z - xy.z * xz.y, xy.z * xz.x - xy.x * xz.z, xy.x * xz.y - xy.y * xz.x).x + Vector(xy.y * xz.z - xy.z * xz.y, xy.z * xz.x - xy.x * xz.z, xy.x * xz.y - xy.y * xz.x).y * Vector(xy.y * xz.z - xy.z * xz.y, xy.z * xz.x - xy.x * xz.z, xy.x * xz.y - xy.y * xz.x).y + Vector(xy.y * xz.z - xy.z * xz.y, xy.z * xz.x - xy.x * xz.z, xy.x * xz.y - xy.y * xz.x).z * Vector(xy.y * xz.z - xy.z * xz.y, xy.z * xz.x - xy.x * xz.z, xy.x * xz.y - xy.y * xz.x).z)  * 0.5f;
    double PyC_Area = sqrt(Vector(py.y * pz.z - py.z * pz.y, py.z * pz.x - py.x * pz.z, py.x * pz.y - py.y * pz.x).x * Vector(py.y * pz.z - py.z * pz.y, py.z * pz.x - py.x * pz.z, py.x * pz.y - py.y * pz.x).x + Vector(py.y * pz.z - py.z * pz.y, py.z * pz.x - py.x * pz.z, py.x * pz.y - py.y * pz.x).y * Vector(py.y * pz.z - py.z * pz.y, py.z * pz.x - py.x * pz.z, py.x * pz.y - py.y * pz.x).y + Vector(py.y * pz.z - py.z * pz.y, py.z * pz.x - py.x * pz.z, py.x * pz.y - py.y * pz.x).z * Vector(py.y * pz.z - py.z * pz.y, py.z * pz.x - py.x * pz.z, py.x * pz.y - py.y * pz.x).z) * 0.5f;
    double PzA_Area = sqrt(Vector(pz.y * px.z - pz.z * px.y, pz.z * px.x - pz.x * px.z, pz.x * px.y - pz.y * px.x).x * Vector(pz.y * px.z - pz.z * px.y, pz.z * px.x - pz.x * px.z, pz.x * px.y - pz.y * px.x).x + Vector(pz.y * px.z - pz.z * px.y, pz.z * px.x - pz.x * px.z, pz.x * px.y - pz.y * px.x).y * Vector(pz.y * px.z - pz.z * px.y, pz.z * px.x - pz.x * px.z, pz.x * px.y - pz.y * px.x).y + Vector(pz.y * px.z - pz.z * px.y, pz.z * px.x - pz.x * px.z, pz.x * px.y - pz.y * px.x).z * Vector(pz.y * px.z - pz.z * px.y, pz.z * px.x - pz.x * px.z, pz.x * px.y - pz.y * px.x).z) * 0.5f;
    double Pxy_Area = sqrt(Vector(px.y * py.z - px.z * py.y, px.z * py.x - px.x * py.z, px.x * py.y - px.y * py.x).x * Vector(px.y * py.z - px.z * py.y, px.z * py.x - px.x * py.z, px.x * py.y - px.y * py.x).x + Vector(px.y * py.z - px.z * py.y, px.z * py.x - px.x * py.z, px.x * py.y - px.y * py.x).y * Vector(px.y * py.z - px.z * py.y, px.z * py.x - px.x * py.z, px.x * py.y - px.y * py.x).y + Vector(px.y * py.z - px.z * py.y, px.z * py.x - px.x * py.z, px.x * py.y - px.y * py.x).z * Vector(px.y * py.z - px.z * py.y, px.z * py.x - px.x * py.z, px.x * py.y - px.y * py.x).z) * 0.5f;
    
    r[0] = PyC_Area / xyC_Area;
    r[1] = PzA_Area / xyC_Area;
    r[2] = Pxy_Area / xyC_Area;

    intersectionPosition = {
            position,
            triangles[twoIndex].v[0].color_diffuse * ratio[0] +
            triangles[twoIndex].v[1].color_diffuse * ratio[1] +
            triangles[twoIndex].v[2].color_diffuse * ratio[2],
           triangles[twoIndex].v[0].color_specular * ratio[0] +
                   triangles[twoIndex].v[1].color_specular * ratio[1] +
                   triangles[twoIndex].v[2].color_specular * ratio[2],
            normalise(triangles[twoIndex].v[0].normal * ratio[0] +
                   triangles[twoIndex].v[1].normal * ratio[1] +
                   triangles[twoIndex].v[2].normal * ratio[2]),
            triangles[twoIndex].v[0].shininess * ratio[0] +
                   triangles[twoIndex].v[1].shininess * ratio[1] +
                   triangles[twoIndex].v[2].shininess * ratio[2]
    };
  }
  Vector currentColor = Vector(0,0,0);
  for (int k = 0; k<lightNumber; k++)
  {
    if (!intersectingShadow(lights[k], intersectionPosition))
    {
      Light &light = lights[k];
      Vector specDiff;
      Vector direction = normalise(light.position-intersectionPosition.position);
      double diff = max((direction.x * intersectionPosition.normal.x + direction.y * intersectionPosition.normal.y + direction.z * intersectionPosition.normal.z),0.0);
      Vector temp = Vector(intersectionPosition.color_diffuse.x * diff, intersectionPosition.color_diffuse.y * diff, intersectionPosition.color_diffuse.z * diff);
      diffuse = Vector(light.color.x * temp.x, light.color.x * temp.y, light.color.z * temp.z);
      Vector view = normalise(-intersectionPosition.position);
      Vector reflect = normalise(-direction - intersectionPosition.normal * 2 * (intersectionPosition.normal.x * -direction.x + intersectionPosition.normal.y * -direction.y + intersectionPosition.normal.z * -direction.z));
      double spec = pow(max((view.x * reflect.x + view.y * reflect.y + view.z * reflect.z),0.0), intersectionPosition.shininess);
      Vector specTemp = Vector(intersectionPosition.color_specular.x * spec, intersectionPosition.color_specular.y * spec, intersectionPosition.color_specular.z * spec);
      specular =  Vector(light.color.x * specTemp.x, light.color.y * specTemp.y, light.color.z * specTemp.z);
  
      specDiff = diffuse + specular;
      
      currentColor =  currentColor + specDiff;
    }
  }
  
  // break logic for recursive reflection
  if (times >= 5)
    return currentColor;
  
  times++;
  // reflect logic
  Vector reflect = normalise(ray.direction - intersectionPosition.normal * 2 * (intersectionPosition.normal.x * ray.direction.x + intersectionPosition.normal.y * ray.direction.y + intersectionPosition.normal.z * ray.direction.z));
  Ray reflectRay = {intersectionPosition.position, reflect};
  Vector reflectColor = getLightingModel(reflectRay, times);
  
  return currentColor * 1/2 +  reflectColor * 0.1;
}

// check if rays are intersecting the shadow region
bool intersectingShadow(const Light &light, const Vertex &point)
{
  Vector direction =normalise(light.position-point.position);
  Ray shadow = {point.position + direction * 5 * 1e-5, direction};
  
  double tempSphere, tempTriangle; int sphereIdx, triangleIdx;
  //check intersections
  checkSphereAndTriangleIntersection(shadow, tempSphere, sphereIdx, tempTriangle, triangleIdx);
  
  // no sphere or triangle intersection
  if (doesIntersectSphere == false && doesIntersectTriangle == false)
    return false;
  
  if ((doesIntersectSphere && doesIntersectTriangle == false) || (doesIntersectSphere && doesIntersectTriangle && tempSphere<tempTriangle))
    positionIntersected =  shadow.origin + shadow.direction * tempSphere;
  else
    positionIntersected =  shadow.origin + shadow.direction * tempTriangle;

  double intersectionPt = sqrt((point.position - positionIntersected).x * (point.position - positionIntersected).x + (point.position - positionIntersected).y * (point.position - positionIntersected).y + (point.position - positionIntersected).z * (point.position - positionIntersected).z);
  double lightPoint;
  lightPoint = sqrt((point.position - light.position).x * (point.position - light.position).x + (point.position - light.position).y * (point.position - light.position).y + (point.position - light.position).z * (point.position - light.position).z);
  
  if (intersectionPt - lightPoint > 1e-5)
    return false;
  
  return true;
}

//save state 

//MODIFY THIS FUNCTION
void draw_scene()
{
  for(unsigned int x=0; x<WIDTH; x++)
  {
    glPointSize(2.0);  
    glBegin(GL_POINTS);
    for(unsigned int y=0; y<HEIGHT; y++)
    {
      if (ANTI_ALIASING == false)
      {
        antiAliasingOff(x,y);
      }
      else if(ANTI_ALIASING == true)
      {
        antiAliasingOn(x,y);
      }
      
    }
    glEnd();
    glFlush();
  }
  printf("Done!\n"); fflush(stdout);
}

void plot_pixel_display(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  glColor3f(((float)r) / 255.0f, ((float)g) / 255.0f, ((float)b) / 255.0f);
  glVertex2i(x,y);
}

void plot_pixel_jpeg(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  buffer[y][x][0] = r;
  buffer[y][x][1] = g;
  buffer[y][x][2] = b;
}

void plot_pixel(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  plot_pixel_display(x,y,r,g,b);
  if(mode == MODE_JPEG)
    plot_pixel_jpeg(x,y,r,g,b);
}

void save_jpg()
{
  printf("Saving JPEG file: %s\n", filename);

  ImageIO img(WIDTH, HEIGHT, 3, &buffer[0][0][0]);
  if (img.save(filename, ImageIO::FORMAT_JPEG) != ImageIO::OK)
    printf("Error in Saving\n");
  else 
    printf("File saved Successfully\n");
}

void parse_check(const char *expected, char *found)
{
  if(strcasecmp(expected,found))
  {
    printf("Expected '%s ' found '%s '\n", expected, found);
    printf("Parse error, abnormal abortion\n");
    exit(0);
  }
}

void parse_doubles(FILE* file, const char *check, Vector &p)
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check(check,str);
  fscanf(file,"%lf %lf %lf",&p.x,&p.y,&p.z);
  printf("%s %lf %lf %lf\n",check,p.x,p.y,p.z);
}

void parse_rad(FILE *file, double *r)
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check("rad:",str);
  fscanf(file,"%lf",r);
  printf("rad: %f\n",*r);
}

void parse_shi(FILE *file, double *shi)
{
  char s[100];
  fscanf(file,"%s",s);
  parse_check("shi:",s);
  fscanf(file,"%lf",shi);
  printf("shi: %f\n",*shi);
}

int loadScene(char *argv)
{
  FILE * file = fopen(argv,"r");
  int number_of_objects;
  char type[50];
  Triangle t;
  Sphere s;
  Light l;
  fscanf(file,"%i", &number_of_objects);

  printf("number of objects: %i\n",number_of_objects);

  parse_doubles(file,"amb:",ambient_light);

  for(int i=0; i<number_of_objects; i++)
  {
    fscanf(file,"%s\n",type);
    printf("%s\n",type);
    if(strcasecmp(type,"triangle")==0)
    {
      printf("found triangle\n");
      for(int j=0;j < 3;j++)
      {
        parse_doubles(file,"pos:",t.v[j].position);
        parse_doubles(file,"nor:",t.v[j].normal);
        parse_doubles(file,"dif:",t.v[j].color_diffuse);
        parse_doubles(file,"spe:",t.v[j].color_specular);
        parse_shi(file,&t.v[j].shininess);
      }

      if(triangleNumber == MAX_TRIANGLES)
      {
        printf("too many triangles, you should increase MAX_TRIANGLES!\n");
        exit(0);
      }
      triangles[triangleNumber++] = t;
    }
    else if(strcasecmp(type,"sphere")==0)
    {
      printf("found sphere\n");

      parse_doubles(file,"pos:",s.position);
      parse_rad(file,&s.radius);
      parse_doubles(file,"dif:",s.color_diffuse);
      parse_doubles(file,"spe:",s.color_specular);
      parse_shi(file,&s.shininess);

      if(sphereNumber == MAX_SPHERES)
      {
        printf("too many spheres, you should increase MAX_SPHERES!\n");
        exit(0);
      }
      spheres[sphereNumber++] = s;
    }
    else if(strcasecmp(type,"light")==0)
    {
      printf("found light\n");
      parse_doubles(file,"pos:",l.position);
      parse_doubles(file,"col:",l.color);

      if(lightNumber == MAX_LIGHTS)
      {
        printf("too many lights, you should increase MAX_LIGHTS!\n");
        exit(0);
      }
      lights[lightNumber++] = l;
    }
    else
    {
      printf("unknown type in scene description:\n%s\n",type);
      exit(0);
    }
  }
  return 0;
}
void display()
{
}

void init()
{
  glMatrixMode(GL_PROJECTION);
  glOrtho(0,WIDTH,0,HEIGHT,1,-1);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glClearColor(0,0,0,0);
  glClear(GL_COLOR_BUFFER_BIT);
}

void idle()
{
  //hack to make it only draw once
  static int once=0;
  if(!once)
  {
    if (MODE_SOFT == true){
    int origin_light_num = lightNumber;
    
    for (int i = 0; i<origin_light_num; i++)
    {
        Vector color = lights[i].color/48;
        Vector center = lights[i].position;
        
        lights[i].color = color;
        for (int j = 0; j<(48-1); j++)
        {
        lights[lightNumber].color = color;
        lights[lightNumber].position = Vector (center.x+distribution(generator), center.y+distribution(generator), center.z+distribution(generator));
        lightNumber++;
        }
    }
    }
    draw_scene();
    if(mode == MODE_JPEG)
      save_jpg();
  }
  once=1;
}

int main(int argc, char ** argv)
{
  if ((argc < 2) || (argc > 3))
  {  
    printf ("Usage: %s <input scenefile> [output jpegname]\n", argv[0]);
    exit(0);
  }
  if(argc == 3)
  {
    mode = MODE_JPEG;
    filename = argv[2];
  }
  else if(argc == 2)
    mode = MODE_DISPLAY;

  glutInit(&argc,argv);
  loadScene(argv[1]);

  glutInitDisplayMode(GLUT_RGBA | GLUT_SINGLE);
  glutInitWindowPosition(0,0);
  glutInitWindowSize(WIDTH,HEIGHT);
  int window = glutCreateWindow("Ray Tracer");
  #ifdef __APPLE__
    // This is needed on recent Mac OS X versions to correctly display the window.
    glutReshapeWindow(WIDTH - 1, HEIGHT - 1);
  #endif
  glutDisplayFunc(display);
  glutIdleFunc(idle);
  init();
  glutMainLoop();
}

void antiAliasingOff(unsigned int x, unsigned int y){
    Vector color;
    double positions[3] = {x_min + (2*x+1)*(w_screen/WIDTH)/2.0f,y_min + (2*y+1)*(w_screen/WIDTH)/2.0f, -1};
    color = getLightingModel({Vector(0,0,0),  normalise(Vector(positions[0],positions[1],positions[2]))}, 1) + ambient_light;
    if (color.x > 1.0)
        color.x = 1.0;
    if (color.y > 1.0)
        color.y = 1.0;
    if (color.z > 1.0)
        color.z = 1.0;
    plot_pixel(x, y, (int)(color.x * 255), (int)(color.y * 255), (int)(color.z * 255));
}

void antiAliasingOn(unsigned int x, unsigned int y){
    Ray rays[4];
    Vector color;
    double x_alias = x_min + (2*x+1)*(w_screen/WIDTH)/2.0f;
    double y_alias = y_min + (2*y+1)*(h_screen/HEIGHT)/2.0f;

    rays[0] = {cameraVec,  normalise(Vector(x_alias-(w_screen/WIDTH)/4.0f,y_alias,-1.0))};
    rays[1] = {cameraVec,  normalise(Vector(x_alias+(w_screen/WIDTH)/4.0f,y_alias,-1.0))};
    rays[2] = {cameraVec,  normalise(Vector(x_alias,y_alias-(h_screen/HEIGHT)/4.0f,-1.0))};
    rays[3] = {cameraVec,  normalise(Vector(x_alias,y_alias+(h_screen/HEIGHT)/4.0f,-1.0))};
    for (int k = 0; k<4; k++)
    {
        color = color + getLightingModel(rays[k], 0);
    }
    color = color/4;
    plot_pixel(x, y, (int)(color.x * 255), (int)(color.y * 255), (int)(color.z * 255));
}

void checkSphereAndTriangleIntersection(const Ray &ray, double &sphereTemp, int &sphereIdx, double &triangleTemp, int &triangleIdx){
    sphereTemp = DBL_MAX_10_EXP; sphereIdx = -1; triangleTemp = DBL_MAX_10_EXP; triangleIdx = -1;
    for (int i = 0; i<sphereNumber; i++)
  {
    double b = 2 * (Vector(ray.origin-spheres[i].position).x * ray.direction.x + Vector(ray.origin-spheres[i].position).y * ray.direction.y + Vector(ray.origin-spheres[i].position).z * ray.direction.z);
    double c = (Vector(ray.origin-spheres[i].position).x * Vector(ray.origin-spheres[i].position).x + Vector(ray.origin-spheres[i].position).y * Vector(ray.origin-spheres[i].position).y + Vector(ray.origin-spheres[i].position).z * Vector(ray.origin-spheres[i].position).z) - pow(spheres[i].radius,2);
    if ((pow(b,2) - 4*1*c) < 0) continue;
    double tmin = min(((-b + sqrt(pow(b,2) - 4*1*c))/2),((-b - sqrt(pow(b,2) - 4*1*c))/2));
    if (tmin > 1e-5 && tmin<sphereTemp)
    {
      sphereIdx = i; sphereTemp = tmin;
    }
  }
  for (int j = 0; j<triangleNumber; j++)
  {
    Vector s1 = Vector(ray.direction.y * (triangles[j].v[2].position - triangles[j].v[0].position).z - ray.direction.z * (triangles[j].v[2].position - triangles[j].v[0].position).y, ray.direction.z * (triangles[j].v[2].position - triangles[j].v[0].position).x - ray.direction.x * (triangles[j].v[2].position - triangles[j].v[0].position).z,
    ray.direction.x * (triangles[j].v[2].position - triangles[j].v[0].position).y - ray.direction.y * (triangles[j].v[2].position - triangles[j].v[0].position).x);
    Vector s2 = Vector((ray.origin - triangles[j].v[0].position).y * (triangles[j].v[1].position - triangles[j].v[0].position).z - (ray.origin - triangles[j].v[0].position).z * (triangles[j].v[1].position - triangles[j].v[0].position).y, (ray.origin - triangles[j].v[0].position).z * (triangles[j].v[1].position - triangles[j].v[0].position).x - (ray.origin - triangles[j].v[0].position).x * (triangles[j].v[1].position - triangles[j].v[0].position).z, (ray.origin - triangles[j].v[0].position).x * (triangles[j].v[1].position - triangles[j].v[0].position).y - (ray.origin - triangles[j].v[0].position).y * (triangles[j].v[1].position - triangles[j].v[0].position).x);
    double k = (s1.x * (triangles[j].v[1].position - triangles[j].v[0].position).x + s1.y * (triangles[j].v[1].position - triangles[j].v[0].position).y + s1.z * (triangles[j].v[1].position - triangles[j].v[0].position).z);

    if (k > -1e-5 && k < 1e-5)
      continue;
    double tmin = (s2.x * (triangles[j].v[2].position - triangles[j].v[0].position).x + s2.y * (triangles[j].v[2].position - triangles[j].v[0].position).y + s2.z * (triangles[j].v[2].position - triangles[j].v[0].position).z) / k;
    double b1 = (s1.x * (ray.origin - triangles[j].v[0].position).x + s1.y * (ray.origin - triangles[j].v[0].position).y + s1.z * (ray.origin - triangles[j].v[0].position).z) / k;
    double b2 = (s2.x * (ray.direction).x + s2.y * (ray.direction).y + s2.z * (ray.direction).z) / k;
    double b3 = 1-b1-b2;
    if (tmin < 1e-5 || !(b1 >= 0 && b1 <= 1) || !(b2 >=0 && b2 <= 1) || !(b3 >= 0 && b3 <= 1))
      continue;
    if (tmin < triangleTemp)
    {
      triangleIdx = j; triangleTemp = tmin;
    }
  }
  doesIntersectSphere = (sphereIdx != -1); doesIntersectTriangle = (triangleIdx != -1);
}
#ifdef __APPLE_CC__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include "Mesh.h"

int gridX = 600;
int gridY = 600;
int gridZ = 600;

const double fovy = 50.;
const double clipNear = .01;
const double clipFar = 1000.;
double x = 0;
double y = 0;
double z = -2.5;

std::string path = "/Users/rohansawhney/Desktop/developer/C++/deformation/bunny.obj";
Mesh mesh;
bool dragging = false;
bool success = true;
int handleIndex = -1;
int iterations = 3;

void printInstructions()
{
    std::cerr << "Drag mouse to deform\n"
              << "→/←: increase/decrease iterations\n"
              << "↑/↓: move in/out\n"
              << "w/s: move up/down\n"
              << "a/d: move left/right\n"
              << "r: reload mesh\n"
              << "escape: exit program\n"
              << std::endl;
}

void init()
{
    glClearColor(0.0, 0.0, 0.0, 0.0);
    glEnable(GL_DEPTH_TEST);
}

void setHandleAndAnchors()
{
    double maxY = -INFINITY;
    for (VertexIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++) {
        if (v->position.y() < -0.4) v->anchor = true;
        if (maxY < v->position.y()) {
            maxY = v->position.y();
            handleIndex = v->index;
        }
    }
    
    mesh.vertices[handleIndex].handle = true;
}

void draw()
{
    glColor4f(0.0, 0.0, 1.0, 0.5);
    glLineWidth(1.0);
    glBegin(GL_LINES);
    for (EdgeCIter e = mesh.edges.begin(); e != mesh.edges.end(); e ++) {
        Eigen::Vector3d a = e->he->vertex->position;
        Eigen::Vector3d b = e->he->flip->vertex->position;
            
        glVertex3d(a.x(), a.y(), a.z());
        glVertex3d(b.x(), b.y(), b.z());
    }
    
    glEnd();
    
    for (VertexIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++) {
        if (v->handle) {
            glColor4f(0.0, 1.0, 0.0, 0.5);
            glPointSize(4.0);
            glBegin(GL_POINTS);
            glVertex3d(v->position.x(), v->position.y(), v->position.z());
            glEnd();
        }
        
        if (v->anchor) {
            glColor4f(1.0, 0.0, 0.0, 0.5);
            glPointSize(2.0);
            glBegin(GL_POINTS);
            glVertex3d(v->position.x(), v->position.y(), v->position.z());
            glEnd();
        }
    }
}

void display()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    
    GLint viewport[4];
    glGetIntegerv(GL_VIEWPORT, viewport);
    double aspect = (double)viewport[2] / (double)viewport[3];
    gluPerspective(fovy, aspect, clipNear, clipFar);
    
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    
    gluLookAt(0, 0, z, x, y, 0, 0, 1, 0);
    
    if (success) {
        draw();
    }

    glutSwapBuffers();
}

void worldSpaceCoords(double& x, double& y, double z)
{
    GLdouble modelMatrix[16];
    GLdouble projMatrix[16];
    GLint viewport[4];
    
    glGetIntegerv(GL_VIEWPORT, viewport);
    glGetDoublev(GL_MODELVIEW_MATRIX, modelMatrix);
    glGetDoublev(GL_PROJECTION_MATRIX, projMatrix);
    
    GLdouble nx, ny, nz;
    gluUnProject(x, viewport[1]+viewport[3]-y, 0, modelMatrix, projMatrix, viewport, &nx, &ny, &nz);
    
    GLdouble fx, fy, fz;
    gluUnProject(x, viewport[1]+viewport[3]-y, 1, modelMatrix, projMatrix, viewport, &fx, &fy, &fz);
    
    if(nz == fz) return;
    
    GLfloat t = (nz - z) / (nz - fz);
    x = nx + (fx - nx) * t;
    y = ny + (fy - ny) * t;
}

void mouseMove(int x, int y)
{
    if (dragging) {
        VertexIter v = mesh.vertices.begin() + handleIndex;
        
        double wx = x, wy = y;
        worldSpaceCoords(wx, wy, v->position.z());
        
        v->position.x() = wx;
        v->position.y() = wy;
        mesh.deform(iterations);
        
        glutPostRedisplay();
    }
}

void mousePressed(int button, int state, int x, int y)
{
    if (button == GLUT_LEFT_BUTTON) {
        if (state == GLUT_DOWN) {
            dragging = true;
            
        } else {
            dragging = false;
        }
    }
}

void keyboard(unsigned char key, int x0, int y0)
{
    switch (key) {
        case 27 :
            exit(0);
        case 'a':
            x -= 0.03;
            break;
        case 'd':
            x += 0.03;
            break;
        case 'w':
            y += 0.03;
            break;
        case 's':
            y -= 0.03;
            break;
        case 'r':
            mesh.read(path);
            setHandleAndAnchors();
            break;
    }
    
    glutPostRedisplay();
}

void special(int i, int x0, int y0)
{
    switch (i) {
        case GLUT_KEY_UP:
            z += 0.03;
            break;
        case GLUT_KEY_DOWN:
            z -= 0.03;
            break;
        case GLUT_KEY_LEFT:
            if (iterations > 1) iterations --;
            break;
        case GLUT_KEY_RIGHT:
            iterations ++;
            break;
    }
    
    std::string title = "Deformation, iterations: " + std::to_string(iterations);
    glutSetWindowTitle(title.c_str());
    
    glutPostRedisplay();
}

int main(int argc, char** argv) {
    
    success = mesh.read(path);
    setHandleAndAnchors();

    printInstructions();
    glutInitWindowSize(gridX, gridY);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
    glutInit(&argc, argv);
    std::string title = "Deformation, iterations: " + std::to_string(iterations);
    glutCreateWindow(title.c_str());
    init();
    glutDisplayFunc(display);
    glutMouseFunc(mousePressed);
    glutMotionFunc(mouseMove);
    glutKeyboardFunc(keyboard);
    glutSpecialFunc(special);
    glutMainLoop();
    
    return 0;
}

// original js code from here : https://github.com/royvanrijn/mandelbulb.js/blob/master/mandelbulb.html
#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;
#include <time.h>

#define PI 3.1415926535897932384626433832795f

#define IMAGE_TGA 1
#define IMAGE_BMP 2
const int fileFormat = IMAGE_BMP;

bool imageSave(const char* name, const int w, const int h, const unsigned char* rgb, const int fileFormat = IMAGE_BMP) {
    ofstream imageFile(name, ios_base::binary);
    if (!imageFile) return false;
    // header
    if (fileFormat == IMAGE_TGA) {
        // specific IMAGE_TGA handling code
        imageFile.put(0).put(0);
        imageFile.put(2);                  /* uncompressed RGB */

        imageFile.put(0).put(0);
        imageFile.put(0).put(0);
        imageFile.put(0);

        imageFile.put(0).put(0);           /* X origin */
        imageFile.put(0).put(0);           /* y origin */

        imageFile.put((w & 0x00FF)).put((w & 0xFF00) / 256);
        imageFile.put((h & 0x00FF)).put((h & 0xFF00) / 256);
        imageFile.put(24);                 /* 24 bit bitmap */
        imageFile.put(0);
    }
    else if (fileFormat == IMAGE_BMP) {
        /* header 14 bytes */
        int header1 = 14, header2 = 40, headerSize = header1 + header2;
        // int v = myScene.sizex * myScene.sizey * 3 + headerSize;  // file size
        int v = w * h * 3 + headerSize;  // file size
        imageFile.put('B').put('M');
        imageFile.put(v & 0xFF).put((v >> 8) & 0xFF).put((v >> 16) & 0xFF).put((v >> 24) & 0xFF);
        imageFile.put(0).put(0);
        imageFile.put(0).put(0);
        imageFile.put(headerSize).put(0).put(0).put(0);
        /* DIB header 40 bytes */
        imageFile.put(header2).put(0).put(0).put(0);
        v = w; imageFile.put(v & 0xFF).put((v >> 8) & 0xFF).put((v >> 16) & 0xFF).put((v >> 24) & 0xFF);
        v = h; imageFile.put(v & 0xFF).put((v >> 8) & 0xFF).put((v >> 16) & 0xFF).put((v >> 24) & 0xFF);
        imageFile.put(1).put(0);    // nb color planes
        imageFile.put(24).put(0);   // nb bits per pixel
        imageFile.put(0).put(0).put(0).put(0);  // compression
        imageFile.put(0).put(0).put(0).put(0);  // image size
        imageFile.put((char)0xc4).put(0x0e).put(0).put(0);  // horizontal resolution
        imageFile.put((char)0xc4).put(0x0e).put(0).put(0);  // vertical resolution
        imageFile.put(0).put(0).put(0).put(0);  // nb of colors
        imageFile.put(0).put(0).put(0).put(0);  // nb of important colors
    }
    // data
    if (fileFormat == IMAGE_TGA) {
        for (int i = 0; i < w * h * 3; i++) imageFile.put(rgb[i]);
    }
    else if (fileFormat == IMAGE_BMP) {
        for (int i = h - 1; i >= 0; i--) {
            for (int j = 0; j < w; j++) {
                imageFile.put(rgb[(i * w + j) * 3 + 2]).put(rgb[(i * w + j) * 3 + 1]).put(rgb[(i * w + j) * 3]);
            }
        }
    }
    // close
    imageFile.close();
    return true;
}
void die(const char* str) {
    fprintf(stderr, "%s\n", str);
    exit(1);
}
typedef struct{
    int w, h;
} SCENE;
SCENE myScene;
typedef struct {
    double x, y, z;
} TRIPLET;

unsigned char* rgb;

int currenty = 0;
int cHeight, cWidth;
double MAX_ITER = 100.0;
double DEPTH_OF_FIELD = 2.5;
double eyeDistanceFromNearField = 2.2;
double eyeDistanceFromNearFieldDelta = 0.0;

double halfPixel;
double pixel;

double lightAngle = 140.0;
double viewAngle = 150.0;

double smallStep = 0.01;
double bigStep = 0.02;

int Iterations = 20;
double Power = 8.0;
double PowerDelta = 0.0;

TRIPLET lightLocation, lightDirection, nearFieldLocation, viewDirection, reverseDirection, eyeLocation, pixelLocation;
TRIPLET vNUL;

/**
 * Here we change the camera position and light(s)
 */
double lightAngleDelta = 2.0, viewAngleDelta = 2.0;
void animateCamera() {
    if (lightAngleDelta != 0.0) { lightAngle += lightAngleDelta; if (lightAngle > 360.0) lightAngle -= 360.0; printf("lightAngle=%lf\n", lightAngle); }
    if (viewAngleDelta != 0.0) { viewAngle += viewAngleDelta; if (viewAngle > 360.0) viewAngle -= 360.0; printf("viewAngle=%lf\n", viewAngle); }
    if (eyeDistanceFromNearFieldDelta != 0.0) { eyeDistanceFromNearField += eyeDistanceFromNearFieldDelta; printf("eyeDistanceFromNearField=%lf\n", eyeDistanceFromNearField); }
    if (PowerDelta != 0.0){ Power += PowerDelta; printf("Power=%lf\n", Power); }
    // Power = 8 + 2 * cos(0.1*tick);

    // lightDirection.y = cos(0.1*tick);
}

__inline double min(const double a, const double b) { return a < b ? a : b; }
__inline double max(const double a, const double b) { return a > b ? a : b; }

__inline double clamp(const double n, double nmin, double nmax) {
    return max(nmin, min(n, nmax));
}

TRIPLET *clampVec(TRIPLET *v1, const double min, const double max) {
    v1->x = clamp(v1->x, min, max);
    v1->y = clamp(v1->y, min, max);
    v1->z = clamp(v1->z, min, max);
    return v1;
}
double dotProduct(const TRIPLET v1, const TRIPLET v2) {
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

double toRad(const double r) {
    return r * PI / 180.0;
}

double saturate(const double n) {
    // return clampVec(n, 0.0, 1.0);
    return clamp(n, 0.0, 1.0);
}

__inline double length(const TRIPLET z) {
    return sqrt(z.x * z.x + z.y * z.y + z.z * z.z);
}

__inline TRIPLET *scalarMultiply(TRIPLET *a, const double amount) {
    a->x *= amount;
    a->y *= amount;
    a->z *= amount;
    return a;
}

__inline TRIPLET *normalize(TRIPLET *a) {
    return scalarMultiply(a, 1 / length(*a));
}


TRIPLET *add(TRIPLET *v1, const TRIPLET v2) {
    v1->x += v2.x;
    v1->y += v2.y;
    v1->z += v2.z;
    return v1;
}

TRIPLET *subtract(TRIPLET *v1, const TRIPLET v2) {
    v1->x -= v2.x;
    v1->y -= v2.y;
    v1->z -= v2.z;
    return v1;
}

TRIPLET *setTo(TRIPLET *v1, const TRIPLET v2) {
    v1->x = v2.x;
    v1->y = v2.y;
    v1->z = v2.z;
    return v1;
}

TRIPLET *turnOrthogonal(TRIPLET *v1) {
    double inverse = 1.0 / sqrt(v1->x * v1->x + v1->z * v1->z);
    double oldX = v1->x;
    v1->x = -inverse * v1->z;
    v1->z = inverse * oldX;
    return v1;
}

TRIPLET *crossProduct(TRIPLET *v1, const TRIPLET v2) {
    double oldX = v1->x;
    double oldY = v1->y;
    v1->x = v2.y * v1->z - v2.z * oldY;
    v1->y = v2.z * oldX - v2.x * v1->z;
    v1->z = v2.x * oldY - v2.y * oldX;
    return v1;
}


void TRprint(TRIPLET v) { printf("{x:%lf, y:%lf, z:%lf}\n", v.x, v.y, v.z); }
void zero(TRIPLET *v) { v->x = v->y = v->z = 0.0; }

void setupScene() {
    // z.x = z.y = z.z = 0.0;

    // zero(&lightLocation); zero(&lightDirection); zero(&nearFieldLocation); zero(&viewDirection); zero(&reverseDirection); zero(&eyeLocation); zero(&pixelLocation);
    zero(&vNUL);    // vNUL.x = vNUL.y = vNUL.z = 0.0;

    double rad = toRad(lightAngle);
    double lightX = ((cos(rad) * DEPTH_OF_FIELD / 2));
    double lightZ = ((sin(rad) * DEPTH_OF_FIELD / 2));

    lightLocation.x = lightX;
    lightLocation.y = (DEPTH_OF_FIELD / 2);
    lightLocation.z = lightZ;

    normalize(subtract(setTo(&lightDirection, vNUL), lightLocation));

    double viewRad = toRad(viewAngle);
    double viewX = ((cos(viewRad) * DEPTH_OF_FIELD / 2));
    double viewZ = ((sin(viewRad) * DEPTH_OF_FIELD / 2));

    nearFieldLocation.x = viewX;
    nearFieldLocation.y = 0.0;
    nearFieldLocation.z = viewZ;

    normalize(subtract(setTo(&viewDirection, vNUL), nearFieldLocation));

    scalarMultiply(setTo(&reverseDirection, viewDirection), eyeDistanceFromNearField);
    subtract(setTo(&eyeLocation, nearFieldLocation), reverseDirection);
}
double mandelbulb(TRIPLET pos) {
    // Power += 0.00001;
    TRIPLET z;
    setTo(&z, pos);
    double dr = 1.0;
    double r = 0.0;
    for (int i = 0; i < Iterations; i++) {
        r = length(z);
        if (r > DEPTH_OF_FIELD) break;

        double theta = acos(z.z / r);
        double phi = atan2(z.y, z.x);
        dr = pow(r, Power - 1.0) * Power * dr + 1.0;
        double zr = pow(r, Power);
        theta *= Power;
        phi *= Power;
        double sinTheta = sin(theta);
        z.x = sinTheta * cos(phi);
        z.y = sin(phi) * sinTheta;
        z.z = cos(theta);
        add(scalarMultiply(&z, zr), pos);
    }
    return 0.5 * log(r) * r / dr;
}
double map(TRIPLET z) {
    return mandelbulb(z);
    // return sphere(z, NUL, 1);
}

/**
 * In this method we calculate the 'soft' shadows
 * From: http://www.iquilezles.org/www/articles/rmshadows/rmshadows.htm
 */
double shadow(const double mint, const double maxt, const double k, const TRIPLET rayLocation) {
    double res = 1.0;
    TRIPLET rd, ro;
    for (double t = mint; t < maxt; ) {
        scalarMultiply(setTo(&rd, lightDirection), t);
        subtract(setTo(&ro, rayLocation), rd);
        double h = map(ro);
        if (h < 0.001) {
            return 0.0;
        }
        res = min(res, k * h / t);
        t += h;
    }
    return res;
}


/**
 * The main draw function for a scanline
 * Make sure setupScene is called first after adjusting the camera and/or light
 */
int shade = 1, shadowOn = 1;
unsigned char *draw(unsigned char* imageData, const int y) {
    const int cHalfWidth = cWidth / 2;
    const int ny = y - cHeight / 2;

    TRIPLET temp;
    // zero(&temp);
    TRIPLET rayLocation, tempViewDirectionX1, tempViewDirectionX2, tempViewDirectionY;
    TRIPLET rayDirection, normal, halfway;

    scalarMultiply(crossProduct(turnOrthogonal(setTo(&tempViewDirectionY, viewDirection)), viewDirection), ny * pixel);
    turnOrthogonal(setTo(&tempViewDirectionX1, viewDirection));

    int pixels = 0; //  3 * (y * cWidth);
    for (int x = 0; x < cWidth; x++) {

        int nx = x - cHalfWidth;

        setTo(&pixelLocation, nearFieldLocation);

        scalarMultiply(setTo(&tempViewDirectionX2, tempViewDirectionX1), nx * pixel);
        add(&pixelLocation, tempViewDirectionX2);
        add(&pixelLocation, tempViewDirectionY);

        setTo(&rayLocation, pixelLocation);

        normalize(subtract(setTo(&rayDirection, rayLocation), eyeLocation));

        double distanceFromCamera = 0.0;
        double d = map(rayLocation);

        int iterations = 0;
        for (; iterations < MAX_ITER; iterations++) {

            if (d < halfPixel) break;

            // Increase rayLocation with direction and d:
            // add(&rayLocation, scalarMultiply(&rayDirection, d));
            scalarMultiply(&rayDirection, d);
            add(&rayLocation, rayDirection);

            // And reset ray direction:
            normalize(&rayDirection);

            // Move the pixel location:
            // distanceFromCamera = length(subtract(setTo(&temp, nearFieldLocation), rayLocation));
            setTo(&temp, nearFieldLocation);
            subtract(&temp, rayLocation);
            distanceFromCamera = length(temp);

            if (distanceFromCamera > DEPTH_OF_FIELD) break;

            d = map(rayLocation);
        }

        if (distanceFromCamera < DEPTH_OF_FIELD && distanceFromCamera > 0) {

            double red, green, blue;
            if (shade == 0) {
                red = 255.0 * (double)iterations / (double)MAX_ITER;
                green = blue = 0.0;
            }
            else {
                rayLocation.x -= smallStep;
                double locationMinX = map(rayLocation);
                rayLocation.x += bigStep;
                double locationPlusX = map(rayLocation);
                rayLocation.x -= smallStep;

                rayLocation.y -= smallStep;
                double locationMinY = map(rayLocation);
                rayLocation.y += bigStep;
                double locationPlusY = map(rayLocation);
                rayLocation.y -= smallStep;

                rayLocation.z -= smallStep;
                double locationMinZ = map(rayLocation);
                rayLocation.z += bigStep;
                double locationPlusZ = map(rayLocation);
                rayLocation.z -= smallStep;

                // Calculate the normal:
                normal.x = (locationMinX - locationPlusX);
                normal.y = (locationMinY - locationPlusY);
                normal.z = (locationMinZ - locationPlusZ);
                normalize(&normal);

                // Calculate the ambient light:
                double dotNL = dotProduct(lightDirection, normal);
                double diff = saturate(dotNL);

                // Calculate specular light:
                normalize(add(setTo(&halfway, rayDirection), lightDirection));

                double dotNH = dotProduct(halfway, normal);
                double spec = pow(saturate(dotNH), 35);

                double shad, brightness;
                if (shadowOn) {
                    shad = shadow(1.0, DEPTH_OF_FIELD, 16.0, rayLocation) + 0.25;
                    brightness = (10.0 + (200.0 + spec * 45.0) * shad * diff) / 270.0;
                }
                else brightness = (10.0 + (200.0 + spec * 45.0) * diff) / 270.0;

                red = 10 + (380 * brightness);
                green = 10 + (280 * brightness);
                blue = (180 * brightness);

                red = clamp(red, 0, 255.0);
                green = clamp(green, 0, 255.0);
                blue = clamp(blue, 0, 255.0);
            }

            imageData[pixels] = (unsigned char)red;
            imageData[pixels + 1] = (unsigned char)green;
            imageData[pixels + 2] = (unsigned char)blue;
            // printf("rgb=[%d %d %d] ", imageData[pixels], imageData[pixels+1], imageData[pixels+2]);
        }
        else {
            imageData[pixels] = 30/*155*/ + (unsigned char)clamp(iterations * 1.5, 0.0, 100.0);
            imageData[pixels + 1] = 30/*205*/; // + (unsigned char)clamp(iterations * 1.5, 0.0, 50.0);
            imageData[pixels + 2] = 30/*255*/;
        }
        pixels += 3;
    }
    return imageData;
}
char fout[1024];
int tick = 0;
void animate() {
    tick ++ ;
    if (currenty == 0) {
        animateCamera();
        setupScene();
    }
#pragma omp parallel for
    for (int icurrenty = 0; icurrenty < cHeight; icurrenty++) {
        // printf("draw line %d...\n", currenty);
        draw(rgb + icurrenty * (3 * cWidth), icurrenty);
    }

    char fout_anim[1024];
    sprintf(fout_anim, "%s_%04d.bmp", fout, tick);
    fprintf(stderr, "saving %s...\n", fout_anim);
    imageSave(fout_anim, cWidth, cHeight, rgb, IMAGE_BMP);

    currenty = 0;
}
clock_t  t1, t2;
void load_ini() {
    FILE* f = fopen("mandelbulb.ini", "r");
    char buf[512], w[80], kw[80];
    while (fgets(buf, 512, f)) {
        if (buf[0] == '\n' || buf[0] == '\r' || buf[0] == '\0' || buf[0] == '#') continue;
        sscanf(buf, "%s", &w);
        strcpy(kw, "viewAngleDelta"); if (!strncmp(w, kw, strlen(kw))) { sscanf(buf + strlen(kw) + 1, "%lf", &viewAngleDelta); continue; }
        strcpy(kw, "viewAngle"); if (!strncmp(w, kw, strlen(kw))) { sscanf(buf + strlen(kw) + 1, "%lf", &viewAngle); continue; }
        strcpy(kw, "lightAngleDelta"); if (!strncmp(w, kw, strlen(kw))) { sscanf(buf + strlen(kw) + 1, "%lf", &lightAngleDelta); continue; }
        strcpy(kw, "lightAngle"); if (!strncmp(w, kw, strlen(kw))) { sscanf(buf + strlen(kw) + 1, "%lf", &lightAngle); continue; }
        strcpy(kw, "MAX_ITER"); if (!strncmp(w, kw, strlen(kw))) { sscanf(buf + strlen(kw) + 1, "%lf", &MAX_ITER); continue; }
        strcpy(kw, "eyeDistanceFromNearFieldDelta"); if (!strncmp(w, kw, strlen(kw))) { sscanf(buf + strlen(kw) + 1, "%lf", &eyeDistanceFromNearFieldDelta); continue; }
        strcpy(kw, "eyeDistanceFromNearField"); if (!strncmp(w, kw, strlen(kw))) { sscanf(buf + strlen(kw) + 1, "%lf", &eyeDistanceFromNearField); continue; }
        
        strcpy(kw, "DEPTH_OF_FIELD"); if (!strncmp(w, kw, strlen(kw))) { sscanf(buf + strlen(kw) + 1, "%lf", &DEPTH_OF_FIELD); continue; }
        strcpy(kw, "PowerDelta"); if (!strncmp(w, kw, strlen(kw))) { sscanf(buf + strlen(kw) + 1, "%lf", &PowerDelta); continue; }
        strcpy(kw, "Power"); if (!strncmp(w, kw, strlen(kw))) { sscanf(buf + strlen(kw) + 1, "%lf", &Power); continue; }
        strcpy(kw, "shade"); if (!strncmp(w, kw, strlen(kw))) { sscanf(buf + strlen(kw) + 1, "%d", &shade); if (shade != 0 && shade != 1) die("shade"); continue; }
        strcpy(kw, "shadowOn"); if (!strncmp(w, kw, strlen(kw))) { sscanf(buf + strlen(kw) + 1, "%d", &shadowOn); if (shadowOn != 0 && shadowOn != 1) die("shadowOn"); continue; }
        fprintf(stderr, "ini file : option '%s' not known.\n", w);
    }
    // printf("lightAngleDelta=%lf\n", lightAngleDelta); die("OK");
}
int main(int argc, char* argv[]) {
    if (argc < 3){
        fprintf(stderr, "Usage : %s <size> <output> [nb]\n", argv[0]);
        return -1;
    }
    load_ini();
    if (false) {
        TRIPLET z;
        TRprint(z);
        z.x = 3.14;
        TRprint(z);
        zero(&z);
        TRprint(z);
        die("OK");
    }
    sscanf(argv[1], "%d", &cWidth);
    // cWidth = 128;
    cHeight = (int)((double)cWidth / ( 1920.0 / 1080.0));

    strcpy(fout, argv[2]);
    
    int nb = 1; if (argc > 3) sscanf(argv[3], "%d", &nb);

    printf("Create %d images %s or size %dx%d.\n", nb, fout, cWidth, cHeight);

    // create rgb buffer
    rgb = (unsigned char*)malloc(cWidth * cHeight * 3);
    if (rgb == (unsigned char*)NULL) die("Can't alloc memory");

    pixel = (DEPTH_OF_FIELD) / ((cHeight + cWidth) / 2);
    halfPixel = pixel / 2;

    long dtTotal = 0;
    for (int i = 0; i < nb; i++) {
        t1 = clock();
        animate();
        // Power += 0.00001;

        t2 = clock();
        long dt = t2 - t1;
        dtTotal += dt;
        printf("%d  (%.3lf sec.)%c", i+1, ((double)dt) / CLOCKS_PER_SEC, 13);
    }
    printf("Finished. %d images generated in %.3lf sec (avg=%.3lf).\n", nb, ((double)dtTotal) / CLOCKS_PER_SEC, ((double)dtTotal) / CLOCKS_PER_SEC / (double)nb);
}

#ifndef _SIMULATION_H_
#define _SIMULATION_H_

#include <stdio.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <signal.h>
#include <windows.h>
#include <sys/time.h>
#define NameI(x) {#x, &x, N_I, sizeof (x) / sizeof (int)}
#define NameR(x) {#x, &x, N_R, sizeof (x) / sizeof (real)}
#define NP_I ((int *) (nameList[k].vPtr) + j)
#define NP_R ((real *) (nameList[k].vPtr) + j)
#define NDIM 2
#define CHAR_MINUS '-'
#define CHAR_ZERO '0'
#define FALSE 0
#define TRUE 1
#define ReadF(x) fread (&x, sizeof (x), 1, fp)
#define ReadFN(x, n) fread (x, sizeof (x[0]), n, fp)
#define WriteF(x) fwrite (&x, sizeof (x), 1, fp)
#define WriteFN(x, n) fwrite (x, sizeof (x[0]), n, fp)

#define VSCopy(v2, s1, v1) \
(v2).x = (s1) * (v1).x, \
(v2).y = (s1) * (v1).y
#define VProd(v) ((v).x * (v).y)
//here
#define VCSum(v) ((v).x + (v).y)

#define VMul(v1, v2, v3) \
(v1).x = (v2).x * (v3).x, \
(v1).y = (v2).y * (v3).y
#define VAdd(v1, v2, v3) \
(v1).x = (v2).x + (v3).x, \
(v1).y = (v2).y + (v3).y

#define VDiv(v1, v2, v3) \
(v1).x = (v2).x / (v3).x 

#define VScale(v, s) \
(v).x *= s, \
(v).y *= s
#define VVAdd(v1, v2) VAdd (v1, v1, v2)

#define AllocMem(a, n, t) a = (t *) malloc ((n) * sizeof (t))

#define VSub(v1, v2, v3) \
(v1).x = (v2).x - (v3).x, \
(v1).y = (v2).y - (v3).y

#define VDot(v1, v2) \
((v1).x * (v2).x + (v1).y * (v2).y)

#define VSAdd(v1, v2, s3, v3) \
(v1).x = (v2).x + (s3) * (v3).x, \
(v1).y = (v2).y + (s3) * (v3).y

#define VSet(v, sx, sy) \
(v).x = sx, \
(v).y = sy

#define VSetAll(v, s) VSet(v, s, s)

#define VZero(v) VSetAll(v, 0)

#define VVSAdd(v1, s2, v2) VSAdd(v1, v1, s2, v2)

#define VLenSq(v) VDot(v, v)

#define VWrap(v, t) \
if (v.t >= 0.5 * region.t) v.t -= region.t; \
else if (v.t < -0.5 * region.t) v.t += region.t

#define VWrapAll(v) \
{VWrap(v, x); \
VWrap(v, y);}

#define Sqr(x) ((x) * (x))

#define Cube(x) ((x) * (x) * (x))

#define DO_MOL for (n = 0; n < nMol; n++)

#define PropZero(v) \
v.sum = 0., \
v.sum2 = 0.
#define PropAccum(v) \
v.sum += v.val,\
v.sum2 += Sqr (v.val)

#define PropAvg(v, n) \
v.sum /= n, \
v.sum2 = sqrt (Max (v.sum2 / n - Sqr (v.sum), 0.))
#define PropEst(v) v.sum, v.sum2

#define Min(x1, x2) (((x1) < (x2)) ? (x1) : (x2))

#define Max(x1, x2) (((x1) > (x2)) ? (x1) : (x2))
#define VLen(v) sqrt (VDot (v, v))
int countVel, limitVel, sizeHistVel, stepVel;

#define MAX_STRING_LENGTH 100
#define IADD 453806245
#define IMUL 314159269
#define MASK 2147483647
#define SCALE 0.4656612873e-9
int randSeedP = 17;

enum {ERR_NONE, ERR_BOND_SNAPPED, ERR_CHECKPT_READ, ERR_CHECKPT_WRITE,
ERR_COPY_BUFF_FULL, ERR_EMPTY_EVPOOL, ERR_MSG_BUFF_FULL,
ERR_OUTSIDE_REGION, ERR_SNAP_READ, ERR_SNAP_WRITE,
ERR_SUBDIV_UNFIN, ERR_TOO_MANY_CELLS, ERR_TOO_MANY_COPIES,
ERR_TOO_MANY_LAYERS, ERR_TOO_MANY_LEVELS, ERR_TOO_MANY_MOLS,
ERR_TOO_MANY_MOVES, ERR_TOO_MANY_NEBRS, ERR_TOO_MANY_REPLICAS};
char *errorMsg[] = {"", "bond snapped", "read checkpoint data",
"write checkpoint data", "copy buffer full", "empty event pool",
"message buffer full", "outside region", "read snap data",
"write snap data", "subdivision unfinished", "too many cells",
"too many copied mols", "too many layers", "too many levels",
"too many mols", "too many moved mols", "too many neighbors",
"too many replicas"};


int doCheckpoint = FALSE;


enum {FL_CHECKA, FL_CHECKB, FL_CKLAST, FL_SNAP};
char *fileNameR[] = {"xxnnchecka.data", "xxnncheckb.data",
"xxnncklast.data", "xxnnsnap.data"}, fileName[5][20];
char *progId = "md";
int runId;
int newRun;
int recordSnap;


typedef struct {
    double x, y, z;
} Vector3D;

typedef struct {
    char* pigmentColor;
    double ambient;
    double diffuse;
    double specular;
    double reflection;
    double refraction;
} Finish;

typedef struct {
    Vector3D position;
    double radius;
    Finish finish;
} Sphere;

typedef struct {
    char* fileName;
    FILE* file;
} POVRayWriter;

POVRayWriter* writer;

typedef double real;
typedef struct {
    real val, sum, sum2;
} Prop;

typedef struct {
    real x, y, z;
} VecR;

typedef struct {
    VecR r, rv, ra;
} Mol;

typedef struct {
    int x, y;
} VecI;

typedef enum {N_I, N_R} VType;

typedef struct {
    char *vName;
    void *vPtr;
    VType vType;
    int vLen, vStatus;
} NameList;

Mol *mol;
VecR region, vSum;
VecI initUcell;
Prop kinEnergy, pressure, totEnergy;
real deltaT, density, rCut, temperature, timeNow, uSum, velMag, virSum, vvSum, *histVel, rangeVel;
int moreCycles, nMol, stepAvg, stepCount, stepEquil, stepLimit;

NameList nameList[] = {
    NameR (deltaT),
    NameR (density),
    NameI (initUcell),
    NameI (stepAvg),
    NameI (stepEquil),
    NameI (stepLimit),
    NameR (temperature)
};

// Function declarations
int GetNameList(int argc, char **argv);
void PrintNameList(FILE *file);
void SetParams();
void SetupJob();
void SingleStep();
void ComputeForces();
#endif/*_SIMULATION_H_*/
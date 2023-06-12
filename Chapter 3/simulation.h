#ifndef _SIMULATION_H_
#define _SIMULATION_H_

#include <stdio.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <signal.h>
#include <sys/time.h>
#include <windows.h>

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
(v2).y = (s1) * (v1).y, \
(v2).y = (s1) * (v1).z

#define VProd(v) ((v).x * (v).y * (v).z)

#define VCSum(v) ((v).x + (v).y + (v).z)

#define VMul(v1, v2, v3) \
(v1).x = (v2).x * (v3).x, \
(v1).y = (v2).y * (v3).y, \
(v1).z = (v2).z * (v3).z

#define VAdd(v1, v2, v3) \
(v1).x = (v2).x + (v3).x, \
(v1).y = (v2).y + (v3).y, \
(v1).z = (v2).z + (v3).z

#define VDiv(v1, v2, v3) \
(v1).x = (v2).x / (v3).x, \
(v1).y = (v2).y / (v3).y,\
(v1).z = (v2).z / (v3).z

#define VScale(v, s) \
(v).x *= s, \
(v).y *= s, \
(v).z *= s

#define VVAdd(v1, v2) VAdd (v1, v1, v2)

#define AllocMem(a, n, t) a = (t *) malloc ((n) * sizeof (t))

#define VSub(v1, v2, v3) \
(v1).x = (v2).x - (v3).x, \
(v1).y = (v2).y - (v3).y, \
(v1).z = (v2).z - (v3).z

#define VDot(v1, v2) \
((v1).x * (v2).x + (v1).y * (v2).y + (v1).z * (v2).z)

#define VSAdd(v1, v2, s3, v3) \
(v1).x = (v2).x + (s3) * (v3).x, \
(v1).y = (v2).y + (s3) * (v3).y, \
(v1).z = (v2).z + (s3) * (v3).z

#define VSet(v, sx, sy, sz) \
(v).x = sx, \
(v).y = sy, \
(v).z = sz

#define VSetAll(v, s) VSet(v, s, s, s)

#define VZero(v) VSetAll(v, 0)

#define VVSAdd(v1, s2, v2) VSAdd(v1, v1, s2, v2)
#define MAX_STRING_LENGTH 100
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

#define IADD 453806245
#define IMUL 314159269
#define MASK 2147483647
#define SCALE 0.4656612873e-9
int randSeedP = 17;
int stepAdjustTemp;

#define N_OFFSET 14

#define OFFSET_VALS \
{{0,0,0}, {1,0,0}, {1,1,0}, {0,1,0}, {-1,1,0}, {0,0,1}, \
{1,0,1}, {1,1,1}, {0,1,1}, {-1,1,1}, {-1,0,1}, \
{-1,-1,1}, {0,-1,1}, {1,-1,1}}

#define VCellWrap(t) \
if (m2v.t >= cells.t) { \
    m2v.t = 0; \
    shift.t = region.t; \
} else if (m2v.t < 0) { \
    m2v.t = cells.t - 1; \
    shift.t = -region.t; \
}

#define VLinear(p, s) \
(((p).z * (s).y + (p).y) * (s).x + (p).x)

#define VCellWrapAll() \
{ \
    VCellWrap(x); \
    VCellWrap(y); \
    VCellWrap(z); \
}
#define VVSub(v1, v2) VSub (v1, v1, v2)
#define DO_CELL(j, m) \
for (j = cellList[m]; j >= 0; j = cellList[j])
#define PCR4(r, ro, v, a, a1, a2, t) \
r.t = ro.t + deltaT * v.t + wr * (cr[0] * a.t + \
cr[1] * a1.t + cr[2] * a2.t)
#define PCV4(r, ro, v, a, a1, a2, t) \
v.t = (r.t - ro.t) / deltaT + wv * (cv[0] * a.t + \ 
cv[1] * a1.t + cv[2] * a2.t)
#define PR(t) \
PCR4 (mol[n].r, mol[n].r, mol[n].rv, mol[n].ra, \
mol[n].ra1, mol[n].ra2, t)
#define PRV(t) \
PCV4 (mol[n].r, mol[n].ro, mol[n].rv, mol[n].ra, \ 
mol[n].ra1, mol[n].ra2, t)
#define CR(t) \
PCR4 (mol[n].r, mol[n].ro, mol[n].rvo, mol[n].ra, \
mol[n].ra1, mol[n].ra2, t)
#define CRV(t) \ 
PCV4 (mol[n].r, mol[n].ro, mol[n].rv, mol[n].ra, \
mol[n].ra1, mol[n].ra2, t)


#define VComp(v, k) \
*((k == 0) ? &(v).x : ((k == 1) ? &(v).y : &(v).z))

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
int stepInitlzTemp;

enum {FL_CHECKA, FL_CHECKB, FL_CKLAST, FL_SNAP};
char *fileNameR[] = {"xxnnchecka.data", "xxnncheckb.data",
"xxnncklast.data", "xxnnsnap.data"}, fileName[5][20];
char *progId = "md";
int runId;
int newRun;
int recordSnap;



typedef double real;
typedef struct {
    real val, sum, sum2;
} Prop;

typedef struct {
    real x, y, z;
} VecR;

typedef struct {
    real u1, u2, u3, u4;
} Quat;

typedef struct {
    VecR r;     // position
    VecR rv;    // velocity
    VecR ra;    // acceleration
    VecR ra1;   // acceleration at time step (t - Δt)
    VecR ra2;   // acceleration at time step (t - 2Δt)
    VecR ro;    // position at time step (t - Δt)
    VecR rvo;   // velocity at time step (t - Δt)
    Quat q, qv, qa, qa1, qa2, qo, qvo;
    VecR torq;
} Mol; 

typedef struct {
    int x, y, z;
} VecI;


VecI cells;
int *cellList;

typedef enum {N_I, N_R} VType;

typedef struct {
    char *vName;
    void *vPtr;
    VType vType;
    int vLen, vStatus;
} NameList;

/**
 * @brief Struct defining an edge.
 */
typedef struct {
    int f[2]; ///< Array of face indices connected to the edge.
    int v[2]; ///< Array of vertex indices connected to the edge.
    int stat; ///< Status of the edge.
} Edge;

/**
 * @brief Struct defining a face.
 */
typedef struct {
    real dist; ///< Distance of the face.
    int fPtr; ///< Pointer to another face.
    int stat; ///< Status of the face.
    int vFar; ///< Index of the farthest vertex.
} Face;

/**
 * @brief Struct defining a list of faces.
 */
typedef struct {
    int e; ///< Index of an edge.
    int link; ///< Index of the next entry in the list.
    int v; ///< Index of a vertex.
} Flist;

/**
 * @brief Struct defining a vertex.
 */
typedef struct {
    VecR pos;        // Position of the vertex.
    real distSq;     // Square distance of the vertex.
    int e[3];        // Array of edge indices connected to the vertex.
    int stat;        // Status of the vertex.
    Flist* faces;    // List of faces connected to the vertex.
} Vert;



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

Mol *mol;
VecR region, vSum;
VecI initUcell;
Prop kinEnergy, pressure, totEnergy;
real *valTrajDev, pertTrajDev, deltaT, density, rCut, temperature, timeNow, uSum, velMag, virSum, vvSum, *histVel, rangeVel, dispHi, rNebrShell, kinEnInitSum;
int countTrajDev, limitTrajDev, stepTrajDev, moreCycles, nMol, stepAvg, stepCount, stepEquil, stepLimit, *nebrTab, nebrNow, nebrTabFac, nebrTabLen, nebrTabMax, stepInitlzTemp;
NameList nameList[] = {
    NameR (deltaT),
    NameR (density),
    NameI (initUcell),
    NameI (stepAvg),
    NameI (stepEquil),
    NameI (stepLimit),
    NameI (stepInitlzTemp),
    NameR (temperature),
    NameI (nebrTabFac),
    NameI (stepInitlzTemp),
    NameI (stepAdjustTemp),
    NameI (limitTrajDev),
    NameR (pertTrajDev),
    NameI (stepTrajDev),
    NameR (rNebrShell)
};

// Function declarations
int GetNameList(int argc, char **argv);
void PrintNameList(FILE *file);
void SetParams();
void SetupJob();
void SingleStep();
void ComputeForces();
void BuildNebrList();
#endif/*_SIMULATION_H_*/
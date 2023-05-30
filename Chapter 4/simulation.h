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
#define MAX_EDGE 200
#define MAX_FACE 50
#define MAX_FLIST 500
#define MAX_ITEM 50
#define MAX_VERT 200

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
(v1).x = (v2).x / (v3).x

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

#define VShift(v, t) \
    if (v.t >= 0.5 * region.t) shift.t -= region.t; \
    else if (v.t < -0.5 * region.t) shift.t += region.t
#define VShiftAll(v) \
    {VShift (v, x); \ 
    VShift (v, y); \
    VShift (v, z);}
/**
 * \def VCross(v1, v2, v3)
 * \brief Compute the cross product of two vectors.
 *
 * This macro calculates the cross product of two vectors: v1 = v2 x v3.
 * It assigns the resulting components of the cross product to v1.
 *
 * \param v1 The output vector storing the cross product.
 * \param v2 The first input vector.
 * \param v3 The second input vector.
 */
#define VCross(v1, v2, v3) \
    (v1).x = (v2).y * (v3).z - (v2).z * (v3).y, \
    (v1).y = (v2).z * (v3).x - (v2).x * (v3).z, \
    (v1).z = (v2).x * (v3).y - (v2).y * (v3).x


/**
 * \def VInterp(v1, s2, v2, v3)
 * \brief Interpolates between two vectors.
 *
 * This macro calculates the interpolation between two vectors v2 and v3, using the weight factor s2.
 * The interpolated vector is stored in v1.
 *
 * \param v1 The interpolated vector.
 * \param s2 The weight factor for v2.
 * \param v2 The first vector.
 * \param v3 The second vector.
 */
#define VInterp(v1, s2, v2, v3) \
    VSSAdd(v1, s2, v2, 1. - (s2), v3)
/**
 * \def VSSAdd(v1, s2, v2, s3, v3)
 * \brief Adds two scaled vectors and stores the result in v1.
 *
 * This macro calculates the sum of two scaled vectors v2 and v3, using the weight factors s2 and s3.
 * The result is stored in v1.
 *
 * \param v1 The resulting vector.
 * \param s2 The weight factor for v2.
 * \param v2 The first vector.
 * \param s3 The weight factor for v3.
 * \param v3 The second vector.
 */
#define VSSAdd(v1, s2, v2, s3, v3) \
    (v1).x = (s2) * (v2).x + (s3) * (v3).x, \
    (v1).y = (s2) * (v2).y + (s3) * (v3).y, \
    (v1).z = (s2) * (v2).z + (s3) * (v3).z

#define SCALE_FAC 32767.
#define VToLin(a, n, v) \
    a[(n) + 0] = (v).x, \
    a[(n) + 1] = (v).y, \
    a[(n) + 2] = (v).z
#define VFromLin(v, a, n) \
    VSet (v, a[(n) + 0], a[(n) + 1], a[(n) + 2])

/**
 * \brief Set all components of a vector to the same value.
 */
#define VSetAll(v, s) \
    VSet(v, s, s, s)

/**
 * \brief Add a constant to each component of a vector.
 */
#define VAddCon(v1, v2, s) \
    (v1).x = (v2).x + (s), \
    (v1).y = (v2).y + (s), \
    (v1).z = (v2).z + (s)

/**
 * \brief Calculate the product of vector components.
 */
#define VProd(v) \
    ((v).x * (v).y * (v).z)

/**
 * \brief Compare if vector v1 is greater than or equal to vector v2 component-wise.
 */
#define VGe(v1, v2) \
    ((v1).x >= (v2).x && (v1).y >= (v2).y && (v1).z >= (v2).z)

/**
 * \brief Compare if vector v1 is less than vector v2 component-wise.
 */
#define VLt(v1, v2) \
    ((v1).x < (v2).x && (v1).y < (v2).y && (v1).z < (v2).z)

/**
 * \brief Compute the linear index of a 3D point in a grid of size s.
 */
#define VLinear(p, s) \
    (((p).z * (s).y + (p).y) * (s).x + (p).x)

/**
 * \brief Compute the sum of vector components.
 */
#define VCSum(v) \
    ((v).x + (v).y + (v).z)

/**
 * \brief Access a specific component of a vector.
 */
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
    VecR r, rv, ra, ra1, ra2, ro, rvo;
    Quat q, qv, qa, qa1, qa2, qo, qvo;
    VecR torq;
    int inCell;
    int inClust;   ///< Flag indicating whether the molecule is in a cluster.
} Mol;

/**
 * \brief Structure representing a cluster.
 */
typedef struct {
    int head;   ///< Index of the head molecule in the cluster.
    int next;   ///< Index of the next molecule in the cluster.
    int size;   ///< Size of the cluster.
} Clust; 



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
/**
 * \brief Struct representing an edge.
 */
typedef struct {
    int f[2], v[2], stat;
} Edge;

/**
 * \brief Struct representing a face.
 */
typedef struct {
    real dist;
    int fPtr, stat, vFar;
} Face;

/**
 * \brief Struct representing an entry in the face list.
 */
typedef struct {
    int e, link, v;
} Flist;

/**
 * \brief Struct representing a vertex.
 */
typedef struct {
    VecR pos;
    real distSq;
    int e[3], stat;
} Vert;

Mol *mol;
VecR region, vSum;
VecI initUcell, cells;
int *cellList;
Prop kinEnergy, pressure, totEnergy;
real deltaT, density, rCut, temperature, timeNow, uSum, velMag, virSum, vvSum, *histVel, rangeVel, dispHi, rNebrShell, kinEnInitSum, *histRdf, rangeRdf,  latticeCorr, *distSq, cellRatio, eulerSum, fParamS, fracPolyVol, rangeLim, regionVol, timeNow, vDistSqMax;
int moreCycles, stepSnap, nMol, bigSize, nCellEdge, nClust, nSingle, stepAvg, rClust, stepCount, stepEquil, stepLimit, *nebrTab, nebrNow, nebrTabFac, nebrTabLen, nebrTabMax, stepInitlzTemp, countRdf, limitRdf, sizeHistRdf, stepRdf,
 *cellList, *eCut, *eDel, *eNew, *fCut, *fDel, *siteSeq, *testSites, *vDel, blockNum, blockSize, curSite, eLast, eLastP, fLast, fListLast, nCell, neCut, neDel, neNew, nfCut, nfDel, nMol, nTestSites, nvDel, runId, siteA, siteB, stepCount, vLast;
Edge *edge;
Face *face;
Flist *flist;
Vert *vert;
VecR *r, fParamV, region;
VecI cells;
Prop polyGeom[4], polyArea, polyVol,cSize;
FILE *fp;


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
    NameR (rNebrShell),
    NameI (limitRdf),
    NameR (rangeRdf),
    NameI (sizeHistRdf),
    NameI (stepRdf),
    NameI (stepSnap)
};
// Function declarations
int GetNameList(int argc, char **argv);
void PrintNameList(FILE *file);
void SetParams();
void SetupJob();
void SingleStep();
void ComputeForces();
void BuildNebrList();
void Sort(real *a, int *seq, int n);
#endif/*_SIMULATION_H_*/
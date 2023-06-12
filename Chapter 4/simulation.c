#include "simulation.h"

// Function definitions
void ErrExit (int code)
{
    printf ("Error: %s\n", errorMsg[code]);
    exit (0);
}

real RandR ()
{
    randSeedP = (randSeedP * IMUL + IADD) & MASK; 
    return (randSeedP * SCALE);
}

void VRand(VecR *p)
{
    real s, x, y;
    s = 2.;

    while (s > 1.) {
        x = 2. * RandR() - 1.;
        y = 2. * RandR() - 1.;
        s = Sqr(x) + Sqr(y);
    }

    p->z = 1. - 2. * s;
    s = 2. * sqrt(1. - s);
    p->x = s * x;
    p->y = s * y;
}

int GetNameList(int argc, char **argv)
{
    int j, k, match, ok;
    char buff[80], *token, currentDir[MAX_PATH];
    FILE *fp;
    
    strcpy (buff, argv[0]);
    strcat (buff, ".in");
    
    if (GetCurrentDirectory(MAX_PATH, currentDir) != 0) {
        strcat(currentDir, "\\simulation.in");
    } else {
        perror("GetCurrentDirectory() error");
        return 1;
    }
    
    if ((fp = fopen(currentDir, "r")) == NULL){
        printf("Failed to open file: %s\n", buff);
        return (0);
    }
    
    for (k = 0; k < sizeof(nameList) / sizeof(NameList); k++)
        nameList[k].vStatus = 0;
    ok = 1;
    while (1)
    {
        fgets(buff, 80, fp);
        if (feof(fp))
            break;
        token = strtok(buff, " \t\n");
        if (!token)
            break;
        match = 0;
        for (k = 0; k < sizeof(nameList) / sizeof(NameList); k++)
        {
            if (strcmp(token, nameList[k].vName) == 0)
            {
                match = 1;
                if (nameList[k].vStatus == 0)
                {
                    nameList[k].vStatus = 1;
                    for (j = 0; j < nameList[k].vLen; j++)
                    {
                        token = strtok(NULL, ", \t\n");
                        if (token)
                        {
                            switch (nameList[k].vType)
                            {
                            case N_I:
                                *NP_I = atol(token);
                                break;
                            case N_R:
                                *NP_R = atof(token);
                                break;
                            }
                        }
                        else
                        {
                            nameList[k].vStatus = 2;
                            ok = 0;
                        }
                    }
                    token = strtok(NULL, ", \t\n");
                    if (token)
                    {
                        nameList[k].vStatus = 3;
                        ok = 0;
                    }
                    break;
                }
                else
                {
                    nameList[k].vStatus = 4;
                    ok = 0;
                }
            }
        }
        if (!match)
            ok = 0;
    }
    fclose(fp);
    for (k = 0; k < sizeof(nameList) / sizeof(NameList); k++)
    {
        if (nameList[k].vStatus != 1)
            ok = 0;
    }
    return (ok);
}


void PrintNameList(FILE *fp)
{
    int j, k;
    fprintf(fp, "NameList -- data\n");
    for (k = 0; k < sizeof(nameList) / sizeof(NameList); k++)
    {
        fprintf(fp, "%s\t", nameList[k].vName);
        if (strlen(nameList[k].vName) < 8)
            fprintf(fp, "\t");
        if (nameList[k].vStatus > 0)
        {
            for (j = 0; j < nameList[k].vLen; j++)
            {
                switch (nameList[k].vType)
                {
                case N_I:
                    fprintf(fp, "%d ", *NP_I);
                    break;
                case N_R:
                    fprintf(fp, "%#g ", *NP_R);
                    break;
                }
            }
        }
        switch (nameList[k].vStatus)
        {
        case 0:
            fprintf(fp, "** no data");
            break;
        case 1:
            break;
        case 2:
            fprintf(fp, "** missing data");
            break;
        case 3:
            fprintf(fp, "** extra data");
            break;
        case 4:
            fprintf(fp, "** multiply defined");
            break;
        }
        fprintf(fp, "\n");
    }
    fprintf(fp, "----\n");
}


void PrintVelDist(FILE *fp)
{
    real vBin;
    int n;
    
    printf("vdist (%.3f)\n", timeNow);
    
    for (n = 0; n < sizeHistVel; n++)
    {
        vBin = (n + 0.5) * rangeVel / sizeHistVel;
        fprintf(fp, "%8.3f %8.3f\n", vBin, histVel[n]);
    }
}

void EvalVelDist()
{
    real deltaV, histSum;
    int j, n;

    if (countVel == 0)
    {
        for (j = 0; j < sizeHistVel; j++)
            histVel[j] = 0.;
    }

    deltaV = rangeVel / sizeHistVel;

    DO_MOL
    {
        j = VLen(mol[n].rv) / deltaV;
        ++histVel[Min(j, sizeHistVel - 1)];
    }

    ++countVel;

    if (countVel == limitVel)
    {
        histSum = 0.;

        for (j = 0; j < sizeHistVel; j++)
            histSum += histVel[j];

        for (j = 0; j < sizeHistVel; j++)
            histVel[j] /= histSum;

        PrintVelDist(stdout);

        countVel = 0;
    }
}

void SetParams()
{
    // Calculate the cutoff distance for the potential function
    rCut = pow(2.0, 1.0 / 6.0);
    // Copy the values of 1.0 / sqrt(density) into the region array
    VSCopy(region, 1.0 / sqrt(density), initUcell);

    // Calculate the total number of molecules in the system
    nMol = VProd(initUcell);

    VSCopy (cells, 1. / rCut, region);

    // Calculate the magnitude of the velocity
    velMag = sqrt(NDIM * (1.0 - 1.0 / nMol) * temperature);

    VSCopy (cells, 1. / (rCut + rNebrShell), region);
    nebrTabMax = nebrTabFac * nMol;
    kinEnInitSum = 0;
}

void LeapfrogStep(int part)
{
    int n;
    if (part == 1) {
        DO_MOL {
            VVSAdd(mol[n].rv, 0.5 * deltaT, mol[n].ra);
            VVSAdd(mol[n].r, deltaT, mol[n].rv);
        }
    } else {
        DO_MOL VVSAdd(mol[n].rv, 0.5 * deltaT, mol[n].ra);
    }
}

void ApplyBoundaryCond()
{
    int n;
    DO_MOL VWrapAll(mol[n].r);
}
void EvalProps()
{
    real vvMax;
    real vv;
    int n;
    VZero(vSum);
    vvSum = 0;
    vvMax = 0.;

    DO_MOL {
        VVAdd(vSum, mol[n].rv);
        vv = VLenSq(mol[n].rv);
        vvSum += vv;
        vvMax = Max (vvMax, vv);
    }

    kinEnergy.val = 0.5 * vvSum / nMol;
    totEnergy.val = kinEnergy.val + uSum / nMol;
    pressure.val = density * (vvSum + virSum) / (nMol * NDIM);
    dispHi += sqrt (vvMax) * deltaT;
    if (dispHi > 0.5 * rNebrShell) nebrNow = 1;
}

void AccumProps(int icode)
{
    if (icode == 0) {
        PropZero(totEnergy);
        PropZero(kinEnergy);
        PropZero(pressure);
    } else if (icode == 1) {
        PropAccum(totEnergy);
        PropAccum(kinEnergy);
        PropAccum(pressure);
    } else if (icode == 2) {
        PropAvg(totEnergy, stepAvg);
        PropAvg(kinEnergy, stepAvg);
        PropAvg(pressure, stepAvg);
    }
}

void PrintSummary(FILE *fp)
{
    fprintf(fp,
            "%5d %8.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f\n",
            stepCount, timeNow, VCSum(vSum) / nMol, PropEst(totEnergy),
            PropEst(kinEnergy), PropEst(pressure));
}

void PredictorStep()
{
    real cr[] = {19., -10., 3.};
    real cv[] = {27., -22., 7.};
    real div = 24.;
    real wr, wv;
    int n;

    wr = Sqr(deltaT) / div;
    wv = deltaT / div;

    DO_MOL {
        mol[n].ro = mol[n].r;
        mol[n].rvo = mol[n].rv;

        PR(x);
        PRV(x);
        PR(y);
        PRV(y);
        PR(z);
        PRV(z);

        mol[n].ra2 = mol[n].ra1;
        mol[n].ra1 = mol[n].ra;
    }
}

void CorrectorStep()
{
    real cr[] = {3., 10., -1.};
    real cv[] = {7., 6., -1.};
    real div = 24.;
    real wr, wv;
    int n;

    wr = Sqr(deltaT) / div;
    wv = deltaT / div;

    DO_MOL {
        CR(x);
        CRV(x);
        CR(y);
        CRV(y);
        CR(z);
        CRV(z);
    }
}

void AdjustInitTemp()
{
    real vFac;
    int n;

    kinEnInitSum += kinEnergy.val;

    if (stepCount % stepInitlzTemp == 0) {
        kinEnInitSum /= stepInitlzTemp;
        vFac = velMag / sqrt(2. * kinEnInitSum);

        DO_MOL
        {
            VScale(mol[n].rv, vFac);
        }

        kinEnInitSum = 0.;
    }
}

void SingleStep () {
    ++ stepCount;
    timeNow = stepCount * deltaT;
    // LeapfrogStep (1);
    PredictorStep ();
    ApplyBoundaryCond ();
    if (nebrNow) {
        nebrNow = 0;
        dispHi = 0.;
        BuildNebrList ();
        }
    ComputeForces ();
    // LeapfrogStep (2);
    CorrectorStep ();
    EvalProps ();
    AccumProps (1);
    if (stepCount % stepAvg == 0) {
        PutConfig ();
        AccumProps (2);
        if (stepCount < stepEquil) AdjustInitTemp ();
        PrintSummary (stdout);
        AccumProps (0);
    }
    if (stepCount >= stepEquil &&
    (stepCount - stepEquil) %stepRdf == 0) EvalRdf ();
}

void InitCoords()
{
    VecR c, gap;
    int j, n, nx, ny, nz;

    VDiv(gap, region, initUcell);

    // An obvious way of reducing equilibration time is to base the initial state on the final state of a previous run.
    n = 0;
    for (nz = 0; nz < initUcell.z; nz++) {
        for (ny = 0; ny < initUcell.y; ny++) {
            for (nx = 0; nx < initUcell.x; nx++) {
                VSet(c, nx + 0.25, ny + 0.25, nz + 0.25);
                VMul(c, c, gap);
                
                VVSAdd(c, -0.5, region);
                for (j = 0; j < 4; j++) {
                    mol[n].r = c;
                    if (j != 3) {
                        if (j != 0) mol[n].r.x += 0.5 * gap.x;
                        if (j != 1) mol[n].r.y += 0.5 * gap.y;
                        if (j != 2) mol[n].r.z += 0.5 * gap.z;
                    }
                    ++n;
                }
            }
        }
    }
}


void InitVels()
{
    int n;
    VZero(vSum);

    DO_MOL {
        VRand(&mol[n].rv);
        VScale(mol[n].rv, velMag);
        VVAdd(vSum, mol[n].rv);
    }

    DO_MOL
    {
        VVSAdd(mol[n].rv, -1. / nMol, vSum);
    }
}

void InitAccels()
{
    int n;
    DO_MOL
    {
        VZero(mol[n].ra);
    }
}
void AllocArrays ()
{
        /**
     * @brief Allocates memory for mol array.
     */
    AllocMem (mol, nMol, Mol);
    AllocMem (histVel, sizeHistVel, real);
    AllocMem (nebrTab, 2 * nebrTabMax, int);
    AllocMem (cellList, VProd (cells) + nMol, int);
    AllocMem (histRdf, sizeHistRdf, real);
    //AllocMem (clust, nMol, Clust);
    //nCellEdge = region.x / rClust;
    //AllocMem (cellList, Cube (nCellEdge) + nMol, int);
  
    /**
     * @brief Allocates memory for distSq array.
     */
    AllocMem(distSq, nMol, real);

    /**
     * @brief Allocates memory for siteSeq array.
     */
    AllocMem(siteSeq, nMol, int);

    /**
     * @brief Allocates memory for testSites array.
     */
    AllocMem(testSites, nMol, int);

    /**
     * @brief Allocates memory for cellList array.
     */
    AllocMem(cellList, VProd(cells) + nMol, int);

    /**
     * @brief Allocates memory for edge array.
     */
    AllocMem(edge, MAX_EDGE, Edge);

    /**
     * @brief Allocates memory for face array.
     */
    AllocMem(face, MAX_FACE, Face);

    /**
     * @brief Allocates memory for flist array.
     */
    AllocMem(flist, MAX_FLIST, Flist);

    /**
     * @brief Allocates memory for vert array.
     */
    AllocMem(vert, MAX_VERT, Vert);

    /**
     * @brief Allocates memory for eCut array.
     */
    AllocMem(eCut, MAX_ITEM, int);

    /**
     * @brief Allocates memory for eDel array.
     */
    AllocMem(eDel, MAX_ITEM, int);

    /**
     * @brief Allocates memory for eNew array.
     */
    AllocMem(eNew, MAX_ITEM, int);

    /**
     * @brief Allocates memory for fCut array.
     */
    AllocMem(fCut, MAX_ITEM, int);

    /**
     * @brief Allocates memory for fDel array.
     */
    AllocMem(fDel, MAX_ITEM, int);

    /**
     * @brief Allocates memory for vDel array.
     */
    AllocMem(vDel, MAX_ITEM, int);
}

void SetupJob () {
    AllocArrays ();
    stepCount = 0;
    countVel = 0;
    countRdf = 0;
    InitCoords ();
    InitVels ();
    InitAccels ();
    AccumProps (0);
    nebrNow = 1;
}

void ComputeForces()
{
    VecR dr;
    real fcVal, rr, rrCut, rri, rri3;
    int j1, j2, n;
    
    rrCut = Sqr(rCut);
    DO_MOL VZero(mol[n].ra);
    uSum = 0.;
    virSum = 0.;
    
    for (j1 = 0; j1 < nMol - 1; j1++)
    {
        for (j2 = j1 + 1; j2 < nMol; j2++)
        {
            VSub(dr, mol[j1].r, mol[j2].r);
            VWrapAll(dr);
            rr = VLenSq(dr);
            
            if (rr < rrCut)
            {
                rri = 1. / rr;
                rri3 = Cube(rri);
                fcVal = 48. * rri3 * (rri3 - 0.5) * rri;
                VVSAdd(mol[j1].ra, fcVal, dr);
                VVSAdd(mol[j2].ra, -fcVal, dr);
                uSum += 4. * rri3 * (rri3 - 1.) + 1.;
                virSum += fcVal * rr;
            }
        }
    }
}
// Predictor and Corrector Methods


// Function to create necessary files and set up filenames for checkpointing and snapshot records
void SetupFiles()
{
    FILE *fp;
    int k;

    // Initialize filenames for checkpointing and snapshot records
    for (k = 0; k < sizeof(fileNameR); k++)
    {
        strcpy(fileName[k], fileNameR[k]);
        fileName[k][0] = progId[0];
        fileName[k][1] = progId[1];
        fileName[k][2] = runId / 10 + CHAR_ZERO;
        fileName[k][3] = runId % 10 + CHAR_ZERO;
    }

    // Handle checkpointing and snapshot files

    if (!doCheckpoint)
    {
        newRun = 1;
    }
    else if ((fp = fopen(fileName[FL_CKLAST], "r")) != 0)
    {
        newRun = 0;
        fclose(fp);
    }
    else
    {
        newRun = 1;
        fp = fopen(fileName[FL_CHECKA], "w");
        fclose(fp);
        fp = fopen(fileName[FL_CHECKB], "w");
        fclose(fp);
        fp = fopen(fileName[FL_CKLAST], "w");
        fputc(CHAR_ZERO + FL_CHECKA, fp);
        fclose(fp);
    }

    // Create snapshot file if it's a new run
    if (newRun && recordSnap)
    {
        fp = fopen(fileName[FL_SNAP], "w");
        fclose(fp);
    }
}

// Function to write checkpoint data to disk
void PutCheckpoint()
{
    int fOk, fVal;
    FILE *fp;

    fOk = 0;
    if ((fp = fopen(fileName[FL_CKLAST], "r+")) != 0)
    {
        fVal = FL_CHECKA + FL_CHECKB - (fgetc(fp) - CHAR_ZERO);
        rewind(fp);
        fputc(CHAR_ZERO + fVal, fp);
        fclose(fp);
        fOk = 1;
    }

    if (fOk && (fp = fopen(fileName[fVal], "w")) != 0)
    {
        WriteF(kinEnergy);
        WriteF(stepCount);
        WriteF(timeNow);
        WriteF(totEnergy);
        WriteFN(mol, nMol);
        if (ferror(fp))
            fOk = 0;
        fclose(fp);
    }
    else
        fOk = 0;

    if (!fOk)
        ErrExit(ERR_CHECKPT_WRITE);
}
void BuildNebrList()
{
    VecR dr, invWid, rs, shift;
    VecI cc, m1v, m2v, vOff[] = OFFSET_VALS;
    real rrNebr;
    int c, j1, j2, m1, m1x, m1y, m1z, m2, n, offset;

    rrNebr = Sqr(rCut + rNebrShell);
    VDiv(invWid, cells, region);

    for (n = nMol; n < nMol + VProd(cells); n++)
        cellList[n] = -1;

    DO_MOL {
        VSAdd(rs, mol[n].r, 0.5, region);
        VMul(cc, rs, invWid);
        c = VLinear(cc, cells) + nMol;
        cellList[n] = cellList[c];
        cellList[c] = n;
    }

    nebrTabLen = 0;

    for (m1z = 0; m1z < cells.z; m1z++) {
        for (m1y = 0; m1y < cells.y; m1y++) {
            for (m1x = 0; m1x < cells.x; m1x++) {
                VSet(m1v, m1x, m1y, m1z);
                m1 = VLinear(m1v, cells) + nMol;

                for (offset = 0; offset < N_OFFSET; offset++) {
                    VAdd(m2v, m1v, vOff[offset]);
                    VZero(shift);
                    VCellWrapAll();
                    m2 = VLinear(m2v, cells) + nMol;

                    DO_CELL(j1, m1) {
                        DO_CELL(j2, m2) {
                            if (m1 != m2 || j2 < j1) {
                                VSub(dr, mol[j1].r, mol[j2].r);
                                VVSub(dr, shift);

                                if (VLenSq(dr) < rrNebr) {
                                    if (nebrTabLen >= nebrTabMax)
                                        ErrExit(ERR_TOO_MANY_NEBRS);

                                    nebrTab[2 * nebrTabLen] = j1;
                                    nebrTab[2 * nebrTabLen + 1] = j2;
                                    ++nebrTabLen;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

// Function to read checkpoint data from disk
void GetCheckpoint()
{
    int fOk, fVal;
    FILE *fp;

    fOk = 0;
    if ((fp = fopen(fileName[FL_CKLAST], "r")) != 0)
    {
        fVal = fgetc(fp) - CHAR_ZERO;
        fclose(fp);
        fOk = 1;
    }

    if (fOk && (fp = fopen(fileName[fVal], "r")) != 0)
    {
        ReadF(kinEnergy);
        ReadF(stepCount);
        ReadF(timeNow);
        ReadF(totEnergy);
        ReadFN(mol, nMol);
        if (ferror(fp))
            fOk = 0;
        fclose(fp);
    }
    else
        fOk = 0;

    if (!fOk)
        ErrExit(ERR_CHECKPT_READ);
}

// PovRay Generator
POVRayWriter* createPOVRayWriter(const char* fileName) {
    POVRayWriter* writer = (POVRayWriter*)malloc(sizeof(POVRayWriter));
    writer->fileName = (char*)malloc(MAX_STRING_LENGTH * sizeof(char));
    strcpy(writer->fileName, fileName);
    writer->file = fopen(fileName, "w");
    return writer;
}

void destroyPOVRayWriter(POVRayWriter* writer) {
    if (writer != NULL) {
        if (writer->file != NULL) {
            fclose(writer->file);
        }
        free(writer->fileName);
        free(writer);
    }
}

void writeHeader(POVRayWriter* writer) {
    fprintf(writer->file, "#include \"colors.inc\"\n");
    fprintf(writer->file, "#include \"textures.inc\"\n");
    fprintf(writer->file, "#include \"glass.inc\"\n");
    fprintf(writer->file, "#include \"metals.inc\"\n");
    fprintf(writer->file, "#include \"woods.inc\"\n\n");
}

void writeCamera(POVRayWriter* writer, Vector3D position, Vector3D lookAt, Vector3D up) {
    fprintf(writer->file, "camera {\n");
    fprintf(writer->file, "  location <%f, %f, %f>\n", position.x, position.y, position.z);
    fprintf(writer->file, "  look_at <%f, %f, %f>\n", lookAt.x, lookAt.y, lookAt.z);
    fprintf(writer->file, "  up <%f, %f, %f>\n", up.x, up.y, up.z);
    fprintf(writer->file, "}\n\n");
}

void writeSphere(POVRayWriter* writer, Sphere sphere) {
    fprintf(writer->file, "sphere {\n");
    fprintf(writer->file, "  <%d, %d, %d>, %f\n", sphere.position.x, sphere.position.y, sphere.position.z, sphere.radius);
    fprintf(writer->file, "  texture {\n");
    fprintf(writer->file, "    pigment { %s }\n", sphere.finish.pigmentColor);
    fprintf(writer->file, "    finish {\n");
    fprintf(writer->file, "      ambient %f\n", sphere.finish.ambient);
    fprintf(writer->file, "      diffuse %f\n", sphere.finish.diffuse);
    fprintf(writer->file, "      specular %f\n", sphere.finish.specular);
    fprintf(writer->file, "      reflection %f\n", sphere.finish.reflection);
    fprintf(writer->file, "      refraction %f\n", sphere.finish.refraction);
    fprintf(writer->file, "    }\n");
    fprintf(writer->file, "  }\n");
    fprintf(writer->file, "}\n\n");
}

int init_pov(char * filename) {
    printf("initializing pov scene...\n");
    writer = createPOVRayWriter(filename);
    writeHeader(writer);

    Vector3D cameraPosition = { 0, 0, -20 };
    Vector3D lookAt = { 0, 0, 0 };
    Vector3D up = { 0, 1, 0 };
    writeCamera(writer, cameraPosition, lookAt, up);
    printf("Done.\n");
    return 0;
}

int write_pov(double x, double y, double z){
    
    Vector3D position = { x, y, z};
    double radius = 0.4;

    Sphere sphere = { position, radius, { "Gray", 0.1, 0.6, 0.4, 0.5, 0.0 } };
    writeSphere(writer, sphere);
    return 0;
}

void runPovEngine(const char* sceneFile) {
    printf("Rendering pov scene...\n");
    char command[100], currentDir[MAX_PATH];
    if (GetCurrentDirectory(MAX_PATH, currentDir) != 0) {
        snprintf(command, sizeof(command), "cd \"%s\" && pvengine %s", currentDir, sceneFile);
    } else {
        perror("Get Directory error");
    }
    system(command);
}

void scaleData(int* x, int* y, int* z) {
    // Define the scaling factors for each axis
    double scaleX = 0.00000001;
    double scaleY = 0.00000001;
    double scaleZ = 0.00000001;

    // Scale the coordinates using the scaling factors
    *x *= scaleX;
    *y *= scaleY;
    *z *= scaleZ;
}

//Chapter four Methods

/**

*@brief Calculate and evaluate radial distribution function (RDF).

*This function calculates and evaluates the radial distribution function (RDF) for a system of particles.

*The RDF quantifies the probability of finding a particle at a given distance from a reference particle.

*The RDF is calculated based on the positions of the particles in the system.

*@note The RDF calculation assumes the system is in equilibrium.
*/
void EvalRdf()
{
    VecR dr;
    real deltaR, normFac, rr;
    int j1, j2, n;

    if (countRdf == 0) {
        for (n = 0; n < sizeHistRdf; n++) {
            histRdf[n] = 0.;
        }
    }

    deltaR = rangeRdf / sizeHistRdf;
    for (j1 = 0; j1 < nMol - 1; j1++) {
        for (j2 = j1 + 1; j2 < nMol; j2++) {
            VSub(dr, mol[j1].r, mol[j2].r);
            VWrapAll(dr);
            rr = VLenSq(dr);

            if (rr < Sqr(rangeRdf)) {
                n = sqrt(rr) / deltaR;
                ++histRdf[n];
            }
        }
    }

    ++countRdf;

    if (countRdf == limitRdf) {
        normFac = VProd(region) / (2. * M_PI * Cube(deltaR) * Sqr(nMol) * countRdf);

        for (n = 0; n < sizeHistRdf; n++) {
            histRdf[n] *= normFac / Sqr(n - 0.5);
        }

        PrintRdf(stdout);
        countRdf = 0;
    }
}

/**
 * @brief Prints the radial distribution function (RDF) to a file.
 *
 * This function prints the RDF values to the specified file pointer.
 *
 * @param fp The file pointer to which the RDF will be printed.
 */
void PrintRdf(FILE* fp)
{
    real rb; // Radial bin position
    int n;

    fprintf(fp, "rdf\n"); // Print a header indicating RDF

    for (n = 0; n < sizeHistRdf; n++) {
        rb = (n + 0.5) * rangeRdf / sizeHistRdf; // Calculate the radial position
        fprintf(fp, "%8.4f %8.4f\n", rb, histRdf[n]); // Print the radial position and RDF value
    }
}

/**
 * \brief Evaluate lattice correlation.
 *
 * This function calculates the lattice correlation based on the given molecular positions.
 * It computes the real and imaginary parts of the correlation using the wavevector and molecular positions.
 * The final lattice correlation value is the square root of the sum of squares of the real and imaginary parts,
 * divided by the total number of molecules.
 */
void EvalLatticeCorr()
{
    VecR kVec;
    real si, sr, t;
    int n;

    // Calculate the wavevector components
    kVec.x = 2. * M_PI * initUcell.x / region.x;
    kVec.y = -kVec.x;
    kVec.z = kVec.x;

    sr = 0.;
    si = 0.;

    // Iterate over all molecules
    DO_MOL {
        // Calculate the dot product between wavevector and molecular position
        t = VDot(kVec, mol[n].r);

        // Accumulate real and imaginary parts of the correlation
        sr += cos(t);
        si += sin(t);
    }

    // Calculate the lattice correlation as the square root of the sum of squares of real and imaginary parts,
    // divided by the total number of molecules
    latticeCorr = sqrt(Sqr(sr) + Sqr(si)) / nMol;
}

/**
 * @brief Analyzes Voronoi polygons.
 *
 * This function analyzes Voronoi polygons by performing various operations
 * such as sorting, initialization, bisecting planes, processing deleted
 * vertices, processing cut edges, processing cut faces, processing new vertices,
 * removing old polygons, finding distance vertices, and calculating polygon geometry and size.
 */
void AnalVorPoly()
{
    int nf;

    // Sort the squared distances and site sequence of test sites
    Sort(distSq, siteSeq, nTestSites);

    // Initialize Voronoi polygons
    InitVorPoly();

    // Iterate over the test sites
    for (curSite = 0; curSite < nTestSites; curSite++) {
        if (distSq[siteSeq[curSite]] >= 4. * vDistSqMax)
            break;

        siteB = testSites[siteSeq[curSite]];

        nvDel = 0;
        neNew = 0;
        neDel = 0;
        neCut = 0;
        nfDel = 0;
        nfCut = 0;

        // Bisect the plane
        BisectPlane();

        // Process deleted vertices
        if (nvDel > 0)
            ProcDelVerts();

        // Process cut edges
        if (neCut > 0)
            ProcCutEdges();

        // Process cut faces
        if (nfCut > 0)
            ProcCutFaces();

        // Process new vertices
        if (neNew > 0)
            ProcNewVerts();

        // Process new face
        if (nfCut > 0)
            ProcNewFace();

        // Remove old polygons
        RemoveOld();

        // Find distance vertices
        if (nfCut > 0)
            FindDistVerts();
    }

    // Check if all faces have been properly subdivided
    for (nf = 0; nf < 4; nf++)
        if (face[nf].stat != 0)
            ErrExit(ERR_SUBDIV_UNFIN);

    // Calculate polygon geometry
    PolyGeometry();

    // Calculate polygon size
    PolySize();
}
/**
 * \brief Find test sites around a given atom.
 *
 * \param na Index of the atom.
 */
void FindTestSites(int na)
{
    VecR dr;
    VecI cn;
    int c, cx, cy, cz, i, ofx, ofy, ofz;

    // Calculate the cell indices for the given atom
    cx = mol[na].inCell % cells.x;
    cy = (mol[na].inCell / cells.x) % cells.y;
    cz = mol[na].inCell / (cells.x * cells.y);

    // Reset the number of test sites
    nTestSites = 0;

    // Iterate over neighboring cells in a 3x3x3 cube
    for (ofz = -1; ofz <= 1; ofz++) {
        // Calculate the current z-index of the neighboring cell
        cn.z = (cz + ofz + cells.z) % cells.z;

        for (ofy = -1; ofy <= 1; ofy++) {
            // Calculate the current y-index of the neighboring cell
            cn.y = (cy + ofy + cells.y) % cells.y;

            for (ofx = -1; ofx <= 1; ofx++) {
                // Calculate the current x-index of the neighboring cell
                cn.x = (cx + ofx + cells.x) % cells.x;

                // Calculate the linear cell index
                c = VLinear(cn, cells) + nMol;

                // Iterate over atoms in the current cell
                DO_CELL(i, c) {
                    // Calculate the displacement vector between the atom and the current atom in the cell
                    VSub(dr, mol[na].r, mol[i].r);
                    // Wrap the displacement vector across periodic boundaries
                    VWrapAll(dr);

                    // Store the index and squared distance of the test site
                    testSites[nTestSites] = i;
                    distSq[nTestSites] = VLenSq(dr);
                    ++nTestSites;
                }
            }
        }
    }
}

/**
 * @brief Initialize Voronoi polyhedron.
 */
void InitVorPoly()
{
    VecR w, vPosI[] = {{-1., -1., -1.}, {1., -1., -1.}, {0., 2., -1.}, {0., 0., 3.}};
    real r2, r6;
    int m, n, ne, nf, nv, s,
        vValI[] = {0,2,5,0,1,4,1,2,3,3,4,5},
        eFacesI[] = {0,3,0,1,0,2,1,2,1,3,2,3},
        eVertsI[] = {0,1,1,2,0,2,2,3,1,3,0,3},
        eI[] = {0,1,2,4,3,1,2,3,5,5,4,0},
        vI[] = {0,1,2,1,3,2,0,2,3,0,3,1};

    r2 = sqrt(2.) * rangeLim;
    r6 = sqrt(6.) * rangeLim;
    siteA = testSites[siteSeq[0]];
    eLast = 5;
    fLast = 3;
    vLast = 3;
    m = 0;

    for (nv = 0; nv <= vLast; nv++) {
        vert[nv].pos = mol[siteA].r;
        VSet(w, r6 / 3., r2 / 3., rangeLim / 3.);
        VMul(w, w, vPosI[nv]);
        VVAdd(vert[nv].pos, w);
        vert[nv].distSq = Sqr(rangeLim);
        vert[nv].stat = 2;

        for (n = 0; n < 3; n++) {
            vert[nv].e[n] = vValI[m];
            ++m;
        }
    }

    vDistSqMax = vert[0].distSq;

    for (ne = 0; ne <= eLast; ne++) {
        edge[ne].v[0] = eVertsI[2 * ne];
        edge[ne].f[0] = eFacesI[2 * ne];
        edge[ne].v[1] = eVertsI[2 * ne + 1];
        edge[ne].f[1] = eFacesI[2 * ne + 1];
        edge[ne].stat = 3;
    }

    for (s = 0; s < MAX_FLIST - 1; s++) {
        flist[s].link = s + 1;
    }

    s = 0;

    for (nf = 0; nf <= fLast; nf++) {
        face[nf].vFar = vI[s];
        face[nf].stat = 3;
        face[nf].fPtr = s;

        for (n = 0; n < 3; n++) {
            flist[s].v = vI[s];
            flist[s].e = eI[s];
            ++s;
        }

        flist[s - 1].link = face[nf].fPtr;
    }

    fListLast = s - 1;
}

/**
 * \brief Bisects a plane defined by two points in 3D space.
 */
void BisectPlane()
{
    VecR dr, shift;
    real d1, d2, d3;
    int nv;

    d1 = 0.;
    fParamS = 0.;

    // Calculate the vector between siteB and siteA
    VSub(fParamV, mol[siteB].r, mol[siteA].r);
    VZero(shift);
    VShiftAll(fParamV);
    VVAdd(fParamV, shift);

    // Calculate dot product and add shift to the vector
    d1 = VDot(fParamV, mol[siteA].r);
    VAdd(dr, mol[siteB].r, shift);

    // Calculate the value of fParamS
    fParamS = 0.5 * (VLenSq(dr) - VLenSq(mol[siteA].r));

    // Iterate over vertices and check if they lie on the plane
    for (nv = 0; nv <= vLast; nv++) {
        if (vert[nv].stat != 0) {
            d2 = VDot(fParamV, vert[nv].pos);

            // Check if the vertex is on the plane
            if (d1 != d2) {
                d3 = (fParamS - d1) / (d2 - d1);

                // Check if the vertex is between siteA and siteB
                if (d3 > 0. && d3 < 1.) {
                    vDel[nvDel] = nv;
                    ++nvDel;
                    vert[nv].stat = 1;
                }
            }
        }
    }
}

/**
 * @brief Process deleted vertices.
 */
void ProcDelVerts()
{
    int e, m, n, nv;

    // Iterate over the deleted vertices.
    for (nv = 0; nv < nvDel; nv++)
    {
        // Iterate over the edges connected to the vertex.
        for (m = 0; m < 3; m++)
        {
            e = vert[vDel[nv]].e[m];
            --edge[e].stat;

            // Check the status of the edge.
            if (edge[e].stat == 2)
            {
                eCut[neCut] = e;
                ++neCut;
            }
            else
            {
                eDel[neDel] = e;
                ++neDel;
            }

            // Iterate over the faces connected to the edge.
            for (n = 0; n < 2; n++)
            {
                // Check the status of the face.
                if (face[edge[e].f[n]].stat == 3)
                {
                    fCut[nfCut] = edge[e].f[n];
                    ++nfCut;
                    face[edge[e].f[n]].stat = 2;
                }
            }
        }
    }
}

/**
 * @brief Process the cut edges.
 */
void ProcCutEdges()
{
    VecR dr;
    real d, dt1, dt2;
    int nd, ne, vt1, vt2;

    for (ne = 0; ne < neCut; ne++) {
        // Check if the edge status is 2
        if (edge[eCut[ne]].stat == 2) {
            // Update the edge status to 3
            edge[eCut[ne]].stat = 3;
            vt1 = edge[eCut[ne]].v[0]; // Get the first vertex of the edge
            vt2 = edge[eCut[ne]].v[1]; // Get the second vertex of the edge

            dt1 = VDot(fParamV, vert[vt1].pos); // Calculate dot product 1
            dt2 = VDot(fParamV, vert[vt2].pos); // Calculate dot product 2

            // Determine the neighbor vertex index based on status
            if (vert[vt1].stat == 1)
                nd = 0;
            else if (vert[vt2].stat == 1)
                nd = 1;

            ++vLast;
            vert[vLast].stat = 2;
            vert[vLast].distSq = 0.;

            // Interpolate the position of the new vertex
            d = (fParamS - dt1) / (dt2 - dt1);
            VInterp(vert[vLast].pos, d, vert[vt2].pos, vert[vt1].pos);

            VSub(dr, vert[vLast].pos, mol[siteA].r);
            vert[vLast].distSq = VLenSq(dr);

            edge[eCut[ne]].v[nd] = vLast;
            vert[vLast].e[0] = eCut[ne];
            vert[vLast].e[1] = 0;
            vert[vLast].e[2] = 0;
        }
    }
}

/**
 * \brief Process cut faces.
 *
 * This function processes cut faces by updating their status and pointers.
 */
void ProcCutFaces()
{
    int faceGone, nf, s, s1, s2, s3, s4, v1, v2, vDelCount;

    eLastP = eLast;
    ++fLast;

    /**
     * \par Iterate over the cut faces
     */
    for (nf = 0; nf < nfCut; nf++) {

        s = face[fCut[nf]].fPtr;
        faceGone = 0;

        /**
         * \par Traverse face pointers until a valid vertex or the initial face is reached
         */
        while (vert[flist[s].v].stat != 2 && !faceGone) {
            s = flist[s].link;
            if (s == face[fCut[nf]].fPtr)
                faceGone = 1;
        }

        /**
         * \par Update face status and pointers
         */
        if (faceGone) {
            fDel[nfDel] = fCut[nf];
            face[fCut[nf]].stat = 1;
            ++nfDel;
        } else {
            face[fCut[nf]].stat = 3;
            face[fCut[nf]].fPtr = s;

            s1 = s;
            s2 = flist[s1].link;

            /**
             * \par Find the last valid vertex and count deleted vertices
             */
            for (; vert[flist[s2].v].stat == 2; s2 = flist[s1].link)
                s1 = s2;

            vDelCount = 1;

            for (s3 = s2, s4 = flist[s3].link; vert[flist[s4].v].stat != 2; s4 = flist[s3].link) {
                ++vDelCount;
                s3 = s4;
            }

            v1 = edge[flist[s1].e].v[0] + edge[flist[s1].e].v[1] - flist[s1].v;
            v2 = edge[flist[s3].e].v[0] + edge[flist[s3].e].v[1] - flist[s4].v;

            ++eLast;
            flist[s3].v = v2;

            /**
             * \par Handle the case when a single vertex is deleted
             */
            if (vDelCount == 1) {
                ++fListLast;
                s = fListLast;
                flist[s1].link = s;
                flist[s].link = s2;
                flist[s].v = v1;
                flist[s].e = eLast;
            } else {
                flist[s2].v = v1;
                flist[s2].e = eLast;
                if (vDelCount > 2)
                    flist[s2].link = s3;
            }

            edge[eLast].v[0] = v1;
            edge[eLast].v[1] = v2;
            edge[eLast].f[0] = fCut[nf];
            edge[eLast].f[1] = fLast;
            edge[eLast].stat = 2;

            eNew[neNew] = eLast;
            ++neNew;
        }
    }
}
/**
 * \brief Process new vertices.
 */
void ProcNewVerts()
{
    int ne, v;

    for (ne = 0; ne < neNew; ne++)
    {
        if (eNew[ne] > eLastP)
        {
            v = edge[eNew[ne]].v[0];

            if (vert[v].e[1] == 0)
                vert[v].e[1] = eNew[ne];
            else
                vert[v].e[2] = eNew[ne];

            v = edge[eNew[ne]].v[1];

            if (vert[v].e[1] == 0)
                vert[v].e[1] = eNew[ne];
            else
                vert[v].e[2] = eNew[ne];
        }
    }
}

/**
 * \brief Process new face.
 */
void ProcNewFace()
{
    int e, n, ne, v;

    for (n = 0; n < neNew; n++)
    {
        ++fListLast;

        if (n == 0)
        {
            e = eNew[0];
            face[fLast].fPtr = fListLast;
            v = edge[e].v[0];
        }
        else
        {
            ne = 1;

            for (e = eNew[ne]; edge[e].v[0] != v && edge[e].v[1] != v || edge[e].stat == 3; e = eNew[ne])
                ++ne;
        }

        flist[fListLast].v = v;
        v = edge[e].v[0] + edge[e].v[1] - v;
        flist[fListLast].e = e;
        edge[e].stat = 3;
    }

    face[fLast].stat = 3;
    flist[fListLast].link = face[fLast].fPtr;
    face[fLast].dist = 0.5 * sqrt(distSq[siteSeq[curSite]]);
}
/**
 * @brief Remove old vertices, edges, and faces.
 *
 * This function removes old vertices, edges, and faces by setting their statuses to 0.
 */
void RemoveOld()
{
    int n;

    // Remove old vertices
    for (n = 0; n < nvDel; n++) {
        /**
         * @note The status of the vertex is set to 0, indicating that it is removed.
         */
        vert[vDel[n]].stat = 0;
    }

    // Remove old edges
    for (n = 0; n < neDel; n++) {
        if (edge[eDel[n]].stat == 1) {
            /**
             * @note The status of the edge is set to 0, indicating that it is removed.
             */
            edge[eDel[n]].stat = 0;
        }
    }

    // Remove old faces
    for (n = 0; n < nfDel; n++) {
        /**
         * @note The status of the face is set to 0, indicating that it is removed.
         */
        face[fDel[n]].stat = 0;
    }
}

/**
 * @brief Find the distance between vertices.
 */
void FindDistVerts()
{
    real dd; /**< Distance variable */
    int nf, s; /**< Loop variables */

    fCut[nfCut] = fLast; /**< Set the last face in the cut array */

    /**
     * Iterate over the cut faces.
     */
    for (nf = 0; nf < nfCut + 1; nf++) {
        /**
         * Check if the face is not in a specific state.
         */
        if (face[fCut[nf]].stat != 0) {
            s = face[fCut[nf]].fPtr;
            dd = vert[flist[s].v].distSq; /**< Compute the squared distance */
            face[fCut[nf]].vFar = flist[s].v; /**< Set the far vertex of the face */

            /**
             * Iterate over the linked faces to find the farthest vertex.
             */
            for (s = flist[s].link; s != face[fCut[nf]].fPtr; s = flist[s].link) {
                /**
                 * Check if the current vertex distance is greater than the previous farthest distance.
                 * If true, update the farthest vertex and the distance.
                 */
                if (vert[flist[s].v].distSq > dd) {
                    dd = vert[flist[s].v].distSq;
                    face[fCut[nf]].vFar = flist[s].v;
                }
            }
        }
    }

    vDistSqMax = 0.; /**< Initialize the maximum vertex distance squared */

    /**
     * Iterate over all faces and update the maximum vertex distance squared.
     */
    for (nf = 0; nf <= fLast; nf++) {
        /**
         * Check if the face is not in a specific state and if the vertex distance squared is greater than the current maximum.
         * If true, update the maximum vertex distance squared.
         */
        if (face[nf].stat != 0 && vDistSqMax < vert[face[nf].vFar].distSq) {
            vDistSqMax = vert[face[nf].vFar].distSq;
        }
    }
}

/**
 * @brief Calculate polygon geometry values.
 */
void PolyGeometry()
{
    int n, ne, nf, nv, s;

    // Initialize polygon geometry values to zero.
    for (n = 0; n < 4; n++)
    {
        polyGeom[n].val = 0.;
    }

    // Count the number of vertices with non-zero status.
    for (nv = 0; nv <= vLast; nv++)
    {
        if (vert[nv].stat != 0)
        {
            ++polyGeom[0].val;
        }
    }

    // Count the number of edges with non-zero status.
    for (ne = 0; ne <= eLast; ne++)
    {
        if (edge[ne].stat != 0)
        {
            ++polyGeom[1].val;
        }
    }

    // Count the number of faces and the number of face links with non-zero status.
    for (nf = 0; nf <= fLast; nf++)
    {
        if (face[nf].stat != 0)
        {
            ++polyGeom[2].val;
            ++polyGeom[3].val;

            // Count additional face links.
            for (s = flist[face[nf].fPtr].link; s != face[nf].fPtr; s = flist[s].link)
            {
                ++polyGeom[3].val;
            }
        }
    }

    // Calculate the average number of face links per face.
    polyGeom[3].val /= polyGeom[2].val;
}

/**
 * \brief Calculate the size of a polygon.
 *
 * This function calculates the area and volume of a polygon.
 */
void PolySize()
{
    VecR ca, d1, d2, d3;
    real a;
    int nf, s, v1, v2;

    // Initialize area and volume to zero
    polyArea.val = 0.;
    polyVol.val = 0.;

    // Iterate over each face of the polygon
    for (nf = 0; nf <= fLast; nf++) {
        // Check if the face is active
        if (face[nf].stat != 0) {
            // Get the first vertex of the face
            s = face[nf].fPtr;
            v1 = flist[s].v;

            // Get the second vertex of the face
            s = flist[s].link;
            v2 = flist[s].v;

            // Calculate the first direction vector
            VSub(d1, vert[v2].pos, vert[v1].pos);

            // Initialize the cross product accumulator
            VZero(ca);

            // Calculate the cross product and accumulate
            for (s = flist[s].link; s != face[nf].fPtr; s = flist[s].link) {
                // Get the next vertex
                v2 = flist[s].v;

                // Calculate the direction vector
                VSub(d2, vert[v2].pos, vert[v1].pos);

                // Calculate the cross product
                VCross(d3, d1, d2);

                // Accumulate the cross product
                VVAdd(ca, d3);

                // Update the first direction vector
                d1 = d2;
            }

            // Calculate the area of the face and add it to the total area
            a = VLen(ca);
            polyArea.val += a / 2.;

            // Calculate the volume of the face and add it to the total volume
            polyVol.val += face[nf].dist * a / 6.;
        }
    }
}

/**
 * @brief Set the size of the simulation cells based on the cell ratio and region size.
 *
 * This function sets the size of the simulation cells by copying the cell ratio to the cells array
 * and calculating the total number of cells (nCell) by taking the product of the cell dimensions.
 */
void SetCellSize()
{
    // Copy the cell ratio to the cells array
    VSCopy(cells, cellRatio, region);

    // Calculate the total number of cells
    nCell = VProd(cells);
}
/**
 * \brief Build clusters based on a distance threshold.
 */
// void BuildClusters()
// {
//     // Calculate the squared distance threshold for cluster formation.
//     real rrClust;
//     rrClust = Sqr(rClust);

//     // Iterate over pairs of particles and add bonded pairs within the distance threshold.
//     // TODO: Add details about the specific implementation.
//     if (VLenSq(dr) < rrClust) {
//         AddBondedPair(j1, j2);
//     }
// }

/**
 * \brief Add a bonded pair to the cluster.
 *
 * \param j1 Index of the first particle.
 * \param j2 Index of the second particle.
 */
// void AddBondedPair(int j1, int j2)
// {
//     int cBig, cSmall, m, mp, nc1, nc2;
//     nc1 = mol[j1].inClust;
//     nc2 = mol[j2].inClust;

//     if (nc1 < 0 && nc2 < 0) {
//         mol[j1].inClust = nClust;
//         mol[j2].inClust = nClust;
//         clust[nClust].size = 2;
//         clust[nClust].head = j1;
//         clust[j1].next = j2;
//         // TODO: Provide details about the specific implementation.
//         ...
//         clust[j2].next = -1;
//         ++nClust;
//     } else if (mol[j1].inClust < 0) {
//         mol[j1].inClust = nc2;
//         clust[j1].next = clust[nc2].head;
//         clust[nc2].head = j1;
//         ++clust[nc2].size;
//     } else if (mol[j2].inClust < 0) {
//         mol[j2].inClust = nc1;
//         clust[j2].next = clust[nc1].head;
//         clust[nc1].head = j2;
//         ++clust[nc1].size;
//     } else {
//         if (nc1 != nc2) {
//             cBig = (clust[nc1].size > clust[nc2].size) ? nc1 : nc2;
//             cSmall = nc1 + nc2 - cBig;
//             for (m = clust[cSmall].head; m >= 0; m = clust[m].next) {
//                 mol[m].inClust = cBig;
//                 mp = m;
//             }
//             clust[mp].next = clust[cBig].head;
//             clust[cBig].head = clust[cSmall].head;
//             clust[cBig].size += clust[cSmall].size;
//             clust[cSmall].size = 0;
//         }
//     }
// }

void InitClusters ()
{
    int n;
    DO_MOL mol[n].inClust = -1;
    nClust = 0;
}
/**
 * \brief Compress the clusters to remove empty clusters.
 */
// void CompressClusters()
// {
//     int j, m, nc;
//     nc = 0;

//     // Iterate over the clusters and copy non-empty clusters to new positions in the array.
//     for (j = 0; j < nClust; j++) {
//         if (clust[j].size > 0) {
//             clust[nc].head = clust[j].head;
//             clust[nc].size = clust[j].size;

//             // Update the cluster index for the particles in the cluster.
//             for (m = clust[nc].head; m >= 0; m = clust[m].next) {
//                 mol[m].inClust = nc;
//             }

//             ++nc;
//         }
//     }

//     // Update the total number of clusters.
//     nClust = nc;
// }

// /**
//  * \brief Analyze cluster sizes.
//  */
// void AnalClusterSize()
// {
//     int cBig, nc, ncUse;
//     PropZero(cSize);
//     ncUse = 0;
//     cBig = 0;

//     for (nc = 0; nc < nClust; nc++) {
//         cSize.val = clust[nc].size;

//         if (cSize.val > clust[cBig].size) {
//             cBig = nc;
//         }

//         if (cSize.val > 1) {
//             ++ncUse;
//             PropAccum(cSize);
//         }
//     }

//     bigSize = clust[cBig].size;
//     nSingle = nMol - cSize.sum;

//     if (ncUse > 0) {
//         PropAvg(cSize, ncUse);
//     }
// }

/**
 * \brief Get configuration data from a file.
 * \return 1 if successful, 0 otherwise.
 */
int GetConfig()
{
    VecR w;
    int fOk, n;
    short *rI;

    fOk = 1;

    if (blockNum == -1) {
        if ((fp = fopen(fileName[FL_SNAP], "r")) == 0) {
            fOk = 0;
        }
    } else {
        // 500 18 More about software
        fseek(fp, blockNum * blockSize, 0);
        ++blockNum;
    }

    if (fOk) {
        ReadF(blockSize);

        if (feof(fp)) {
            return (0);
        }

        ReadF(nMol);
        ReadF(region);
        ReadF(stepCount);
        ReadF(timeNow);

        if (blockNum == -1) {
            SetCellSize();
            AllocArrays();
            blockNum = 1;
        }

        AllocMem(rI, NDIM * nMol, short);
        ReadFN(rI, NDIM * n);

        DO_MOL {
            VFromLin(w, rI, NDIM * n);
            VScale(w, 1. / SCALE_FAC);
            VAddCon(w, w, -0.5);
            VMul(mol[n].r, w, region);
        }

        free(rI);

        if (ferror(fp)) {
            fOk = 0;
        }
    }

    if (!fOk) {
        ErrExit(ERR_SNAP_READ);
    }

    return (1);
}

/**
 * \brief Write configuration data to a file.
 */


void PutConfig()
{
    VecR w;
    int blockSize, fOk, n;
    short *rI;
    FILE *fp;

    fOk = 1;

    blockSize = (NDIM + 1) * sizeof(real) + 3 * sizeof(int) + nMol * NDIM * sizeof(short);

    if ((fp = fopen(fileName[FL_SNAP], "a")) != 0) {
        WriteF(blockSize);
        WriteF(nMol);
        WriteF(region);
        WriteF(stepCount);
        WriteF(timeNow);

        AllocMem(rI, NDIM * nMol, short);

        DO_MOL {
            VDiv(w, mol[n].r, region);
            VAddCon(w, w, 0.5);
            VScale(w, SCALE_FAC);
            VToLin(rI, NDIM * n, w);
        }

        WriteFN(rI, NDIM * nMol);
        free(rI);

        if (ferror(fp)) {
            fOk = 0;
        }

        fclose(fp);
    } else {
        fOk = 0;
    }

    if (!fOk) {
        ErrExit(ERR_SNAP_WRITE);
    }
}

/**
 * \brief Subdivide cells and compute forces between atoms.
 */
void SubdivCells()
{
    VecR dr, invWid, rs, shift;
    VecI cc, m1v, m2v, vOff[] = OFFSET_VALS;
    real fcVal, rr, rrCut, rri, rri3, uVal;
    int c, j1, j2, m1, m1x, m1y, m1z, m2, n, offset;

    rrCut = Sqr(rCut);
    VDiv(invWid, cells, region);

    for (n = nMol; n < nMol + VProd(cells); n++) {
        cellList[n] = -1;
    }

    DO_MOL {
        VSAdd(rs, mol[n].r, 0.5, region);
        VMul(cc, rs, invWid);
        c = VLinear(cc, cells) + nMol;
        cellList[n] = cellList[c];
        cellList[c] = n;
    }

    DO_MOL {
        VZero(mol[n].ra);
    }

    uSum = 0.;
    virSum = 0.;

    for (m1z = 0; m1z < cells.z; m1z++) {
        for (m1y = 0; m1y < cells.y; m1y++) {
            for (m1x = 0; m1x < cells.x; m1x++) {
                VSet(m1v, m1x, m1y, m1z);
                m1 = VLinear(m1v, cells) + nMol;

                for (offset = 0; offset < N_OFFSET; offset++) {
                    VAdd(m2v, m1v, vOff[offset]);
                    VZero(shift);
                    VCellWrapAll();

                    m2 = VLinear(m2v, cells) + nMol;

                    DO_CELL(j1, m1) {
                        DO_CELL(j2, m2) {
                            if (m1 != m2 || j2 < j1) {
                                VSub(dr, mol[j1].r, mol[j2].r);
                                VVSub(dr, shift);
                                rr = VLenSq(dr);

                                if (rr < rrCut) {
                                    rri = 1. / rr;
                                    rri3 = Cube(rri);
                                    fcVal = 48. * rri3 * (rri3 - 0.5) * rri;
                                    uVal = 4. * rri3 * (rri3 - 1.) + 1.;

                                    VVSAdd(mol[j1].ra, fcVal, dr);
                                    VVSAdd(mol[j2].ra, -fcVal, dr);

                                    uSum += uVal;
                                    virSum += fcVal * rr;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}
void Sort(real *a, int *seq, int n)
{
    real q;
    int i, ir, ixt, j, k;

    for (j = 0; j < n; j++)
        seq[j] = j;

    if (n > 1) {
        k = n / 2;
        ir = n - 1;

        while (1) {
            if (k > 0) {
                --k;
                ixt = seq[k];
                q = a[ixt];
            } else {
                ixt = seq[ir];
                q = a[ixt];
                seq[ir] = seq[0];
                --ir;

                if (ir == 0) {
                    seq[0] = ixt;
                    break;
                }
            }

            i = k;
            j = 2 * k + 1;

            while (j <= ir) {
                if (j < ir && a[seq[j]] < a[seq[j + 1]])
                    ++j;

                if (q < a[seq[j]]) {
                    seq[i] = seq[j];
                    i = j;
                    j = 2 * j + 1;
                } else {
                    j = ir + 1;
                }
            }

            seq[i] = ixt;
        }
    }
}


/**
 * @brief Main function that runs the program.
 *
 * @param argc The number of command-line arguments.
 * @param argv An array of strings containing the command-line arguments.
 * @return int The exit status of the program.
 */

int main(int argc, char **argv)
{
    int n, na;
    /**
     * @brief The run ID obtained from the command-line argument.
     */
    runId = atoi(argv[1]);

    cellRatio = 0.5;

    /**
     * @brief Set up the necessary files.
     */
    SetupFiles();

    blockNum = -1;

    while (GetConfig())
    {
        regionVol = VProd(region);

        SubdivCells();

        rangeLim = region.x;

        /**
         * @brief Initialize properties to zero.
         */
        PropZero(polyArea);
        PropZero(polyVol);

        for (n = 0; n < 4; n++)
        {
            PropZero(polyGeom[n]);
        }

        for (na = 0; na < nMol; na++)
        {
            FindTestSites(na);
            AnalVorPoly();
            PropAccum(polyArea);
            PropAccum(polyVol);

            for (n = 0; n < 4; n++)
            {
                PropAccum(polyGeom[n]);
            }
        }

        fracPolyVol = polyVol.sum / regionVol;

        PropAvg(polyArea, nMol);
        PropAvg(polyVol, nMol);

        for (n = 0; n < 4; n++)
        {
            PropAvg(polyGeom[n], nMol);
        }

        polyArea.sum /= pow(regionVol, 2. / 3.);
        polyArea.sum2 /= pow(regionVol, 2. / 3.);
        polyVol.sum /= regionVol;
        polyVol.sum2 /= regionVol;

        eulerSum = polyGeom[0].sum + polyGeom[2].sum - polyGeom[1].sum;
        printf("eulerSum: %d", eulerSum);
        // ... (print the results) ...
    }
}

//***Cluster Analysis Main***

// int main (int argc, char **argv)
// {
//     runId = atoi (argv[1]);
//     rClust = atof (argv[2]);
//     SetupFiles ();
//     blockNum = -1;
//     while (GetConfig ()) {
//         InitClusters ();
//         BuildClusters ();
//         CompressClusters ();
//         AnalClusterSize ();
//         printf ("%d %d %d %.1f %.1f\n", nSingle, nClust, bigSize,
//             PropEst (cSize));
//             }
// }
//int main(int argc, char **argv) {
//     GetNameList(argc, argv);
//     PrintNameList(stdout);
//     SetParams();
//     SetupJob();
//     printf("Starting simulation...\n");
//     int moreCycles = 1;
//     while (moreCycles) {
//         SingleStep();
//         if (stepCount >= stepLimit)
//             moreCycles = 0;
//     }

//     int n;
//     // Generate and run the pov scene
//     init_pov("scene.pov");

//     DO_MOL {
//         // scale coordinate data
//         scaleData(&mol[n].r.x, &mol[n].r.y, &mol[n].r.z);
//         write_pov(mol[n].r.x, mol[n].r.y, mol[n].r.z);
//     }
//     destroyPOVRayWriter(writer);
//     runPovEngine("scene.pov");
//     printf("Done\n");
//     return 0;
// }
// 
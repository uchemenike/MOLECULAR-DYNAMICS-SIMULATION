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
        AccumProps (2);
        if (stepCount < stepEquil) AdjustInitTemp ();
        PrintSummary (stdout);
        AccumProps (0);
    }
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
    AllocMem (mol, nMol, Mol);
    AllocMem (histVel, sizeHistVel, real);
    AllocMem (nebrTab, 2 * nebrTabMax, int);
    AllocMem (cellList, VProd (cells) + nMol, int);
}

void SetupJob () {
    AllocArrays ();
    stepCount = 0;
    countVel = 0;
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
    for (k = 0; k < sizeof(fileNameR) / sizeof(fileNameR[0]); k++)
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
int main(int argc, char **argv) {
    GetNameList(argc, argv);
    PrintNameList(stdout);
    SetParams();
    SetupJob();
    printf("Starting simulation...\n");
    int moreCycles = 1;
    while (moreCycles) {
        SingleStep();
        if (stepCount >= stepLimit)
            moreCycles = 0;
    }

    int n;
    // Generate and run the pov scene
    init_pov("scene.pov");

    DO_MOL {
        // scale coordinate data
        scaleData(&mol[n].r.x, &mol[n].r.y, &mol[n].r.z);
        write_pov(mol[n].r.x, mol[n].r.y, mol[n].r.z);
    }
    destroyPOVRayWriter(writer);
    runPovEngine("scene.pov");
    printf("Done\n");
    return 0;
}
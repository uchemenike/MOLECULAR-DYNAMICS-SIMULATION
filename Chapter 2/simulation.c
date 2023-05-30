/**
 * @file simulation.h
 * @brief Header file for the simulation chapter 2.
 */

#include "simulation.h"

/**
 * @brief Prints an error message and exits the program.
 * @param code The error code.
 */
void ErrExit (int code)
{
    printf ("Error: %s\n", errorMsg[code]);
    exit (0);
}

/**
 * @brief Generates a random real number.
 * @return The random real number.
 */
real RandR ()
{
    randSeedP = (randSeedP * IMUL + IADD) & MASK; 
    return (randSeedP * SCALE);
}
/**
 * @brief Generates a random vector with unit length.
 * @param p Pointer to the vector object.
 */

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

/**
 * @brief Gets the name list from the command line arguments or a file.
 * @param argc The number of command line arguments.
 * @param argv Array of command line arguments.
 * @return 1 if successful, 0 otherwise.
 */

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
    printf("%s", buff);
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

/**
 * @brief Prints the name list to a file.
 * @param fp Pointer to the file object.
 */

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

/**
 * @brief Prints the velocity distribution histogram to a file.
 * @param fp Pointer to the file object.
 */

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

/**
 * @brief Evaluates the velocity distribution histogram.
 */

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

/**
 * @brief Sets the simulation parameters.
 */

void SetParams()
{
    // Calculate the cutoff distance for the potential function
    rCut = pow(2.0, 1.0 / 6.0);

    // Copy the values of 1.0 / sqrt(density) into the region array
    VSCopy(region, 1.0 / sqrt(density), initUcell);

    // Calculate the total number of molecules in the system
    nMol = VProd(initUcell);

    // Calculate the magnitude of the velocity
    velMag = sqrt(NDIM * (1.0 - 1.0 / nMol) * temperature);
}

/**
 * @brief Performs a leapfrog step for particle update.
 * @param part The part of the leapfrog step (1 or 2).
 */

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

/**
 * @brief Applies the boundary conditions to the particles.
 */

void ApplyBoundaryCond()
{
    int n;
    DO_MOL VWrapAll(mol[n].r);
}

/**
 * @brief Evaluates the properties of the system (kinetic energy, total energy, pressure).
 */
void EvalProps()
{
    real vv;
    int n;
    VZero(vSum);
    vvSum = 0;

    DO_MOL {
        VVAdd(vSum, mol[n].rv);
        vv = VLenSq(mol[n].rv);
        vvSum += vv;
    }

    kinEnergy.val = 0.5 * vvSum / nMol;
    totEnergy.val = kinEnergy.val + uSum / nMol;
    pressure.val = density * (vvSum + virSum) / (nMol * NDIM);
}

/**
 * @brief Accumulates the properties of the system.
 * @param icode The accumulation code (0 for reset, 1 for accumulate, 2 for average).
 */

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

/**
 * @brief Prints a summary of the simulation properties to a file.
 * @param fp Pointer to the file object.
 */

void PrintSummary(FILE *fp)
{
    fprintf(fp,
            "%5d %8.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f\n",
            stepCount, timeNow, VCSum(vSum) / nMol, PropEst(totEnergy),
            PropEst(kinEnergy), PropEst(pressure));
}

/**
 * @brief Performs a single step of the simulation.
 */

void SingleStep () {
    ++ stepCount;
    timeNow = stepCount * deltaT;
    LeapfrogStep (1);
    ApplyBoundaryCond ();
    ComputeForces ();
    LeapfrogStep (2);
    EvalProps ();
    AccumProps (1);
    if (stepCount % stepAvg == 0) {
        AccumProps (2);
        PrintSummary (stdout);
        AccumProps (0);
    }
}

/**
 * @brief Initializes the coordinates of the particles.
 */
void InitCoords()
{
    VecR c, gap;
    int n, nx, ny;

    VDiv(gap, region, initUcell);
    n = 0;
    for (ny = 0; ny < initUcell.y; ny++) {
        for (nx = 0; nx < initUcell.x; nx++) {
            VSet(c, nx + 0.5, ny + 0.5);
            VMul(c, c, gap);
            VVSAdd(c, -0.5, region);
            mol[n].r = c;
            ++n;
        }
    }
}

/**
 * @brief Initializes the velocities of the particles.
 */
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

/**
 * @brief Initializes the accelerations of the particles.
 */
void InitAccels()
{
    int n;
    DO_MOL
    {
        VZero(mol[n].ra);
    }
}

/**
 * @brief Allocates memory for arrays used in the simulation.
 */
void AllocArrays ()
{
    AllocMem (mol, nMol, Mol);
    AllocMem (histVel, sizeHistVel, real);
}

/**
 * @brief Sets up the simulation job by allocating memory and initializing particles.
 */
void SetupJob () {
    AllocArrays ();
    stepCount = 0;
    countVel = 0;
    InitCoords ();
    InitVels ();
    InitAccels ();
    AccumProps (0);
}

/**
 * @brief Computes the forces between particles.
 */
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


/**
 * @brief Sets up the necessary files and filenames for checkpointing and snapshot records.
 */
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

/**
 * @brief Writes checkpoint data to disk.
 */

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

/**
 * @brief Reads checkpoint data from disk.
 */

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

/**
 * @brief Creates a POVRayWriter object and initializes it.
 * @param fileName The name of the output file.
 * @return The created POVRayWriter object.
 */
POVRayWriter* createPOVRayWriter(const char* fileName) {
    POVRayWriter* writer = (POVRayWriter*)malloc(sizeof(POVRayWriter));
    writer->fileName = (char*)malloc(MAX_STRING_LENGTH * sizeof(char));
    strcpy(writer->fileName, fileName);
    writer->file = fopen(fileName, "w");
    return writer;
}

/**
 * @brief Destroys a POVRayWriter object and frees associated memory.
 * @param writer The POVRayWriter object to destroy.
 */
void destroyPOVRayWriter(POVRayWriter* writer) {
    if (writer != NULL) {
        if (writer->file != NULL) {
            fclose(writer->file);
        }
        free(writer->fileName);
        free(writer);
    }
}
/**
 * @brief Writes the necessary header includes to the POV-Ray scene file.
 * @param writer The POVRayWriter object.
 */
void writeHeader(POVRayWriter* writer) {
    fprintf(writer->file, "#include \"colors.inc\"\n");
    fprintf(writer->file, "#include \"textures.inc\"\n");
    fprintf(writer->file, "#include \"glass.inc\"\n");
    fprintf(writer->file, "#include \"metals.inc\"\n");
    fprintf(writer->file, "#include \"woods.inc\"\n\n");
}

/**
 * @brief Writes the camera properties to the POV-Ray scene file.
 * @param writer The POVRayWriter object.
 * @param position The position of the camera.
 * @param lookAt The point the camera is looking at.
 * @param up The up vector of the camera.
 */

void writeCamera(POVRayWriter* writer, Vector3D position, Vector3D lookAt, Vector3D up) {
    fprintf(writer->file, "camera {\n");
    fprintf(writer->file, "  location <%f, %f, %f>\n", position.x, position.y, position.z);
    fprintf(writer->file, "  look_at <%f, %f, %f>\n", lookAt.x, lookAt.y, lookAt.z);
    fprintf(writer->file, "  up <%f, %f, %f>\n", up.x, up.y, up.z);
    fprintf(writer->file, "}\n\n");
}
/**
 * @brief Writes a sphere object to the POV-Ray scene file.
 * @param writer The POVRayWriter object.
 * @param sphere The sphere object to write.
 */

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
/**
 * @brief Initializes the POV-Ray scene by creating the POVRayWriter object and writing the header.
 * @param filename The name of the output file.
 * @return 0 on success, or an error code on failure.
 */
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

/**
 * @brief Writes a sphere object with the given coordinates to the POV-Ray scene.
 * @param x The x-coordinate of the sphere.
 * @param y The y-coordinate of the sphere.
 * @param z The z-coordinate of the sphere.
 * @return 0 on success, or an error code on failure.
 */

int write_pov(double x, double y, double z){
    
    Vector3D position = { x, y, z};
    double radius = 0.4;

    Sphere sphere = { position, radius, { "Gray", 0.1, 0.6, 0.4, 0.5, 0.0 } };
    writeSphere(writer, sphere);
    return 0;
}
/**
 * @brief Runs the POV-Ray engine to render the POV-Ray scene file.
 * @param sceneFile The name of the POV-Ray scene file.
 */

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

/**
 * @brief Scales the given coordinates using predefined scaling factors.
 * @param x The x-coordinate to scale.
 * @param y The y-coordinate to scale.
 */

void scaleData(int* x, int* y) {
    // Define the scaling factors for each axis
    double scaleX = 0.00000001;
    double scaleY = 0.00000001;
    //double scaleZ = 0.00000001;

    // Scale the coordinates using the scaling factors
    *x *= scaleX;
    *y *= scaleY;
    //*z *= scaleZ;
}

/**
 * @brief The main program entry point.
 * @param argc The number of command-line arguments.
 * @param argv An array of command-line argument strings.
 * @return 0 on success, or an error code on failure.
 */

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

    int i = 0;
    // Generate and run the pov scene
    init_pov("scene.pov");

    while(i < nMol) {
        // scale coordinate data
        scaleData(&mol[i].r.x, &mol[i].r.y);
        write_pov(mol[i].r.x, mol[i].r.y, 0);
        i++;
    }
    destroyPOVRayWriter(writer);
    runPovEngine("scene.pov");
    printf("Done\n");
    return 0;
}
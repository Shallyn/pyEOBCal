/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel: 
* 20xx.xx.xx, LOC
**/

#include "myIO.h"

INT cmd_mkdir(CHAR *folderName)
{
    CHAR cmd[256] = "mkdir ";
    strcat(cmd, folderName);
    if (access(folderName, 0) == -1)
    {
        system(cmd);
    }
    return CEV_SUCCESS;
}

INT read_waveform(REAL8Vector **time, 
                  REAL8Vector **hreal, 
                  REAL8Vector **himag,
                  FILE *file)
{
    const size_t block = 1024;
    REAL8 start;
    REAL8 end;
    double *t = NULL;
    double *hp = NULL;
    double *hc = NULL;
    size_t bufsz = 0;
    size_t n, l;
    char line[LINE_MAX];

    for (l = 0, n = 0; fgets(line, sizeof(line), file); ++l) 
    {
        int c;
        if (*line == '#')
            continue;
        if (n == bufsz) {       /* allocate more memory */
            bufsz += block;
            hp = realloc(hp, bufsz * sizeof(*hp));
            hc = realloc(hc, bufsz * sizeof(*hc));
        }
        c = sscanf(line, "%le %le %le", t + n, hp + n, hc + n);
        if (c != 3) {
            fprintf(stderr, "error: format error on line %zd: %s\n", l, line);
            free(t);
            free(hp);
            free(hc);
            return CEV_FAILURE;
        }
        ++n;
    }
    hp = realloc(hp, n * sizeof(*hp));
    hc = realloc(hc, n * sizeof(*hp));
    REAL8Vector *tVec = CreateREAL8Vector(n);
    REAL8Vector *hrVec = CreateREAL8Vector(n);
    REAL8Vector *hiVec = CreateREAL8Vector(n);
    memcpy(tVec->data, t, n*sizeof(*t));
    memcpy(hrVec->data, hp, n * sizeof(*hp));
    memcpy(hiVec->data, hc, n * sizeof(*hc));
    free(t);
    free(hp);
    free(hc);
    return CEV_SUCCESS;
}

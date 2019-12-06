/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel: 
* 20xx.xx.xx, LOC
**/

#ifndef __INCLUDE_MYIO__
#define __INCLUDE_MYIO__

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include "myDatatypes.h"
#include "myUtils.h"

INT cmd_mkdir(CHAR folderName[]);
INT read_waveform(REAL8Vector **time, 
                  REAL8Vector **hreal, 
                  REAL8Vector **himag,
                  FILE *file);


#endif


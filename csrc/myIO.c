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

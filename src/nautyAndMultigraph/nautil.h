//
// Created by veetoo on 7/12/24.
//

#ifndef NAUTIL_H
#define NAUTIL_H

#include "nauty.h"

#ifdef __cplusplus
extern "C" {
#endif

    int nextelement(set*, int, int );
    int itos(int, char*);
    void nautil_check(int, int, int, int);

#ifdef __cplusplus
}
#endif

#endif // NAUTIL_H
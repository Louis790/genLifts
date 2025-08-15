//
// Created by veetoo on 7/12/24.
//

#ifndef GTNAUTY_H
#define GTNAUTY_H

#include "nauty.h"

#ifdef __cplusplus
extern "C" {
#endif

    void fcanonise_inv_sg(sparsegraph*, int , int , sparsegraph*, char*,
        void (*invarproc)(graph*,int*,int*,int,int,int,permutation*,int, boolean,int,int), int, int, int, boolean);

#ifdef __cplusplus
}
#endif

#endif //GTNAUTY_H

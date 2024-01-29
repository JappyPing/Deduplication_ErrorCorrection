/*
 *  macros.h
 *
 *  Created on: Jun 1, 2015
 *      Author: gonzales
 */

#ifndef MACROS_H_
#define MACROS_H_
#include <iostream>
#include <vector>
#include <algorithm>
#include <set>
#include <map>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <stdint.h>
#include <assert.h>
#include <limits.h>
#include <zlib.h>
#include <math.h>
#include <thread>
#include <mutex>
#include <unordered_map>
using namespace std;

/*the version and name*/
#define PROGRAM_NAME	"ParDRe"
#define PROGRAM_VERSION "2.1.5-PC170509"

typedef unsigned int UINT32;
typedef unsigned short WORD;
typedef unsigned char BYTE;
typedef BYTE VECTOR[16];

#define FILE_FORMAT_FASTQ   1
#define FILE_FORMAT_FASTA   2

#endif /* MACROS_H_ */

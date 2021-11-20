/*
 * MappingValues.h
 *
 *      Author: Julia
 */

#ifndef MAPPINGVALUES_H_
#define MAPPINGVALUES_H_


#endif /* MAPPINGVALUES_H_ */

//This array maps the ascii values to a decimal value
//$(36): 0; @(64): 1; A(65): 2; C(67): 3; G(71): 4; T(84): 5;
const static int mapVals [] = {
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 3, 0, 0, 0, 0, 0,};

//This array maps the ascii values to a decimal value
//$(36): 0; @(64): 1; A(65): 2; C(67): 3; G(71): 4; T(84): 5;
const static int mapVals2 [] = {
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 2, 0, 3, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 5, 0, 0, 0, 0, 0,};

const static int qualScores [] = {
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
	1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
	1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
	1, 1, 1, 1, 1, 1, 1, 1, 1, 1,};

//This is the inverse of the mapVals array
const static char inverseMapVals [] = {'A', 'C', 'G', 'T'};
const static char inverseComplementMapVals [] = {'T', 'G', 'C', 'A'};

#define maskValR 0xffffffffffffffffULL //all 1s
#define appValL 0x8000000000000000ULL  //1 then all 0s
#define megabyte 1048576

//Exit signals
#define help 911
#define bufOverflow 792
#define inputProblem 716
#define runSerial 112
#define workTag 1
#define dieTag 2

//Graph Param
#define numFlags 5

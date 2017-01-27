#ifndef _const_
#define _const_

#define PYLM_INITIAL_D 	0.2
#define PYLM_INITIAL_THETA 2.0
#define PYLM_INITIAL_A 	1.0
#define PYLM_INITIAL_B 	1.0
#define PYLM_INITIAL_ALPHA 1.0
#define PYLM_INITIAL_BETA  1.0

#define BEGIN_OF_SENTENSE 0
#define END_OF_SENTENSE 0
#define BEGIN_OF_POS 0
#define END_OF_POS 0

typedef struct Word {
	int word_id;
	int tag_id;
} Word;

#endif
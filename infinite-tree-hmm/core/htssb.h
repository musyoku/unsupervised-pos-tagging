#ifndef _htssb_
#define _htssb_
#include "node.h"

class HTSSB{
public:
	Node* _root;
	HTSSB(){
		_root = new Node(NULL);
	}
};

#endif
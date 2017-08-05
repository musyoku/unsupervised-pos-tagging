#pragma once

using id = int;

namespace bhmm {
	 class Node;
	typedef struct Word {
		id _id;
		int _state;
	} Word;
}
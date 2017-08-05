#pragma once

using id = int;

namespace ithmm {
	class Node;
	typedef struct Word {
		id _id;
		Node* _state;
	} Word;
}
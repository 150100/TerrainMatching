#pragma once

template <Element>
class RBTree
{
	Element root;
	RBTree leftSubtree;
	RBTree rightSubtree;
	int color;

public:
	RBTree();
	~RBTree();
};


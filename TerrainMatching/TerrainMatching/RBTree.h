#ifndef RBTree_h
#define RBTree_h

template <T>
class RBTree
{
	class RBTreeNode
	{
	public:
		enum Color { RED, BLACK };

		RBTreeNode(T *_v) : value(_v) : color(BLACK) {}

		bool operator<(RBTreeNode &_node) { return *value < *_node->value; }

		T getValue() { return *value; }

		Color color;
		T* value;
	};

	RBTreeNode root;
	RBTree leftSubtree;
	RBTree rightSubtree;
	int color;

public:
	RBTree();
	~RBTree();
};


#endif //RBTREE_H
#ifndef _splay_
#define _splay_
#include <iostream>
#include <cstdlib>
#include <boost/serialization/serialization.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

typedef int KEY_TYPE;

// typedef struct splay{
// 	KEY_TYPE key;
// 	struct splay* lchild;
// 	struct splay* rchild;
// }splay;

namespace Splay {
	template <typename T>
	class Node {
	public:
		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive& archive, unsigned int version)
		{
			static_cast<void>(version); // No use
			archive & key;
			archive & value;
			archive & lchild;
			archive & rchild;
		}
		int sum_num_child(){
			int sum = 0;
			if(lchild){
				sum += lchild->sum_num_child() + 1;
			}
			if(rchild){
				sum += rchild->sum_num_child() + 1;
			}
			return sum;
		}
		// Node(KEY_TYPE _key, T _value){
		// 	key = _key;
		// 	value = _value;
		// }
		KEY_TYPE key;
		T value;
		Node<T>* lchild;
		Node<T>* rchild;
	};

	template <class T>
	class Node<T*> {
	public:
		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive& archive, unsigned int version)
		{
			static_cast<void>(version); // No use
			archive & key;
			archive & value;
			archive & lchild;
			archive & rchild;
		}
		int sum_num_child(){
			int sum = 0;
			if(lchild){
				sum += lchild->sum_num_child() + 1;
			}
			if(rchild){
				sum += rchild->sum_num_child() + 1;
			}
			return sum;
		}
		// Node(KEY_TYPE _key, T _value){
		// 	key = _key;
		// 	value = _value;
		// }
		KEY_TYPE key;
		T* value;
		Node<T*>* lchild;
		Node<T*>* rchild;
	};

	template <typename T>
	class Tree {
	public:
		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive& archive, unsigned int version)
		{
			static_cast<void>(version); // No use
			archive & _p_node;
			archive & _root;
		}
		Node<T>* _p_node;
		Node<T>* _root;
		Tree(){
			_p_node = NULL;
			_root = NULL;
		}
		bool is_empty(){
			if(_root == NULL){
				return true;
			}
			return false;
		}
		/* RR(Y rotates to the right):

		        k2                   k1
		       /  \                 /  \
		      k1   Z     ==>       X   k2
		     / \                      /  \
		    X   Y                    Y    Z
		*/
		inline Node<T>* RR_Rotate(Node<T>* k2)
		{
			Node<T>* k1 = k2->lchild;
			k2->lchild = k1->rchild;
			k1->rchild = k2;
			return k1;
		}

		/* LL(Y rotates to the left):

		        k2                       k1
		       /  \                     /  \
		      X    k1         ==>      k2   Z
		          /  \                /  \
		         Y    Z              X    Y
		 */
		inline Node<T>* LL_Rotate(Node<T>* k2)
		{
			Node<T>* k1 = k2->rchild;
			k2->rchild = k1->lchild;
			k1->lchild = k2;
			return k1;
		}

		/* An implementation of top-down splay tree 
		 If key is in the tree, then the node containing the key will be rotated to root,
		 else the last non-NULL node (on the search path) will be rotated to root.
		 */
		Node<T>* Splay(KEY_TYPE key, Node<T>* target)
		{
			if(!target){
				return NULL;
			}
			Node<T> header;
			/* header.rchild points to L tree; header.lchild points to R Tree */
			header.lchild = header.rchild = NULL;
			Node<T>* LeftTreeMax = &header;
			Node<T>* RightTreeMin = &header;
			
			/* loop until target->lchild == NULL || target->rchild == NULL; then break!
			   (or when find the key, break too.)
			 The zig/zag mode would only happen when cannot find key and will reach
			 null on one side after RR or LL Rotation.
			 */
			while(1){
				if(key < target->key){
					if(!target->lchild){
						break;
					}
					if(key < target->lchild->key){
						target = RR_Rotate(target); 
						if(!target->lchild){
							break;
						}
					}
					/* Link to R Tree */
					RightTreeMin->lchild = target;
					RightTreeMin = RightTreeMin->lchild;
					target = target->lchild;
					RightTreeMin->lchild = NULL;
				}else if(key > target->key){
					if(!target->rchild){
						break;
					}
					if(key > target->rchild->key){
						target = LL_Rotate(target);
						if(!target->rchild){
							break;
						}
					}
					/* Link to L Tree */
					LeftTreeMax->rchild = target;
					LeftTreeMax = LeftTreeMax->rchild;
					target = target->rchild;
					LeftTreeMax->rchild = NULL;
				}else{
					break;
				}
			}
			LeftTreeMax->rchild = target->lchild;
			RightTreeMin->lchild = target->rchild;
			target->lchild = header.rchild;
			target->rchild = header.lchild;
			return target;
		}
		void insert(KEY_TYPE key, T value){
			_root = _insert(key, value, _root);
		}
		Node<T>* _insert(KEY_TYPE key, T value, Node<T>* root)
		{
			if(!_p_node){
				_p_node = new Node<T>;
				if(!_p_node){
					fprintf(stderr, "Out of memory!\n");
					exit(1);
				}
				_p_node->key = key;
				_p_node->value = value;
				_p_node->lchild = _p_node->rchild = NULL;
			}else{
				 // could take advantage of the node remains because of there was duplicate key before.
				_p_node->key = key;
				_p_node->value = value;
			}
			if(!root){
				root = _p_node;
				_p_node = NULL;
				return root;
			}
			root = Splay(key, root);
			/* This is BST that, all keys <= root->key is in root->lchild, all keys > 
			   root->key is in root->rchild. (This BST doesn't allow duplicate keys.) */
			if(key < root->key)
			{
				_p_node->lchild = root->lchild;
				_p_node->rchild = root;
				root->lchild = NULL;
				root = _p_node;
			}
			else if(key > root->key)
			{
				_p_node->rchild = root->rchild;
				_p_node->lchild = root;
				root->rchild = NULL;
				root = _p_node;
			}
			else
				return root;
			_p_node = NULL;
			return root;
		}
		bool delete_key(KEY_TYPE key)
		{
			Node<T>* temp;
			if(!_root){
				return false;
			}
			_root = Splay(key, _root);
			if(key != _root->key){
				// No such node in splay tree
				return false;
			}else{
				if(!_root->lchild){
					temp = _root;
					_root = _root->rchild;
				}else{
					temp = _root;
					/*Note: Since key == _root->key, so after Splay(key, _root->lchild), 
					  the tree we get will have no right child tree. (key > any key in 
					  _root->lchild)*/
					_root = Splay(key, _root->lchild);
					_root->rchild = temp->rchild;
				}
				free(temp);
			}
			return true;
		}
		Node<T>* search(KEY_TYPE key)
		{
			_root =  Splay(key, _root);
			if(_root->key == key){
				return _root;
			}
			return NULL;
		}
	};

	template <class T>
	class Tree<T*> {
	public:
		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive& archive, unsigned int version)
		{
			static_cast<void>(version); // No use
			archive & _p_node;
			archive & _root;
		}
		Node<T*>* _p_node;
		Node<T*>* _root;
		Tree(){
			_p_node = NULL;
			_root = NULL;
		}
		bool is_empty(){
			if(_root == NULL){
				return true;
			}
			return false;
		}
		int size(){
			if(_root == NULL){
				return 0;
			}
			return 1 + _root->sum_num_child();
		}
		/* RR(Y rotates to the right):

		        k2                   k1
		       /  \                 /  \
		      k1   Z     ==>       X   k2
		     / \                      /  \
		    X   Y                    Y    Z
		*/
		inline Node<T*>* RR_Rotate(Node<T*>* k2)
		{
			Node<T*>* k1 = k2->lchild;
			k2->lchild = k1->rchild;
			k1->rchild = k2;
			return k1;
		}

		/* LL(Y rotates to the left):

		        k2                       k1
		       /  \                     /  \
		      X    k1         ==>      k2   Z
		          /  \                /  \
		         Y    Z              X    Y
		 */
		inline Node<T*>* LL_Rotate(Node<T*>* k2)
		{
			Node<T*>* k1 = k2->rchild;
			k2->rchild = k1->lchild;
			k1->lchild = k2;
			return k1;
		}

		/* An implementation of top-down splay tree 
		 If key is in the tree, then the node containing the key will be rotated to root,
		 else the last non-NULL node (on the search path) will be rotated to root.
		 */
		Node<T*>* Splay(KEY_TYPE key, Node<T*>* target)
		{
			if(!target){
				return NULL;
			}
			Node<T*> header;
			/* header.rchild points to L tree; header.lchild points to R Tree */
			header.lchild = header.rchild = NULL;
			Node<T*>* LeftTreeMax = &header;
			Node<T*>* RightTreeMin = &header;
			
			/* loop until target->lchild == NULL || target->rchild == NULL; then break!
			   (or when find the key, break too.)
			 The zig/zag mode would only happen when cannot find key and will reach
			 null on one side after RR or LL Rotation.
			 */
			while(1){
				if(key < target->key){
					if(!target->lchild){
						break;
					}
					if(key < target->lchild->key){
						target = RR_Rotate(target); 
						if(!target->lchild){
							break;
						}
					}
					/* Link to R Tree */
					RightTreeMin->lchild = target;
					RightTreeMin = RightTreeMin->lchild;
					target = target->lchild;
					RightTreeMin->lchild = NULL;
				}else if(key > target->key){
					if(!target->rchild){
						break;
					}
					if(key > target->rchild->key){
						target = LL_Rotate(target);
						if(!target->rchild){
							break;
						}
					}
					/* Link to L Tree */
					LeftTreeMax->rchild = target;
					LeftTreeMax = LeftTreeMax->rchild;
					target = target->rchild;
					LeftTreeMax->rchild = NULL;
				}else{
					break;
				}
			}
			LeftTreeMax->rchild = target->lchild;
			RightTreeMin->lchild = target->rchild;
			target->lchild = header.rchild;
			target->rchild = header.lchild;
			return target;
		}
		void insert(KEY_TYPE key, T* value){
			_root = _insert(key, value, _root);
		}
		Node<T*>* _insert(KEY_TYPE key, T* value, Node<T*>* root)
		{
			if(!_p_node){
				_p_node = new Node<T*>;
				if(!_p_node){
					fprintf(stderr, "Out of memory!\n");
					exit(1);
				}
				_p_node->key = key;
				_p_node->value = value;
				_p_node->lchild = _p_node->rchild = NULL;
			}else{
				// could take advantage of the node remains because of there was duplicate key before.
				_p_node->key = key;
				_p_node->value = value;
			}
			if(!root)
			{
				root = _p_node;
				_p_node = NULL;
				return root;
			}
			root = Splay(key, root);
			/* This is BST that, all keys <= root->key is in root->lchild, all keys > 
			   root->key is in root->rchild. (This BST doesn't allow duplicate keys.) */
			if(key < root->key)
			{
				_p_node->lchild = root->lchild;
				_p_node->rchild = root;
				root->lchild = NULL;
				root = _p_node;
			}
			else if(key > root->key)
			{
				_p_node->rchild = root->rchild;
				_p_node->lchild = root;
				root->rchild = NULL;
				root = _p_node;
			}
			else{
				return root;
			}
			_p_node = NULL;
			return root;
		}
		bool delete_key(KEY_TYPE key)
		{
			Node<T*>* temp;
			if(!_root){
				return false;
			}
			_root = Splay(key, _root);
			if(key != _root->key){
				// No such node in splay tree
				return false;
			}else{
				if(!_root->lchild){
					temp = _root;
					_root = _root->rchild;
				}else{
					temp = _root;
					/*Note: Since key == _root->key, so after Splay(key, _root->lchild), 
					  the tree we get will have no right child tree. (key > any key in 
					  _root->lchild)*/
					_root = Splay(key, _root->lchild);
					_root->rchild = temp->rchild;
				}
				free(temp);
			}
			return true;
		}
		Node<T*>* search(KEY_TYPE key)
		{
			_root =  Splay(key, _root);
			if(_root->key != key){
				return NULL;
			}
			return _root;
		}
	};
}

#endif
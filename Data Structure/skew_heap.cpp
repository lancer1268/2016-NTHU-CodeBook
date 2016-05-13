#ifndef SKEW_HEAP
#define SKEW_HEAP
template<typename T,typename _Compare=std::less<T> >
class skew_heap{
	private:
		struct node{
			T data;
			node *l,*r;
			node(const T&d):data(d),l(0),r(0){}
		}*root;
		int _size;
		_Compare cmp;
		node *merge(node *a,node *b){
			if(!a||!b)return a?a:b;
			if(cmp(a->data,b->data))return merge(b,a);
			node *t=a->r;
			a->r=a->l;
			a->l=merge(b,t);
			return a;
		}
		void _clear(node *&o){
			if(o)_clear(o->l),_clear(o->r),delete o;
		}
	public:
		skew_heap():root(0),_size(0){}
		~skew_heap(){_clear(root);}
		inline void clear(){
			_clear(root);root=0;_size=0;
		}
		inline void join(skew_heap &o){
			root=merge(root,o.root);
			o.root=0;
			_size+=o._size;
			o._size=0;
		}
		inline void swap(skew_heap &o){
			node *t=root;
			root=o.root;
			o.root=t;
			int st=_size;
			_size=o._size;
			o._size=st;
		}
		inline void push(const T&data){
			_size++;
			root=merge(root,new node(data));
		}
		inline void pop(){
			if(_size)_size--;
			node *tmd=merge(root->l,root->r);
			delete root;
			root=tmd;
		}
		inline const T& top(){return root->data;}
		inline int size(){return _size;}
		inline bool empty(){return !_size;}
};
#endif

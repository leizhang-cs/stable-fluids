#include "test.h"
using namespace std;



int main(){
    return 0;
}


Algorithm
void Build(setTri)
    n = setTri.size;
    Dim = largest Dimension;
    sort on Dim according to bounding box position;

    // construction a complete binary tree
    tree[2*n-1];
    for(i = n-1: 2*n-1)
        tree[i] = Tri[i].box;
    end for

    for(i = n-1: 0)
        tree[i/2-1] = Union(tree[i], tree[i-1]);
    end for    

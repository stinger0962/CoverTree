#ifndef _SINGLE_TREE_TRAVERSAL_H_
#define _SINGLE_TREE_TRAVERSAL_H_
#define prune_type unsigned long long


#include <ratio>
#include <numeric>
#include <chrono>
#include <cfloat>
#include <vector>
#include <omp.h>
#include "Clock.hpp"
#include <time.h>


template <typename Tree, typename Rule>
class SingleTreeTraversal {
  public:
    typedef typename Tree::NodeTree Box;
    Tree& tree;
    Rule& rules;
    double  timeSate = 0;
    std::chrono::high_resolution_clock::time_point start;
    /* Convenience constructor */
    SingleTreeTraversal(Tree& t, Rule& r, int&) : tree(t), rules(r) {}
    SingleTreeTraversal(Tree& t, Rule& r,prune_type&) : tree(t), rules(r) {}

    /* Traverse the kd-tree checking hypersphere-hyperrectangle intersections to discard regions of the space */
    void traverse (Box& k, int t) const;
};

/**
 * Traverse the binary tree while discarding regions of space
 * with hypersphere-hyperrectange intersections.
 */

template <typename Tree, typename Rule>
void   
SingleTreeTraversal<Tree, Rule>::traverse (Box& s, int t) const {
  int tree_size =  SingleTreeTraversal<Tree, Rule>::tree.data.size();
  if (!rules.prune_subtree(s, t)) {
    /* Evaluate the leaf node */
    if (s.is_leaf()) {                                              
      rules.base_case(s, t);
     } 
     else {
       bool first_branch = rules.visit(s, t);
       int child_id1 = s.child + first_branch;
       int child_id2 = s.child + !(first_branch);
       if (s.size() > (tree_size / omp_get_num_threads())){
        #pragma omp task  default(shared) 
          traverse (tree.node_data[child_id1], t);
          traverse (tree.node_data[child_id2], t);
        #pragma omp taskwait
       } 
       else {
         traverse (tree.node_data[child_id1], t);
         traverse (tree.node_data[child_id2], t);
       }
     }
  }
  /* puting the center for the pruned branch */
  else { 
    rules.centroid_case(s, t);
  }
}
#endif

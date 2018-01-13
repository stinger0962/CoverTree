#ifndef _DUAL_TREE_TRAVERSAL_H_
#define _DUAL_TREE_TRAVERSAL_H_
#define prune_type unsigned long long

#include <cfloat>
#include <omp.h>
#include <NBodyLower.h>
#include "IRNode.h"
#include "Binary_tree.h"

class Nbody::NBodyLower;

//#include <cilk/cilk.h>
//#include <cilk/cilk_api.h>
 
template <typename Tree, typename Base, typename pruneGenerator>
class DualTreeTraversal {
 public:
  typedef typename Tree::NodeTree Box;
  Tree& src_tree;
  Tree& trg_tree;
  Base& base;
  pruneGenerator prune;
  real_t temp = -1;
  prune_type& num_prunes;
  int task_level;

  /* Convenience constructor */
  DualTreeTraversal(Tree& st, Tree& tt, Base&  N, pruneGenerator& pg, prune_type& num,int levelT = 10)
      : src_tree(st), trg_tree(tt), base(N), prune(pg), num_prunes(num), task_level(levelT){
     
  }

  /* Traverse the tree checking hypersphere-hyperrectangle intersections to discard regions of the space */
  void traverse (Box& s, Box& t) ;
 private:
  /* Implementation helper function */
  void traverse_impl(Box& s, Box& t, unsigned level) ;
};

/**
 * Traverse the binary tree while discarding regions of space
 * with hypersphere-hyperrectangeintersections.
 */
template <typename Tree, typename Base, typename pruneGenerator>
void
DualTreeTraversal<Tree, Base, pruneGenerator>::traverse (Box& s, Box& t)  {
#pragma omp parallel
  {
#pragma omp single nowait
    {
       traverse_impl(s, t);
    }
    
  }
}

template <typename Tree, typename Base, typename pruneGenerator>
void
DualTreeTraversal<Tree, Base, pruneGenerator>::traverse_impl(Box& s, Box& t,
                                             unsigned level = 0)  {

if (!prune.prune_subtree(s, t)) {
 // if(true){
  /* If both source and target nodes are leaves, evaluate the base case */
    if (s.is_leaf() && t.is_leaf()) {
      base.base_case(s,t);
      // base.adjustSetIterationBounds(1 , s.begin() , s.end());
      // base.adjustSetIterationBounds(0 , t.begin() , t.end());
      // base.execute();
       // temp = (base.getOutput())[0][0];
      temp = base.get_temp();
      prune.set_temp(temp);

    }

    /* Recurse down the target node. Recursion order does not matter. */
    else if (s.is_leaf() && !t.is_leaf()) {
      {
         traverse_impl(s, trg_tree.node_data[t.child], level+1);
      }
      {
         traverse_impl(s, trg_tree.node_data[t.child+1], level+1);
      }
    }

    /* Recurse down the source node. Recursion order does matter. */
    else if (!s.is_leaf() && t.is_leaf()) {
      /* Check which branch should be explored first */
     bool first_branch = prune.visit(s,t); 
      /* Traverse the first (manhattan nearest) branch */
      // int child_id = s.child + first_branch;
      // traverse_impl(src_tree.node_data[child_id], t, level+1);

      /* Check if the further branch can be pruned */
      // child_id = s.child + !(first_branch);
      // traverse_impl(src_tree.node_data[child_id], t, level+1);
    }

    else if (level <= 10) {
    //else if (level <= task_level) {
#if 0
      cilk_spawn 
         traverse_impl(s,
                      trg_tree.node_data[t.child], level+1);
         traverse_impl(s,
                      trg_tree.node_data[t.child+1], level+1);
      cilk_sync;
#endif
// #pragma omp task default(shared) untied //if (level <= 5)
      {
         traverse_impl(s,trg_tree.node_data[t.child], level+1);
      }

#pragma omp task default(shared) untied //if (level <= 5)

      {
         traverse_impl(s,trg_tree.node_data[t.child+1], level+1);
      }
    }
    /* We have to recurse down both source and target nodes. */
    else {
      /* Target descent order does not matter, we choose left branch first */
      {
        int trg_child_id = t.child;

        /* Check which branch should be explored first */
        bool first_branch = prune.visit(s, trg_tree.node_data[trg_child_id] );

        /* Traverse the first (manhattan nearest) branch */
        int src_child_id = s.child + first_branch;
        traverse_impl(src_tree.node_data[src_child_id],trg_tree.node_data[trg_child_id], level+1);

        /* Check if the further branch can be pruned */
        src_child_id = s.child + !(first_branch);
        traverse_impl(src_tree.node_data[src_child_id],trg_tree.node_data[trg_child_id], level+1);
      }

      {
        /* Now recurse down the right target branch */
        int trg_child_id = t.child + 1;

        /* Check which branch should be explored first */
        bool first_branch = prune.visit(s, trg_tree.node_data[trg_child_id]);
        /* Traverse the first (manhattan nearest) branch */
        int src_child_id = s.child + first_branch;
        traverse_impl(src_tree.node_data[src_child_id],trg_tree.node_data[trg_child_id], level+1);

        /* Traverse the second (manhatten farthest) branch */
        src_child_id = s.child + !(first_branch);
        traverse_impl(src_tree.node_data[src_child_id],trg_tree.node_data[trg_child_id], level+1);
      }
    }
  } else {

      prune.centroid_case(s, t);  
      base.adjustPartitionCounter(1 , t.begin() , t.end() , s.size());

   }
}

#endif

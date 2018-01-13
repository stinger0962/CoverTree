#include <limits.h>
#include <stdint.h>
#include <iostream>
#include <assert.h>
#include <algorithm>
#include "Cover_tree.h"
#include "utils.h"

using namespace std;


template<typename Boundary>
  int
CTree<Boundary>::compute_min_cover_level(real_t d) {
  int cover_level = ceil(log2(d) ); 
  return cover_level;
}

template<typename Boundary>
  real_t 
CTree<Boundary>::compute_cover_radius_from_level(int l) {
  if (l == MAX_LEVEL) 
    return numeric_limits<real_t>::max(); /* switch to <double> if not work */
  real_t d = pow(2, l);
  return d;
}
/* Determine if a root at level l covers distance d
 */
template<typename Boundary>
  bool
CTree<Boundary>::is_covered(int l, real_t d) {
  return compute_cover_radius_from_level(l) >= d;

}


/* 
 * Construct a cover tree on a dummy node with all the data points
 * Then flatten the tree nodes and permute points
 */
template<typename Boundary>
  int
CTree<Boundary>::build_covertree() {
  /* Construct cover tree 
   * Dummy is not a real node in the cover tree,
   * It stores info of root node, and all points in the tree
   */
  NodeTree dummy(data[0], 0, MAX_LEVEL);
  dummy.add_points(data);
  p_root = construct_tree(&dummy, data[0].size());
  compute_max_dist();

  /* Permute data points */
  int j = 0;
  permute_data(p_root, j);

  /* Flatten nodes */
  p_root->nodeid = 0;
  node_data.push_back(*p_root);
  flatten_tree(p_root);  
  set_child(p_root);  
  
  /* Compute centroid for each tree node */
  // centroids(this->root());

  /* traverse() performs various validation,
   * and store level message in the map
   */
  /* 
  map<int, int> map;
  traverse(p_root, p_root->level(), map);
  for (auto pair: map)
    cout << pair.first << " -> " << pair.second << endl;
  */
  return 1;
}


/* Batch construction of a cover tree
 * which has root at p_node
 * and covers points stored in p_node.points
 */
template<typename Boundary>
typename CTree<Boundary>::NodeTree*
CTree<Boundary>::construct_tree(NodeTree* p_node, int dim) {

  
  /* Retrieve the point and its level */
  const Point& p = p_node->point();
  int index = p_node->point_index_;
  int node_level = p_node->level();
  Boundary b(p.size());
  /* temp used to retrieve points_ right before insertion */
  list<int> temp = p_node->points_;
  /* points_ will be consumed during the construction */
  list<int>& points = p_node->points_;

  // cout << "Depth: " << i << ", Level: " << node_level << endl; 

  pair<real_t, int> farthest_info = find_farthest_point(p, points);  
  real_t max_dist = farthest_info.first;
  
  /* Base case 
   * Stop if points only contains: self point, or duplicate points
   */
  if (points.size() <= 1 || max_dist == 0) {
    NodeTree* p_retval = new NodeTree(p, index, node_level, b);
    p_retval->points_ = temp;     
    return p_retval;
  }

  /* Assign proper level for root */
  if (node_level == MAX_LEVEL)
  	node_level = compute_min_cover_level(max_dist);
  /* Assign level for child, root_level means it is root of a subtree */
  // int root_level = node_level - 1;
  int root_level = compute_min_cover_level(max_dist) - 1;

  /* 
  cout << "NodeTrees: " << node_data.size() << " | ";  
  cout << "Point: " << p << "Level: " << node_level << " | ";
  cout << "#p: " << points.size() << ", " << points.front() << points.back() << " | ";
  cout << "FarPt: " << farthest_info.second << " Dist: " << max_dist;
  cout << " | child level: " << root_level << endl;
  */

  NodeTree n(p, index, node_level); // only a PH to call compute_distance

  /* Find points that can be inserted at root_level */ 
  list<int> root_points = find_roots(n, index, root_level, points);
  
  /* Construct root nodes */
  vector<NodeTree> roots;
  
  /* On stack: these are not real nodes */
  for (int r: root_points) 
    roots.push_back(NodeTree(data[r], r, root_level) );

  /* Assign points to each child node, also update lo and hi vectors*/
  auto lo_hi = assign_points_to_roots(&roots, &points, dim);

  NodeTree* p_retval = new NodeTree(p, index, node_level, b);  

  /* Update lo and hi for node's bounds */
  for (int i = 0; i < dim; i++) {
    p_retval->bounds.lo[i] = lo_hi.first[i];
    p_retval->bounds.hi[i] = lo_hi.second[i];
  }

  /* Build a sub cover tree based on each root, and points assigned to root
   * These roots are child nodes of p_retval
   */
  for (NodeTree& root: roots) {
    /* Skip if root only covers itself */
    if (root.distance_between_points(root.point(), p) == 0 
        && root.points_.size() == 0)
    /* Remove dependence on metric.h */
    /* if (metric.compute_distance(p, root.point() ) == 0 
        && root.points_.size() == 0) */
      continue;
  
    NodeTree* p_child = construct_tree(&root, dim);
    p_retval->add_child(p_child);
  }


  p_retval->points_ = temp;
  return p_retval;
}

 /* Find the point in pts that is farthest to p
  * Return both the found point and its distance to p
  */
template<typename Boundary>
  pair<real_t, int> 
CTree<Boundary>::find_farthest_point(const Point& p, list<int>& pts) {
  NodeTree n(p, 0, 0);
  real_t far_dist = 0;
  auto p_far = pts.begin();
  for (auto it = pts.begin(); it != pts.end(); ++it) {
    real_t cur_dist = n.distance_between_points(p, data[*it]);
    if (cur_dist > far_dist) {
      p_far = it;
      far_dist = cur_dist;
    }
  }  
  return make_pair(far_dist, *p_far);
}

/* Assign points whose indexes are stored in p_source to a node r in p_roots
 * such that r.point_ has the closest distance to p
 * (force nearest ancestor)
 */
template<typename Boundary>
  pair<Vec, Vec> 
CTree<Boundary>::assign_points_to_roots(vector<NodeTree>* p_roots, list<int>* p_source, int dim) {
  Vec lo(dim);
  Vec hi(dim);
  for (int i = 0; i < dim; i++) {
    lo[i] = DBL_MAX;
    hi[i] = -DBL_MAX;
  }
  if (p_source->empty() )
    return make_pair(lo, hi);
  while (!p_source->empty() ) {
    /* Retrieve one index from p_source */
    int p = p_source->front();

    /* Update bounds if necessary */
    for (int i = 0; i < dim; i++) {
      hi[i] = data[p][i] > hi[i] ? data[p][i] : hi[i];
      lo[i] = data[p][i] < lo[i] ? data[p][i] : lo[i]; 
    }
 
    p_source->pop_front();

    /* Iterate through p_roots, find closest distance */
    auto it = p_roots->begin();
    auto last = p_roots->end();
    auto nearest = it;
    
    real_t min_dist = it->distance_between_points(it->point(), data[p]);
    
    for (++it; it != last; ++it) {
      const Point& q = it->point();
      real_t cur_dist = it->distance_between_points(q, data[p]);
      if (cur_dist < min_dist) {
        nearest = it;
        min_dist = cur_dist;
      }
    }

    list<int>& destination = nearest->points_;
    destination.push_back(p);
  }

  return make_pair(lo, hi);
}


 /*
  * Find all points at level l whose parent is p
  * Point r is the latest added child 
  * Next child must satisfy the separation invariant with r
  */
template<typename Boundary>
  list<int> 
CTree<Boundary>::find_roots(NodeTree& n, int p, int l, list<int>* p_points, int r) {
  /* Disqualify points that cannot be inserted as root */
  auto cannot_be_root = [this, &n, r, l] (int q) {
    real_t dist = n.distance_between_points(data[r], data[q]);
    return dist <= compute_cover_radius_from_level(l);
  };
  
  /* Remove unsatisfied points from list */
  p_points->remove_if(cannot_be_root);

  if (p_points->empty() )
    return list<int>();

  pair<real_t, int> farthest_pair = find_farthest_point(data[p], *p_points);
  /* Do nothing if p_points only contains duplicate points
  if (farthest_pair.first == 0)
    return list<Point>(); 
  */

  int s = farthest_pair.second;  
  list<int> roots = find_roots(n, p, l, p_points, s);
  roots.push_back(s);
  
  return roots; 
}


/*
 * For each node in the tree,
 * Compute the max distance to any of its descendants
 * This distance is the radius of its bound
 * Breadth-first-search
 */
template<typename Boundary>
  void
CTree<Boundary>::compute_max_dist() {
  vector<CTree<Boundary>::NodeTree*> travel;
  vector<CTree<Boundary>::NodeTree*> active;

  CTree<Boundary>::NodeTree* cur = p_root;
  
  p_root->maxdist = 0.0;
  travel.push_back(p_root);

  while (travel.size() > 0) {
    cur = travel.back();
    if (cur->maxdist <= 0.0) {
      while (cur->children.size() > 0) {
        active.push_back(cur);
        for (int i = 0; i < cur->children.size(); ++i) {
          cur->children[i]->maxdist = 0.0;
          travel.push_back(cur->children[i]);
        }
        cur = cur->children[cur->children.size()-1];
      }
    }
    else
      active.pop_back();
    
    for (const auto& n: active) {
      n->maxdist = max(n->maxdist, n->distance_between_points( n->point(), cur->point()) );
      n->bounds.update_bounding_box(n->point(), n->maxdist);
    }
    travel.pop_back(); 
  }
}   

/* Flatten data points, add points to data_perm in an order such that 
 * Points covered by the same node will be stored in the consecutive order
 * Update and return each node's beginning index in the data_perm
 * j increments from 0, serving as index in data_perm
 */
template<typename Boundary>
  int
CTree<Boundary>::permute_data(NodeTree* node, int& j) {
  // cout << endl << "@" << j << ": " << *node;
  // cout << "max_dist " << node->maxdist << endl;
  node->num = node->points_.size();
  /* Leaf contains one distinct point, or multiple duplicate points */
  if (node->is_leaf()) {
    int retval = j;
    /* Permute data so that: data_perm[j] = data[index[j] ]   */
    for (int p: node->points_) {
      index.push_back(p); // or index[j] = p
      data_perm[j++] = data[p];
    }
    node->beg = retval;
    return retval; 
  }
  /* Non-leaf node */
  int retval = 0;
  for (int i = 0; i < node->children.size(); i++) {
    if (i == 0) {
      retval = permute_data(node->children[i], j);
      node->beg = retval;
    }
    else
      permute_data(node->children[i], j);
  }
  return retval;
}

/* Flatten nodes, add nodes to node_data in an order such that  
 * child nodes of a node will be stored in the consecutive order
 */
template<typename Boundary>
  void
CTree<Boundary>::flatten_tree(NodeTree* node) {
  
  if (node->is_leaf()) {
    return; 
  }
  for (int i = 0; i < node->children.size(); i++) {
    node->children[i]->nodeid = node_data.size();
    node_data.push_back(*(node->children[i]) );
  }
  for (int i = 0; i < node->children.size(); i++) {
    flatten_tree(node->children[i]);
  }
  return;
}  

/* Set child for each node, both in tree and array(node_data)
 * child is the index of a node's first child node, -1 for a leaf node
 */
template<typename Boundary>
  void
CTree<Boundary>::set_child(NodeTree* node) {
  
  if (node->is_leaf()) {
    node_data[node->nodeid].child = -1;
    node->child = -1;

    return; 
  }
  node_data[node->nodeid].child = node->children[0]->nodeid; 
  node->child = node->children[0]->nodeid;
  for (int i = 0; i < node->children.size(); i++) {
    set_child(node->children[i]);
  }
  return;
}  

/* Various validation can be performed in this traversal */
template<typename Boundary>
  void
CTree<Boundary>::traverse(NodeTree* node, int lvl, map<int, int>& map) {
  // cout << "#" << node->index() << ": " << centers[node->index()] << endl;

  /* 
  cout << "#" << node->index() << " [lo, hi]: ";
  for (int i = 0; i < node->bounds.lo.size(); i++)
    cout << "[" << node->bounds.lo[i] << "," << node->bounds.hi[i] << "], ";
  cout << endl;
  */
  
  /*
  real_t max_dist_correct = 0;
  for (int i = node->begin(); i < node->end(); i++) {
    max_dist_correct = max (max_dist_correct, 
                  node->distance_between_points(node->point(), data_perm[i]));
  }
  if (node->bounds.radius != max_dist_correct)
    cout << "Nodeid: " << node->nodeid
         << ", Bound Radius: " 
         << node->bounds.radius << ", Correct Maxdist: " 
         << max_dist_correct << endl;
  */
  map[lvl]++;
  if (node->is_leaf()) {
    return;
  }
  for (int i = 0; i < node->children.size(); i++) {
    traverse(node->children[i], lvl - 1, map);
  }
}

/* Batch computation of center of mass for all the nodes in the tree
 * Accept tree's root as parameter
 */
template<typename Boundary>
  void
CTree<Boundary>::centroids(NodeTree& s) {
  int dim = data[0].size();
  vector<real_t> sum_of_dim(dim, 0); 
  /* Here we need to calculate centroids for all nodes, so we need true leaf */
  if (s.is_leaf()) {
    for (long i = s.begin(); i < s.end(); i++) {
      for (int d = 0; d < dim; d++) {
        assert (i < data_perm.size());
        assert (d < sum_of_dim.size());
        assert (d < data_perm[i].size());
        sum_of_dim[d] += data_perm[i][d];
      }      
    }
    for (int d = 0; d < dim; d++)
      centers[s.index()][d] = sum_of_dim[d] / s.size();
  }
  else {
    for (int i = 0; i < s.num_child(); i++) {
      centroids(this->node_data[s.child+i]);
      // cerr << centers[s.child+i] << endl;
    }
    for (int d = 0; d < dim; d++) {
      sum_of_dim[d] = 0;  
      for (int i = 0; i < s.num_child(); i++)
        sum_of_dim[d] += centers[s.child+i][d] 
                         * this->node_data[s.child+i].size();
      centers[s.index()][d] = sum_of_dim[d] / s.size();
    }
  }
}






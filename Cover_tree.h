/** 
 * Defines a cover tree
 */

#ifndef INC_COVER_TREE_H
#define INC_COVER_TREE_H

#include <list>
#include <vector>
#include <map>
#include <math.h>
#include "Metric.h"
#include "kNN_vector.h"
#include "range.h"
#include "vector_of_array.hpp"

#if !defined (USE_FLOAT)
typedef double real_t;
#else
typedef float real_t;
#endif

typedef vector_of_array<real_t> Points_t;
typedef typename Points_t::array_type Point;
typedef fixed_vector<real_t> Vec;
typedef vector_of_array<real_t> Mat;

using namespace std;

/* Wrapping class for a cover tree */
template <typename Boundary>
class CTree {

public:
  
  
  friend std::ostream& operator<< (std::ostream& os, const CTree& tree) {
    for (int i = 0; i < tree.nodes(); i++)
      os << tree.node_data[i] << endl;
    return os;
  }  
  class NodeTree;
  /* Largest double-precision float is 2^1024 
   * and the smallest positive value 2^-1074 
   */
private:
  static const int MAX_LEVEL = 1024;
  static const int MIN_LEVEL = -1074;
  
  /* Virtual leaf size 
   * In traversal, when a node.num_points() <= leaf_size, 
   * we treat it as a leaf node */
  // static const int leaf_size = 64;

  /* Compute the lowest level required for a root of cover tree
   * such that the root covers the distance of d
   */
  int compute_min_cover_level(real_t d);
 
  /* Compute the cover distance from a root at level l */
  real_t compute_cover_radius_from_level(int l);

  /* Determine if a root at level l covers the distance d*/
  bool is_covered(int l, real_t d);
  
  /* Find the point in pts that is farthest to p 
   * Return both the found point and its distance to p 
   */
  pair<real_t, int> find_farthest_point(const Point& p, list<int>& pts);

  /* 
   * Find all points that may be inserted as root as level l
   */ 
  list<int> find_roots(NodeTree& n, int p, int l, list<int>* p_points, int r);
  
  /* Wrapper of find_roots */
  list<int> find_roots(NodeTree& n, int p, int l, list<int> pts) {
    list<int> roots = find_roots(n, p, l, &pts, p);
    roots.push_front(p);
    return roots;
  }
  /* Assign each point p in p_source to a node r in p_roots 
   * where r.point_ has the closest distance to p 
   * (force nearest ancestor) 
   */
  pair<Vec, Vec> assign_points_to_roots(vector<NodeTree>* p_roots, list<int>* p_source, int dim);

  /* Batch construction of a cover tree */
  NodeTree* construct_tree(NodeTree* dummy, int i);

  /* Compute maxdist for each node */
  void compute_max_dist();
 
  /* Flatten data points (generate data_perm), and update beg for each node 
   * @j: increments from 0, serving as index in data_perm
   */
  int permute_data(NodeTree* node, int& j);

  /* Flatten nodes, update child for each node */
  void flatten_tree(NodeTree* node);  

  /* Set child for each node, which is the index of its first child node */
  void set_child(NodeTree* node);

  /* For test */
  void traverse(NodeTree* node, int lvl, map<int, int>& map);


public:
  class NodeTree {

  friend class CTree;

  public:
    friend std::ostream& operator<< (std::ostream& os, const NodeTree& node) {
      os << "L: " << node.level() << " P: " << node.point() << endl;
      os << "Beg: " << node.beg << " Num: " << node.num << endl;
      os << "Maxdist: " << node.maxdist << endl;
      os << "Bounds, raidus: " << node.bounds.radius << " center: ";
      for (int i = 0; i < node.bounds.center.size(); i++)
        os << node.bounds.center[i] << " ";
      os << endl << "Nid: " << node.index() << " 1stC: " << node.child << endl;
      for (int i = 0; i < node.children.size(); ++i)
      os << node.children[i]->level() 
          << " " << node.children[i]->point() << "; " << endl;
      // os << ";" << std::endl;
      return os;
    }

      
    NodeTree(Point p, int i, int l) : point_(p), point_index_(i), level_(l) {}
    NodeTree(Point p, int i, int l, Boundary b) : point_(p), point_index_(i), level_(l), bounds(b) { bounds.update_bounding_box(p, 0); }
    // NodeTree() {}

    const Point& point() const { return point_; }
    const int level() const { return level_; }
    
    Boundary bounds; 
    
    /* Related to child nodes */ 
    vector<NodeTree*> children; 
    void add_child(NodeTree* p_child) { children.push_back(p_child); }  
  
    int child; /* Index of first child node in node_data */

    int num_child() {return children.size(); }
 
    bool is_leaf() {return children.size() == 0 || num_points() <= getenv__leaf_size(); }
 
    /* Public interface of flattened data */
    int begin() const {
      return beg;
    }
 
    int end() const {
      return beg + num;
    }

    int size() const {
      return num;
    }

    int index() const {
      return nodeid;
    }

    /* Actual max distance to any of its descendant, 
     * A better upperbound than 2^level 
     * */    
    real_t maxdist;     
 
    /* Public interface for bounds */
    real_t min_distance(const Point& point) {
      return bounds.min_distance(point);
    }
    
    real_t min_distance(const NodeTree& box) {
      return bounds.min_distance(box);
    }
    
    real_t max_distance(const Point& point) {
      return bounds.max_distance(point);
    }
    
    real_t max_distance(const NodeTree& box) {
      return bounds.max_distance(box);
    }

    Range range_distance(const Point& point) {
      return bounds.range_distance(point);
    }
    
    Range range_distance(const NodeTree& box) {
      return bounds.range_distance(box);
    }
    
    real_t distance_between_points(const Point& p1, const Point& p2) {
        return bounds.distance_between_points(p1, p2);
    }

    Vec center() {
      return bounds.center;
    }
   
    Vec min_point(const Point& point) {
      return bounds.min_point(point);
    }    
 
    Vec max_point(const Point& point) {
      return bounds.max_point(point);
    }
    
    Vec min_point(NodeTree& box) {
      return bounds.min_point(box);
    }
    
    Vec max_point(NodeTree& box) {
      return bounds.max_point(box);
    }
    
    inline int num_points() {
      return points_.size();
    }

    
 
  private:
    Point point_;
    int level_;
  
    /* Index of associated point in the original data[] */
    int point_index_;
 
    /* Index of this node in the vector node_data */
    int nodeid;
    
    /* Starting index of points covered by this node in the vector data_perm */
    int beg;
    
    /* Number of points covered by this node */
    int num;
 
    // EuclideanMetric<Point, Point> metric;
    
    /* This collection is useful in the tree construction 
     * It stores indexes of all points in data[] covered by a sub cover tree
     * It is also useful in geenrating data_perm[]
     */
    list<int> points_; 
    
    /* Add all points from ps to points_*/
    void add_points(Points_t& ps) {
      for (int i = 0; i < ps.size(); i++)
        points_.push_back(i);
    }
  };

  /* All nodes, in post-order */  
  vector<NodeTree> node_data;
  
  /* Pointer to the root node */  
  NodeTree* p_root;

  /* Original dataset */
  Points_t& data;
 
  /* Permuted dataset*/
  Points_t& data_perm;

  /* Center of mass of each node */
  Points_t centers;
 
  /* Tree constructor: initialize the datasets in the tree */
  CTree(Points_t& ps, Points_t& perm) : data(ps), data_perm(perm), centers(perm) {
    // leaf_size = getenv__leaf_size();
    // cerr << leaf_size << endl;
  };

  /* Wrapper of the batch construction function */  
  int build_covertree();
 
  /* These indices are used to pair points in data with points in data_perm 
   * data_perm[i] is data[index[i] ]
   */
  vector<long> index;
 
  /* Return the number of nodes in the tree */
  inline long nodes() const {
    return node_data.size();
  }

  /* Return the number of points in the tree */
  inline long points() const {
    return data.size();
  }
  
  /* Return the root node from the tree */
  NodeTree& root() {
    return *p_root;
  }

  /* Batch computation of center of mass for all the nodes in the tree 
   * Accept root node as parameter 
   */
  void centroids(NodeTree& s);
  
};

#include "Cover_tree.cpp"

#endif

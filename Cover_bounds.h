/**
 * Defines a bounding ball for cover tree
 */

#ifndef INC_COVER_BOUNDS_H
#define INC_COVER_BOUNDS_H

#include "fixed_vector.h"
#include "Metric.h"

#include "Cover_tree.h"

class Cover;

typedef CTree<Cover>::NodeTree Box;

/* Hypersphere bounds for a specified point in the cover tree
 * Radius of the sphere is determined by scale of the point 
 */
class Cover {
public:
  int dim;
  real_t radius; 
  Vec center;

  Vec lo;
  Vec hi;  
  Vec hi_i; // not used anymore
  Vec lo_i; // not used anymore

  Cover() {}

  /* Create the bound with a specified dimension */
  Cover(unsigned int d) : dim(d), radius(0), hi_i(d), lo_i(d), 
                          lo(d), hi(d), center(d) {
    for (int i = 0; i < dim; i++) 
      center[i] = 0;
  }
  /* Set up properties of a bound object
   * bound object is initialized during tree construction
   * and updated after maxdist is calculated
   * */
  void update_bounding_box (const Point& data, real_t maxdist);

  /* Computes the distance between two points */
  real_t distance_between_points (const Point& point1, const Point& point2);

  real_t distance_between_points (const Point& point1, const Point& point2) const;

  /* Computes the point-to-point squared distance */
  real_t distance_to_center (const Point& point);
  
  /* Computes minimum box-to-point squared distance */
  real_t min_distance (const Point& point);

  /* Computes minimum box-to-box squared distance */
  real_t min_distance (const Box& box);
   
  /* Computes maximum box-to-point squared distance */
  real_t max_distance (const Point& point);

  /* Computes maximum box-to-box squared distance */
  real_t max_distance (const Box& box);
  
 
  /* Return coordinates of the point on the bound with min distance to a point */
  Vec min_point (const Point& point);

  /* Return coordiantes of the point on the bound with min distance to a box */
  Vec min_point (Box& box);

  /* Return coordinates of the point on the bound with max distance to a point */
  Vec max_point (const Point& point);

  /* Return coordiantes of the point on the bound with max distance to a box */
  Vec max_point (Box& box);

  /* Computes min and max box-to-point squared distance */
  Range range_distance (const Point& point);

  /* Computes min and max box-to-box squared distance */
  Range range_distance (const Box& box);

 /* Return the diameter of he hypersphere */
  real_t width (unsigned int d) { return 2 * radius; }
};


/* Computes the bounds for a point with a specified scale */
void
Cover::update_bounding_box (const Point& data, real_t maxdist) {
  radius = maxdist;
  for (int i = 0; i < dim; i++)
    center[i] = data[i];
}
  
/* Computes the distance between two points */
real_t 
Cover::distance_between_points (const Point& point1, const Point& point2) {
  EuclideanMetric<Point, Point> metric;
  real_t dist = metric.compute_distance(point1, point2);
  return dist;
}

real_t 
Cover::distance_between_points (const Point& point1, const Point& point2) const {
  EuclideanMetric<Point, Point> metric;
  real_t dist = metric.compute_distance(point1, point2);
  return dist;
}

/* Computes the point-to-point squared distance */
real_t
Cover::distance_to_center (const Point& point) { 
  EuclideanMetric<Vec, Point> metric;
  real_t dist = metric.compute_distance(center, point);
  return dist * dist;

}
/* Computes minimum box-to-point squared distance */
real_t
Cover::min_distance (const Point& point) {
  EuclideanMetric<Vec, Point> metric;
  real_t dist = metric.compute_distance(center, point) - radius;
  dist = (dist + fabs(dist) ) / 2;
  return dist * dist;
}

/* Computes minimum box-to-box squared distance */
real_t 
Cover::min_distance (const Box& box) { 
  EuclideanMetric<Vec, Vec> metric;
  real_t dist = metric.compute_distance(center, box.bounds.center) - radius - box.bounds.radius;
  dist = (dist + fabs(dist) ) / 2;
  return dist * dist;
}

/* Computes maximum box-to-point squared distance */
real_t
Cover::max_distance (const Point& point) {
  EuclideanMetric<Vec, Point> metric;
  real_t dist = metric.compute_distance(center, point) + radius;
  return dist * dist;
}

/* Computes maximum box-to-box squared distance */
real_t 
Cover::max_distance (const Box& box) { 
  EuclideanMetric<Vec, Vec> metric;
  real_t dist = metric.compute_distance(center, box.bounds.center) + radius + box.bounds.radius;
  return dist * dist;
}

/* Return coordinates of the point on the bound with min distance to a point */

Vec 
Cover::min_point (const Point& point) {
  Vec p(dim);
  EuclideanMetric<Vec, Point> metric;
  for (int i = 0; i < dim; i++) {
    if (metric.compute_distance(center, point) > radius)
      p[i] = center[i] - radius;
  }
  return p;
}


/* Return coordiantes of the point on the bound with min distance to a box */

Vec 
Cover::min_point (Box& box) {
  Vec p(dim);
  EuclideanMetric<Vec, Vec> metric;
  for (int i = 0; i < dim; i++) {
    if (metric.compute_distance(center, box.center() ) > (radius + box.bounds.radius) )
      p[i] = center[i] - radius - box.bounds.radius;
  }
  return p;
}


/* Return coordinates of the point on the bound with max distance to a point */

Vec 
Cover::max_point (const Point& point) {
  Vec p(dim);
  EuclideanMetric<Vec, Point> metric;
  for (int i = 0; i < dim; i++) {
    if (metric.compute_distance(center, point) > radius)
      p[i] = center[i] + radius;
  }
  return p;
}


/* Return coordiantes of the point on the bound with max distance to a box */

Vec 
Cover::max_point (Box& box) { 
  Vec p(dim);
  EuclideanMetric<Vec, Vec> metric;
  for (int i = 0; i < dim; i++) {
    if (metric.compute_distance(center, box.center() ) > (radius + box.bounds.radius) )
      p[i] = center[i] + radius + box.bounds.radius;
  }
  return p;
}


/* Computes min and max box-to-point squared distance */
Range 
Cover::range_distance (const Point& point) {
  EuclideanMetric<Vec, Point> metric;
  real_t dist = metric.compute_distance(center, point);
  real_t min_dist = pow( ( (dist - radius) + fabs(dist - radius) ) / 2, 2);
  real_t max_dist = pow(dist + radius, 2);
  return Range(min_dist, max_dist);
}

/* Computes min and max box-to-box squared distance */
Range 
Cover::range_distance (const Box& box) {
  EuclideanMetric<Vec, Vec> metric;
  real_t dist = metric.compute_distance(center, box.bounds.center);
  real_t sum_radius = radius + box.bounds.radius;
  real_t min_dist = pow( ( (dist - sum_radius) + fabs(dist - sum_radius) ) / 2, 2);
  real_t max_dist = pow(dist + sum_radius, 2);
  return Range(min_dist, max_dist);
}

#endif

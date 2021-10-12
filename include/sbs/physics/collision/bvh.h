// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2009 Ilya Baran <ibaran@mit.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

/**
 * Modified by Quoc-Minh Ton-That
 */

#ifndef SBS_EIGEN_BVH_MODULE_H
#define SBS_EIGEN_BVH_MODULE_H

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/StdVector>
#include <algorithm>
#include <queue>

namespace Eigen {

/**
  * \defgroup BVH_Module BVH module
  * \brief This module provides generic bounding volume hierarchy algorithms
  * and reference tree implementations.
  *
  *
  * \code
  * #include <unsupported/Eigen/BVH>
  * \endcode
  *
  * A bounding volume hierarchy (BVH) can accelerate many geometric queries.  This module provides a
  generic implementation
  * of the two basic algorithms over a BVH: intersection of a query object against all objects in
  the hierarchy and minimization
  * of a function over the objects in the hierarchy.  It also provides intersection and minimization
  over a cartesian product of
  * two BVH's.  A BVH accelerates intersection by using the fact that if a query object does not
  intersect a volume, then it cannot
  * intersect any object contained in that volume.  Similarly, a BVH accelerates minimization
  because the minimum of a function
  * over a volume is no greater than the minimum of a function over any object contained in it.
  *
  * Some sample queries that can be written in terms of intersection are:
  *   - Determine all points where a ray intersects a triangle mesh
  *   - Given a set of points, determine which are contained in a query sphere
  *   - Given a set of spheres, determine which contain the query point
  *   - Given a set of disks, determine if any is completely contained in a query rectangle
  (represent each 2D disk as a point \f$(x,y,r)\f$
  *     in 3D and represent the rectangle as a pyramid based on the original rectangle and shrinking
  in the \f$r\f$ direction)
  *   - Given a set of points, count how many pairs are \f$d\pm\epsilon\f$ apart (done by looking at
  the cartesian product of the set
  *     of points with itself)
  *
  * Some sample queries that can be written in terms of function minimization over a set of objects
  are:
  *   - Find the intersection between a ray and a triangle mesh closest to the ray origin (function
  is infinite off the ray)
  *   - Given a polyline and a query point, determine the closest point on the polyline to the query
  *   - Find the diameter of a point cloud (done by looking at the cartesian product and using
  negative distance as the function)
  *   - Determine how far two meshes are from colliding (this is also a cartesian product query)
  *
  * This implementation decouples the basic algorithms both from the type of hierarchy (and the
  types of the bounding volumes) and
  * from the particulars of the query.  To enable abstraction from the BVH, the BVH is required to
  implement a generic mechanism
  * for traversal.  To abstract from the query, the query is responsible for keeping track of
  results.
  *
  * To be used in the algorithms, a hierarchy must implement the following traversal mechanism (see
  KdBVH for a sample implementation): \code typedef Volume  //the type of bounding volume typedef
  Object  //the type of object in the hierarchy typedef Index   //a reference to a node in the
  hierarchy--typically an int or a pointer typedef VolumeIterator //an iterator type over node
  children--returns Index typedef ObjectIterator //an iterator over object (leaf) children--returns
  const Object & Index getRootIndex() const //returns the index of the hierarchy root const Volume
  &getVolume(Index index) const //returns the bounding volume of the node at given index void
  getChildren(Index index, VolumeIterator &outVBegin, VolumeIterator &outVEnd, ObjectIterator
  &outOBegin, ObjectIterator &outOEnd) const
      //getChildren takes a node index and makes [outVBegin, outVEnd) range over its node children
      //and [outOBegin, outOEnd) range over its object children
    \endcode
  *
  * To use the hierarchy, call BVIntersect or BVMinimize, passing it a BVH (or two, for cartesian
  product) and a minimizer or intersector.
  * For an intersection query on a single BVH, the intersector encapsulates the query and must
  provide two functions:
  * \code
      bool intersectVolume(const Volume &volume) //returns true if the query intersects the volume
      bool intersectObject(const Object &object) //returns true if the intersection search should
  terminate immediately \endcode
  * The guarantee that BVIntersect provides is that intersectObject will be called on every object
  whose bounding volume
  * intersects the query (but possibly on other objects too) unless the search is terminated
  prematurely.  It is the
  * responsibility of the intersectObject function to keep track of the results in whatever manner
  is appropriate.
  * The cartesian product intersection and the BVMinimize queries are similar--see their individual
  documentation.
  *
  * The following is a simple but complete example for how to use the BVH to accelerate the search
  for a closest red-blue point pair:
  * \include BVH_Example.cpp
  * Output: \verbinclude BVH_Example.out
  */
}

namespace Eigen {

namespace internal {

#ifndef EIGEN_PARSED_BY_DOXYGEN
template <typename BVH, typename Intersector>
bool intersect_helper(const BVH& tree, Intersector& intersector, typename BVH::Index root)
{
    typedef typename BVH::Index Index;
    typedef typename BVH::VolumeIterator VolIter;
    typedef typename BVH::ObjectIterator ObjIter;

    VolIter vBegin = VolIter(), vEnd = VolIter();
    ObjIter oBegin = ObjIter(), oEnd = ObjIter();

    std::vector<Index> todo(1, root);

    while (!todo.empty())
    {
        tree.getChildren(todo.back(), vBegin, vEnd, oBegin, oEnd);
        todo.pop_back();

        for (; vBegin != vEnd; ++vBegin) // go through child volumes
            if (intersector.intersectVolume(tree.getVolume(*vBegin)))
                todo.push_back(*vBegin);

        for (; oBegin != oEnd; ++oBegin) // go through child objects
            if (intersector.intersectObject(*oBegin))
                return true; // intersector said to stop query
    }
    return false;
}
#endif // not EIGEN_PARSED_BY_DOXYGEN

template <typename Volume1, typename Object1, typename Object2, typename Intersector>
struct intersector_helper1
{
    intersector_helper1(const Object2& inStored, Intersector& in)
        : stored(inStored), intersector(in)
    {
    }
    bool intersectVolume(const Volume1& vol)
    {
        return intersector.intersectVolumeObject(vol, stored);
    }
    bool intersectObject(const Object1& obj)
    {
        return intersector.intersectObjectObject(obj, stored);
    }
    Object2 stored;
    Intersector& intersector;

  private:
    intersector_helper1& operator=(const intersector_helper1&);
};

template <typename Volume2, typename Object2, typename Object1, typename Intersector>
struct intersector_helper2
{
    intersector_helper2(const Object1& inStored, Intersector& in)
        : stored(inStored), intersector(in)
    {
    }
    bool intersectVolume(const Volume2& vol)
    {
        return intersector.intersectObjectVolume(stored, vol);
    }
    bool intersectObject(const Object2& obj)
    {
        return intersector.intersectObjectObject(stored, obj);
    }
    Object1 stored;
    Intersector& intersector;

  private:
    intersector_helper2& operator=(const intersector_helper2&);
};

} // end namespace internal

/**  Given a BVH, runs the query encapsulated by \a intersector.
  *  The Intersector type must provide the following members: \code
     bool intersectVolume(const BVH::Volume &volume) //returns true if volume intersects the query
     bool intersectObject(const BVH::Object &object) //returns true if the search should terminate
  immediately \endcode
  */
template <typename BVH, typename Intersector>
void BVIntersect(const BVH& tree, Intersector& intersector)
{
    internal::intersect_helper(tree, intersector, tree.getRootIndex());
}

/**  Given two BVH's, runs the query on their Cartesian product encapsulated by \a intersector.
  *  The Intersector type must provide the following members: \code
     bool intersectVolumeVolume(const BVH1::Volume &v1, const BVH2::Volume &v2) //returns true if
  product of volumes intersects the query bool intersectVolumeObject(const BVH1::Volume &v1, const
  BVH2::Object &o2) //returns true if the volume-object product intersects the query bool
  intersectObjectVolume(const BVH1::Object &o1, const BVH2::Volume &v2) //returns true if the
  volume-object product intersects the query bool intersectObjectObject(const BVH1::Object &o1,
  const BVH2::Object &o2) //returns true if the search should terminate immediately \endcode
  */
template <typename BVH1, typename BVH2, typename Intersector>
void BVIntersect(
    const BVH1& tree1,
    const BVH2& tree2,
    Intersector& intersector) // TODO: tandem descent when it makes sense
{
    typedef typename BVH1::Index Index1;
    typedef typename BVH2::Index Index2;
    typedef internal::intersector_helper1<
        typename BVH1::Volume,
        typename BVH1::Object,
        typename BVH2::Object,
        Intersector>
        Helper1;
    typedef internal::intersector_helper2<
        typename BVH2::Volume,
        typename BVH2::Object,
        typename BVH1::Object,
        Intersector>
        Helper2;
    typedef typename BVH1::VolumeIterator VolIter1;
    typedef typename BVH1::ObjectIterator ObjIter1;
    typedef typename BVH2::VolumeIterator VolIter2;
    typedef typename BVH2::ObjectIterator ObjIter2;

    VolIter1 vBegin1 = VolIter1(), vEnd1 = VolIter1();
    ObjIter1 oBegin1 = ObjIter1(), oEnd1 = ObjIter1();
    VolIter2 vBegin2 = VolIter2(), vEnd2 = VolIter2(), vCur2 = VolIter2();
    ObjIter2 oBegin2 = ObjIter2(), oEnd2 = ObjIter2(), oCur2 = ObjIter2();

    std::vector<std::pair<Index1, Index2>> todo(
        1,
        std::make_pair(tree1.getRootIndex(), tree2.getRootIndex()));

    while (!todo.empty())
    {
        tree1.getChildren(todo.back().first, vBegin1, vEnd1, oBegin1, oEnd1);
        tree2.getChildren(todo.back().second, vBegin2, vEnd2, oBegin2, oEnd2);
        todo.pop_back();

        for (; vBegin1 != vEnd1; ++vBegin1)
        { // go through child volumes of first tree
            const typename BVH1::Volume& vol1 = tree1.getVolume(*vBegin1);
            for (vCur2 = vBegin2; vCur2 != vEnd2; ++vCur2)
            { // go through child volumes of second tree
                if (intersector.intersectVolumeVolume(vol1, tree2.getVolume(*vCur2)))
                    todo.push_back(std::make_pair(*vBegin1, *vCur2));
            }

            for (oCur2 = oBegin2; oCur2 != oEnd2; ++oCur2)
            { // go through child objects of second tree
                Helper1 helper(*oCur2, intersector);
                if (internal::intersect_helper(tree1, helper, *vBegin1))
                    return; // intersector said to stop query
            }
        }

        for (; oBegin1 != oEnd1; ++oBegin1)
        { // go through child objects of first tree
            for (vCur2 = vBegin2; vCur2 != vEnd2; ++vCur2)
            { // go through child volumes of second tree
                Helper2 helper(*oBegin1, intersector);
                if (internal::intersect_helper(tree2, helper, *vCur2))
                    return; // intersector said to stop query
            }

            for (oCur2 = oBegin2; oCur2 != oEnd2; ++oCur2)
            { // go through child objects of second tree
                if (intersector.intersectObjectObject(*oBegin1, *oCur2))
                    return; // intersector said to stop query
            }
        }
    }
}

namespace internal {

#ifndef EIGEN_PARSED_BY_DOXYGEN
template <typename BVH, typename Minimizer>
typename Minimizer::Scalar minimize_helper(
    const BVH& tree,
    Minimizer& minimizer,
    typename BVH::Index root,
    typename Minimizer::Scalar minimum)
{
    typedef typename Minimizer::Scalar Scalar;
    typedef typename BVH::Index Index;
    typedef std::pair<Scalar, Index> QueueElement; // first element is priority
    typedef typename BVH::VolumeIterator VolIter;
    typedef typename BVH::ObjectIterator ObjIter;

    VolIter vBegin = VolIter(), vEnd = VolIter();
    ObjIter oBegin = ObjIter(), oEnd = ObjIter();
    std::priority_queue<QueueElement, std::vector<QueueElement>, std::greater<QueueElement>>
        todo; // smallest is at the top

    todo.push(std::make_pair(Scalar(), root));

    while (!todo.empty())
    {
        tree.getChildren(todo.top().second, vBegin, vEnd, oBegin, oEnd);
        todo.pop();

        for (; oBegin != oEnd; ++oBegin) // go through child objects
            minimum = (std::min)(minimum, minimizer.minimumOnObject(*oBegin));

        for (; vBegin != vEnd; ++vBegin)
        { // go through child volumes
            Scalar val = minimizer.minimumOnVolume(tree.getVolume(*vBegin));
            if (val < minimum)
                todo.push(std::make_pair(val, *vBegin));
        }
    }

    return minimum;
}
#endif // not EIGEN_PARSED_BY_DOXYGEN

template <typename Volume1, typename Object1, typename Object2, typename Minimizer>
struct minimizer_helper1
{
    typedef typename Minimizer::Scalar Scalar;
    minimizer_helper1(const Object2& inStored, Minimizer& m) : stored(inStored), minimizer(m) {}
    Scalar minimumOnVolume(const Volume1& vol)
    {
        return minimizer.minimumOnVolumeObject(vol, stored);
    }
    Scalar minimumOnObject(const Object1& obj)
    {
        return minimizer.minimumOnObjectObject(obj, stored);
    }
    Object2 stored;
    Minimizer& minimizer;

  private:
    minimizer_helper1& operator=(const minimizer_helper1&);
};

template <typename Volume2, typename Object2, typename Object1, typename Minimizer>
struct minimizer_helper2
{
    typedef typename Minimizer::Scalar Scalar;
    minimizer_helper2(const Object1& inStored, Minimizer& m) : stored(inStored), minimizer(m) {}
    Scalar minimumOnVolume(const Volume2& vol)
    {
        return minimizer.minimumOnObjectVolume(stored, vol);
    }
    Scalar minimumOnObject(const Object2& obj)
    {
        return minimizer.minimumOnObjectObject(stored, obj);
    }
    Object1 stored;
    Minimizer& minimizer;

  private:
    minimizer_helper2& operator=(const minimizer_helper2&);
};

} // end namespace internal

/**  Given a BVH, runs the query encapsulated by \a minimizer.
  *  \returns the minimum value.
  *  The Minimizer type must provide the following members: \code
     typedef Scalar //the numeric type of what is being minimized--not necessarily the Scalar type
  of the BVH (if it has one) Scalar minimumOnVolume(const BVH::Volume &volume) Scalar
  minimumOnObject(const BVH::Object &object) \endcode
  */
template <typename BVH, typename Minimizer>
typename Minimizer::Scalar BVMinimize(const BVH& tree, Minimizer& minimizer)
{
    return internal::minimize_helper(
        tree,
        minimizer,
        tree.getRootIndex(),
        (std::numeric_limits<typename Minimizer::Scalar>::max)());
}

/**  Given two BVH's, runs the query on their cartesian product encapsulated by \a minimizer.
  *  \returns the minimum value.
  *  The Minimizer type must provide the following members: \code
     typedef Scalar //the numeric type of what is being minimized--not necessarily the Scalar type
  of the BVH (if it has one) Scalar minimumOnVolumeVolume(const BVH1::Volume &v1, const BVH2::Volume
  &v2) Scalar minimumOnVolumeObject(const BVH1::Volume &v1, const BVH2::Object &o2) Scalar
  minimumOnObjectVolume(const BVH1::Object &o1, const BVH2::Volume &v2) Scalar
  minimumOnObjectObject(const BVH1::Object &o1, const BVH2::Object &o2) \endcode
  */
template <typename BVH1, typename BVH2, typename Minimizer>
typename Minimizer::Scalar BVMinimize(const BVH1& tree1, const BVH2& tree2, Minimizer& minimizer)
{
    typedef typename Minimizer::Scalar Scalar;
    typedef typename BVH1::Index Index1;
    typedef typename BVH2::Index Index2;
    typedef internal::minimizer_helper1<
        typename BVH1::Volume,
        typename BVH1::Object,
        typename BVH2::Object,
        Minimizer>
        Helper1;
    typedef internal::minimizer_helper2<
        typename BVH2::Volume,
        typename BVH2::Object,
        typename BVH1::Object,
        Minimizer>
        Helper2;
    typedef std::pair<Scalar, std::pair<Index1, Index2>> QueueElement; // first element is priority
    typedef typename BVH1::VolumeIterator VolIter1;
    typedef typename BVH1::ObjectIterator ObjIter1;
    typedef typename BVH2::VolumeIterator VolIter2;
    typedef typename BVH2::ObjectIterator ObjIter2;

    VolIter1 vBegin1 = VolIter1(), vEnd1 = VolIter1();
    ObjIter1 oBegin1 = ObjIter1(), oEnd1 = ObjIter1();
    VolIter2 vBegin2 = VolIter2(), vEnd2 = VolIter2(), vCur2 = VolIter2();
    ObjIter2 oBegin2 = ObjIter2(), oEnd2 = ObjIter2(), oCur2 = ObjIter2();
    std::priority_queue<QueueElement, std::vector<QueueElement>, std::greater<QueueElement>>
        todo; // smallest is at the top

    Scalar minimum = (std::numeric_limits<Scalar>::max)();
    todo.push(std::make_pair(Scalar(), std::make_pair(tree1.getRootIndex(), tree2.getRootIndex())));

    while (!todo.empty())
    {
        tree1.getChildren(todo.top().second.first, vBegin1, vEnd1, oBegin1, oEnd1);
        tree2.getChildren(todo.top().second.second, vBegin2, vEnd2, oBegin2, oEnd2);
        todo.pop();

        for (; oBegin1 != oEnd1; ++oBegin1)
        { // go through child objects of first tree
            for (oCur2 = oBegin2; oCur2 != oEnd2; ++oCur2)
            { // go through child objects of second tree
                minimum = (std::min)(minimum, minimizer.minimumOnObjectObject(*oBegin1, *oCur2));
            }

            for (vCur2 = vBegin2; vCur2 != vEnd2; ++vCur2)
            { // go through child volumes of second tree
                Helper2 helper(*oBegin1, minimizer);
                minimum =
                    (std::min)(minimum, internal::minimize_helper(tree2, helper, *vCur2, minimum));
            }
        }

        for (; vBegin1 != vEnd1; ++vBegin1)
        { // go through child volumes of first tree
            const typename BVH1::Volume& vol1 = tree1.getVolume(*vBegin1);

            for (oCur2 = oBegin2; oCur2 != oEnd2; ++oCur2)
            { // go through child objects of second tree
                Helper1 helper(*oCur2, minimizer);
                minimum = (std::min)(
                    minimum,
                    internal::minimize_helper(tree1, helper, *vBegin1, minimum));
            }

            for (vCur2 = vBegin2; vCur2 != vEnd2; ++vCur2)
            { // go through child volumes of second tree
                Scalar val = minimizer.minimumOnVolumeVolume(vol1, tree2.getVolume(*vCur2));
                if (val < minimum)
                    todo.push(std::make_pair(val, std::make_pair(*vBegin1, *vCur2)));
            }
        }
    }
    return minimum;
}

} // end namespace Eigen

namespace Eigen {

namespace internal {

// internal pair class for the BVH--used instead of std::pair because of alignment
template <typename Scalar, int Dim>
struct vector_int_pair
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW_IF_VECTORIZABLE_FIXED_SIZE(Scalar, Dim)
    typedef Matrix<Scalar, Dim, 1> VectorType;

    vector_int_pair(const VectorType& v, int i) : first(v), second(i) {}

    VectorType first;
    int second;
};

// these templates help the tree initializer get the bounding boxes either from a provided
// iterator range or using bounding_box in a unified way
template <typename ObjectList, typename VolumeList, typename BoxIter>
struct get_boxes_helper
{
    void
    operator()(const ObjectList& objects, BoxIter boxBegin, BoxIter boxEnd, VolumeList& outBoxes)
    {
        outBoxes.insert(outBoxes.end(), boxBegin, boxEnd);
        eigen_assert(outBoxes.size() == objects.size());
        EIGEN_ONLY_USED_FOR_DEBUG(objects);
    }
};

template <typename ObjectList, typename VolumeList>
struct get_boxes_helper<ObjectList, VolumeList, int>
{
    void operator()(const ObjectList& objects, int, int, VolumeList& outBoxes)
    {
        outBoxes.reserve(objects.size());
        for (int i = 0; i < (int)objects.size(); ++i)
            outBoxes.push_back(bounding_box(objects[i]));
    }
};

} // end namespace internal

/** \class KdBVH
 *  \brief A simple bounding volume hierarchy based on AlignedBox
 *
 *  \param _Scalar The underlying scalar type of the bounding boxes
 *  \param _Dim The dimension of the space in which the hierarchy lives
 *  \param _Object The object type that lives in the hierarchy.  It must have value semantics.
 * Either bounding_box(_Object) must be defined and return an AlignedBox<_Scalar, _Dim> or bounding
 * boxes must be provided to the tree initializer.
 *
 *  This class provides a simple (as opposed to optimized) implementation of a bounding volume
 * hierarchy analogous to a Kd-tree. Given a sequence of objects, it computes their bounding boxes,
 * constructs a Kd-tree of their centers and builds a BVH with the structure of that Kd-tree.  When
 * the elements of the tree are too expensive to be copied around, it is useful for _Object to be a
 * pointer.
 */
template <typename _Scalar, int _Dim, typename _Object>
class KdBVH
{
  public:
    enum { Dim = _Dim };
    typedef _Object Object;
    typedef std::vector<Object, aligned_allocator<Object>> ObjectList;
    typedef _Scalar Scalar;
    typedef AlignedBox<Scalar, Dim> Volume;
    typedef std::vector<Volume, aligned_allocator<Volume>> VolumeList;
    typedef int Index;
    typedef const int* VolumeIterator; // the iterators are just pointers into the tree's vectors
    typedef const Object* ObjectIterator;

    KdBVH() {}

    /** Given an iterator range over \a Object references, constructs the BVH.  Requires that
     * bounding_box(Object) return a Volume. */
    template <typename Iter>
    KdBVH(Iter begin, Iter end)
    {
        init(begin, end, 0, 0);
    } // int is recognized by init as not being an iterator type

    /** Given an iterator range over \a Object references and an iterator range over their bounding
     * boxes, constructs the BVH */
    template <typename OIter, typename BIter>
    KdBVH(OIter begin, OIter end, BIter boxBegin, BIter boxEnd)
    {
        init(begin, end, boxBegin, boxEnd);
    }

    /** Given an iterator range over \a Object references, constructs the BVH, overwriting whatever
     * is in there currently. Requires that bounding_box(Object) return a Volume. */
    template <typename Iter>
    void init(Iter begin, Iter end)
    {
        init(begin, end, 0, 0);
    }

    /** Given an iterator range over \a Object references and an iterator range over their bounding
     * boxes, constructs the BVH, overwriting whatever is in there currently. */
    template <typename OIter, typename BIter>
    void init(OIter begin, OIter end, BIter boxBegin, BIter boxEnd)
    {
        objects.clear();
        boxes.clear();
        children.clear();

        objects.insert(objects.end(), begin, end);
        int n = static_cast<int>(objects.size());

        if (n < 2)
            return; // if we have at most one object, we don't need any internal nodes

        VolumeList objBoxes;
        VIPairList objCenters;

        // compute the bounding boxes depending on BIter type
        internal::get_boxes_helper<ObjectList, VolumeList, BIter>()(
            objects,
            boxBegin,
            boxEnd,
            objBoxes);

        objCenters.reserve(n);
        boxes.reserve(n - 1);
        children.reserve(2 * n - 2);

        for (int i = 0; i < n; ++i)
            objCenters.push_back(VIPair(objBoxes[i].center(), i));

        build(objCenters, 0, n, objBoxes, 0); // the recursive part of the algorithm

        ObjectList tmp(n);
        tmp.swap(objects);
        for (int i = 0; i < n; ++i)
            objects[i] = tmp[objCenters[i].second];
    }

    /** \returns the index of the root of the hierarchy */
    inline Index getRootIndex() const { return (int)boxes.size() - 1; }

    /** Given an \a index of a node, on exit, \a outVBegin and \a outVEnd range over the indices of
     * the volume children of the node and \a outOBegin and \a outOEnd range over the object
     * children of the node */
    EIGEN_STRONG_INLINE void getChildren(
        Index index,
        VolumeIterator& outVBegin,
        VolumeIterator& outVEnd,
        ObjectIterator& outOBegin,
        ObjectIterator& outOEnd) const
    { // inlining this function should open lots of optimization opportunities to the compiler
        if (index < 0)
        {
            outVBegin = outVEnd;
            if (!objects.empty())
                outOBegin = &(objects[0]);
            outOEnd =
                outOBegin +
                objects.size(); // output all objects--necessary when the tree has only one object
            return;
        }

        int numBoxes = static_cast<int>(boxes.size());

        int idx = index * 2;
        if (children[idx + 1] < numBoxes)
        { // second index is always bigger
            outVBegin = &(children[idx]);
            outVEnd   = outVBegin + 2;
            outOBegin = outOEnd;
        }
        else if (children[idx] >= numBoxes)
        { // if both children are objects
            outVBegin = outVEnd;
            outOBegin = &(objects[children[idx] - numBoxes]);
            outOEnd   = outOBegin + 2;
        }
        else
        { // if the first child is a volume and the second is an object
            outVBegin = &(children[idx]);
            outVEnd   = outVBegin + 1;
            outOBegin = &(objects[children[idx + 1] - numBoxes]);
            outOEnd   = outOBegin + 1;
        }
    }

    /** \returns the bounding box of the node at \a index */
    inline const Volume& getVolume(Index index) const { return boxes[index]; }

    inline Volume& getVolume(Index index) { return boxes[index]; }

    inline std::vector<int> const& getChildren() const { return children; }

    inline VolumeList const& getBoxes() const { return boxes; }
    inline ObjectList const& getObjects() const { return objects; }

  private:
    typedef internal::vector_int_pair<Scalar, Dim> VIPair;
    typedef std::vector<VIPair, aligned_allocator<VIPair>> VIPairList;
    typedef Matrix<Scalar, Dim, 1> VectorType;
    struct VectorComparator // compares vectors, or more specifically, VIPairs along a particular
                            // dimension
    {
        VectorComparator(int inDim) : dim(inDim) {}
        inline bool operator()(const VIPair& v1, const VIPair& v2) const
        {
            return v1.first[dim] < v2.first[dim];
        }
        int dim;
    };

    // Build the part of the tree between objects[from] and objects[to] (not including objects[to]).
    // This routine partitions the objCenters in [from, to) along the dimension dim, recursively
    // constructs the two halves, and adds their parent node.  TODO: a cache-friendlier layout
    void build(VIPairList& objCenters, int from, int to, const VolumeList& objBoxes, int dim)
    {
        eigen_assert(to - from > 1);
        if (to - from == 2)
        {
            boxes.push_back(
                objBoxes[objCenters[from].second].merged(objBoxes[objCenters[from + 1].second]));
            children.push_back(
                from + (int)objects.size() - 1); // there are objects.size() - 1 tree nodes
            children.push_back(from + (int)objects.size());
        }
        else if (to - from == 3)
        {
            int mid = from + 2;
            std::nth_element(
                objCenters.begin() + from,
                objCenters.begin() + mid,
                objCenters.begin() + to,
                VectorComparator(dim)); // partition
            build(objCenters, from, mid, objBoxes, (dim + 1) % Dim);
            int idx1 = (int)boxes.size() - 1;
            boxes.push_back(boxes[idx1].merged(objBoxes[objCenters[mid].second]));
            children.push_back(idx1);
            children.push_back(mid + (int)objects.size() - 1);
        }
        else
        {
            int mid = from + (to - from) / 2;
            nth_element(
                objCenters.begin() + from,
                objCenters.begin() + mid,
                objCenters.begin() + to,
                VectorComparator(dim)); // partition
            build(objCenters, from, mid, objBoxes, (dim + 1) % Dim);
            int idx1 = (int)boxes.size() - 1;
            build(objCenters, mid, to, objBoxes, (dim + 1) % Dim);
            int idx2 = (int)boxes.size() - 1;
            boxes.push_back(boxes[idx1].merged(boxes[idx2]));
            children.push_back(idx1);
            children.push_back(idx2);
        }
    }

    std::vector<int> children; // children of x are children[2x] and children[2x+1], indices bigger
                               // than boxes.size() index into objects.
    VolumeList boxes;
    ObjectList objects;
};

} // end namespace Eigen

#endif // SBS_EIGEN_BVH_MODULE_H

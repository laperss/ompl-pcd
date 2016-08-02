#ifndef __PCD_CELL_H_INCLUDED__
#define __PCD_CELL_H_INCLUDED__ 
//#include "PCD.h"
#include <deque>
#include <vector>
#include <iostream>
#include "ompl/control/SpaceInformation.h"
#include "ompl/geometric/planners/PlannerIncludes.h"

using namespace std;

namespace ob = ompl::base;
#define PCD_DBL_EPS 1.0e-12
#define PCD_MAX_NUM_SAMPLES_PER_CELL 400

class PCD;
class Cell;
typedef std::vector<Cell*> PCD_Graph;

struct PathSegment {
PathSegment(Cell& c): cell(&c) {}
Cell* cell;
const ob::State* p;
const ob::State* q;
};

typedef deque<PathSegment> CellPath;

class Cell{
public:
    typedef boost::uint64_t BoundaryFlagT;
    enum    CellType {POSS_FREE,  POSS_FULL,  MIXED};
    typedef unsigned int Cell_ID;
    enum    RedundancyCheck { DO_CHECK, NO_CHECK };
    enum {REG_UNSPEC = 1, REG_START = 2, REG_GOAL = 3,BRIDGE_FACTOR = 6};

    // data
    vector<const ob::State*> samples_in_cell_;
    vector<double>           lower_;
    vector<double>           upper_;
    vector<double>           range_;
    vector<double>           centroid_;
    vector<Cell*>            neighbors_;
    Cell*                    parent_;
    Cell*                    lo_child_;
    Cell*                    up_child_;

    CellType                 type_;
    unsigned int             region_type_;

    // used for A* search
    unsigned int             visited_;
    unsigned int             closed_;
    double                   init_dist_;
    double                   goal_dist_;
    Cell*                    previous_node_;  
    static unsigned int      num_samples;
    
    
    static ob::SpaceInformationPtr si_;
    static void setSI(const ob::SpaceInformationPtr&);

    // Constructor 
    Cell(const ompl::base::SpaceInformationPtr& si);
    Cell(const vector<double>& min_vals, 
	 const vector<double>& max_vals,
	 Cell&                 parent_cell);

    bool strictlyContains(const ob::State* sample) const;
    bool contains(const ob::State* sample) const;
    bool addSample(const ob::State* sample, RedundancyCheck check = DO_CHECK);
    ob::State*  getSample(ompl::RNG &rng_); 
    bool isOnBoundary(const ob::State* state, unsigned int& coord_indx) const;
    bool isEmpty() const;
    bool SampleExists(const ob::State& sample) const;
    bool isBoundaryCell() const;
    bool isLeftBoundary(unsigned int i) const;
    bool isRightBoundary(unsigned int i) const;
    bool isLeftAndRightBoundary(unsigned int i) const;
    CellType getType() const;
    Cell_ID getID() const;  
    void strip();
    void removeLeftBoundary(unsigned int i);
    void removeRightBoundary(unsigned int i);
    void print(std::ostream& os) const;

    // static functions
    static PCD_Graph createGraph(const ob::State*   start,
				 const ob::State*   goal,
				 vector<double>&    min_vals,
				 vector<double>&    max_vals);
    static double cellDistance(Cell& cell1, Cell& cell2);
    static void destroyGraph(PCD_Graph& graph);
    static void makeNeighbors(Cell& a, Cell& b);
    static bool isFreeBridgeCell(const Cell& cell);
    static void getCommonCenter(PathSegment& a, PathSegment& b);
    static void getPath(CellPath& cell_path);
    static bool cellsSharePlane(const Cell& a, const Cell& b, unsigned int i);
    static unsigned int getNumCells();
    static unsigned int getNumSamples();
    void splitCell(const ob::State*   state,
		   vector<bool>       valid_directions,
		   PCD_Graph&         cell_graph);
private:
    BoundaryFlagT         left_boundary_;
    BoundaryFlagT         right_boundary_;

    static unsigned int   num_cells;
    static unsigned int   next_cell_id;
    unsigned int          id_;

    Cell(const vector<double>& min_vals,
	 const vector<double>& max_vals);
};


struct cellGreater {bool operator()(const Cell* p, const Cell* q) const;};
bool areEqual(double a, double b);
bool overlap(double lo1, double hi1, double lo2, double hi2);
bool disjoint(double lo1, double hi1, double lo2, double hi2);
bool areAdjacent(const Cell& a, const Cell& b);

// INLINE FUNCTIONS
inline bool Cell::isEmpty() const
{
    return samples_in_cell_.empty();
}

inline Cell::CellType Cell::getType() const
{
    return type_;
}

inline Cell::Cell_ID Cell::getID() const
{
    return id_;
}

inline bool Cell::isBoundaryCell() const
{
    return left_boundary_ != 0 || right_boundary_ != 0;
}

inline bool Cell::isLeftBoundary(unsigned int i) const
{
    assert(i < 64);
    return (left_boundary_ & (UINT64_C(1) << i)) != 0;
}

inline bool Cell::isRightBoundary(unsigned int i) const
{
    assert(i < 64);
    return (right_boundary_ & (UINT64_C(1) << i)) != 0;
}

inline bool Cell::isLeftAndRightBoundary(unsigned int i) const
{
    return isLeftBoundary(i) && isRightBoundary(i);
}

inline void Cell::removeLeftBoundary(unsigned int i)
{
  assert(i < 64);
  left_boundary_ &= ~(UINT64_C(1) << i);
  assert(!isLeftBoundary(i));
  return;
}


inline void Cell::removeRightBoundary(unsigned int i)
{
  assert(i < 64);
  right_boundary_ &= ~(UINT64_C(1) << i);
  assert(!isRightBoundary(i));
  return;
}

inline bool Cell::cellsSharePlane(const Cell& a, 
				  const Cell& b, unsigned int i)
{
    return  areEqual(a.lower_[i], b.upper_[i])         ||
	areEqual(a.upper_[i], b.lower_[i])             ||
	(a.isLeftBoundary(i)  && b.isRightBoundary(i)) ||
	(a.isRightBoundary(i) && b.isLeftBoundary(i));
}

inline bool cellGreater::operator()(const Cell* p, const Cell* q) const
{ 
    return (p->init_dist_ + p->goal_dist_) > (q->init_dist_ + q->goal_dist_); 
} 

inline bool areEqual(double a, double b) 
{
    return fabs(a - b) <= PCD_DBL_EPS;
}

inline bool overlap(double lo1, double hi1, double lo2, double hi2) 
{
    return !disjoint(lo1, hi1, lo2, hi2);
}

inline bool disjoint(double lo1, double hi1, double lo2, double hi2)
{
    assert((lo1 <= hi1) && (lo2 <= hi2));
    return (hi1 <= lo2) || (lo1 >= hi2);
}

template<typename T>
inline void clear(std::vector<T>& v)
{
    std::vector<T>().swap(v);
    return;
    }
#endif 

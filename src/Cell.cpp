#include "../Cell.h"
#include "../PCD.h"
#include <iostream>
#include "ompl/base/SpaceInformation.h"

using namespace std;
namespace ob = ompl::base;


#define MIN(a,b) (((a) > (b))? (b) : (a))
#define MAX(a,b) (((a) > (b))? (a) : (b))


ob::SpaceInformationPtr Cell::si_;
unsigned int Cell::num_cells    = 0;
unsigned int Cell::next_cell_id = 0;
unsigned int Cell::num_samples  = 0;


void Cell::setSI(const ob::SpaceInformationPtr& si){
	si_ = si;
}


// --------------------------- CONSTRUCTORS ----------------------------------

Cell::Cell(const vector<double>& min_vals,
           const vector<double>& max_vals,
           Cell&                 parent_cell):
    lower_(min_vals),
    upper_(max_vals),
    centroid_(),
    type_(parent_cell.type_),
    region_type_(REG_UNSPEC),
    neighbors_(),
    parent_(&parent_cell),
    lo_child_(0),
    up_child_(0),
    visited_(0),
    closed_(0),
    init_dist_(INFINITY),
    goal_dist_(-1.0),
    previous_node_(0),
    id_(next_cell_id++),
    left_boundary_(parent_cell.left_boundary_),
    right_boundary_(parent_cell.right_boundary_)
{
    assert(!min_vals.empty());
    assert(min_vals.size() == max_vals.size());
    unsigned int i = min_vals.size();
    centroid_.resize(i);
    range_.resize(i);
    for ( ; (i)-- != 0; ) {
	range_[i] = max_vals[i]-min_vals[i];
	assert(min_vals[i] < max_vals[i]);
	centroid_[i] = 0.5 * (min_vals[i] + max_vals[i]);
    }

     ++num_cells;
}

Cell::Cell(const vector<double>& min_vals, 
	   const vector<double>& max_vals):
//define boundaries
    lower_(min_vals),
    upper_(max_vals),
    region_type_(REG_UNSPEC),
    parent_(0),
    lo_child_(0),
    up_child_(0),
    visited_(0),
    closed_(0),
    type_(POSS_FREE),
    id_(next_cell_id++),
    goal_dist_(-1.0),
    left_boundary_(UINT64_C(0)),
    right_boundary_(UINT64_C(0))
{
    assert(!min_vals.empty());
    assert(min_vals.size() == max_vals.size());
    ++num_cells;

    unsigned int i = lower_.size();
    centroid_.resize(i);
    range_.resize(i);
    for ( ; (i)-- != 0; ) {
	range_[i] = max_vals[i]-min_vals[i];
	centroid_[i] = 0.5 * (lower_[i] + upper_[i]);
	assert(lower_[i] <= upper_[i]);
    }
}

PCD_Graph Cell::createGraph(const ob::State* start,
			    const ob::State* goal,
			    vector<double>&  min_vals,
			    vector<double>&  max_vals)
{
    Cell* const root_cell = new Cell(min_vals, max_vals);
    root_cell->addSample(start, NO_CHECK);
    root_cell->addSample(goal, NO_CHECK);
    return PCD_Graph(1, root_cell);
}

void Cell::destroyGraph(PCD_Graph& graph)
{
    unsigned int i = graph.size();
    for ( ; (i)-- != 0; ) 
    {
	delete graph[i];
    }
    clear(graph);
    assert(graph.empty() && graph.capacity() == 0);
}

double Cell::cellDistance(Cell& cell1, Cell& cell2)
{
    vector<double> v1 = cell1.centroid_;
    vector<double> v2 = cell1.centroid_;
    double dist          = 0.0;
    for (unsigned int i = 0; i < v1.size(); ++i) 
    {
	dist += (v1[i]-v2[i])<0? -(v1[i]-v2[i]) : (v1[i]-v2[i]);
    }
    return dist;
}

void Cell::strip()
{
    clear(samples_in_cell_);
    clear(neighbors_);
    clear(centroid_); 
    return;
}


// ========================= SAMPLE FUNCTIONS ================================

bool Cell::addSample(const ob::State* sample, RedundancyCheck check)
{
  ++num_samples;
  // copy and add state
  ob::State* p_sample =   si_->getStateSpace()->cloneState(sample);
  samples_in_cell_.push_back(p_sample);
  return p_sample != 0;
}

ob::State*  Cell::getSample(ompl::RNG &rng_) 	
{
    ob::ScopedState<> state(si_->getStateSpace());
    for (unsigned int i = 0; i < lower_.size(); ++i ) {
        state[i] = lower_[i] + rng_.uniform01()*(upper_[i] - lower_[i]);
    }
    ob::State* state2 =  si_->cloneState(state.get());
    return state2;
}

bool Cell::strictlyContains(const ob::State* config) const
{
    ob::ScopedState<> scoped_state(si_->getStateSpace());
    scoped_state = config;
    
    int dim = scoped_state.getSpace()->getDimension();
    for (int i=0; i<dim; ++i ) 
    {
	if (scoped_state[i] <= lower_[i] || scoped_state[i] >= upper_[i]) 
	    return false;
    }
    return true;
}

bool Cell::contains(const ob::State* config) const
{
    ob::ScopedState<> scoped_state(si_->getStateSpace());
    scoped_state = config;
    int dim = scoped_state.getSpace()->getDimension();
    for (int i=0; i<dim; ++i ) 
    {
	if ((scoped_state[i] > upper_[i])) 
	    return false;
	if ((scoped_state[i] < lower_[i])) 
	    return false;
    }
    return true;
}

bool areAdjacent(const Cell& a, const Cell& b)
{
    assert(&a != &b);
  
    const unsigned int dim = a.lower_.size();
    unsigned int i         = dim;
  
    if (a.isBoundaryCell() || b.isBoundaryCell()) 
    {
	for ( ; (i)-- != 0; ) 
	{
	    if (Cell::cellsSharePlane(a, b, i)) 
	    {
		bool adjacent  = true;
		unsigned int j = dim;
		for ( ; (j)-- != 0; ) 
		{
		    if ((i != j) && disjoint(a.lower_[j], a.upper_[j], b.lower_[j], b.upper_[j])) 
		    {
			adjacent = false;
			break;
		    }
		}
		if (adjacent) 
		    return true;
	    }
	}
    } else 
    {
	for ( ; (i)-- != 0; ) 
	{
	    if (areEqual(a.lower_[i], b.upper_[i]) || areEqual(a.upper_[i], b.lower_[i])) 
	    {
		bool adjacent  = true;
		unsigned int j = dim;
		for ( ; (j)-- != 0; ) {
		    if ((i != j) && disjoint(a.lower_[j], a.upper_[j], b.lower_[j], b.upper_[j])) 
		    {
			adjacent = false;
			break;
		    }
		}
		if (adjacent) 
		    return true;
	    }
	}
    }
    return false;
}

void Cell::makeNeighbors(Cell& a, Cell& b)
{
    assert(&a != &b);
    assert(areAdjacent(a, b));
    assert(areAdjacent(b, a));
    assert(a.type_ != Cell::MIXED);
    assert(b.type_ != Cell::MIXED);
    assert(find(a.neighbors_.begin(), a.neighbors_.end(), &b) == a.neighbors_.end());
    assert(find(b.neighbors_.begin(), b.neighbors_.end(), &a) == b.neighbors_.end());
  
    a.neighbors_.push_back(&b);
    b.neighbors_.push_back(&a);
}

void Cell::getCommonCenter(PathSegment& a, PathSegment& b)
{
    assert(a.cell != b.cell);
    assert(areAdjacent(*a.cell, *b.cell));
    assert(areAdjacent(*b.cell, *a.cell));
  
    unsigned int dim = a.cell->lower_.size();
    vector<double> aq;
    vector<double> bp;
    for (int i=0; i<dim; ++i) 
    {
	const double a_low = a.cell->lower_[i];
	const double a_upp = a.cell->upper_[i];
	const double b_low = b.cell->lower_[i];
	const double b_upp = b.cell->upper_[i];
	if (overlap(a_low, a_upp, b_low, b_upp)) 
	{
	    aq.push_back(0.5 * (MIN(a_upp, b_upp) 
				+ MAX(a_low, b_low)));
	    bp.push_back(0.5 * (MIN(a_upp, b_upp) 
				+ MAX(a_low, b_low)));
	} 
	else if (areEqual(a_upp, b_low)) 
	{
	    aq.push_back(b_low);
	    bp.push_back(b_low);
	} else if (areEqual(a_low, b_upp)) 
	{
	    aq.push_back(a_low);
	    bp.push_back(a_low);
	} else if (a.cell->isLeftBoundary(i) && b.cell->isRightBoundary(i)) 
	{
	    aq.push_back(a_low);
	    bp.push_back(b_upp);
	} else if (b.cell->isLeftBoundary(i) && a.cell->isRightBoundary(i)) 
	{
	    aq.push_back( a_upp);
	    bp.push_back(b_low);
	} 
    }
    ob::ScopedState<> aq_ss(si_->getStateSpace());
    ob::ScopedState<> bp_ss(si_->getStateSpace());
    for (int i = 0; i<dim; ++i) 
    {
	aq_ss[i] = aq[i];
	bp_ss[i] = bp[i];
    }
    ob::State* aq_state = si_->cloneState(aq_ss.get());
    ob::State* bp_state = si_->cloneState(bp_ss.get());
 
    a.q =  aq_state;
    b.p =  bp_state;
    return;
}

void Cell::getPath(CellPath& cell_path)
{
    const unsigned int n = cell_path.size();
    for (unsigned int i = 0, j = 1; j < n; i = j++) 
	getCommonCenter(cell_path[i], cell_path[j]);
}

unsigned int Cell::getNumCells()
{
    return num_cells;
}

unsigned int Cell::getNumSamples()
{
  return num_samples;
}

void Cell::print(ostream& os) const
{
  const unsigned dim = lower_.size();
  os << "------------------\n";
  os << "cell id: " << getID() << "\n"
    << "("<<upper_[0]<<", "<<upper_[1]<<"), ("<<lower_[0]<<", "<<lower_[1]<<")\n"
      //<< "dim: "     << dim     << "\n"
     << "type: ";
  if (type_ == POSS_FREE) {
    os << "POSS_FREE\n";
  } else if (type_ == POSS_FULL) {
    os << "POSS_FULL\n";
  } else {
    os << "MIXED\n";
  }
  os << "region: ";
  if (region_type_ == REG_UNSPEC) {
    os << " unspecified\n";
  } else if (region_type_ == REG_START) {
    os << " start region\n";
  } else if (region_type_ == REG_GOAL) {
    os << " goal region\n";
  } else {
    os << region_type_ << ", BUG! (invalid region type)\n";
  }
  
  os << "num neighbors: "    << neighbors_.size()       << "\n"
     << "num samples:   "    << samples_in_cell_.size() << "\n";
  os << "neigbors: (";
  for(int i=0; i<neighbors_.size();++i)
      os << neighbors_[i]->getID() << " ";

  os << ")\n------------------\n";
}

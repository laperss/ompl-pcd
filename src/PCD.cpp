#include "../PCD.h"
#include "../Cell.h"
#include "../SplitAdvisor.h"

#include <iostream>
#include <ompl/base/spaces/RealVectorStateSpace.h>

#define DEF_MAX_NEW_FREE_CELLS 10
#define DEF_MAX_NUM_ITER 100000

// ======================== CONSTRUCTORS ============================
og::PCD::PCD(const ob::SpaceInformationPtr &si, vector<Cell*>& cell_division) : 
    ob::Planner(si, "PCD"),
    astar_timer_(0),
    pcd_graph_(0),
    max_num_iter_(DEF_MAX_NUM_ITER),

    max_new_free_cells_(DEF_MAX_NEW_FREE_CELLS),
    split_advisor_(0),
    max_step_size_(1.0 *3.14/180),
    cell_division_(cell_division)
{
    Cell::setSI(si);
	
    // Planner specs
    Planner::specs_.multithreaded = false;
    Planner::specs_.approximateSolutions = false; //For now!
    Planner::specs_.optimizingPaths = false;
    Planner::specs_.directed = false;
    Planner::specs_.provingSolutionNonExistence = false;
}

og::PCD::PCD(const ob::SpaceInformationPtr &si) : 
    ob::Planner(si, "PCD"),
    astar_timer_(0),
    pcd_graph_(0),
    max_num_iter_(DEF_MAX_NUM_ITER),

    max_new_free_cells_(DEF_MAX_NEW_FREE_CELLS),
    split_advisor_(0),
    max_step_size_(1.0 *3.14/180),
    cell_division_(0)
{
    Cell::setSI(si);  // makes it possible to access si_ in cell
    Planner::specs_.multithreaded = false;
    Planner::specs_.approximateSolutions = false; 
    Planner::specs_.optimizingPaths = false;
    Planner::specs_.directed = false;
    Planner::specs_.provingSolutionNonExistence = false;
}

// free memory in destructor
og::PCD::~PCD()
{
    freeMemory();
}

void ompl::geometric::PCD::freeMemory()
{
}

// =============================== PROBLEM SETUP ========================================
void og::PCD::setProblemDefinition(const ob::ProblemDefinitionPtr &pdef)
{
    Planner::setProblemDefinition(pdef);
}

void og::PCD::setup()
{
    Planner::setup();
    // THIS IS NOW DEPENDENT ON THE STATE SPACE AND MUST BE CHANGED!
    const ob::RealVectorBounds & bounds = si_->getStateSpace()->as<ob::RealVectorStateSpace>()->getBounds () ;
    
    min_bound_ = bounds.low ;
    max_bound_ = bounds.high ;
    split_directions_.assign(2, true);
    Cell::destroyGraph(pcd_graph_);
}

bool og::PCD::setSplitAdvisor(const og::SplitAdvisor& advisor)
{
    og::SplitAdvisor* const tmp = new og::SplitAdvisor(advisor);
    if (split_advisor_ != 0) 
	delete split_advisor_;
    split_advisor_ = tmp;
    return true;
}

void og::PCD::clear()
{
    Planner::clear();
    sampler_.reset();
}

void og::PCD::getPlannerData(ob::PlannerData &data) const
{
    Planner::getPlannerData(data);
}


// ============================= SOLUTION METHODS =================================

ob::PlannerStatus og::PCD::solve(const ob::PlannerTerminationCondition &ptc)
{
    const unsigned int free_cells_sample_period = 5000;
    unsigned int next_free_cells_sampling       = free_cells_sample_period;
    CellPath cell_path;
    bool success            = false;
    bool approximate        = false;
    unsigned int num_iter   = 0;

    start_time = ompl::time::now();
    start_ptr =   pis_.nextStart();
    goal_ptr  =   pis_.nextGoal();

    assert(isSatisfied(start_ptr));
    assert(isSatisfied(goal_ptr));
    
    pcd_graph_ = Cell::createGraph(start_ptr,goal_ptr,min_bound_,max_bound_);
 
    checkValidity();
    while ((num_iter < max_num_iter_) && !success) 
    {
	++num_iter;
	cout << "\n* find cell path\n";
	if (findCellPath(cell_path)) 
	{
	    cout << "* get cell path\n";
	    Cell::getPath(cell_path);
	    cout << "(";
	    for(int i=0; i<cell_path.size();++i)
		cout << cell_path[i].cell->getID()<< ", ";
	    cout << ")\n";
	    if (checkPath(cell_path)) 
	    {
		success = true;
	    }

	} else {
	    if (pcd_graph_.size() >= next_free_cells_sampling) 
	    {
		next_free_cells_sampling += free_cells_sample_period;
		cout << "* sample free cell\n";
		sampleFreeCells();
	    } else		
	    {
		cout << "* sample occupied cell\n";
		sampleOccCells();		
	    }	
	}
	cout << "\n";
    }
    if (success == true)
    {
	og::PathGeometric *path = new og::PathGeometric(si_);
	path->append(cell_path[0].p);
	for (int i= 0 ; i <  cell_path.size(); ++i)
	{
	    path->append(cell_path[i].q);
	}
	pdef_->addSolutionPath(ob::PathPtr(path));
    }
    getAllCells(cell_division_, pcd_graph_.front());
    return ob::PlannerStatus(success , approximate);
}

bool og::PCD::findCellPath(CellPath& cell_path)
{
    ++astar_timer_;
    cell_path.clear();
    assert(!pcd_graph_.empty());

    Cell& start_cell = getNonMixedCell(*pcd_graph_.front(), start_ptr);
    Cell& goal_cell  = getNonMixedCell(*pcd_graph_.front(), goal_ptr);

    if (&start_cell == &goal_cell) 
    {
	assert(start_cell.type_ == Cell::POSS_FREE);    
	PathSegment segment(start_cell);
	segment.p = start_ptr;
	segment.q = goal_ptr;
	cell_path.push_front(segment);
	return true;
    }

    priority_queue<Cell*, vector<Cell*>, cellGreater> open_cells;
    deque<Cell*> new_cells;

    assert(start_cell.contains(start_ptr));
    assert(goal_cell.contains(goal_ptr));
    assert(start_cell.type_ == Cell::POSS_FREE);
    assert(goal_cell.type_  == Cell::POSS_FREE);
    assert(&start_cell != &goal_cell);
  
    // initiate A* search:
    Cell* current_node     = &start_cell;
    current_node->init_dist_   = 0.0;
    current_node->goal_dist_ = cellDistance(start_cell, goal_cell);
    current_node->previous_node_ = 0;	
    current_node->visited_ = astar_timer_;
    open_cells.push(current_node);
    // start A* search
    for (;;)
    {
	while (!open_cells.empty() && open_cells.top()->closed_ == astar_timer_)
	{
	    open_cells.pop();
	}
	if (open_cells.empty()) 
	    break; 
	current_node = open_cells.top();
	open_cells.pop();
	new_cells.clear();
	current_node->closed_ = astar_timer_; 	
	unsigned int i = current_node->neighbors_.size();

	for ( ; (i)-- != 0; )
	{
	    Cell& neighbor = *current_node->neighbors_[i];
	    if ((neighbor.type_ == Cell::POSS_FREE) && (neighbor.closed_ != astar_timer_)) 
	    {
		const double path_dist = current_node->init_dist_ + 
		    cellDistance(*current_node, neighbor); 
		if (neighbor.visited_ != astar_timer_ || path_dist < neighbor.init_dist_) 
		{
		    new_cells.push_front(&neighbor);
		    neighbor.init_dist_     = path_dist;
		    neighbor.previous_node_ = current_node;
		    if (neighbor.visited_  != astar_timer_) 
		    {
			neighbor.goal_dist_ = cellDistance(neighbor, goal_cell);;
		    }
		}
	    }
	}

	const unsigned int n = new_cells.size();
	for (i = 0; i < n; ++i) 
	{
	    current_node = new_cells[i];
	    if (current_node == &goal_cell)  
	    {
		cell_path.push_back(PathSegment(*current_node));
		while (current_node->previous_node_ != 0) 
		{   
		    current_node = current_node->previous_node_;
		    assert(current_node->type_ == Cell::POSS_FREE);
		    cell_path.push_front(PathSegment(*current_node));
		}
		cell_path.front().p = start_ptr;
		cell_path.back().q  = goal_ptr;
		return true;
	    }
	    current_node->visited_ = astar_timer_;
	    open_cells.push(current_node);
	}
    }
    return false;
}

double og::PCD::cellDistance(Cell& cell1, Cell& cell2)
{
    vector<double> c1 = cell1.centroid_;
    vector<double> c2 = cell2.centroid_;
    double dist = distance(c1, c2);
    return dist;
}

double og::PCD::distance(vector<double> v1, vector<double> v2)
{
    double dist          = 0.0;
    const unsigned int n = v1.size();
    for (unsigned int i = 0; i < n; ++i) 
    {
	dist += (v1[i]-v2[i])<0? -(v1[i]-v2[i]) : (v1[i]-v2[i]);
    }
    return dist;
}

double og::PCD::distance(const ob::State* s1, const ob::State* s2, int i) const
{
    double dist;
    ob::ScopedState<> state1(si_->getStateSpace());
    ob::ScopedState<> state2(si_->getStateSpace());
    state1 = s1;
    state2 = s2;
    dist = (state1[i]-state2[i] > 0)? (state1[i]-state2[i]):-(state1[i]-state2[i]);
    return dist;
}

Cell& og::PCD::getNonMixedCell(Cell& cell, const ob::State* state)
{
    Cell* leaf = &cell;
    while (leaf->type_ == Cell::MIXED) 
    {
	assert(leaf->lo_child_ != 0 && leaf->up_child_ != 0); 
	assert((leaf->lo_child_->contains(state) ^ 
		leaf->up_child_->strictlyContains(state)) || 
	       (leaf->up_child_->contains(state) ^ 
		leaf->lo_child_->strictlyContains(state))); 

	leaf = leaf->lo_child_->contains(state) ? 
	    leaf->lo_child_ : 
	    leaf->up_child_; 
    }
    return *leaf;
}





// ======================= CHECK PATH FUNCTIONS ===============================

bool og::PCD::checkPath(CellPath& cell_path)
{
    if (!parallelCheck(cell_path)) 
    	return false;
    const unsigned int n = cell_path.size();

    for (unsigned int i = 1; i < n; ++i)  
    {
	const ob::State* state = cell_path[i].p;
	if (!isSatisfied(state)) 
	{
	    Cell& cell = *cell_path[i].cell;
	    splitCell(cell, state, split_directions_, pcd_graph_);
	    sampleOccCell(getNonMixedCell(cell, state), 100);
	    return false;
	}
    }
    return serialCheck(cell_path);
}

bool og::PCD::parallelCheck(CellPath& cell_path)
{
    static const double keys[] = {0.5, 0.25, 0.75, 0.125, 0.375, 0.625, 0.875, 0.0625};   
    const unsigned int NUM_CHECKS = (sizeof(keys) / sizeof((keys)[0]));
    const unsigned int n          = cell_path.size();

    ob::State* state =  si_->getStateSpace()->allocState();

    for (unsigned int j = 0; j < NUM_CHECKS; ++j) 
    {
	const double t = keys[j];
	unsigned int i = n;
	// we test the path backwards because it gives shorter planning times
	for  ( ; (i)-- != 0; )  
	{
	    Cell& cell = *cell_path[i].cell;
	    si_->getStateSpace()->interpolate(cell_path[i].p, cell_path[i].q, t, state);      
	    if (isSatisfied(state)) 
	    {
		cell.addSample(state);
	    } else 
	    {
	     	splitCell(cell, state, split_directions_, pcd_graph_);
	     	sampleOccCell(getNonMixedCell(cell, state));
		return false;
	    }
	}
    }
    return true;
}

bool og::PCD::serialCheck(CellPath& cell_path)
{
    // loop over each segment
    unsigned int i = cell_path.size();
    for  ( ; (i)-- != 0; )  
    {
	if (!isSegmentOK(cell_path[i].p, cell_path[i].q, *cell_path[i].cell)) 
	    return false;
    } 
    return true;
}

bool og::PCD::isSegmentOK(const ob::State* from,
			  const ob::State* to,
			  Cell&     cell)
{
    static ob::State* conf_lo;
    static ob::State* conf_hi;

    assert(isSatisfied(from));
    assert(isSatisfied(to));

    assert(cell.contains(from));
    assert(cell.contains(to));
  
    const unsigned int max_samples_to_add = 10;
    unsigned int num_samples_added        = 0;
  
    deque<pair<double, double> > intervals;
    intervals.push_back(make_pair(0.0, 1.0));

    ob::State* conf_low =   si_->getStateSpace()->allocState();
    ob::State* conf_high =  si_->getStateSpace()->allocState();
    si_->getStateSpace()->copyState(conf_low,from);
    si_->getStateSpace()->copyState(conf_high,to);

    double t_low = 0.0;
    double t_high = 1.0;
  
    while (!intervals.empty()) 
    {
	intervals.pop_front();
	if (si_->distance (conf_low, conf_high) >= max_step_size_)
	{
	    const double t_mid = 0.5 * (t_low + t_high);
	    si_->getStateSpace()->interpolate(from,to,t_mid,conf_high);   
	    if (!isSatisfied(conf_high)) 
	    {
		assert(cell.getType() == Cell::POSS_FREE);
		splitCell(cell, conf_high, split_directions_, pcd_graph_);
		sampleOccCell(getNonMixedCell(cell, conf_high), 100);
		return false;
	    } else if (num_samples_added < max_samples_to_add) 
	    {
		cell.addSample(conf_high);
		++num_samples_added;
	    }
	    intervals.push_back(make_pair(t_low, t_mid));
	    intervals.push_back(make_pair(t_mid, t_high));
	}
	if (!intervals.empty()) 
	{
	    t_low = intervals.front().first;  
	    t_high = intervals.front().second;
	    si_->getStateSpace()->interpolate(from, to, t_low, conf_low);
	    si_->getStateSpace()->interpolate(from, to, t_high, conf_high);
	}
    }
  
    return true;
}


// =============================== CELL FUNCTIONS ===================================

void  og::PCD::splitCell(Cell&              cell,
			 const ob::State*   state,
			 vector<bool>       valid_directions,
			 PCD_Graph&         cell_graph)
{
    cout << "* split cell "<< cell.getID() <<"\n";

    const unsigned int  dim = si_->getStateDimension();
    // preconditions

    assert(valid_directions.size() == dim);
    assert(find(valid_directions.begin(), valid_directions.end(), true) != valid_directions.end());
    assert(cell.parent_ == 0 || cell.parent_->type_ == Cell::MIXED);
    assert(cell.lo_child_ == 0 && cell.up_child_ == 0);  
    assert(cell.type_ != Cell::MIXED);
    assert(cell.contains(state));
    assert((cell.type_ == Cell::POSS_FREE && !isSatisfied(state)) ||
	   (cell.type_ == Cell::POSS_FULL &&  isSatisfied(state)));
    assert((cell.contains(start_ptr)   && cell.type_ == Cell::POSS_FREE) || !cell.contains(start_ptr));
    assert((cell.contains(goal_ptr) && cell.type_ == Cell::POSS_FREE) || !cell.contains(goal_ptr));
  
    const bool midpoint_split                = false;
    const Cell::CellType orig_cell_type  = cell.type_;
    const Cell::CellType split_cell_type = orig_cell_type == Cell::POSS_FREE ?
	Cell::POSS_FULL                   :
	Cell::POSS_FREE;

    Cell* split_cell = &cell;   // the cell currently being split
    unsigned int num_iter = 0;


    while (!split_cell->isEmpty()) {
	++num_iter;
        assert(split_cell->contains(state)); // contains the colliding sample
        assert(split_cell->type_ != Cell::MIXED);
	
        // find the sample in the cell closest to config
        const ob::State* nearest  =  findNearestSample(state, split_cell->samples_in_cell_); 
        unsigned int split_dir = 0;

        findSplitDirection(state, nearest, valid_directions, *split_cell, split_dir);
	double split_t = 0.5;

	// find state to divide, can this be done in a nicer way?
	ob::State* mid_state = si_->getStateSpace()->allocState();
	si_->getStateSpace()->interpolate(state,nearest, 0.5, mid_state);
	ob::ScopedState<> split_state(si_->getStateSpace());
	split_state = mid_state;

        const double split_coord = split_state[split_dir];

	Cell*  lo_child = new Cell(split_cell->lower_,
				   split_cell->upper_,
				   *split_cell);
    
	Cell*  up_child = new Cell(split_cell->lower_,
				   split_cell->upper_,
				   *split_cell);
	cell_graph.push_back(lo_child);
	cell_graph.push_back(up_child);

    	// update the bounds of the cells
	lo_child->upper_[split_dir]    = split_coord;
	lo_child->centroid_[split_dir] = 0.5 * (lo_child->lower_[split_dir] + split_coord);
	up_child->lower_[split_dir]    = split_coord;
	up_child->centroid_[split_dir] = 0.5 * (up_child->upper_[split_dir] + split_coord);

	// adjust the boundary conditions if necessary
	lo_child->removeRightBoundary(split_dir);
	up_child->removeLeftBoundary(split_dir);
	assert(areAdjacent(*lo_child, *up_child));
	assert(areAdjacent(*up_child, *lo_child));

	// set the child relations
	split_cell->lo_child_ = lo_child;
	split_cell->up_child_ = up_child;
	assert(split_cell->contains(state)); // contains the colliding sample

	assert((lo_child->contains(state) ^ up_child->strictlyContains(state)) ||
	       (up_child->contains(state) ^ lo_child->strictlyContains(state)));
	if (lo_child->contains(state)) {
	    lo_child->type_ = split_cell_type;
	    up_child->type_ = orig_cell_type;
	} else {
	    up_child->type_ = split_cell_type;
	    lo_child->type_ = orig_cell_type;
	}

	vector<const ob::State*>::iterator iter        = split_cell->samples_in_cell_.begin();
	const vector<const ob::State*>::iterator itend =       split_cell->samples_in_cell_.end();
	while (iter != itend) {
	    assert(split_cell->contains(*iter));
	    assert((lo_child->contains(*iter) ^ up_child->strictlyContains(*iter)) ||
		   (up_child->contains(*iter) ^ lo_child->strictlyContains(*iter)));
      
	    if (lo_child->contains(*iter)) {
		lo_child->samples_in_cell_.push_back(*iter);
	    } else {
		up_child->samples_in_cell_.push_back(*iter);
	    }
	    ++iter;
	}
    
	split_cell->type_ = Cell::MIXED;

	// Reassign the neighbors
	Cell::makeNeighbors(*lo_child, *up_child); 
    
	unsigned int i = split_cell->neighbors_.size();
	for ( ; (i)-- != 0; ) {
	    Cell& neighbor = *split_cell->neighbors_[i];
	    assert(areAdjacent(*split_cell, neighbor));
	    assert(find(neighbor.neighbors_.begin(),
			neighbor.neighbors_.end(),
			split_cell) != neighbor.neighbors_.end());
	    unsigned int j = neighbor.neighbors_.size();
	    for ( ; (j)-- != 0; ) {
		if (neighbor.neighbors_[j] == split_cell) {
		    neighbor.neighbors_.erase(neighbor.neighbors_.begin() + j);
		    break;
		}
	    }
	    assert(find(neighbor.neighbors_.begin(),
			neighbor.neighbors_.end(),
			split_cell) == neighbor.neighbors_.end());
	    assert(areAdjacent(*lo_child, neighbor) || areAdjacent(*up_child, neighbor));
      
	    if (areAdjacent(*lo_child, neighbor)) {
		Cell::makeNeighbors(*lo_child, neighbor);
	    }
	    if (areAdjacent(*up_child, neighbor)) {
		Cell::makeNeighbors(*up_child, neighbor);
	    }
	}
	// A MIXED cell is hereafter only used as a graph component, thus we can
	// remove all dynamically allocated memory associated with it.
	split_cell->strip();
	for (i = 0; i < lo_child->neighbors_.size(); ++i) {
	    assert(lo_child->neighbors_[i]->type_ != Cell::MIXED);
	}
	for (i = 0; i < up_child->neighbors_.size(); ++i) {
	    assert(up_child->neighbors_[i]->type_ != Cell::MIXED);
	}    
	// determine the next cell to split
	split_cell = (lo_child->type_ == split_cell_type)? lo_child : up_child;
    }
  
    // post conditions
    assert(split_cell->type_ != orig_cell_type);
    assert(split_cell->type_ != Cell::MIXED);
    assert(split_cell->parent_ == 0 || split_cell->parent_->type_ == Cell::MIXED);
    assert(split_cell->contains(state));
    assert(split_cell->isEmpty());
    split_cell->addSample(state);
}

bool og::PCD::findSplitDirection(const ob::State*    state1,
				 const ob::State*    state2,
				 vector<bool>&       valid_directions,
				 const Cell&         cell,
				 unsigned int&       split_coord_indx) const
{

    int  dim = si_->getStateDimension();
    bool split_found          = false;
    unsigned int num_failures = 0;
    double max_diff           = 0.0;
   
    while (!split_found && num_failures < 2) 
    {
	for (unsigned int i = 0; i < dim; ++i) 
	{
	    if (valid_directions[i]) 
	    {
		// min value
		const double epsilon = (((1.0e-8) > ( (cell.upper_[i] - cell.lower_[i]) / 1e3))? 
					( (cell.upper_[i] - cell.lower_[i]) / 1e3) : 
					(1.0e-8));
		double 	diff = distance(state1, state2, i) ;
		if (diff > epsilon) 
		{
		    // normalize with the range of this DOF
   		    diff /= ( max_bound_[i] - min_bound_[i]);
   		    if (diff > max_diff) {
   			max_diff         = diff;
   			split_coord_indx = i;
   			split_found      = true;
   		    }
   		}
	    }
	}
	if (!split_found) 
	{
	    ++num_failures;
            // set directions to valid and try again
	    valid_directions.assign(dim, true);
	}
    }
    return split_found;
}


void og::PCD::getAllCells(vector<Cell*>& cells, Cell* node)
{
    if (node->lo_child_ == 0 && node->up_child_ == 0) 
    {
	cells.push_back(node);
    }else if (node->lo_child_ != 0 || node->up_child_ != 0)
    {
	if(node->lo_child_ != 0)
	    getAllCells(cells, node->lo_child_);
	if(node->up_child_ != 0)
	    getAllCells(cells, node->up_child_);
    }

}


// ============================= SAMPLE FUNCITONS ====================================

const ob::State* og::PCD::findNearestSample(const ob::State*  state,
					    const vector<const ob::State*> samples)
{
    assert(!samples.empty());
    double min_dist       = 1.0e38;// infinity in copp...
    const ob::State* nearest =  si_->getStateSpace()->allocState();
    unsigned int i = samples.size();
    for ( ; (i)-- != 0; ) 
    {
	const ob::State* test_config = samples[i];
	double 	diff = si_->distance ( test_config, state) ;
	if (diff < min_dist) 
	{
	    min_dist = diff;
	    nearest  = test_config;
	}	
    }	
    return nearest;
}

bool og::PCD::isSatisfied(const ob::State* state)	
{
    const bool satisfied = si_->isValid(state); 
    if (split_advisor_ != 0) 
    {
	// only use the advisor for configurations that are colliding
	if (!satisfied) 
	{
	    // NOT IMPLEMENTED YET
	    //split_advisor_->GetSplitDirections(split_directions_);
	    split_directions_.assign(split_directions_.size(), true);
	} else {
	    split_directions_.assign(split_directions_.size(), true);
	}
    }
    return satisfied;
}

void og::PCD::sampleOccCells()
{
    vector<Cell*> obstCells;
    obstCells.reserve(pcd_graph_.size());
    // get all occupied cells and reset region_type (used for identification of
    // the start region and and the goal region)
    unsigned int i = pcd_graph_.size();

    for ( ; (i)-- != 0; ) 
    {
	Cell& cell = *pcd_graph_[i];
	if (cell.getType() == Cell::POSS_FULL) 
	{
	    obstCells.push_back(&cell);
	}

	// when sampling, the previous A*-search was unsuccessful and 
	// therefore the complete start region is closed
	cell.region_type_ = cell.closed_ == astar_timer_ ?
	    Cell::REG_START        :
	    Cell::REG_UNSPEC;
    }
    // shuffle the obstacle cells
    random_shuffle(obstCells.begin(), obstCells.end());

    // identify goal region
    const ob::State*  goal     = goal_ptr;
    Cell* currCell     = &getNonMixedCell(*pcd_graph_.front(), goal);
    currCell->region_type_ = Cell::REG_GOAL;
    deque<Cell*> tempCells;
    tempCells.push_back(currCell);
  
    while (!tempCells.empty()) 
    {
	currCell = tempCells.front();
	tempCells.pop_front();
	i = currCell->neighbors_.size();
	for ( ; (i)-- != 0; ) 
	{
	    Cell& cell = *currCell->neighbors_[i];

	    if ((cell.getType()    == Cell::POSS_FREE) &&
		(cell.region_type_ == Cell::REG_UNSPEC)) 
	    {
		tempCells.push_back(&cell);
		cell.region_type_ = Cell::REG_GOAL;
	    }
	}
    }

    vector<Cell*> bridgeCells;
    bridgeCells.reserve(obstCells.size());
    vector<Cell*>::iterator cellIt   = obstCells.begin();
    vector<Cell*>::iterator cellsEnd = obstCells.end();
    while (cellIt != cellsEnd) 
    {
	unsigned int prod = 1;
	i = (*cellIt)->neighbors_.size();
	for ( ; (i)-- != 0; ) 
	{
	    prod *= (*cellIt)->neighbors_[i]->region_type_;
	    if (prod % Cell::BRIDGE_FACTOR == 0) 
	    {
		bridgeCells.push_back(*cellIt);
		break;
	    }
	}
	++cellIt;
    }
    random_shuffle(bridgeCells.begin(), bridgeCells.end());

    // Use a random number to determine whether we will use bridge-cell sampling or not.
    const unsigned int numObstCells      = obstCells.size();
    const unsigned int numBridgeCells    = bridgeCells.size();
    const unsigned int max_bridge_trials = (rng_.uniform01() < 0.8)? 10 : 0; // 30

    ob::State* state;
    unsigned int num_new_cells  = 0;
    unsigned int num_trials     = 0;
    bool obst_cell_found        = true;

    // sample until a free bridge cell connects the start and goal regions
    while (num_new_cells == 0 && num_trials < max_bridge_trials && obst_cell_found) 
    {
	obst_cell_found = false;
	++num_trials;

	for (i = 0; i < numBridgeCells && num_new_cells <= max_new_free_cells_; ++i) 
	{
	    Cell& cell = *bridgeCells[i];

	    if ((cell.getType() == Cell::POSS_FULL) &&
		(distance(cell.lower_, cell.upper_) > 1.0e-14)) 
	    {
		obst_cell_found = true;
		state = cell.getSample(rng_); 
		if (isSatisfied(state)) 
		{
		    splitCell(cell, state, split_directions_, pcd_graph_);
		    ++num_new_cells;
		    return;
		}
	    } else 
	    {
		cell.addSample(state, Cell::NO_CHECK);
	    }
	}
    }
    obst_cell_found = true;
    while (num_new_cells == 0 && obst_cell_found) 
    {
    	obst_cell_found = false;
	
    	for (i = 0; i < numObstCells && num_new_cells < max_new_free_cells_; ++i) 
	{
    	    Cell& cell = *obstCells[i];
    	    if ((cell.getType() == Cell::POSS_FULL) &&	
    		(distance(cell.lower_, cell.upper_) > 1.0e-14)) 
	    {
    		obst_cell_found = true;
    		state = cell.getSample(rng_); 
    		if (isSatisfied(state)) 
		{
    		    splitCell(cell,  state, split_directions_, pcd_graph_); 
    		    ++num_new_cells;
    		} else 
		{
    		    cell.addSample(state, Cell::NO_CHECK);
    		}
    	    }
    	}
    }
    return;
}

void og::PCD::sampleOccCell(Cell& occ_cell, unsigned int max_num_tries)
{
    unsigned int num_tries = 0;
    ob::State* sample;
    while (num_tries < max_num_tries) 
    {
	// sample from occ_cell
	sample = occ_cell.getSample(rng_);
	if (isSatisfied(sample)) 
	{
	    splitCell(occ_cell, sample, split_directions_, pcd_graph_);
	    return;
	} else 
	{
	    if (!occ_cell.addSample(sample, Cell::NO_CHECK)) 
	    {
		return; // the cell is already full with samples
	    }
	}
	++num_tries;
    }
    return;
}

void og::PCD::sampleFreeCells()
{
    unsigned int num_new_cells = 0;
    bool free_cell_found       = true;
    ob::State* config;
  
    while (num_new_cells == 0 && free_cell_found) 
    {
	free_cell_found = false;
	for (unsigned int i = 0; i < pcd_graph_.size() && num_new_cells < max_new_free_cells_; ++i)
	{
	    Cell& cell = *pcd_graph_[i];
	    if (cell.getType() == Cell::POSS_FREE) 
	    {
		free_cell_found = true;
		if (distance(cell.lower_, cell.upper_) > 1.0e-14) 
		{
		    config = cell.getSample(rng_);
		    if (!isSatisfied(config)) 
		    {
			splitCell(cell, config, split_directions_, pcd_graph_);
			++num_new_cells;
		    } else 
		    {
			cell.addSample(config, Cell::NO_CHECK);
		    }
		}
	    }
	}
    }
    return;
}

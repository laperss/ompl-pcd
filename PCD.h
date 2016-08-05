#ifndef __PCD_H_INCLUDED__
#define __PCD_H_INCLUDED__ 
#include <vector>
#include <deque>
#include <queue>

#include "ompl/geometric/planners/PlannerIncludes.h"
#include "Cell.h"
#include <ompl/base/Planner.h>

using namespace std;
namespace ob = ompl::base;
namespace og = ompl::geometric;

namespace ompl
{
    namespace geometric
    {
	class PCD: public ob::Planner
	{
	public:	
	    PCD(const ob::SpaceInformationPtr &si);
	    virtual void setProblemDefinition(const ob::ProblemDefinitionPtr &pdef);
	    virtual ob::PlannerStatus solve(const ob::PlannerTerminationCondition &ptc);
	    virtual void clear(void);
	    virtual void setup(void);
	    virtual void getPlannerData(ob::PlannerData &data) const;
	    // return the number of collision checks performed
	    int    getNChecks();
	    // fron a vector of states, find the one closest to a given sample
	    const ob::State* findNearestSample(const ob::State*                     state,
					       const std::vector<const ob::State*>  samples);
	    // cell_division and all_cell_divisions_ are used to plot the solution process
	    // and can be removed for final use. 
	    vector<Cell*> cell_division_;
	    vector<vector<Cell*>> all_cell_divisions_;
	    
	    // The bounds of the entire problem
	    struct SpaceBounds {
	        vector<double> min;
	        vector<double> max;
	    } ;
	private:
	    // Performs A* search to find cell path from start to goal
	    bool findCellPath(CellPath& cell_path);
	    // The preferred split directions of a cell
	    std::vector<bool> split_directions_;
	    // The largest step size we consider when searching for obstacles
	    double                   max_step_size_;
	    // Maximum number of new cells
	    unsigned int             max_new_free_cells_;
	    // Number of performed collision checks;
	    int                      collision_checks_;
	    // Check if a path is collision free
	    bool checkPath(CellPath& cell_path);
	    // First part of checkPath, all cells checked at 
	    bool parallelCheck(CellPath& cell_path);
	    // Check if each segment is collision free
	    bool serialCheck(CellPath& cell_path);
	    bool isSegmentOK(const ob::State*   from,
			     const ob::State*   to,
			     Cell*              cell);
	    // Sample among all occupied cells
	    void  sampleOccCells();
	    // Sample from a given occupied cell
	    void  sampleOccCell(Cell* occ_cell, unsigned int max_num_tries = 10);
	    // Sample from a possibly free cell
	    void  sampleFreeCells();
	    // Return free or occupied cell containing given state 
	    Cell* getNonMixedCell(Cell* cell, const ob::State* state);
	    // Get the spatial limits of the problem
	    void  getSpaceLimits(ob::StateSpacePtr stateSpace, SpaceBounds& bounds);
	    // Returns a vector containing free or occupied cells
	    void  getAllCells(vector<Cell*>& cells, Cell* node);
	    // Split a cell into free/possibly occupied cells
	    void  splitCell(Cell*             cell,
			    const ob::State*  config,
			    std::vector<bool> valid_directions,
			    PCD_Graph&        cell_graph);
	    // Get a suitable split direction
	    bool findSplitDirection(const ob::State*      conf1,
				    const ob::State*      conf2,
				    std::vector<bool>&    valid_directions,
				    const Cell*           cell,
				    unsigned int&         split_coord_indx) const;
	    // Check if a sample collides with an obstacle
	    bool   isSatisfied(const ob::State* state);
	    // Calculate the eucledian distance between two vectors
	    double distance(vector<double> v1, vector<double> v2);
	    // Calculate the distance in dimension i
	    double distance(const ob::State* s1, const ob::State* s2, int i) const;

	protected:
	    RNG                           rng_;
	    ob::StateSamplerPtr           sampler_;
	    PCD_Graph                     pcd_graph_; 
	    unsigned int                  astar_timer_;    
	    int                           max_num_it;
	    const  ob::State*             start_ptr;
	    const  ob::State*             goal_ptr;
	    double                        goalBias_;
	    double                        minValidPathFraction_;
	    SpaceBounds                   spaceBounds_;

	};

    } // namespace geometric
} // namespace ompl
#endif

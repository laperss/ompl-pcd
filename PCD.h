#ifndef __PCD_H_INCLUDED__
#define __PCD_H_INCLUDED__ 
#include <vector>
#include <deque>
#include <queue>

#include "ompl/geometric/planners/PlannerIncludes.h"
#include "Cell.h"


using namespace std;
namespace ob = ompl::base;
namespace og = ompl::geometric;

namespace ompl
{
    namespace geometric
    {
	// FORWARD DECLARATION
	// split advisor is not used at this point
	class SplitAdvisor;

	class PCD: public ob::Planner
	{
	public:	
	    PCD(const ob::SpaceInformationPtr &si, vector<Cell*>& cell_division);
	    PCD(const ob::SpaceInformationPtr &si);
	    void splitCell(Cell& cell);
	    virtual ~PCD();
	    virtual void setProblemDefinition(const ob::ProblemDefinitionPtr &pdef);
	    virtual      ob::PlannerStatus solve(const ob::PlannerTerminationCondition &ptc);
	    virtual void clear(void);
	    virtual void setup(void);
	    virtual void getPlannerData(ob::PlannerData &data) const;
	    void freeMemory();
	    void sample_path();

	    bool setSplitAdvisor(const SplitAdvisor& advisor);
	    const ob::State* findNearestSample(const ob::State*                     state,
					       const std::vector<const ob::State*>  samples);
	    vector<Cell*>   cell_division_;

	private:
	    bool                           findCellPath(CellPath& cell_path);
	    SplitAdvisor*                  split_advisor_;
	    std::vector<bool>              split_directions_;
	    std::vector<double>            range_vals_;
	    double                         max_step_size_;
	    unsigned int                   max_new_free_cells_;
	    unsigned int                   max_num_iter_;

	    bool checkPath(CellPath& cell_path);
	    bool parallelCheck(CellPath& cell_path);
	    bool serialCheck(CellPath& cell_path);
	    bool isSegmentOK(const ob::State*   from,
			     const ob::State*   to,
			     Cell&              cell);

	    void  sampleOccCells();
	    void  sampleOccCell(Cell& occ_cell, unsigned int max_num_tries = 10);
	    Cell& getNonMixedCell(Cell& cell, const ob::State* state);
	    void  sampleFreeCells();

	    void  getAllCells(vector<Cell*>& cells, Cell* node);
	    void  splitCell(Cell&             cell,
			    const ob::State*  config,
			    std::vector<bool> valid_directions,
			    PCD_Graph&        cell_graph);

	    bool findSplitDirection(const ob::State*      conf1,
				    const ob::State*      conf2,
				    std::vector<bool>&    valid_directions,
				    const Cell&           cell,
				    unsigned int&         split_coord_indx) const;

	    bool   isSatisfied(const ob::State* state);
	    double cellDistance(Cell& cell1, Cell& cell2);
	    double distance(vector<double> v1, vector<double> v2);
	    double distance(const ob::State* s1, const ob::State* s2, int i) const;

	protected:
	    RNG                           rng_;
	    ob::StateSamplerPtr           sampler_;
	    std::vector<double>           max_bound_;
	    std::vector<double>           min_bound_;
	    PCD_Graph                     pcd_graph_; // structure containing all (current?) cells
	    unsigned int                  astar_timer_;    
	    int                           max_num_it;
	    boost::posix_time::ptime      start_time;
	    const  ob::State*             start_ptr;
	    const  ob::State*             goal_ptr;
	    double                        goalBias_;
	    double                        minValidPathFraction_;
	};

    } // namespace geometric
} // namespace ompl
#endif
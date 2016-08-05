#include "../Cell.h"
#include "../PCD.h"
#include <ompl/config.h>
#include <ompl/base/SpaceInformation.h>
#include <ompl/base/spaces/RealVectorStateSpace.h>
#include <iostream>
#include <fstream>
#include <ctime>

namespace ob = ompl::base;
namespace og = ompl::geometric;

//Random blocks obstables
vector<vector<double>> obst1{{-0.8,-0.2},{0.6, 0.4}};
vector<vector<double>> obst2{{-0.5, -0.5},{-0.3,-0.2}};
vector<vector<double>> obst3{{0.3, 0.4},{0.4,0.9}};
vector<vector<double>> obst4{{0.7, -0.8},{0.9,-0.4}};
vector<vector<double>> obst5{{0.3, -1},{0.5,-0.8}};
vector<vector<double>> obst6{{-1, 0.7},{-0.8,1}};
vector<vector<double>> obst7{{-0.7, 0.8},{-0.5,0.95}};
vector<vector<double>> obst8{{0.0, -0.9},{0.05,-0.3}};
vector<vector<double>> obst9{{-0.95, -0.9},{-0.7, -0.8}};
vector<vector<vector<double>>> obstacles{obst1,obst2,obst3,obst4,obst5,obst6,obst7,obst8,obst9};

bool isStateValid(const ob::State *state){
    const ob::RealVectorStateSpace::StateType *rv = state->as<ob::RealVectorStateSpace::StateType>();
    for (int i = 0; i < obstacles.size();++i)
    {
	vector<vector<double>> obstacle = obstacles[i]; 
	if(rv->values[0]>obstacle[0][0] && 
	   rv->values[0]<obstacle[1][0] && 
	   rv->values[1]>obstacle[0][1] && 
	   rv->values[1]<obstacle[1][1])
	    return false;
    }
    return true;
}


int main()
{
    ob::StateSpacePtr space(new ob::RealVectorStateSpace(2));
    ob::RealVectorBounds bounds(2);
    bounds.setLow(-1);
    bounds.setHigh(1);
    space->as<ob::RealVectorStateSpace>()->setBounds(bounds);
    ob::SpaceInformationPtr si(new ob::SpaceInformation(space));

    si->setStateValidityChecker(std::bind(&isStateValid,std::placeholders::_1));
    ob::ScopedState<> start(space);
    start.random();
    while(!isStateValid(start->as<ob::State>())){
    	start.random();
    }
    ob::ScopedState<> goal(space);
    goal.random();
    while(!isStateValid(goal->as<ob::State>())){
    	goal.random();
    }
    ob::ProblemDefinitionPtr pdef(new ob::ProblemDefinition(si)); 
    pdef->setStartAndGoalStates(start, goal);

    og::PCD* pcd = new og::PCD(si);
    ob::PlannerPtr planner(pcd);  
    planner->setProblemDefinition(pdef);
    planner->setup();

    ob::PlannerStatus solved = planner->solve(1.0);
    if (solved)
    {
	ofstream cell_data;
        time_t now;
	now = time(NULL);
	char filename[45];
	filename[0] = '\0';
	strftime(filename, 45, "./plot_cells/data/cells-%y%m%d%H%M%S.txt", gmtime(&now));
	cell_data.open (filename);
	cell_data << "obstacles: " <<obstacles.size() << "\n";
	for(int i = 0; i < obstacles.size();++i)
	{
	    vector<vector<double>> obstacle = obstacles[i]; 
	    cell_data << obstacle[0][0] << " " 
		      << obstacle[0][1] << " "
		      << obstacle[1][0] << " "
		      << obstacle[1][1]<<"\n"; 
	}
    	ob::PathPtr path = pdef->getSolutionPath();
        path->print(cell_data);
	cell_data << "collision_checks: " << pcd->getNChecks() <<"\n";
	vector<vector<Cell*>> all_cell_divisions = pcd->all_cell_divisions_;
	for (unsigned int i = 0; i< all_cell_divisions.size();++i)
	{
	    cell_data << "Iteration: " << i+1 << "\n";
	    vector<Cell*> cell_division = all_cell_divisions[i];
	    for(unsigned int id = 0; id < cell_division.size(); ++id)
	    {
		cell_division[id]->print(cell_data);
		vector<const ob::State*> samples = cell_division[id]->samples_in_cell_;
		cell_data << "samples: \n";
		for(int i = 0; i<samples.size(); ++i)
		    si->printState(samples[i],cell_data);
	    }
	}
	cell_data.close();
    }
    return 0;
}


//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "InclusionRefinement.h"

#include "libmesh/implicit_system.h"
#include "libmesh/periodic_boundary.h"
#include "libmesh/dof_map.h"

// registerMooseObject("MooseApp", InclusionRefinement);

template <>
InputParameters
	validParams<InclusionRefinement>()
{
    
	InputParameters params = validParams<MeshModifier>();
	params.addClassDescription("Changes the subdomain ID of elements either (XOR) inside or outside "
		"the specified box to the specified ID.");    
	params.addRequiredParam<std::string>("fractureMeshModifier","fractureMeshModifier");
	params.addParam< bool >("doBoundaryRefinement","doBoundaryRefinement");
	params.addParam< std::vector<int> >("refinements","refinements");
	params.addParam<std::string>("outputFileName","outputFileName");
    
    
	return params;
}

InclusionRefinement::InclusionRefinement(const InputParameters & parameters) :
MeshModifier(parameters),
_meshModifierName(getParam<std::string>("fractureMeshModifier")),
_myMeshModifier(_app.getMeshModifier(_meshModifierName.c_str())),
_inclusionsMeshModifier( dynamic_cast<InclusionsMeshModifier const &>(_myMeshModifier) ),
_hasRefinementVector( isParamValid("refinements") ),
_hasOutputFileName( isParamValid("outputFileName") )
{
	if (_hasRefinementVector)
	{
		_refinements=getParam<std::vector<int> >("refinements");
		if ( _refinements.size()%2 != 0 )
		{
			_console<<"The size of refinements has to be even, exiting...\n";
			exit(1);
		}
	}
	else
		_console<<"Maybe print mesh\n";

	if (_hasOutputFileName)
	{
		_outputFileName=getParam<std::string >("outputFileName");
	}
	if ( isParamValid("doBoundaryRefinement") )
	{
		_doBoundary=getParam<bool>("doBoundaryRefinement");
	}
	else
		_doBoundary=false;	
}

void InclusionRefinement::modify()
{		
	_meshBase=&_mesh_ptr->getMesh();
	//Parallel::Communicator const & _pp_comm( _mesh_ptr->getMesh().comm() );
	distributedMesh=dynamic_cast<UnstructuredMesh *>(_meshBase);

	// TEST FOR PERIODIC // ONLY 2D	
    Real xmin=_mesh_ptr->getMinInDimension(0);
	Real xmax=_mesh_ptr->getMaxInDimension(0);
    Real ymin=_mesh_ptr->getMinInDimension(1);
	Real ymax=_mesh_ptr->getMaxInDimension(1);
		
	EquationSystems equation_systems (distributedMesh[0]);
	ImplicitSystem & system = equation_systems.add_system<ImplicitSystem> ("system");
	system.add_variable ("variable", FIRST);
	
	DofMap & dof_map = system.get_dof_map();
    PeriodicBoundary horz(RealVectorValue(xmax-xmin, 0., 0.));
	horz.myboundary = 3;
	horz.pairedboundary = 1;
	dof_map.add_periodic_boundary(horz);
	
	PeriodicBoundary vert(RealVectorValue(0., ymax-ymin, 0.));
	vert.myboundary = 0;
	vert.pairedboundary = 2;
	dof_map.add_periodic_boundary(vert);
	equation_systems.init ();	
	//END TEST
	meshRefinement=new MeshRefinement(distributedMesh[0]);
	meshRefinement[0].set_periodic_boundaries_ptr(dof_map.get_periodic_boundaries());

	if (_hasRefinementVector)
		doRefine(_refinements);

	if (_hasOutputFileName)
	{
		_console<<"\nWriting mesh file\n";
		distributedMesh->write(_outputFileName);
		_console<<"Done!\n";
	}
}

void InclusionRefinement::doAMR()
{
	MeshBase::const_element_iterator     el=distributedMesh->active_elements_begin();
	MeshBase::const_element_iterator end_el=distributedMesh->active_elements_end  ();

	std::cout<<"  Looping over elements\n";
	for ( ; el != end_el ; ++el)
	{
		Elem * elem = *el;
		if (_doBoundary)
		{
			if ( _inclusionsMeshModifier.IsOnBoundary(elem[0]) )
			{
				elem[0].set_refinement_flag(Elem::REFINE);
			}
		}
		else
		{
			if ( _inclusionsMeshModifier.DoesIntersect(elem[0]) )
				elem[0].set_refinement_flag(Elem::REFINE);
		}
	}

	std::cout<<"  Actual refinement\n";
	meshRefinement->refine_elements();
	std::cout<<"  Done!\n";
}

void InclusionRefinement::doUMR(int i)
{
	meshRefinement->uniformly_refine(i);
}


void InclusionRefinement::doRefine(std::vector<int> const & refinements)
{
	std::cout<<"We are going to perform this sequence of refinements\n";
	for (int i=0; i<refinements.size(); ++++i)
	{
		std::cout<<refinements.at(i)<<" AMR and ";
		std::cout<<refinements.at(i+1)<<" UMR\n";
	}

	for (int i=0;i<refinements.size(); ++++i)
	{
		std::cout<<"\nAdaptivty loop " <<i/2<<": "<<refinements.at(i)<<" AMR and "<<refinements.at(i+1)<<" UMR"<<std::endl;
		for (int j=0; j<refinements.at(i); ++j)
		{
			std::cout<<" Doing AMR "<<j+1<<std::endl;
			doAMR();
		}

		std::cout<<" Doing "<<refinements.at(i+1)<<" UMR"<<std::endl;
		if (refinements.at(i+1)>0)
			doUMR(refinements.at(i+1));
	}

}

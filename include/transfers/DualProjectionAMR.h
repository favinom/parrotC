/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#ifndef DUALPROJECTIONAMR_H
#define DUALPROJECTIONAMR_H

// MOOSE includes
#include "MultiAppTransfer.h"

#include "libmesh/petsc_matrix.h"
#include "libmesh/petsc_vector.h"

// Forward declarations
class DualProjectionAMR;

template<>
InputParameters validParams<DualProjectionAMR>();

/**
 * Copy the value to the target domain from the nearest node in the source domain.
 */
class DualProjectionAMR : public MultiAppTransfer
{
public:
    DualProjectionAMR(const InputParameters & parameters);
    
    virtual void initialSetup() override;
    
    virtual void execute() override;
    
    virtual ~DualProjectionAMR();
    
protected:
    /**
     * Return the nearest node to the point p.
     * @param p The point you want to find the nearest node to.
     * @param distance This will hold the distance between the returned node and p
     * @param nodes_begin - iterator to the beginning of the node list
     * @param nodes_end - iterator to the end of the node list
     * @return The Node closest to point p.
     */
    bool is_inside();
    void constructInterpolationMatrix();
    // this method is needed just because there is a bug in the nodal kernel:
    // entries in the hanging nodes are not distributed to the closest entries as it should be
    // We construct this matrix which is responsible for this, maybe in future it has to be removed
    void constructReorganizationMatrix();
    
    std::vector<VariableName> _to_var_name;
    std::vector<VariableName> _from_var_name;
    
    bool const _add;
    bool const _dual;
    bool const _has_hanging_nodes;
    
    Real _xmin,_xmax,_ymin,_ymax,_zmin,_zmax;
    Point _centroid;
    
    
    
    FEProblemBase * _from_problem;
    FEProblemBase * _to_problem;
    
    MeshBase const * _from_mesh;
    MeshBase const * _to_mesh;
    
    MeshBase const * _fine_mesh;
    MeshBase const * _coarse_mesh;
    
    PetscMatrix<Number> * _proj_matrix;
    PetscMatrix<Number> * _reorganization_matrix;
    
    PetscVector<Number> * _from_vector;
    PetscVector<Number> * _to_vector;

    Parallel::Communicator * _my_communicator;
    Parallel::Communicator const * _from_communicator;
    Parallel::Communicator const * _to_communicator;

    unsigned _total_from_nodes;
    unsigned _total_to_nodes;
    
    unsigned _local_from_nodes;
    unsigned _local_to_nodes;
    
};

#endif /* MULTIAPPINTERPOLATIONTRANSFER_H */

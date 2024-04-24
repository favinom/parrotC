#pragma once

// MOOSE includes
#include "MeshModifier.h"
#include "MooseMesh.h"
#include "BoxInclusion.h"

// Forward Declarations
class InclusionsMeshModifier;

template <>
InputParameters validParams<InclusionsMeshModifier>();

/**
 * A user object that runs over all the nodes and does an aggregation
 * step to compute a single value.
 */
class InclusionsMeshModifier : public MeshModifier
{
public:
    InclusionsMeshModifier(const InputParameters & parameters);
    
    void                      modify       ()                                              override;
    bool                      isInside     (RealVectorValue const & point, Real bound=0.0) const;
    std::vector<unsigned int> whichIsInside(RealVectorValue const & point, Real bound=0.0) const;
	bool                      DoesIntersect(Elem const & elem)                             const;
	bool                      IsOnBoundary (Elem const & elem)                             const;
    //
    // std::vector<int> whichIsInside(RealVectorValue const & point,Real bound=0.0) const ;
    //
    // virtual bool isInsideRegion(RealVectorValue const & point, int region, Real & bound) const = 0;
    //
    unsigned int get_fn()       const { return _input_fn; };
    unsigned int get_total_fn() const { return _total_fn; };

    
protected:
    
    int _input_fn;
	int _total_fn;
	
	int _nx_min,_nx_max,_ny_min,_ny_max,_nz_min,_nz_max;
	Real _Lx,_Ly,_Lz;
	
    std::vector<Real> _fx;
    std::vector<Real> _fy;
    std::vector<Real> _fz;
    std::vector<Real> _fa1;
    std::vector<Real> _fa2;
    std::vector<Real> _fa3;
    std::vector<Real> _fd1;
    std::vector<Real> _fd2;
    std::vector<Real> _fd3;
	
	int _fractureDim;
	
	std::vector<BoxInclusion> _boxInclusions;
};

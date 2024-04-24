#pragma once

#include<vector>

#include "libmesh/vector_value.h"
#include "libmesh/tensor_value.h"
#include "libmesh/elem.h"
//#include "libmesh/point.h"

typedef libMesh::RealVectorValue RealVectorValue;
typedef libMesh::RealTensorValue RealTensorValue;
typedef libMesh::Real Real;
typedef libMesh::Elem Elem;
typedef libMesh::Point Point;

class BoxInclusion
{
public:
    BoxInclusion();
	void Initialize   (RealVectorValue center, RealVectorValue angles, RealVectorValue axisSizes, int dim);
	bool IsInside     (RealVectorValue const & point, Real bound=0.0) const;
    bool IsInside     (Elem const & elem, Real bound=0.0)             const;
	bool IsOnBoundary (Elem const & elem)                             const;
	bool DoesIntersect(Elem const & elem)                             const;
	void Print()                                                      const;
	
protected:
	
	bool            _initialized;
	int             _dim;
    RealVectorValue _center;
    RealVectorValue _angles;
	RealVectorValue _axisSizes;
	
	RealTensorValue _rotation;
	std::vector<RealVectorValue> _axis;
	// Here we store the constant of the line defining the axis/plane?
	std::vector<Real> _d;

    std::vector<RealVectorValue> _vertices;
    
    void ComputeNormalsFromAngles();
    bool IsOnBoundary2D(Elem const & elem) const;
    bool IsOnBoundary3D(Elem const & /*elem*/) const;
    
};

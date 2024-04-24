#include "BoxInclusion.h"

#include <cassert>

#include "geohelpers.h"

BoxInclusion::BoxInclusion()
{
	_initialized=false;
	_center=RealVectorValue(0.0,0.0,0.0);
	_angles=RealVectorValue(0.0,0.0,0.0);
	_axisSizes=RealVectorValue(0.0,0.0,0.0);
}

void BoxInclusion::Initialize(RealVectorValue center, RealVectorValue angles, RealVectorValue axisSizes, int dim)
{
	_dim      =dim;
	for (int i=0; i<_dim; ++i)
	{
		_center(i)   =center   (i);
		_angles(i)   =angles   (i);
		_axisSizes(i)=axisSizes(i);
		
	}
	_axis.resize(_dim);
	ComputeNormalsFromAngles();
	
	_d.resize(_dim);
	for (int d=0; d<_dim; ++d)
	{
		_d.at(d)=_axis.at(d)*_center;
	}
	
	if (_dim==2)
		_vertices.resize(4);
	else if (_dim==3)
		_vertices.resize(8);
	
	int counter=0;
	if (_dim==2)
	{
		for (int jj=-1; jj<=1; ++++jj)
			for (int ii=-1; ii<=1; ++++ii)
		{
			RealVectorValue & t=_vertices.at(counter);
			t=RealVectorValue( ii*0.5*_axisSizes(0),jj*0.5*_axisSizes(1),0.0 );
			t=_rotation*t;
			t=t+_center;
			counter++;
		}
	}
	else if (_dim==3)
	{
		for (int kk=-1; kk<=1; ++++kk)
			for (int jj=-1; jj<=1; ++++jj)
				for (int ii=-1; ii<=1; ++++ii)
		{
			RealVectorValue & t=_vertices.at(counter);
			t=RealVectorValue( ii*0.5*_axisSizes(0),jj*0.5*_axisSizes(1),kk*0.5*_axisSizes(2) );
			t=_rotation*t;
			t=t+_center;
			counter++;
		}
	}
	
	_initialized=true;
}

void BoxInclusion::ComputeNormalsFromAngles()
{
	RealTensorValue R1;
	RealTensorValue R2;
	RealTensorValue R3;
    
	R1(0,0)=std::cos(_angles(0));
	R1(0,1)=-std::sin(_angles(0));
	R1(0,2)=0.0;
	R1(1,0)=std::sin(_angles(0));
	R1(1,1)=std::cos(_angles(0));
	R1(1,2)=0.0;
	R1(2,0)=0.0;
	R1(2,1)=0.0;
	R1(2,2)=1.0;
    
	R2(0,0)=std::cos(_angles(1));
	R2(0,1)=0.0;
	R2(0,2)=-std::sin(_angles(1));
	R2(1,0)=0.0;
	R2(1,1)=1.0;
	R2(1,2)=0.0;
	R2(2,0)=std::sin(_angles(1));
	R2(2,1)=0.0;
	R2(2,2)=std::cos(_angles(1));
    
	R3(0,0)=1.0;
	R3(0,1)=0.0;
	R3(0,2)=0.0;
	R3(1,0)=0.0;
	R3(1,1)=std::cos(_angles(2));
	R3(1,2)=-std::sin(_angles(2));
	R3(2,0)=0.0;
	R3(2,1)=std::sin(_angles(2));
	R3(2,2)=std::cos(_angles(2));
    
	_rotation=R1*R2*R3;
	
	for (int d=0; d<_dim; ++d)    
		for (int i=0; i<_dim; ++i)
			_axis.at(d)(i)=_rotation(i,d);
	
}

bool BoxInclusion::IsInside(RealVectorValue const & point, Real bound) const
{
	RealVectorValue localbounds;
	for (int d=0; d<_dim; ++d)
	{
		localbounds(d)=_axisSizes(d)/2.0+bound;
	}
	for (int d=0; d<_dim; ++d)
	{
		Real temp=std::fabs( _axis.at(d)*point-_d.at(d) );
		if (temp>localbounds(d))
		{
			return false;
		}
	}
	return true;
}

bool BoxInclusion::IsInside(Elem const & elem, Real bound) const
{
	assert( elem.n_vertices()==elem.n_nodes () );
	
	for (unsigned int v=0;v<elem.n_vertices(); ++v)
	{
		Point const & point=elem.point(v);
		RealVectorValue rvv(point);
		
		if (! IsInside(rvv, bound))
		{
			return false;
		}
	}
	
	return true;
	
}

bool BoxInclusion::IsOnBoundary(Elem const & elem) const
{
	if (_dim==2)
	{
		return IsOnBoundary2D(elem);
	}
	if (_dim==3)
	{
		return IsOnBoundary3D(elem);
	}
	assert(false);
	return false;
}

bool BoxInclusion::IsOnBoundary2D(Elem const & elem) const
{
	assert( elem.n_vertices()==elem.n_nodes () );
	
	RealVectorValue cmax(-1e9,-1e9,-1e9);
	RealVectorValue cmin( 1e9, 1e9, 1e9);

	for (unsigned int v=0; v<elem.n_vertices(); ++v)
	{
		Point const & point=elem.point(v);
		RealVectorValue rvv(point);
		for (int j=0; j<3; ++j)
		{
			cmax(j)=std::max( cmax(j),rvv(j) );
			cmin(j)=std::min( cmin(j),rvv(j) );
		}
	}

	RealVectorValue const & p0=_vertices.at(0);
	RealVectorValue const & p1=_vertices.at(1);
	RealVectorValue const & p2=_vertices.at(2);
	RealVectorValue const & p3=_vertices.at(3);
	
	if ( doesEdgeIntersectElement_2D(p0, p1, cmin, cmax) )
		return true;
	if ( doesEdgeIntersectElement_2D(p0, p2, cmin, cmax) )
		return true;
	if ( doesEdgeIntersectElement_2D(p1, p3, cmin, cmax) )
		return true;
	if ( doesEdgeIntersectElement_2D(p2, p3, cmin, cmax) )
		return true;
	
	return false;

}

bool BoxInclusion::IsOnBoundary3D(Elem const & /*elem*/) const
{
	assert(false);
	return false;
}

bool BoxInclusion::DoesIntersect(Elem const & elem) const
{
	if ( IsInside(elem) || IsOnBoundary(elem) )
	{
		return true;
	}
	else
	{
		return false;
	}
}

void BoxInclusion::Print()  const
{
	std::cout<<"Center\n"<<_center<<std::endl;
	std::cout<<"Angles\n"<<_angles<<std::endl;
	std::cout<<"AxisSizes\n"<<_axisSizes<<std::endl;
	std::cout<<"Rotation\n"<<_rotation<<std::endl;
	
	std::cout<<"Axis\n";
	for (int i=0; i<_axis.size(); ++i)
	{
		std::cout<<_axis.at(i)<<std::endl;
	}
	
	std::cout<<"d\n";
	for (int i=0; i<_d.size(); ++i)
	{
		std::cout<<_d.at(i)<<std::endl;
	}

	std::cout<<"Vertices\n";
	for (int i=0; i<_vertices.size(); ++i)
	{
		std::cout<<_vertices.at(i)<<std::endl;
	}
	

}

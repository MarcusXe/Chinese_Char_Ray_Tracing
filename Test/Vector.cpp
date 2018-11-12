#include "stdafx.h"
#include "Vector.h"
#include "math.h"

CVector::CVector()//z������
{
	x = 0.0;
	y = 0.0;
	z = 1.0;
}

CVector::~CVector()
{

}

CVector::CVector(const CPi3 &p)
{
	x=p.x;
	y=p.y;
	z=p.z;	
}
CVector::CVector(const CPi3 &p0,const CPi3 &p1)
{
	x=p1.x-p0.x;
	y=p1.y-p0.y;
	z=p1.z-p0.z;
}
CVector::CVector(double x, double y, double z)
{
	this->x = x;
	this->y = y;
	this->z = z;	
}

CVector::CVector(const CP3 &p)
{
	x = p.x;
	y = p.y;
	z = p.z;
}

CVector::CVector(const CP3 &p0, const CP3 &p1)
{
	x = p1.x - p0.x;
	y = p1.y - p0.y;
	z = p1.z - p0.z;
}

double CVector::Magnitude()//ʸ����ģ
{
	return sqrt(x * x + y * y + z * z);
}

CVector CVector::Normalize()//��һ��Ϊ��λʸ��
{
	CVector result;
	double Mag = Magnitude();
	if(fabs(Mag) < 1e-6)
		Mag = 1.0;
	result.x = x / Mag;
	result.y = y / Mag;
	result.z = z / Mag;
	return result;
}

CVector operator + (const CVector &v0, const CVector &v1)//ʸ���ĺ�
{
	CVector result;
	result.x = v0.x + v1.x;
	result.y = v0.y + v1.y;
	result.z = v0.z + v1.z;
	return result;
}

CVector operator - (const CVector &v0, const CVector &v1)//ʸ���Ĳ�
{
	CVector result;
	result.x = v0.x - v1.x;
	result.y = v0.y - v1.y;
	result.z = v0.z - v1.z;
	return result;
}

CVector operator * (const CVector &v, double scalar)//ʸ���볣���Ļ�
{
	CVector result;
	result.x = v.x * scalar;
	result.y = v.y * scalar;
	result.z = v.z * scalar;
	return result;
}

CVector operator * (double scalar, const CVector &v)//������ʸ���Ļ�
{
	CVector result;
	result.x = v.x * scalar;
	result.y = v.y * scalar;
	result.z = v.z * scalar;
	return result;
}

CVector operator / (const CVector &v, double scalar)//ʸ������
{
	if(fabs(scalar) < 1e-6)
		scalar = 1.0;
	CVector result;
	result.x = v.x / scalar;
	result.y = v.y / scalar;
	result.z = v.z / scalar;
	return result;
}

CVector operator += (CVector &v0, CVector &v1)//"+="���������
{
	v0.x += v1.x;
	v0.y += v1.y;
	v0.z += v1.z;
	return v0;
}

CVector operator -= (CVector &v0, CVector &v1)//"-="���������
{
	v0.x -= v1.x;
	v0.y -= v1.y;
	v0.z -= v1.z;
	return v0;
}

CVector operator *= (CVector &v0, CVector &v1)//"*="���������
{
	v0.x *= v1.x;
	v0.y *=  v1.y;
	v0.z *=  v1.z;
	return v0;
}

CVector operator /= (CVector &v, double scalar)//"/="���������
{
	v.x /= scalar;
	v.y /= scalar;
	v.z /= scalar;
	return v;
}

double DotProduct(const CVector &v0, const CVector &v1)//ʸ���ĵ��
{
	return(v0.x * v1.x + v0.y * v1.y + v0.z * v1.z);
}

CVector CrossProduct(const CVector &v0, const CVector &v1)//ʸ���Ĳ��
{
	CVector result;
	result.x = v0.y * v1.z - v0.z * v1.y;
	result.y = v0.z * v1.x - v0.x * v1.z;
	result.z = v0.x * v1.y - v0.y * v1.x;
	return result;
}
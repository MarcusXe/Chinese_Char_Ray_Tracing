#pragma once
#include "P3.h"
#include "Pi3.h"
class CVector
{
public:
	CVector();
	virtual ~CVector();
	CVector(double x, double y, double z);
	CVector(const CP3 &p);
	CVector(const CPi3 &);
	CVector(const CP3 &p0, const CP3 &p1);
	CVector(const CPi3 &,const CPi3 &);	
	double Magnitude();//ʸ����ģ
	CVector Normalize();//��λʸ��
	friend CVector operator+ (const CVector &v0, const CVector &v1);//���������
	friend CVector operator- (const CVector &v0, const CVector &v1);
	friend CVector operator* (const CVector &v, double scalar);
	friend CVector operator* (double scalar, const CVector &v);
	friend CVector operator/ (const CVector &v,double scalar);
	friend CVector operator+= (CVector &v0, CVector &v1);
	friend CVector operator-= (CVector &v0, CVector &v1);
	friend CVector operator*= (CVector &v0, CVector &v1);
	friend CVector operator/= (CVector &v, double scalar);
	friend double  DotProduct(const CVector &v0, const CVector &v1);  //ʸ�����
	friend CVector CrossProduct(const CVector &v0, const CVector &v1);//ʸ�����
public:
	double x,y,z;
};


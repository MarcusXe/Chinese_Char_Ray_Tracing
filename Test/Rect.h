#pragma once
/*
	��Ƭ��
	���ԣ�
		��ķ�ʸ��
		���D
		���x��y��z����ȡֵ��Χ
*/
#include "primitive.h"

class Rect :
	public CPrimitive
{
public:
	Rect(void);
	~Rect(void);
	Rect(CP3 positionP, CMaterial pMaterial, double, double, double, double, double, double);
	bool GetInterPoint(CRay Ray, CInterPoint &InPoint); //��ȡ ֱ�� �� ��Ƭ�Ľ���
public:
	double X_min, X_max;
	double Y_min, Y_max;
	double Z_min, Z_max;
};


#pragma once
/*
	����
	����:
		��ĵ�����
		��İ뾶
*/

#include "Primitive.h"
class CSphere:
	public CPrimitive
{
public:
	CSphere();
	virtual ~CSphere();
	CSphere(double r, CP3 positionP,CMaterial pMaterial);
    bool GetInterPoint(CRay Ray,CInterPoint &InPoint); //��ȡ ֱ�� �� ��Ľ���
};


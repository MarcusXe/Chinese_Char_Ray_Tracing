#pragma once
/*
	������
	����:
		�����λ��
		����Ĳ���
*/
#include "P3.h"
#include"InterPoint.h"
#include"Material.h"
#include"Ray.h"
class CPrimitive
{
public:
	CPrimitive(void);
	CPrimitive(double r, CP3 positionP,CMaterial pMaterial);
	~CPrimitive(void);
	virtual bool GetInterPoint(CRay Ray,CInterPoint &InPoint); //��ȡ ֱ�� �� ����Ľ��� 
public:
	CP3 positionP;//���򷽳��е������ģ����淽���е���������
	double r;//���򷽳��е����뾶�����淽���е���D
	CMaterial pMaterial;//ÿ������Ĳ���
};


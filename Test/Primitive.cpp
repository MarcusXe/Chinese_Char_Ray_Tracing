#include "stdafx.h"
#include "Primitive.h"


CPrimitive::CPrimitive(void)
{
}
CPrimitive::CPrimitive(double r, CP3 positionP,CMaterial pMaterial)
{
	this->r = r;
	this->positionP = positionP;
	this->pMaterial = pMaterial;
} 

CPrimitive::~CPrimitive(void)
{
}
bool CPrimitive::GetInterPoint(CRay Ray,CInterPoint &InPoint) //��ȡ ֱ�� �� ����Ľ��� 
{
	return true;
}
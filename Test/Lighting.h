#pragma once
#include "Vector.h"
#include "Light.h"
#include "Material.h"
#include "Ray.h"
class CLighting
{
public:
	CLighting();
	CLighting(int);
	virtual ~CLighting();
	void SetLightNumber(int);//���ù�Դ����
	CRGB Lighting(CP3,CP3,CVector,CMaterial *,int i);//�������	
public:
	int LightNum;//��Դ����
	CLight *Light;//��Դ����
	CRGB Ambient;//������
};


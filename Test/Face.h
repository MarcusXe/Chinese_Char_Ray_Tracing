#pragma once
#include "Vector.h"
class CFace
{
public:
	CFace();
	virtual ~CFace();
	void SetNum(int ptNum);//������Ķ�����
	void SetFaceNormal(CP3 pt0, CP3 pt1, CP3 pt2);//����С�淨ʸ��
public:
	int ptNum; //��Ķ�����
	int* ptI;//��Ķ�������
	CVector fNormal; //С��ķ�ʸ��
};


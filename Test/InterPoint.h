#pragma once
/*
	����������Ľ�����
	���ԣ�
		����t
		��������
		����ķ�ʸ��
		���㴦�Ĳ���
*/
#include"Vector.h"
#include"Material.h"
class CInterPoint   
{
public:
	CInterPoint();
	virtual ~CInterPoint();

	double t;//���ܲ���T
	CP3 IntersectionPoint; // ֱ��������Ľ�������
    CVector Nformal;//����ķ�ʸ��
	int type;
	CMaterial pMaterial;//���㴦�Ĳ���
};


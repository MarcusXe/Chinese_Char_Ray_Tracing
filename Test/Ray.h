#pragma once
/*
	������
	����:
		���ߵ��������
		���ߵ�ʸ��
*/


#include"Vector.h"
#include"RGB.h"
class CRay  
{
public:
	CRay();
	virtual ~CRay();
    CRay::CRay(CP3 origin, CVector dir);//��ȡ���ߵ���ʼ��͵�λ������ʸ��
    CP3 GetPoint(double t);//��ȡ����������Ľ���
	void Normalized();//��λ������
public:
	double alpha;//����������������
    CP3 origin; //���ߵ����
	CVector dir;   //���ߵķ���
};


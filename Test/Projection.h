#pragma once
/*
	͸��ͶӰ��
*/

#include "P3.h"
#include "Pi3.h"
class CProjection  
{
public:
	CProjection(void);
	~CProjection(void);
	void InitalParameter();
	CP2 PerspectiveProjection(CP3 WorldP);//��ά͸��ͶӰ
	CP2 OrthogonalProjection(CP3 WorldP);//����ͶӰ
	CP2 CavalierProjection(CP3 WorldP);//б�Ȳ�ͶӰ
	CP2 CabinetProjection(CP3 WorldP);//б����ͶӰ
	CPi3 PerspectiveProjection3(CP3 WorldP);//��ά͸��ͶӰ
public:
	double k[9];//͸��ͶӰ����
	double R, Theta, Phi, d;//�ӵ����û�����ϵ�е�������
	double Near, Far;//Զ��������
	CP3 ViewPoint;
};


#pragma once
#include "P3.h"

class CLight  
{
public:
	CLight();
	virtual ~CLight();
	void SetDiffuse(CRGB);//���ù�Դ���������
	void SetSpecular(CRGB);//���ù�Դ�ľ��淴���
	void SetPosition(double,double,double);//���ù�Դ��ֱ������λ��
	void SetGlobal(double,double,double);//���ù�Դ��������λ��
	void SetCoef(double,double,double);//���ù�ǿ��˥��ϵ��
	void SetOnOff(bool);//���ù�Դ����״̬	
	void GlobalToXYZ();//������ת��Ϊֱ������
public:
	CRGB L_Diffuse;//�����������ɫ	
	CRGB L_Specular;//��ľ���߹���ɫ
	CP3  L_Position;//��Դ��λ��
	double L_R,L_Phi,L_Theta;//��Դ������
	double L_C0;//����˥��ϵ��
	double L_C1;//����˥��ϵ��
	double L_C2;//����˥��ϵ��
	bool L_OnOff;//��Դ����
	bool b_Specular;
	bool b_Diffuse;
	bool b_Ambient;
};


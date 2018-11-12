#include "stdafx.h"
#include "Lighting.h"
#include "math.h"
#define  PI  3.14159265
#define  MIN(a,b)  ((a<b)?(a):(b))
#define  MAX(a,b)  ((a>b)?(a):(b))


CLighting::CLighting()
{
	LightNum=1;
	Light=new CLight[LightNum];
	Ambient=CRGB(77,77,77);//������㶨����
}

CLighting::~CLighting()
{
	if(Light)
	{
		delete []Light;
		Light=NULL;
	}
}
void CLighting::SetLightNumber(int lnum)
{
	if(Light)
		delete []Light;
	LightNum=lnum;
	Light=new CLight[lnum];
}

CLighting::CLighting(int lnum)
{
	LightNum=lnum;
	Light=new CLight[lnum];
	Ambient=CRGB(77,77,77);	
}


/*
	��ֲ����յ�RGB���ֲ�����ֻ������������Լ����淴�������
		ViewPoint���ӵ�
		Point������
		Normal�����㷨��
		PMaterial������
		i����Դ�ĸ���
*/
CRGB CLighting::Lighting(CP3 ViewPoint,CP3 Point,CVector Normal,CMaterial *pMaterial,int i)
{	
	CRGB LastC=pMaterial->M_Emit;//��������ɢɫΪ��ʼֵ	
	for(int i=0;i<LightNum;i++)//���Թ�Դ
	{	
		if(Light[i].L_OnOff)
		{		
			CRGB InitC;
			InitC.red=0,InitC.green=0,InitC.blue=0;
			CVector VL(Point,Light[i].L_Position);//VLΪָ���Դ��ʸ��
			double d=VL.Magnitude();//dΪ�⴫���ľ��룬���ڹ�ʸ��VL��ģ
			VL=VL.Normalize();//��λ����ʸ��
			CVector VN=Normal;
			VN=VN.Normalize();//��λ����ʸ��			
			//��1���������������
			double CosTheta=MAX(DotProduct(VL,VN),0);			
			InitC.red+=Light[i].L_Diffuse.red*pMaterial->M_Diffuse.red*CosTheta;
			InitC.green+=Light[i].L_Diffuse.green*pMaterial->M_Diffuse.green*CosTheta;
			InitC.blue+=Light[i].L_Diffuse.blue*pMaterial->M_Diffuse.blue*CosTheta;
			//��2�������뾵�淴���
			CVector VV(Point,ViewPoint);//VVΪ��ʸ��
			VV=VV.Normalize();//��λ����ʸ��
			CVector VH=(VL+VV)/(VL+VV).Magnitude();//VHΪƽ��ʸ��	
			double nHN=pow(MAX(DotProduct(VH,VN),0),pMaterial->M_n);
			InitC.red+=Light[i].L_Specular.red*pMaterial->M_Specular.red*nHN;
			InitC.green+=Light[i].L_Specular.green*pMaterial->M_Specular.green*nHN;
			InitC.blue+=Light[i].L_Specular.blue*pMaterial->M_Specular.blue*nHN;	
			//��3������ǿ˥��
			double c0=Light[i].L_C0;//c0Ϊ����˥������
			double c1=Light[i].L_C1;//c1����˥������
			double c2=Light[i].L_C2;//c2����˥������			
			double f=(1.0/(c0+c1*d+c2*d*d));//��ǿ˥������
			f=MIN(1.0,f);		
			LastC+=InitC*f;		
		}
		else
			LastC+=Point.c;//����������ɫ		
	}
	//��4�������뻷����
	LastC+=Ambient*pMaterial->M_Ambient;
	//��5������ɫ��һ����[0,1]����
	LastC.Normalize();		
	//��6�������������㶥��Ĺ�ǿ��ɫ
	return LastC;
}



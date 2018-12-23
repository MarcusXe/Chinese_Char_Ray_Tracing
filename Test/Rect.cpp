#include "stdafx.h"
#include "Rect.h"
#include "stdafx.h"
#include "math.h"
#define Round(d) int(floor(d+0.5))


Rect::Rect(void)
{
}

Rect::Rect(CP3 positionP, CMaterial pMaterial, double X_min, double X_max, double Y_min, double Y_max, double Z_min, double Z_max)
{
	this->r = 150;
	this->positionP = positionP;
	this->pMaterial = pMaterial;
	this->X_min = X_min;
	this->X_max = X_max;
	this->Y_min = Y_min;
	this->Y_max = Y_max;
	this->Z_min = Z_min;
	this->Z_max = Z_max;
}
Rect::~Rect(void)
{
}
bool Rect::GetInterPoint(CRay ray, CInterPoint &InPoint) //��ȡ ֱ�� �� ��Ľ���
{
	CInterPoint InPoint1;
	bool sign = false;
	//���淨ʸ������ߵĽ���
	double mid = DotProduct(positionP, ray.dir);
	if (mid != 0)//���������н���
	{
		//t=-(N*O+d)/(N*V)
		double ans = -(DotProduct(positionP, ray.origin) + r) / mid;
		//if (ans > 0.00001)
		//{
			InPoint1.t = ans;
			InPoint1.IntersectionPoint = ray.GetPoint(InPoint1.t);   //��������
			InPoint1.Nformal = CVector(positionP);//����ķ�ʸ��
			InPoint1.type = 1;
		//}
		if (Round(InPoint1.IntersectionPoint.x) >= X_min && Round(InPoint1.IntersectionPoint.x) <= X_max
			&& Round(InPoint1.IntersectionPoint.y) >= Y_min && Round(InPoint1.IntersectionPoint.y) <= Y_max
			&& Round(InPoint1.IntersectionPoint.z) >= Z_min && Round(InPoint1.IntersectionPoint.z) <= Z_max)
		{	
				InPoint = InPoint1;
				InPoint.pMaterial = this->pMaterial;
				InPoint.pMaterial.SetAmbient(CRGB(0.0, 55.0, 0.0));//���ʶԻ�����ķ�����
				InPoint.pMaterial.SetDiffuse(CRGB(0.0, 55.0, 0.0));//���ʶԻ�������������ķ��������
				InPoint.pMaterial.SetSpecular(CRGB(0.0, 55.0, 0.0));//���ʶԾ��淴���ķ�����
				InPoint.pMaterial.SetEmit(CRGB(0.0, 55.0, 0.0));//��������ɢ����ɫ
				InPoint.pMaterial.M_n = 50.0;//�߹�ָ��
				InPoint.pMaterial.sigma = 0;
				return true;
		}
	}
	return false;
}


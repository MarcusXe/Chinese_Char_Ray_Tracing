#include "stdafx.h"
#include "Rect.h"
#include "stdafx.h"
#include "math.h"
//��������
#define Round(d) int(floor(d+0.5))


Rect::Rect(void)
{
}

Rect::Rect(CMaterial pMaterial, double X_min, double Z_min, double X_max, double Z_max, double Y_min, double Y_max)
{
	CP3 p0(0, 0, 0);
	CVector c1(X_min, 0, Z_min);
	CVector c2(X_max, 0, Z_max);
	double c4 = (p0.x - X_min)*(p0.x - X_min)/* + (p0.y - Y_min)*(p0.y - Y_min)*/ + (p0.z - Z_min)*(p0.z - Z_min);
	double c3 =sqrt((X_max-X_min)*(X_max - X_min)/*+ (Y_max - Y_min)*(Y_max - Y_min)*/+ (Z_max - Z_min)*(Z_max - Z_min));
	double c5=DotProduct(CVector((p0.x - X_min),0, (p0.z - Z_min)), (c1 - c2))/c3;
	//�浽Զ��ľ���
	double tt = sqrt(c4 - c5 * c5);//sqrt(c4-c5*c5);
	//�浽ԭ���"�ؾ�",������ķ�������˵����λ��λ�÷�����
	this->r = sqrt(c4 - c5 * c5);
	if (X_min <= X_max && Z_min <= Z_max)
		this->r = -tt;
	if (Z_min >= Z_max && X_min <= X_max)
		this->r = tt;
	if (Z_min <= Z_max && X_min >= X_max)
		this->r = -tt;
	if (Z_min > Z_max&&X_min > X_max)
		this->r = -tt;
	//��ȡ����ķ���ʸ��
	this->positionP =CP3(-CVector(Z_min - Z_max, 0, X_max - X_min).Normalize().x,0,-CVector(Z_min - Z_max, 0, X_max - X_min).Normalize().z);
	this->pMaterial = pMaterial;
	this->X_min = X_min>X_max?X_max:X_min;
	this->X_max = X_min < X_max ? X_max : X_min;
	this->Y_min =Y_min;
	this->Y_max = Y_max;
	this->Z_min = Z_min>Z_max?Z_max:Z_min;
	this->Z_max = Z_min < Z_max ? Z_max : Z_min;
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
		double ans =-(DotProduct(positionP, ray.origin) + r) / mid;
		if (ans > 0.00001)
		{
			InPoint1.t = ans;
			InPoint1.IntersectionPoint = ray.GetPoint(InPoint1.t);   //��������
			InPoint1.Nformal = CVector(positionP);//����ķ�ʸ��
			InPoint1.type = 1;
		}
		if (Round(InPoint1.IntersectionPoint.x) >= X_min && Round(InPoint1.IntersectionPoint.x) <= X_max
			&& Round(InPoint1.IntersectionPoint.y) >= Y_min && Round(InPoint1.IntersectionPoint.y) <= Y_max
			&& Round(InPoint1.IntersectionPoint.z) >= Z_min && Round(InPoint1.IntersectionPoint.z) <= Z_max)
		{	
				InPoint = InPoint1;
				InPoint.pMaterial = this->pMaterial;
				InPoint.pMaterial.SetAmbient(CRGB(55.0, 0, 100.0));//���ʶԻ�����ķ�����
				InPoint.pMaterial.SetDiffuse(CRGB(55.0, 0, 100.0));//���ʶԻ�������������ķ��������
				InPoint.pMaterial.SetSpecular(CRGB(55.0, 0, 100.0));//���ʶԾ��淴���ķ�����
				InPoint.pMaterial.SetEmit(CRGB(55.0, 0, 100.0));//��������ɢ����ɫ
				InPoint.pMaterial.M_n = 50.0;//�߹�ָ��
				InPoint.pMaterial.sigma = 0;
				return true;
		}
	}
	return false;
}


#include "stdafx.h"
#include "Rect.h"
#include "stdafx.h"
#include "math.h"
//四舍五入
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
	//面到远点的距离
	double tt = sqrt(c4 - c5 * c5);//sqrt(c4-c5*c5);
	//面到原点的"截距",根据面的方向和两端点相对位置位置分正负
	this->r = sqrt(c4 - c5 * c5);
	if (X_min <= X_max && Z_min <= Z_max)
		this->r = -tt;
	if (Z_min >= Z_max && X_min <= X_max)
		this->r = tt;
	if (Z_min <= Z_max && X_min >= X_max)
		this->r = -tt;
	if (Z_min > Z_max&&X_min > X_max)
		this->r = -tt;
	//获取侧面的方向矢量
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
bool Rect::GetInterPoint(CRay ray, CInterPoint &InPoint) //获取 直线 与 面的交点
{
	CInterPoint InPoint1;
	bool sign = false;
	//曲面法矢量与光线的交点
	double mid = DotProduct(positionP, ray.dir);
	if (mid != 0)//光线与面有交点
	{
		//t=-(N*O+d)/(N*V)
		double ans =-(DotProduct(positionP, ray.origin) + r) / mid;
		if (ans > 0.00001)
		{
			InPoint1.t = ans;
			InPoint1.IntersectionPoint = ray.GetPoint(InPoint1.t);   //交点坐标
			InPoint1.Nformal = CVector(positionP);//交点的法矢量
			InPoint1.type = 1;
		}
		if (Round(InPoint1.IntersectionPoint.x) >= X_min && Round(InPoint1.IntersectionPoint.x) <= X_max
			&& Round(InPoint1.IntersectionPoint.y) >= Y_min && Round(InPoint1.IntersectionPoint.y) <= Y_max
			&& Round(InPoint1.IntersectionPoint.z) >= Z_min && Round(InPoint1.IntersectionPoint.z) <= Z_max)
		{	
				InPoint = InPoint1;
				InPoint.pMaterial = this->pMaterial;
				InPoint.pMaterial.SetAmbient(CRGB(55.0, 0, 100.0));//材质对环境光的反射率
				InPoint.pMaterial.SetDiffuse(CRGB(55.0, 0, 100.0));//材质对环境光和漫反射光的反射率相等
				InPoint.pMaterial.SetSpecular(CRGB(55.0, 0, 100.0));//材质对镜面反射光的反射率
				InPoint.pMaterial.SetEmit(CRGB(55.0, 0, 100.0));//材质自身发散的颜色
				InPoint.pMaterial.M_n = 50.0;//高光指数
				InPoint.pMaterial.sigma = 0;
				return true;
		}
	}
	return false;
}


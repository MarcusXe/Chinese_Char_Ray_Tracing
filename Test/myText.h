#include "Rect.h"
#pragma once
#define Round(d) int(floor(d+0.5))
#define CREATE_FONT(nHeight, sName) CreateFont(\
    nHeight,\
    0,\
    0,\
    0,\
    FW_NORMAL,\
    FALSE,\
    FALSE,\
    0,\
    DEFAULT_CHARSET,\
    OUT_DEFAULT_PRECIS,\
    CLIP_DEFAULT_PRECIS,\
    CLEARTYPE_QUALITY,\
    DEFAULT_PITCH | FF_SWISS,\
    sName\
    )
#include "primitive.h"
/******定义结构体用于活性边表AET和新边表NET***********************************/
typedef struct XET
{
	float x;
	float dx, ymax;
	XET* next;
} AET, NET;

/******定义点结构体point******************************************************/
typedef struct Point
{
	double x;
	double y;
}point;
//闭合多边形
typedef struct Poly
{
	point *p;
	int len;
}poly;
typedef struct Zi
{
	poly ploy[20];
	int count;
}zi;
class myText :
	public CPrimitive
{
private:
	CDC *dc;
	LPVOID lpv = nullptr;
	UINT m_nTTPOLYDataLen;
	TTPOLYCURVE m_stTTPOLYGONHEADERVector;
public:
	//字
	zi Czi;
	wchar_t nChar;
	//大小
	UINT size = 200;
	//宽度
	int w;
	//方向
	Rect* rect[80];
	CVector dit;
	double X_min, X_max;
	double Y_min, Y_max;
	double Z_min, Z_max;
	double Zi_X_min=50000, Zi_X_max=-5000;
	double Zi_Y_min=50000, Zi_Y_max=-5000;
	int mycount=0;

public:
	//myText(CDC &dc,wchar_t CZi, int Size, int w, CP3 Position, CVector dit, CMaterial pMaterial);
	bool GetInterPoint(CRay ray, CInterPoint & InPoint);
	myText(CDC & Dc, wchar_t CZi, UINT Size, int W, CVector Dit, double r, CP3 Position, CMaterial Material, double X_min, double X_max, double Y_min, double Y_max, double Z_min, double Z_max);
	//myText(CDC & Dc, double r, wchar_t CZi, int Size, int W, CP3 Position, CVector Dit, CMaterial Materialdouble, double X_min, double X_max, double Y_min, double Y_max, double Z_min, double Z_max);
	bool pointInPoly(point p, poly Cpoly);
	bool pointInZi(point p, zi Czi);
	bool pointInPolyBian(point p, poly Cpoly);
	bool pointInZiBian(point p, zi Czi);
	//void PolyFill(HDC & hdc, point *polypoint, DWORD POINTNUM);
	int FIXEDToInt(FIXED & fixed);
	void initZi();
	//bool GetInterPoint(CPaintDC & dc);
	myText();
	~myText();
};


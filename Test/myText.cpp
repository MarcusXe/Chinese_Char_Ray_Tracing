#include "stdafx.h"
#include "myText.h"
#include "test.h"
#include "math.h"
#define Round(d) int(floor(d+0.5))

bool myText::GetInterPoint(CRay ray, CInterPoint &InPoint) //��ȡ ֱ�� �� ��Ľ���
{
	CInterPoint InPoint1;
	//bool sign = false;
	//���淨ʸ������ߵĽ���
	double mid = DotProduct(positionP, ray.dir);
	if (mid != 0 || 1)
	{
		//t=-(N*O+d)/(N*V)
		double ans = -(DotProduct(positionP, ray.origin) + r) / mid;

		if (ans > 0.00001)
		{
			InPoint1.t = ans;
			InPoint1.IntersectionPoint = ray.GetPoint(InPoint1.t);  //��������
			InPoint1.Nformal = CVector(positionP);//����ķ�ʸ��
			InPoint1.type = 1;
		}
		point p;
		float tx = InPoint1.IntersectionPoint.x - X_min;
		float ty = InPoint1.IntersectionPoint.y - Y_min;
		float tz = InPoint1.IntersectionPoint.z - Z_min;
		p.x = tx+ DotProduct(positionP, CVector(1.0,0.0,0.0))*ty/2;
		p.y = tz +DotProduct(positionP, CVector(0.0, 0.0, 1.0))*ty/2;
		if (Round(p.x) <= Zi_X_max - Zi_X_min && Round(p.y) <= Zi_Y_max - Zi_Y_min)
			if (pointInZi(p, Czi))
			{
				InPoint = InPoint1;
				InPoint.pMaterial = this->pMaterial;
				{
					InPoint.pMaterial.SetAmbient(CRGB(255.0, 255.0, 0.0));//���ʶԻ�����ķ�����
					InPoint.pMaterial.SetDiffuse(CRGB(125.0, 255.0, 0.0));//���ʶԻ�������������ķ��������
					InPoint.pMaterial.SetSpecular(CRGB(125.0, 255.0, 0.0));//���ʶԾ��淴���ķ�����
					InPoint.pMaterial.SetEmit(CRGB(125.0, 255.0, 0.0));//��������ɢ����ɫ
					InPoint.pMaterial.M_n = 50.0;//�߹�ָ��
					InPoint.pMaterial.sigma = 0;

				}
				return true;
			}

	}
	return false;
}

//
//bool CalPlaneLineIntersectPoint(CVector planeVector, CPi3 planePoint, CVector lineVector, CPi3 linePoint, CPi3 p1, CPi3 p2,CPi3 &Result)
//{
//	bool isIn=false;
//	float vp1, vp2, vp3, n1, n2, n3, v1, v2, v3, m1, m2, m3, t, vpt;
//	vp1 = planeVector.x;
//	vp2 = planeVector.y;
//	vp3 = planeVector.z;
//	n1 = planePoint.x;
//	n2 = planePoint.y;
//	n3 = planePoint.z;
//	v1 = lineVector.x;
//	v2 = lineVector.y;
//	v3 = lineVector.z;
//	m1 = linePoint.x;
//	m2 = linePoint.y;
//	m3 = linePoint.z;
//	vpt = v1 * vp1 + v2 * vp2 + v3 * vp3;
//	//�����ж�ֱ���Ƿ���ƽ��ƽ��
//	if (vpt!=0)
//	{
//		t = ((n1 - m1) * vp1 + (n2 - m2) * vp2 + (n3 - m3) * vp3) / vpt;
//		Result.x = m1 + v1 * t;
//		Result.y = m2 + v2 * t;
//		Result.z = m3 + v3 * t;
//		if ((Result.x - p1.x)*(Result.x - p2.x) + (Result.y - p1.y)*(Result.y - p2.y) > 0)
//			isIn = false;
//	}
//	return isIn;
//}

myText::myText(CDC &Dc, wchar_t CZi, UINT Size, int W, CVector Dit, double r, CP3 Position, CMaterial Material, double X_min, double X_max, double Y_min, double Y_max, double Z_min, double Z_max) :dc(&Dc) {
	this->dit = Dit;
	this->r = r;
	this->positionP = Position;
	this->pMaterial = Material;
	this->X_min = X_min;
	this->X_max = X_max;
	this->Y_min = Y_min;
	this->Y_max = Y_max;
	this->Z_min = Z_min;
	this->Z_max = Z_max;
	this->nChar = CZi;
	this->size = Size;
	this->w = W;
	initZi();
	//	int MaxY = 0;

}
//bool myText::GetInterPoint(CRay ray, CInterPoint &InPoint) {
//	CInterPoint InPoint1;
//	bool isIn = false;
//	CVector planeVertor=;
//	CalPlaneLineIntersectPoint()
//
//}

myText::~myText() {

}
//�жϵ��ڶ���α���
bool myText::pointInPolyBian(point p, poly Cpoly) {
	int px = p.x,
		py = p.y;


	for (int i = 1, j = i - 1, l = Cpoly.len; i <= l; i++, j = i - 1) {
		int sx = Cpoly.p[j].x,
			sy = Cpoly.p[j].y,
			tx = Cpoly.p[i].x,
			ty = Cpoly.p[i].y;

		// �������ζ����غ�
		if ((sx == px && sy == py) || (tx == px && ty == py)) {
			return 1;
		}

		// �ж��߶����˵��Ƿ�����������
		if ((sy < py && ty >= py) || (sy >= py && ty < py)) {
			// �߶��������� Y ������ͬ�ĵ�� X ����
			double x = sx * 1.0 + ((double)py - (double)sy) * ((double)tx - (double)sx) / ((double)ty - (double)sy)*1.0;
			// ���ڶ���εı���
			if (Round(x) == px) {
				return true;
			}
		}
	}
	// ���ߴ�������α߽�Ĵ���Ϊ����ʱ���ڶ������
	return false;
}


bool myText::pointInPoly(point p, poly Cpoly) {
	int px = p.x,
		py = p.y,
		flag = false;

	for (int i = 1, j = i - 1, l = Cpoly.len; i <= l; i++, j = i - 1) {
		int sx = Cpoly.p[j].x -Zi_X_min,
			sy = Cpoly.p[j].y - Zi_Y_min,
			tx = Cpoly.p[i].x - Zi_X_min,
			ty = Cpoly.p[i].y - Zi_Y_min;

		// �������ζ����غ�
		if ((sx == px && sy == py) || (tx == px && ty == py)) {
			return 1;
		}

		// �ж��߶����˵��Ƿ�����������
		if ((sy < py && ty >= py) || (sy >= py && ty < py)) {
			// �߶��������� Y ������ͬ�ĵ�� X ����
			double x = sx * 1.0 + ((double)py - (double)sy) * ((double)tx - (double)sx) / ((double)ty - (double)sy)*1.0;

			// ���ڶ���εı���
			if (x == px) {
				return 1;
			}

			// ���ߴ�������εı߽�
			if (x > px) {
				flag = !flag;
			}
		}
	}
	// ���ߴ�������α߽�Ĵ���Ϊ����ʱ���ڶ������
	return flag;
}


/*�жϵ�������*/
bool myText::pointInZi(point p, zi Czi) {
	bool flag = false;
	for (int i = 0; i <Czi.count; i++)
		if (pointInPoly(p, Czi.ploy[i]))
			flag = !flag;
	return flag;
}
/*�жϵ����ֱ���*/
bool myText::pointInZiBian(point p, zi Czi) {
	for (int i = 0; i < Czi.count; i++)
		if (pointInPolyBian(p, Czi.ploy[i]))
			return true;
	return false;
}

/*�������*/
//void myText::PolyFill(HDC &hdc, point *polypoint, DWORD POINTNUM)
//{
//	/******������ߵ��y����(ɨ�赽�˽���)****************************************/
//	int MaxY = 0;
//	int i;
//	for (i = 0; i < POINTNUM; i++)
//		if (polypoint[i].y > MaxY)
//			MaxY = polypoint[i].y;
//
//	/*******��ʼ��AET��***********************************************************/
//	AET *pAET = new AET;
//	pAET->next = NULL;
//
//	/******��ʼ��NET��************************************************************/
//	NET *pNET[1024];
//	for (i = 0; i <= MaxY; i++)
//	{
//		pNET[i] = new NET;
//		pNET[i]->next = NULL;
//	}
//	/******ɨ�貢����NET��*********************************************************/
//	for (i = 0; i <= MaxY; i++)
//	{
//		for (int j = 0; j < POINTNUM; j++)
//			if (polypoint[j].y == i)
//			{
//				//һ�����ǰ���һ�����γ�һ���߶Σ�������ĵ�Ҳ�γ��߶�
//				if (polypoint[(j - 1 + POINTNUM) % POINTNUM].y > polypoint[j].y)
//				{
//					NET *p = new NET;
//					p->x = polypoint[j].x;
//					p->ymax = polypoint[(j - 1 + POINTNUM) % POINTNUM].y;
//					p->dx = (polypoint[(j - 1 + POINTNUM) % POINTNUM].x - polypoint[j].x) / (polypoint[(j - 1 + POINTNUM) % POINTNUM].y - polypoint[j].y);
//					p->next = pNET[i]->next;
//					pNET[i]->next = p;
//				}
//				if (polypoint[(j + 1 + POINTNUM) % POINTNUM].y > polypoint[j].y)
//				{
//					NET *p = new NET;
//					p->x = polypoint[j].x;
//					p->ymax = polypoint[(j + 1 + POINTNUM) % POINTNUM].y;
//					p->dx = (polypoint[(j + 1 + POINTNUM) % POINTNUM].x - polypoint[j].x) / (polypoint[(j + 1 + POINTNUM) % POINTNUM].y - polypoint[j].y);
//					p->next = pNET[i]->next;
//					pNET[i]->next = p;
//				}
//			}
//	}
//	/******���������»��Ա߱�AET*****************************************************/
//	for (i = 0; i <= MaxY; i++)
//	{
//		//�����µĽ���x,����AET
//		NET *p = pAET->next;
//		while (p)
//		{
//			p->x = p->x + p->dx;
//			p = p->next;
//		}
//		//���º���AET������*************************************************************/
//		//�ϱ�����,���ٿ��ٿռ�
//		AET *tq = pAET;
//		p = pAET->next;
//		tq->next = NULL;
//		while (p)
//		{
//			while (tq->next && p->x >= tq->next->x)
//				tq = tq->next;
//			NET *s = p->next;
//			p->next = tq->next;
//			tq->next = p;
//			p = s;
//			tq = pAET;
//		}
//		//(�Ľ��㷨)�ȴ�AET����ɾ��ymax==i�Ľ��****************************************/
//		AET *q = pAET;
//		p = q->next;
//		while (p)
//		{
//			if (p->ymax == i)
//			{
//				q->next = p->next;
//				delete p;
//				p = q->next;
//			}
//			else
//			{
//				q = q->next;
//				p = q->next;
//			}
//		}
//		//��NET�е��µ����AET,���ò��뷨��Xֵ��������**********************************/
//		p = pNET[i]->next;
//		q = pAET;
//		while (p)
//		{
//			while (q->next && p->x >= q->next->x)
//				q = q->next;
//			NET *s = p->next;
//			p->next = q->next;
//			q->next = p;
//			p = s;
//			q = pAET;
//		}
//		/******��������ɫ***************************************************************/
//
//		p = pAET->next;
//		while (p && p->next)
//		{
//			for (float j = p->x; j <= p->next->x; j++)
//				SetPixel(hdc, static_cast<int>(j), i, RGB(50, 50, 50));
//			p = p->next->next;//���Ƕ˵����
//		}
//	}
//}


/*����ת����ƽ������*/
int myText::FIXEDToInt(FIXED& fixed)
{
	if (fixed.fract >= 0x8000)
		return(fixed.value + 1);
	else
		return(fixed.value);
}

void myText::initZi()
{
	//HDC hDC = dc.GetSafeHdc();
	CFont font;

	//�����
	poly Cpoly;
	//������������θ���
	Czi.count = 0;
	VERIFY(font.CreateFont(size, 0, 0, 0, FW_NORMAL, FALSE, FALSE, 0, ANSI_CHARSET, OUT_DEFAULT_PRECIS, CLIP_DEFAULT_PRECIS, DEFAULT_QUALITY, DEFAULT_PITCH | FF_SWISS, _T("����")));
	CFont *pOldFont = dc->SelectObject(&font);
	//���岢��ʼ���任����
	MAT2 mat2;
	memset(&mat2, 0, sizeof(mat2));
	mat2.eM11.value = 1;
	mat2.eM22.value = 1;

	GLYPHMETRICS metrics;    //�����ַ������Ϣ
	DWORD dwDataSize = 0;	//��ʼ���ַ����ݻ�������С

	//ͨ������GetGlyphOutline()ȷ���洢�ַ��ṹ�Ŀռ�
	dwDataSize = dc->GetGlyphOutline(nChar, GGO_NATIVE, &metrics, 0, NULL, &mat2);
	if ((dwDataSize != 0) && (dwDataSize != GDI_ERROR))
	{
		//���������ַ����ݻ�������С
		LPBYTE pPixels = new BYTE[dwDataSize];
		ASSERT(pPixels != NULL);
		TTPOLYGONHEADER *pHeader = (TTPOLYGONHEADER*)pPixels;
		//��ȡ�����ױ���������
		dwDataSize = dc->GetGlyphOutline(nChar, GGO_NATIVE, &metrics, dwDataSize, pPixels, &mat2);

		while (dwDataSize > 0)
		{
			//�����ַ���������㣬ת������
			int xOld = FIXEDToInt(pHeader->pfxStart.x);
			int yOld = FIXEDToInt(pHeader->pfxStart.y);
			point ploypoint[500];
			//����TTF����ṹ��ȡ�ַ�����
			//::MoveToEx(hDC, iXpos + xOld, iYpos - yOld, NULL);
			ploypoint[0].x = xOld;
			ploypoint[0].y =-yOld;
			TTPOLYCURVE *pCurrentCurve = (TTPOLYCURVE*)(pHeader + 1);
			int remainByte = pHeader->cb - sizeof(TTPOLYGONHEADER);
			int currentSize = 0;
			//�����ѭ��
			while (remainByte > 0)
			{
				CPoint lpPoint[1000];

				//CPoint bezi[2];
				int index;
				//����������ѭ��
				for (index = 0; index < pCurrentCurve->cpfx; ++index)
				{
					lpPoint[index].x = FIXEDToInt(pCurrentCurve->apfx[index].x);
					lpPoint[index].y = -FIXEDToInt(pCurrentCurve->apfx[index].y);

				}
				//������ʱ
				switch (pCurrentCurve->wType)
				{
				case TT_PRIM_LINE:
				case TT_PRIM_QSPLINE:
					for (index = 0; index < pCurrentCurve->cpfx; index++)
					{
						//::LineTo(hDC, lpPoint[index].x, lpPoint[index].y);
						ploypoint[currentSize + index + 1].x = lpPoint[index].x;
						ploypoint[currentSize + index + 1].y = lpPoint[index].y;
					}
					currentSize += pCurrentCurve->cpfx;
					break;
				default:
					break;
				}

				//PolyScan(hDC, ploypoint, pCurrentCurve->cpfx);
				int count = sizeof(TTPOLYCURVE) + (pCurrentCurve->cpfx - 1) * sizeof(POINTFX);
				pCurrentCurve = (TTPOLYCURVE*)((char*)pCurrentCurve + count);

				remainByte -= count;
			}
			//ploypoint[currentSize+1] = ploypoint[0];
			Cpoly.p = (point *)malloc(sizeof(point)*(currentSize + 1));
			memcpy(Cpoly.p, ploypoint, sizeof(point)*(currentSize + 1));
			Cpoly.len = currentSize;
			Czi.ploy[Czi.count] = Cpoly;
			Czi.count++;

			for (int i = 0; i < currentSize + 1; i++) {
				if (ploypoint[i].x > Zi_X_max)
					Zi_X_max = ploypoint[i].x;
				if (ploypoint[i].y > Zi_Y_max)
					Zi_Y_max = ploypoint[i].y;
				if (ploypoint[i].y < Zi_Y_min)
					Zi_Y_min = ploypoint[i].y;
				if (ploypoint[i].x < Zi_X_min)
					Zi_X_min = ploypoint[i].x;
			}
			//PolyScan(hDC, ploypoint, currentSize+2);
			//::LineTo(hDC, iXpos + xOld, iYpos - yOld);
			dwDataSize -= pHeader->cb;
			pHeader = (TTPOLYGONHEADER*)((char*)pHeader + pHeader->cb);
		}

		//for (int i = 0; i < 500; i++)
		//	for (int j = 0; j < 500; j++)
		//		if (pointInZi(point{ i,j }, Czi))
		//		{
		//			//::SetPixel(hDC, i, j, RGB(0, 80, 0));
		//		}
		delete[] pPixels;
		//free(Cpoly.p);
	}
}
myText::myText()
{
}


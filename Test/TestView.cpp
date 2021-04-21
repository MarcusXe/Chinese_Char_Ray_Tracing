
// TestView.cpp : CTestView ���ʵ��
//

#include "stdafx.h"
// SHARED_HANDLERS ������ʵ��Ԥ��������ͼ������ɸѡ�������
// ATL ��Ŀ�н��ж��壬�����������Ŀ�����ĵ����롣
#ifndef SHARED_HANDLERS
#include "Test.h"
#endif

#include "TestDoc.h"
#include "TestView.h"
#include "math.h"
#include "myText.h"
#include "Rect.h"
#define  PI 3.1415926
#define  MIN(a,b)  ((a<b)?(a):(b))
#define  MAX(a,b)  ((a>b)?(a):(b))
#define Round(d) int(floor(d+0.5))


#ifdef _DEBUG
#define new DEBUG_NEW
#endif


// CTestView

IMPLEMENT_DYNCREATE(CTestView, CView)

BEGIN_MESSAGE_MAP(CTestView, CView)
	// ��׼��ӡ����
	ON_COMMAND(ID_FILE_PRINT, &CView::OnFilePrint)
	ON_COMMAND(ID_FILE_PRINT_DIRECT, &CView::OnFilePrint)
	ON_COMMAND(ID_FILE_PRINT_PREVIEW, &CView::OnFilePrintPreview)
END_MESSAGE_MAP()

// CTestView ����/����

CTestView::CTestView()
{
	// TODO: �ڴ˴���ӹ������
	Positionlight[0] = CP3(0,0,-200);//���ù�Դλ������
	Positionlight[1] = CP3(-0,200,200);//���ù�Դλ������
	nLightCount =2;
	 
 	pLight=new CLighting(nLightCount);//һά��Դ��̬����
	pLight->Light[0].SetPosition(Positionlight[0].x,Positionlight[0].y,Positionlight[0].z);//���ù�Դλ������
	pLight->Light[0].SetPosition(Positionlight[1].x,Positionlight[1].y,Positionlight[1].z);//���ù�Դλ������	
	for(int i=0;i<nLightCount;i++)
	{
		pLight->Light[0].L_Diffuse = CRGB(255.0, 255.0, 255.0); //��Դ����������ɫ	
		pLight->Light[0].L_Specular = CRGB(255.0, 255.0, 255.0);//��Դ����߹���ɫ
		pLight->Light[0].L_C0 = 1.0;//����˥��ϵ��
		pLight->Light[0].L_C1 = 0.0000001;//����˥��ϵ��
		pLight->Light[0].L_C2 = 0.00000001;//����˥��ϵ��
		pLight->Light[0].L_OnOff = true;//��Դ����	

	}

	pMaterial[0].SetAmbient(CRGB(0.1,0.1,0.1));//���ʶԻ�����ķ�����
	pMaterial[0].SetDiffuse(CRGB(0.2,0.2,0.2));//���ʶԻ�������������ķ��������
	pMaterial[0].SetSpecular(CRGB(0.99,0.99,0.99));//���ʶԾ��淴���ķ�����
	pMaterial[0].SetEmit(CRGB(0.0,0.0,0.0));//��������ɢ����ɫ
	pMaterial[0].M_n=100.0;//�߹�ָ��
	pMaterial[0].sigma = 1.1;

	pMaterial[1].SetAmbient(CRGB(0.0,0.1,0.0));//���ʶԻ�����ķ�����
	pMaterial[1].SetDiffuse(CRGB(0.0,0.508,0.0));//���ʶԻ�������������ķ��������
	pMaterial[1].SetSpecular(CRGB(0.0,0.508,0.0));//���ʶԾ��淴���ķ�����
	pMaterial[1].SetEmit(CRGB(0,0,0));//��������ɢ����ɫ
	pMaterial[1].M_n=50.0;//�߹�ָ��
	pMaterial[1].sigma = 0;

	pMaterial[2].SetAmbient(CRGB(0.1,0.0,0.0));//���ʶԻ�����ķ�����
	pMaterial[2].SetDiffuse(CRGB(1.0,0.0,0.0));//���ʶԻ�������������ķ��������
	pMaterial[2].SetSpecular(CRGB(1.0,0.2,0.2));//���ʶԾ��淴���ķ�����
	pMaterial[2].SetEmit(CRGB(1.0,0,0));//��������ɢ����ɫ
	pMaterial[2].M_n=50.0;//�߹�ָ��
	pMaterial[2].sigma = 0;

	pMaterial[3].SetAmbient(CRGB(0.1,0.1,0.0));//���ʶԻ�����ķ�����
	pMaterial[3].SetDiffuse(CRGB(1.0,1.0,0.0));//���ʶԻ�������������ķ��������
	pMaterial[3].SetSpecular(CRGB(1.0,0.1,0.1));//���ʶԾ��淴���ķ�����
	pMaterial[3].SetEmit(CRGB(0.5,0.5,0));//��������ɢ����ɫ
	pMaterial[3].M_n=20.0;//�߹�ָ��
	pMaterial[3].sigma = 0;

	is_pppp = true;
}

CTestView::~CTestView()
{
	delete pLight;
	for(int i=0; i<7; i++)
	{
		delete objects[i];
	}
}

BOOL CTestView::PreCreateWindow(CREATESTRUCT& cs)
{
	// TODO: �ڴ˴�ͨ���޸�
	//  CREATESTRUCT cs ���޸Ĵ��������ʽ

	return CView::PreCreateWindow(cs);
}

// CTestView ����

void CTestView::OnDraw(CDC* pDC)
{
	CTestDoc* pDoc = GetDocument();
	ASSERT_VALID(pDoc);
	if (!pDoc)
		return;

	// TODO: �ڴ˴�Ϊ����������ӻ��ƴ���
	DoubleBuffer(pDC);
}


// CTestView ��ӡ

BOOL CTestView::OnPreparePrinting(CPrintInfo* pInfo)
{
	// Ĭ��׼��
	return DoPreparePrinting(pInfo);
}

void CTestView::OnBeginPrinting(CDC* /*pDC*/, CPrintInfo* /*pInfo*/)
{
	// TODO: ��Ӷ���Ĵ�ӡǰ���еĳ�ʼ������
}

void CTestView::OnEndPrinting(CDC* /*pDC*/, CPrintInfo* /*pInfo*/)
{
	// TODO: ��Ӵ�ӡ����е��������
}


// CTestView ���

#ifdef _DEBUG
void CTestView::AssertValid() const
{
	CView::AssertValid();
}

void CTestView::Dump(CDumpContext& dc) const
{
	CView::Dump(dc);
}

CTestDoc* CTestView::GetDocument() const // �ǵ��԰汾��������
{
	ASSERT(m_pDocument->IsKindOf(RUNTIME_CLASS(CTestDoc)));
	return (CTestDoc*)m_pDocument;
}
#endif //_DEBUG


// CTestView ��Ϣ�������
void CTestView::DoubleBuffer(CDC* pDC)//˫����
{
	CRect rect;//�ͻ�������
	GetClientRect(&rect);//��ÿͻ����Ĵ�С
	Width = rect.Width();
    Height = rect.Height();
	pDC->SetMapMode(MM_ANISOTROPIC);//pDC�Զ�������ϵ
	pDC->SetWindowExt(rect.Width(), rect.Height());//���ô��ڷ�Χ
	pDC->SetViewportExt(rect.Width(), -rect.Height());//x��ˮƽ���ң�y�ᴹֱ����
	pDC->SetViewportOrg(rect.Width() / 2,rect.Height() / 2);//�ͻ�������Ϊԭ��
	CDC memDC;//�ڴ�DC	
	memDC.CreateCompatibleDC(pDC);//��������ʾpDC���ݵ�memDC
	CBitmap NewBitmap, *pOldBitmap;//�ڴ��г���ͼ�����ʱλͼ 
	NewBitmap.CreateCompatibleBitmap(pDC, rect.Width(), rect.Height());//��������λͼ 
	pOldBitmap = memDC.SelectObject(&NewBitmap);//������λͼѡ��memDC
	rect.OffsetRect(-rect.Width() / 2, -rect.Height() / 2);
	memDC.SetMapMode(MM_ANISOTROPIC);//memDC�Զ�������ϵ
	memDC.SetWindowExt(rect.Width(), rect.Height());
	memDC.SetViewportExt(rect.Width(), -rect.Height());
	memDC.SetViewportOrg(rect.Width() / 2,rect.Height() / 2);
	//pDC->FillSolidRect(rect,pDC->GetBkColor());//��ԭ���������ͻ����������Ǻ�ɫ
	pDC->FillSolidRect(rect,RGB(0,0,0));

	DrawObject(pDC);//���ƶ���
	
//	pDC->BitBlt(rect.left, rect.top, rect.Width(), rect.Height(), &memDC, -rect.Width() / 2, -rect.Height() / 2, SRCCOPY);//���ڴ�memDC�е�λͼ��������ʾpDC��
	memDC.SelectObject(pOldBitmap);//�ָ�λͼ
	NewBitmap.DeleteObject();//ɾ��λͼ
	memDC.DeleteDC();//ɾ��memDC
}

void CTestView::DrawObject(CDC* pDC)
{
	objects[5] = new myText(*pDC, L'��', (UINT)150, 10, CVector(0, 0, 0), 200, CP3(0.0,1.0,0.0), pMaterial[0],0, 400, -200, -200,5, 600);
	objects[1] = new myText(*pDC, L'��', (UINT)150, 10, CVector(0, 0, 0), 200, CP3(0.0, 1.0, 0.0), pMaterial[0], 0, 400, -200, -200, 150, 600);
	objects[0] = new CPlane(300,CP3(0.0,1.0,0.0),pMaterial[0],-400,400,-300,-300,-1000,200);//�����ǽ
	//objects[0] = new myText(*pDC, L'3', (UINT)150,50, CVector(0,0,0), 0, CP3(0.0,0.707,0.707), pMaterial[0], 0, 300, 0, 300, 0, 300);
	objects[2] = new CSphere(70,CP3(20,-110, 80),pMaterial[0]);//��ȡ��İ뾶��λ�ã�����
	//objects[3] = new Rect(pMaterial[0], -100, 50, 200, -250, -50, 50);
	objects[3] = new myText(*pDC, L'X', (UINT)150, 50, CVector(0, 0, 0), 200, CP3(0.0, 1.0, 0.0), pMaterial[0], 150, 400, -200, -200, 5, 600);
	//objects[6] = new myText(*pDC, L'3', (UINT)150, 10, CVector(0, 0, 0), 200, CP3(0.0, 1.0, 0.0), pMaterial[0], 10, 400, -200, -200, 150, 600);
	objects[4] = new CSphere(100,CP3(-20,-70, 00),pMaterial[1]);//��ȡ��İ뾶��λ�ã�����
	//objects[5] = new CSphere(50,CP3(-150,-150,-50),pMaterial[2]);//��ȡ��İ뾶��λ�ã�����
	//objects[1] = new CSphere(90,CP3(150,-120, -500),pMaterial[3]);//��ȡ��İ뾶��λ�ã�����
	nObjectCount =6;

	for (int i = -150; i< Width/2+200; i++)//��������
	{
		for (int j = -Height/2-200; j< Height/2+200; j++)
		{
			is_pppp = true;
			CVector EyeToScreenP(projection.ViewPoint,CP3(i,j,-200)); //����Ļ�� ĳһ�ӵ�  �������ĵ� ʸ��	
			CRay ray(projection.ViewPoint,EyeToScreenP); //���ߵ� ��㣬 ����ʸ��
			ray.alpha = -1;
			ray.Normalized();  //��λ����ʸ��

			CRGB color = Trace(ray,4);   //��������ʸ��		

			CP2 p = projection.PerspectiveProjection(ppppp);//͸��ͶӰ
			color.Normalize();
			pDC->SetPixelV(Round(p.x-200),Round(p.y+300),RGB(color.red,color.green,color.blue));
		}
	}
}
CRGB CTestView::Trace(CRay ray, int count) // ���ӵ� ������ ��һ�� ʸ��  �� �ݹ����
{
	CRGB ret(0.0,0.8,0.6); //��ʼ ��ɫ ����ɫ
	if(!count)//�ݹ����Ϊcount
	{
		return ret;
	}
	CInterPoint min_hit; 
	for (int i = 0; i<nObjectCount; i++)// ѭ���������壬�ҵ� ��������������ཻ ��һ��
	{
		CInterPoint hit;
		if (objects[i]->GetInterPoint(ray,hit) && hit.t > 0.00001 &&hit.t < min_hit.t)
		{
			min_hit = hit; //�ҵ�����Ľ��㱣��
		}
	}

	if (min_hit.t == 100000)  //û���ҵ�
	{
		return ret;
	}
	else
	{
		if (is_pppp)//����ÿһ�θտ�ʼ�Ľ������꣬����͸��ͶӰ
		{
			ppppp = min_hit.IntersectionPoint;
		}
		is_pppp = FALSE;
		return Shade(ray, min_hit, count);  //�ҵ��˹���������ĵ�һ������min_hit
	}
}
CRGB CTestView::Shade(CRay ray, CInterPoint  hit, int count)//��������ʸ��  ���� �ݹ����
{
	CRGB ret(0.0,0.2,0.0);
	CP3 point = hit.IntersectionPoint;  //����������Ľ���
	for (int i = 0; i < nLightCount; i++)
	{
		CRay shadow_ray(point, CVector(point,pLight->Light[i].L_Position)); //���㵽 ��Դ�� ����
		shadow_ray.Normalized();//����λ��
	                                             
		CInterPoint min_hit;
		for (int j = 0; j<nObjectCount; j++)//�ж��� ����� ��Դ֮���Ƿ��������赲
			{
				CInterPoint hit1;
				if (objects[j]->GetInterPoint(shadow_ray, hit1) && hit1.t > 0.00001 && hit1.t < min_hit.t)
				{
					min_hit = hit1;	

				}
			}
		//�������Դ�ľ���
		double dist = sqrt(pow(pLight->Light[i].L_Position.x-point.x,2)+ 
			pow(pLight->Light[i].L_Position.y-point.y,2)+ 
			pow(pLight->Light[i].L_Position.z-point.z,2) );


		if (min_hit.t >= dist)
		{
			 ret += pLight->Lighting(ray.origin ,point, hit.Nformal,&hit.pMaterial, i);
		}
	}

///////////	
	////����ݹ�
	CRGB TotalRGB ;
	CRGB s = Trace(CalculateReflection(ray, hit), count - 1); //����
	CRGB t(0.0,0.0,0.0);
	t = Trace(CalculateRefraction(ray, hit), count - 1);
	if(hit.pMaterial.sigma != 0)//����͸��
	{
		TotalRGB = 0.8 * t + (1 - 0.5 - 0.08) * ret + 0.08 * s;
	}
	else
	{
		TotalRGB = 0.02 * t + 0.280 * s + 0.7 * ret ;
	}
	return TotalRGB;
}     

CRay CTestView::CalculateReflection(CRay in, CInterPoint hit) //��֪������� �ͽ��� ���� ����
{
	CRay ret;
	CVector VL(-in.dir.x, -in.dir.y, -in.dir.z);//��λ��֮����������ʸ��
	CVector VN(hit.Nformal.x, hit.Nformal.y, hit.Nformal.z);//���㷨����ʸ��
	VN = VN.Normalize();//��λ��
	CVector R = 2.0 * VN * fabs(DotProduct(VN, VL)) - VL;

	ret.dir.x = R.x;
	ret.dir.y = R.y;
	ret.dir.z = R.z;//�¹��ߵķ���ʸ��
	ret.origin = hit.IntersectionPoint;//�¹��ߵ���ʼ��
	ret.Normalized();

	return ret;
}
CRay CTestView::CalculateRefraction(CRay in, CInterPoint hit) //����
{
	CRay ret;
	hit.Nformal = hit.Nformal.Normalize();//���㷨����ʸ����λ��

	double angle1 = -(in.dir.x*hit.Nformal.x + in.dir.y*hit.Nformal.y + in.dir.z*hit.Nformal.z);
	double angle2 = sqrt(1-((1-angle1*angle1)/(hit.pMaterial.sigma*hit.pMaterial.sigma)));
	
	double tmp = angle1/hit.pMaterial.sigma - angle2;

	ret.dir.x = tmp*hit.Nformal.x + in.dir.x/hit.pMaterial.sigma;//�¹��ߵķ���ʸ��
	ret.dir.y = tmp*hit.Nformal.y + in.dir.y/hit.pMaterial.sigma;
	ret.dir.z = tmp*hit.Nformal.z + in.dir.z/hit.pMaterial.sigma;

	ret.origin = hit.IntersectionPoint;//�¹��ߵ���ʼ��

	ret.Normalized();
	return ret;
} 
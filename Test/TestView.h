
// TestView.h : CTestView ��Ľӿ�
//

#pragma once
#include "Sphere.h"
#include "ZBuffer.h"
#include "Lighting.h"
#include "Material.h"
#include "Projection.h"
#include "Plane.h"

class CTestView : public CView
{
protected: // �������л�����
	CTestView();
	DECLARE_DYNCREATE(CTestView)

// ����
public:
	CTestDoc* GetDocument() const;

// ����
public:
	void DoubleBuffer(CDC* pDC);//˫����
	void DrawObject(CDC* pDC);//��������

    CRay CalculateReflection(CRay in, CInterPoint hit); //��֪������� �ͽ��� ���� ����
    CRGB Shade(CRay ray, CInterPoint  hit, int count);//��������ʸ��  ���� �ݹ����
    CRGB Trace(CRay ray, int count); // ���ӵ� ������ ��һ�� ʸ��  �� �ݹ����
    CRay CalculateRefraction(CRay in, CInterPoint hit);

// ��д
public:
	virtual void OnDraw(CDC* pDC);  // ��д�Ի��Ƹ���ͼ
	virtual BOOL PreCreateWindow(CREATESTRUCT& cs);
protected:
	virtual BOOL OnPreparePrinting(CPrintInfo* pInfo);
	virtual void OnBeginPrinting(CDC* pDC, CPrintInfo* pInfo);
	virtual void OnEndPrinting(CDC* pDC, CPrintInfo* pInfo);

// ʵ��
public:
	virtual ~CTestView();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

protected:
	int	 LightNum;//��Դ����
	CLighting *pLight;//���ջ���
	CMaterial pMaterial[4];//�������
    CProjection  projection;
	int Width;
	int Height;
	CP3 ViewPoint;
	CPrimitive *objects[7];
	int nObjectCount;
	CP3 Positionlight[2];
	int nLightCount;
	CP3 ppppp;
	bool is_pppp;

// ���ɵ���Ϣӳ�亯��
protected:
	DECLARE_MESSAGE_MAP()
};

#ifndef _DEBUG  // TestView.cpp �еĵ��԰汾
inline CTestDoc* CTestView::GetDocument() const
   { return reinterpret_cast<CTestDoc*>(m_pDocument); }
#endif


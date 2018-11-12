#pragma once
#include "Bucket.h"//Ͱͷ�ļ�
#include "Vector.h"//ʸ��ͷ�ļ�
class CZBuffer
{
public:
	CZBuffer();
    virtual ~CZBuffer();
	void CreateBucket();//����Ͱ
	void CreateEdge();//�߱�
	void AddEt(CAET *);//�ϲ�ET��
	void ETOrder();
	void Gouraud(CDC *);//���
	void InitDeepBuffer(int,int,double);//��ʼ����Ȼ�����
    CRGB Interpolation(double,double,double,CRGB,CRGB);//���Բ�ֵ
	void SetPoint(CPi3 *,int);
	void ClearMemory();//�����ڴ�
	void DeleteAETChain(CAET* pAET);//ɾ���߱�
protected:
	int PNum;//�������
	CPi3 *P;//��������
	CAET *pHeadE,*pCurrentE,*pEdge;//��Ч�߱���ָ��
	CBucket *pCurrentB,*pHeadB;
	double **zBuffer;//����ȳ���
	int Width,Height;//����������
};


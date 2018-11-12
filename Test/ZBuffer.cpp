#include "stdafx.h"
#include "ZBuffer.h"
#include "math.h"
#define Round(d) int(floor(d+0.5))//��������궨��

CZBuffer::CZBuffer()
{
	P=NULL;
	pHeadE=NULL;
	pCurrentB=NULL;
	pEdge=NULL;
	pCurrentE=NULL;
	pHeadB=NULL;
	zBuffer=NULL;
}

CZBuffer::~CZBuffer()
{
	if(P!=NULL)
	{
		delete []P;
		P=NULL;
	}
	for(int i=0;i<Width;i++)
	{
		delete[] zBuffer[i];
		zBuffer[i]=NULL;
	}
	if(zBuffer!=NULL)
	{
		delete zBuffer;
		zBuffer=NULL;
	}	
	ClearMemory();
}

void CZBuffer::SetPoint(CPi3 p[],int m)
{
	if(P!=NULL)
	{
		delete []P;
		P=NULL;
	}
	P=new CPi3[m];
    for(int i=0;i<m;i++)
		P[i]=p[i];	
	PNum=m;
}

void CZBuffer::CreateBucket()//����Ͱ��
{
	int yMin,yMax;
	yMin=yMax=Round(P[0].y);
	for(int i=1;i<PNum;i++)//���Ҷ���������ǵ���С�����ɨ����
	{
		if(P[i].y<yMin)
			yMin=Round(P[i].y);//ɨ���ߵ���Сֵ
		if(P[i].y>yMax)
			yMax=Round(P[i].y);//ɨ���ߵ����ֵ
	}
	for(int y=yMin;y<=yMax;y++)
	{
		if(yMin==y)//����Ͱͷ���
		{
			pHeadB=new CBucket;//����Ͱ��ͷ���
			pCurrentB=pHeadB;//pCurrentBΪCBucket��ǰ���ָ��
			pCurrentB->ScanLine=yMin;
			pCurrentB->pET=NULL;//û�����ӱ߱�
			pCurrentB->pNext=NULL;
		}
		else//����ɨ����
		{
			pCurrentB->pNext=new CBucket;//����Ͱ���������
			pCurrentB=pCurrentB->pNext;
			pCurrentB->ScanLine=y;
			pCurrentB->pET=NULL;
			pCurrentB->pNext=NULL;
		}
	}
}

void CZBuffer::CreateEdge()//�����߱�
{
	for(int i=0;i<PNum;i++)
	{
		pCurrentB=pHeadB;
	//	int j=(i+1)%PNum;//�ߵĵڶ������㣬P[i]��P[j]���ɱ�
	   int j = i+1;
	   if(i == PNum-1) j=0;
		if(P[i].y<P[j].y)//�ߵ��յ������
		{
			pEdge=new CAET;
			pEdge->x=P[i].x;//����ET����ֵ
			pEdge->yMax=Round(P[j].y);
			pEdge->k=(P[j].x-P[i].x)/(P[j].y-P[i].y);//����1/k
			pEdge->ps=P[i];//�󶨶������ɫ
			pEdge->pe=P[j];
			pEdge->pNext=NULL;
			while(pCurrentB->ScanLine!=P[i].y)//��Ͱ��Ѱ�Ҹñߵ�yMin
			{
				pCurrentB=pCurrentB->pNext;//�Ƶ�yMin���ڵ�Ͱ���
			}		
		}
		if(P[j].y<P[i].y)//�ߵ��յ������
		{
			pEdge=new CAET;
			pEdge->x=P[j].x;
			pEdge->yMax=Round(P[i].y);
			pEdge->k=(P[i].x-P[j].x)/(P[i].y-P[j].y);
			pEdge->ps=P[i];
			pEdge->pe=P[j];
			pEdge->pNext=NULL;
			while(pCurrentB->ScanLine!=P[j].y)
			{
				pCurrentB=pCurrentB->pNext;
			}
		}
		if(int(P[j].y)!=P[i].y)
		{
			pCurrentE=pCurrentB->pET;
			if(pCurrentE==NULL)
			{
				pCurrentE=pEdge;
				pCurrentB->pET=pCurrentE;
			}
			else
			{
				while(pCurrentE->pNext!=NULL)
				{
					pCurrentE=pCurrentE->pNext;
				}
				pCurrentE->pNext=pEdge;
			}
		}
	}
}

void CZBuffer::Gouraud(CDC *pDC)//�������
{
	double	CurDeep=0.0;//��ǰɨ���ߵ����
	double	DeepStep=0.0;//��ǰɨ��������x��������Ȳ���
	double	A,B,C,D;//ƽ�淽��Ax+By+Cz��D=0��ϵ��
	CVector V01(P[0],P[1]),V02(P[0],P[2]);
	CVector VN=CrossProduct(V01,V02);
	A=VN.x;B=VN.y;C=VN.z;
	D=-A*P[0].x-B*P[0].y-C*P[0].z;
	DeepStep=-A/C;//����ɨ������Ȳ�������
	CAET *pT1,*pT2;
	pT2 = NULL;
	pHeadE=NULL;	
	for(pCurrentB=pHeadB;pCurrentB!=NULL;pCurrentB=pCurrentB->pNext)
	{
		for(pCurrentE=pCurrentB->pET;pCurrentE!=NULL;pCurrentE=pCurrentE->pNext)
		{
			pEdge=new CAET;
			pEdge->x=pCurrentE->x;
			pEdge->yMax=pCurrentE->yMax;
			pEdge->k=pCurrentE->k;
			pEdge->ps=pCurrentE->ps;
			pEdge->pe=pCurrentE->pe;
			pEdge->pNext=NULL;
			AddEt(pEdge);
		}
		ETOrder();	
		pT1=pHeadE;
		if(pT1==NULL)
			return;
		while(pCurrentB->ScanLine>=pT1->yMax)//�±��Ͽ�
		{
			CAET * pAETTEmp=pT1;			
			pT1=pT1->pNext;
			delete pAETTEmp;
			pHeadE=pT1;
			if(pHeadE==NULL)
				return;
		}
		if(pT1->pNext!=NULL)
		{
			pT2=pT1;
			pT1=pT2->pNext;
		}
		while(pT1!=NULL)
		{
			if(pCurrentB->ScanLine>=pT1->yMax)//�±��Ͽ�
			{
				CAET* pAETTemp =pT1;
				pT2->pNext=pT1->pNext;				
				pT1=pT2->pNext;
				delete pAETTemp;
			}
			else
			{
				pT2=pT1;
				pT1=pT2->pNext;
			}
		}
		CRGB ca,cb,cf;//ca��cb��������������ɫ��cf����������������ɫ
		ca=Interpolation(pCurrentB->ScanLine,pHeadE->ps.y,pHeadE->pe.y,pHeadE->ps.c,pHeadE->pe.c);
		cb=Interpolation(pCurrentB->ScanLine,pHeadE->pNext->ps.y,pHeadE->pNext->pe.y,pHeadE->pNext->ps.c,pHeadE->pNext->pe.c);
		BOOL bInFlag=FALSE;//����������Ա�־����ʼֵΪ�ٱ�ʾ�����ⲿ
		double xleft,xright;//ɨ���ߺ���Ч���ཻ����������յ�����
		for(pT1=pHeadE;pT1!=NULL;pT1=pT1->pNext)
		{
			if(FALSE==bInFlag)
			{
				xleft=pT1->x;
				CurDeep=-(xleft*A+pCurrentB->ScanLine*B+D)/C;//z=-(Ax+By-D)/C
				bInFlag=TRUE;
			}
			else
			{
				xright=pT1->x;
				for(double x=xleft;x<xright;x++)//����ҿ�
				{
					cf=Interpolation(x,xleft,xright,ca,cb);
					if(CurDeep<=zBuffer[Round(x)+Width/2][pCurrentB->ScanLine+Height/2])//�����ǰ����������С��֡��������ԭ����������
					{
						zBuffer[Round(x)+Width/2][pCurrentB->ScanLine+Height/2]=CurDeep;//ʹ�õ�ǰ���������ȸ�����Ȼ�����
						pDC->SetPixelV(Round(x),pCurrentB->ScanLine,RGB(cf.red,cf.green,cf.blue));//���Ƶ�ǰ������
					}
					CurDeep+=DeepStep;
				}
				bInFlag=FALSE;
			}
		}
		for(pT1=pHeadE;pT1!=NULL;pT1=pT1->pNext)//�ߵ�������
			pT1->x=pT1->x+pT1->k;
	}
}

void CZBuffer::AddEt(CAET *pNewEdge)//�ϲ�ET��
{
	CAET *pCE;
	pCE=pHeadE;
	if(pCE==NULL)
	{
		pHeadE=pNewEdge;
		pCE=pHeadE;
	}
	else
	{
		while(pCE->pNext!=NULL)
		{
			pCE=pCE->pNext;
		}
		pCE->pNext=pNewEdge;
	}
}

void CZBuffer::ETOrder()//�߱���ð�������㷨
{
	CAET *pT1,*pT2;
	int Count=1;
	pT1=pHeadE;
	if(pT1==NULL)
		return;
	if(pT1->pNext==NULL)//�����ET��û������ET��
		return;//Ͱ���ֻ��һ���ߣ�����Ҫ����
	while(pT1->pNext!=NULL)//ͳ�Ʊ߽��ĸ���
	{
		Count++;
		pT1=pT1->pNext;
	}
	for(int i=0;i<Count-1;i++)//ð������
	{
		CAET **pPre=&pHeadE;
		pT1=pHeadE;
		for(int j=0;j<Count-1-i;j++)
		{
			pT2=pT1->pNext;
		
			if ((pT1->x>pT2->x)||((pT1->x==pT2->x)&&(pT1->k>pT2->k)))
			{
				pT1->pNext=pT2->pNext;
				pT2->pNext=pT1;
				*pPre=pT2;
			pPre=&(pT2->pNext);//����λ��Ϊ�´α���׼��
			}
			else
			{
				pPre=&(pT1->pNext);
				pT1=pT1->pNext;
			}
		}
	}
}

CRGB CZBuffer::Interpolation(double t,double t1,double t2,CRGB clr1,CRGB clr2)//��ɫ���Բ�ֵ
{
	CRGB color;
	color=(t-t2)/(t1-t2)*clr1+(t-t1)/(t2-t1)*clr2;
	return color;
}

void CZBuffer::InitDeepBuffer(int Width,int Height,double Depth)//��ʼ����Ȼ���
{
	this->Width=Width,this->Height=Height;
	zBuffer=new double *[Width];
	for(int i=0;i<Width;i++)
		zBuffer[i]=new double[Height];
	for(int i=0;i<Width;i++)//��ʼ����Ȼ���
		for(int j=0;j<Height;j++)
			zBuffer[i][j]=Depth;
}

void CZBuffer::ClearMemory()
{
	DeleteAETChain(pHeadE);
	CBucket *pBucket=pHeadB;
	while (pBucket !=NULL)//���ÿһ��Ͱ
	{
		CBucket *pBucketTemp=pBucket->pNext;
		DeleteAETChain(pBucket->pET);
		delete pBucket;
		pBucket=pBucketTemp;
	}
	pHeadB=NULL;
	pHeadE=NULL;
}

void CZBuffer::DeleteAETChain(CAET* pAET)
{
	while (pAET!=NULL)
	{
		CAET *pAETTemp=pAET->pNext;
		delete pAET;
		pAET=pAETTemp;
	}
}
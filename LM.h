
#pragma once
#include <vector>
#include <string>
#include "nr.h"

using namespace std;

static float tmp;
#define CHANGE(x,y) {tmp=(x);(x)=(y);(y)=tmp;}

struct hist{
	short value;
	int freq;
};



//���f���֐��̏����l����p
struct funcdata{
	int height;
	int width;
	int position;
};

class LM
{
public:
	void curveFitting(string patientName);//�t�B�b�e�B���O����֐�
	short getThreshold();//Kang��̎�@�ɂ�邵�����l(LT)��Ԗ�

	float testVal;
	float *funcParam; //�֐��̃p�����[�^
	int funcNum; //�K�E�X�֐��̌�

	//�R���X�g���N�^
	LM(void);
	~LM(void);

private:
	//�Ȑ��t�B�b�e�B���O�ŗ��p
	vector<hist> dataRead(string filename);
	vector<funcdata> hillclimbing(vector<hist> histvec);
	float calcGauss(float x, float a[], int fnum);

	//Numerical Recipe�̂���
	void mrqmin(float x[], float y[], float sig[], int ndata, float a[], int ia[], int ma, float **covar, float **alpha, float *chisq, void (*funcs)(float, float [], float *, float [], int), float *alamda);
	void covsrt(float **covar, int ma, int lista[], int mfit);
	void gaussj(float **a, int n , float **b, int m);
	void mrqcof(float x[], float y[], float sig[], int ndata, float a[], int ia[], int ma, float **alpha, float beta[], float *chisq, void (*funcs)(float, float [], float *, float [], int));
	void covstr(float **covar, int ma, int lista[], int mfit);

};
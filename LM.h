
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



//モデル関数の初期値推定用
struct funcdata{
	int height;
	int width;
	int position;
};

class LM
{
public:
	void curveFitting(string patientName);//フィッティング制御関数
	short getThreshold();//Kangらの手法によるしきい値(LT)を返戻

	float testVal;
	float *funcParam; //関数のパラメータ
	int funcNum; //ガウス関数の個数

	//コンストラクタ
	LM(void);
	~LM(void);

private:
	//曲線フィッティングで利用
	vector<hist> dataRead(string filename);
	vector<funcdata> hillclimbing(vector<hist> histvec);
	float calcGauss(float x, float a[], int fnum);

	//Numerical Recipeのもの
	void mrqmin(float x[], float y[], float sig[], int ndata, float a[], int ia[], int ma, float **covar, float **alpha, float *chisq, void (*funcs)(float, float [], float *, float [], int), float *alamda);
	void covsrt(float **covar, int ma, int lista[], int mfit);
	void gaussj(float **a, int n , float **b, int m);
	void mrqcof(float x[], float y[], float sig[], int ndata, float a[], int ia[], int ma, float **alpha, float beta[], float *chisq, void (*funcs)(float, float [], float *, float [], int));
	void covstr(float **covar, int ma, int lista[], int mfit);

};
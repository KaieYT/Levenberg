#include "LM.h"
#include <float.h>
#include <math.h>
#include <tgmath.h>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

//new description

//Levenberg-Marquardtñ@ Numerical Recipes in C pp.507~
/*
x[], y[]:ÉfÅ[É^ì_Å@sig[] : y[]ÇÃïWèÄïŒç∑ ÇÉpÉâÉÅÅ[É^a[]Çä‹ÇﬁîÒê¸å`ä÷êîÇ≈ìñÇƒÇÕÇﬂÇÈÅD
îzóÒlista[]Ç…ì¸ÇÍÇΩÉpÉâÉÅÅ[É^î‘çÜÇÃì‡ÅCç≈èâÇÃmfitå¬Ç™é¿ç€ÇÃìñÇƒÇÕÇﬂÇ…ÇÊÇ¡ÇƒåàíËÇ≥ÇÍÇÈÉpÉâÉÅÅ[É^Ç…ëäìñÇ∑ÇÈÅD
écÇËÇÃma-mfitå¬ÇÕÅCì¸óÕéûÇÃílÇ…å≈íËÇ≥ÇÍÇÈÅDÉãÅ[É`ÉìÇ™ï‘Ç∑ÇÃÇÕÅCçXêVÇ≥ÇÍÇΩmaå¬ÇÃÉpÉâÉÅÅ[É^aÇ∆x^2=chisqÇÃílÇ≈Ç†ÇÈÅD
îzóÒcovar[1~ma][1~ma]Ç∆alpha[1~ma][1~ma]ÇÃÇ§Çø[1~mfit][1~mfit]ÇÃïîï™ÇÕÅCÇŸÇ∆ÇÒÇ«ÇÃåJÇËï‘ÇµéûÇ…çÏã∆óÃàÊÇ∆ÇµÇƒóòópÇ≥ÇÍÇÈÅD
óòópé“íËã`ä÷êîfuncsÇÕìñÇƒÇÕÇﬂÇÈä÷êîyfitÅCÇ®ÇÊÇ—ÇªÇÃÉpÉâÉÅÅ[É^aÇ…Ç¬Ç¢ÇƒÇÃì±ä÷êîdyda[1..ma]ÇÃxÇ≈ÇÃílÇåvéZÇ∑ÇÈÇ‡ÇÃÇ∆Ç∑ÇÈÅD
1âÒñ⁄ÇÃåƒÇ—èoÇµÇ≈ÇÕÅCÉpÉâÉÅÅ[É^aÇÃèâä˙êÑë™ílÇó^Ç¶ÅCèâä˙âªÇÃÇΩÇﬂÇÃalamda < 0Ç∆ÇµÇƒåƒÇ—ñﬂÇ∑ÅDÇ±ÇÃÉãÅ[É`ÉìÇ™é˚ë©ÇµÇΩÇ∆Ç´ÅC
covarÇ…ã§ï™éUçsóÒÅCalphaÇ…ã»ó¶çsóÒÇ™ì¸ÇÈÅD
*/
/*
void mrqmin(float x[], float y[], float sig[], int ndata, float a[], int ma, int lista[], int mfit, float **covar, float **alpha, float *chisq, void (*funcs)(float, float [], float *, float [], int), float *alamda){
	void covsrt(float **covar, int ma, int lista[], int mfit);
	void gaussj(float **a, int n , float **b, int m);
	void mrqcof(float x[], float y[], float sig[], int ndata, float a[], int ma, int lista[], int mfit, float **alpha, float beda[], float *chisq, void (*funcs)(float, float *, float *, float *, int));
	int k, kk, j, ihit;
	static float *da, *atry, **oneda, *beta, ochisq;

	if(*alamda < 0.0){ //èâä˙âª
		oneda = matrix(1, mfit, 1, 1); //Ç±ÇÍÇÁÇÃïœêîÇÕç≈å„ÇÃåƒÇ—èoÇµ(alambda
		atry = Yvector(1, ma);
		da = Yvector(1, ma);
		beta = Yvector(1, ma);
		kk = mfit + 1;

		for(int j = 1 ; j <= ma ; j++){ //listaÇÕê≥ÇµÇ¢åWêîÇÃèáóÒÇ©ÅH
			ihit = 0;
			for(k = 1 ; k <= mfit ; k++){
				if(lista[k] == j)
					ihit++;
			if(ihit == 0)
				lista[kk++] = j;
			else if(ihit > 1)
				nrerror("Bad lista permutation in mrqmin-1");
		}
			if(kk != ma+1)
				nrerror("Bad lista permutation in meqmin-2");
			*alamda = 0.001;
			mrqcof(x, y, sig, ndata, a, ma, lista, mfit, alpha, beta, chisq, funcs);
			ochisq = (*chisq);
		}
		for(int j = 1; j <=mfit ; j++){ //ê¸å`âªìñÇƒÇÕÇﬂçsóÒÇÃëŒäpóvëfÇëùÇ‚Ç∑
			for(int k = 1 ; k <= mfit ; k++)
				covar[j][k] = alpha[j][k];
			covar[j][j] = alpha[j][j] * (1.0 + (*alamda));
			oneda[j][1] = beta[j];
		}
		gaussj(covar, mfit, oneda,1); //çsóÒï˚íˆéÆÇâÇ≠
		for(int j = 1 ; j <= mfit ; j++)
			da[j] = oneda[j][1];
		if(*alamda == 0.0){ //é˚ë©ÇµÇΩÇÁalamda = 0Ç≈ã§ï™éUçsóÒÇãÅÇﬂÇÈ
			covsrt(covar, ma, lista, mfit);
			free_vector(beta, 1, ma);
			free_vector(da, 1, ma);
			free_vector(atry ,1, ma);
			free_matrix(oneda, 1, mfit, 1, 1);
			return;
		}
		for(int j = 1 ; j <= ma ; j++)
			atry[j] = a[j];
		for(int j = 1 ; j <= mfit ; j++) //ééçsÇÕê¨å˜ÇµÇΩÇ©ÅH
			atry[lista[j]] += da[j];
		mrqcof(x, y, sig, ndata, a, ma, lista, mfit, alpha, beta, chisq, funcs);
		if(*chisq < ochisq){ //ê¨å˜Å®êVÇµÇ¢âÇçÃóp
			*alamda *= 0.1;
			ochisq = ( *chisq);
			for(int j = 1; j <= mfit ; j++){
				for(int k = 1 ; k <= mfit ; k++)
					alpha[j][k] = covar[j][k];
				beta[j] = da[j];
				a[lista[j]] = atry[lista[j]];
			}
		}else{
			*alamda += 10.0; //é∏îsÅ®alamdaÇëùâ¡ÇµÇƒñﬂÇÈ
			*chisq = ochisq;
		}
		return;
	}
}
*/

/*
void mrqcof(float x[], float y[], float sig[], int ndata, float a[], int ma, int lista[], int mfit, float **alpha, float beta[], float *chisq, void (*funcs)(float, float *, float *, float *, int)){
	int k, j , i;
	float ymod, wt, sig2i, dy, *dyda;

	dyda = Yvector(1, ma);
	for(int j = 1 ; j <= mfit ; j++){
		for(k = 1 ; k <= j ; k++)
			alpha[j][k] = 0;
		beta[j] = 0.0;
	}

	*chisq = 0.0;
	for(int i = 1 ; i <= ndata ; i++){
		(*funcs)(x[i], a, &ymod, dyda, ma);
		sig2i = 1.0 / (sig[i] * sig[i]);
		dy = y[i] - ymod;
		for(int j = 1 ; j <= mfit ; j++){
			wt = dyda[lista[j]] * sig2i;
			for(int k = 1 ; k <= j ; k++)
				alpha[k][j] += wt * dyda[lista[k]];
			beta[j] += dy * wt;
	}
		(*chisq)  = dy * dy * sig2i;
	}
	for(int j = 2 ; j <= mfit ; j++)
		for(int k = 1 ; k <= j-1 ; k++)
			alpha[k][j] = alpha[j][k];
	free_vector(dyda, 1, ma);
}
*/

//Numerical Recipe in C 2nd editionî≈ mrqmin
void LM::mrqmin(float x[], float y[], float sig[], int ndata, float a[], int ia[], int ma, float **covar, 
float **alpha, float *chisq, void (*funcs)(float, float [], float *, float [], int), float *alamda){
	/*
	void covsrt(float **covar, int ma, int ia[], int mfit);
	void gaussj(float **a, int n, float **b, int m);
	void mrqcof(float x[], float y[], float sig[], int ndata, float a[], int ia[], int ma, float **alpha, float beta[], float *chisq, void (*funcs)(float, float [], float *, float [], int));
	*/
	int j,k,l;
	static int mfit;
	static float ochisq,*atry,*beta,*da,**oneda;
	if (*alamda < 0.0) { //Initialization.
		atry=Yvector(1,ma);
		beta=Yvector(1,ma);
		da=Yvector(1,ma);
		for (mfit=0,j=1;j<=ma;j++)
			if (ia[j]) mfit++;
		oneda=matrix(1,mfit,1,1);
		*alamda=0.001;
		mrqcof(x,y,sig,ndata,a,ia,ma,alpha,beta,chisq,funcs);
		ochisq=(*chisq);
		for (j=1;j<=ma;j++) atry[j]=a[j];
	}
	for (j=1;j<=mfit;j++) { //Alter linearized fitting matrix, by augmenting diagonal elements.
			for (k=1;k<=mfit;k++) covar[j][k]=alpha[j][k];
		covar[j][j]=alpha[j][j]*(1.0+(*alamda));
		oneda[j][1]=beta[j];
	}
	gaussj(covar,mfit,oneda,1); //Matrix  solution.
		for (j=1;j<=mfit;j++) da[j]=oneda[j][1];
	if (*alamda == 0.0) { //Once converged, evaluate covariance matrix.
		covsrt(covar,ma,ia,mfit);
		covsrt(alpha,ma,ia,mfit); //Spread  out alpha to  its  full  size  too.
		free_matrix(oneda,1,mfit,1,1);
		free_vector(da,1,ma);
		free_vector(beta,1,ma);
		free_vector(atry,1,ma);
		return;
	}
	for (j=0,l=1;l<=ma;l++) //Did  the  trial  succeed?
		if (ia[l]) atry[l]=a[l]+da[++j];
	mrqcof(x,y,sig,ndata,atry,ia,ma,covar,da,chisq,funcs);
	if (*chisq < ochisq) { //Success, accept the  new solution.
			*alamda *= 0.1;
		ochisq=(*chisq);
		for (j=1;j<=mfit;j++) {
			for (k=1;k<=mfit;k++) alpha[j][k]=covar[j][k];
			beta[j]=da[j];
		}
		for (l=1;l<=ma;l++) a[l]=atry[l];
	} else { //Failure,  increase alamda and return.
			*alamda *= 10.0;
		*chisq=ochisq;
	}
}

//2nd editionî≈ mrqcof
void LM::mrqcof(float x[], float y[], float sig[], int ndata, float a[], int ia[], int ma, float **alpha, float beta[], float *chisq, void (*funcs)(float, float [], float *, float [], int)) //Used by mrqmin to evaluate the linearized tting matrix alpha, and vector beta as in (15.5.8), and  calculate 2.
{
	int i,j,k,l,m,mfit=0;
	float ymod,wt,sig2i,dy,*dyda;
	dyda=Yvector(1,ma);
	for (j=1;j<=ma;j++)
		if (ia[j]) mfit++;
	for (j=1;j<=mfit;j++) { //Initialize  (symmetric) alpha, beta.
		for (k=1;k<=j;k++) alpha[j][k]=0.0;
		beta[j]=0.0;
	}
	*chisq=0.0;
	for (i=1;i<=ndata;i++) { //Summation loop over all  data.
		(*funcs)(x[i],a,&ymod,dyda,ma);
		sig2i=1.0/(sig[i]*sig[i]);
		dy=y[i]-ymod;
		for (j=0,l=1;l<=ma;l++) {
			if (ia[l]) {
				wt=dyda[l]*sig2i;
				for (j++,k=0,m=1;m<=l;m++)
					if (ia[m]) alpha[j][++k] += wt*dyda[m];
				beta[j] += dy*wt;
			}
		}
	*chisq += dy*dy*sig2i; //And  find chi^2
	}
	for (j=2;j<=mfit;j++) //Fill  in  the  symmetric  side.
	for (k=1;k<j;k++) alpha[k][j]=alpha[j][k];
	free_vector(dyda,1,ma);
}

//na/3å¬ÇÃï°êîÇÃÉKÉEÉXä÷êîÇÃòaÇåvéZ
//a[1..na]éüå≥ÅCdyda[1..na]éüå≥
//1ä÷êîÇ…Ç¬Ç´a[i] : çÇÇ≥ a[i+1] : íÜêSà íu a[i+2] : ïù ÇïKóvÇ∆Ç∑ÇÈ
//
void fgauss(float x, float a[], float *y, float dyda[], int na){
	int i;
	float fac, ex, arg;

	*y = 0.0;
	
	for(int i = 1 ; i <= na-1 ; i += 3){
		arg = ( x-a[i+1]) / a[i+2];
		ex = exp(-arg * arg);
		fac = a[i] * ex * 2.0 * arg;
		*y += a[i] * ex;
		dyda[i] = ex;
		dyda[i+1] = fac / a[i+2];
		dyda[i+2] = fac * arg / a[i+2];
	}
}

/*
void covstr(float **covar, int ma, int lista[], int mfit){
	float swap;

	for(int j=1 ; j <ma ; j++)
		for(int i = j+1 ; j <= mfit ; j++)
			covar[i][j] = 0.0;
	for(int i=1 ; i < mfit ; i++)
		for(int j=i+1 ; j<=mfit ; j++){
			if(lista[j] > lista[i])
				covar[lista[j]][lista[i]] = covar[i][j];
			else
				covar[lista[i]][lista[j]] = covar[j][i];
		}
		swap = covar[1][1];

		for(int j=1 ; j<=ma ; j++){
			covar[1][j] = covar[j][j];
			covar[j][j] = 0.0;
		}
		covar[lista[1]][lista[1]] = swap;
		for(int j=2 ; j<=mfit ; j++)
			covar[lista[j]][lista[j]] = covar[1][j];
		for(int j=2 ; j<=ma ; j++)
			for(int i=1 ; i <=j-1 ; i++)
				covar[i][j] = covar[j][i];
}
*/
//2nd editionî≈
void LM::covsrt(float **covar, int ma, int ia[], int mfit)//Expand in storage the covariance matrix covar, so as to take into account parameters that are being  held fixed.  (For  the  latter,  return  zero  covariances.)
{
	int i,j,k;
	float swap;
	for (i=mfit+1;i<=ma;i++)
		for (j=1;j<=i;j++) covar[i][j]=covar[j][i]=0.0;
	k=mfit;
	for (j=ma;j>=1;j--) {
		if (ia[j]) {
			for (i=1;i<=ma;i++) CHANGE(covar[i][k],covar[i][j])
			for (i=1;i<=ma;i++) CHANGE(covar[k][i],covar[j][i])
			k--;
		}
	}
}

void LM::gaussj(float **a, int n, float **b, int m){
	int *indxc,*indxr,*ipiv;
	int i,icol,irow,j,k,l,ll;
	float big,dum,pivinv,temp;
	indxc=ivector(1,n);
	indxr=ivector(1,n);
	ipiv=ivector(1,n);
	for (j=1;j<=n;j++) ipiv[j]=0;
	for (i=1;i<=n;i++) {
		big=0.0;
		for (j=1;j<=n;j++)
			if (ipiv[j] != 1)
				for (k=1;k<=n;k++) {
					if (ipiv[k] == 0) {
						if (fabs(a[j][k]) >= big) {
							big=fabs(a[j][k]);
							irow=j;
							icol=k;
						}
					} else if (ipiv[k] > 1) nrerror("gaussj: Singular Matrix-1");
				}
				++(ipiv[icol]);
				if (irow != icol) {
					for (l = 1;l <= n;l++) 
						CHANGE(a[irow][l], a[icol][l])
					for (l = 1;l <= m;l++) 
						CHANGE(b[irow][l], b[icol][l])
				}
				indxr[i]=irow;

				indxc[i]=icol;
				if (a[icol][icol] == 0.0) 
					nrerror("gaussj: Singular Matrix-2");
				pivinv=1.0 / a[icol][icol];
				a[icol][icol]=1.0;
				for (l=1;l<=n;l++) a[icol][l] *= pivinv;
				for (l=1;l<=m;l++) b[icol][l] *= pivinv;
				for (ll=1;ll<=n;ll++)
					if (ll != icol) {
						dum=a[ll][icol];
						a[ll][icol]=0.0;
						for (l=1;l<=n;l++) a[ll][l] -= a[icol][l]*dum;
						for (l=1;l<=m;l++) b[ll][l] -= b[icol][l]*dum;
					}
	}

	for (l=n;l>=1;l--) {
		if (indxr[l] != indxc[l])
			for (k=1;k<=n;k++)
				CHANGE(a[k][indxr[l]],a[k][indxc[l]]);
	}
	free_ivector(ipiv,1,n);
	free_ivector(indxr,1,n);
	free_ivector(indxc,1,n);
}

//CSVÉfÅ[É^ì«Ç›çûÇ›
vector<hist> LM::dataRead(string patientName){
	string str;

	//filename =  "C:\\Users\\YASUDA\\Desktop\\hist_NCCHE_3DABD_1_3_sample.csv";
	/*
	string filename;
	filename = SOURCEDIR;
	filename += "\\";
	filename += "hist_";
	filename += patientName;
	filename += ".csv";
     */
    string filename = "//Users//t_yasuda//Desktop//hist_NCCHE_3DABD_1_1.csv";

	stringstream ss;
	int hoge, huga;

	vector<hist> histvec;

	ifstream ifs(filename);
	if(!ifs){
        std::cout << filename << endl;
        std::cout << "Error : Input data file is not found" << endl;
        exit(1);
	}

	while(getline(ifs, str)){
		string token;
		stringstream ss, ss2;

		istringstream stream(str);

		while(getline(stream, token, ',')){
			ss << token;
			ss >> hoge;
		}
		ss2 << token;
		ss2 >> huga;
		hist h;
		h.value = hoge;
		h.freq= huga;
		histvec.push_back(h);
		//printf("num : %d, data : %d\n", hoge, huga);
	}
	

	return histvec;
}

//ã…íl(ã…ëÂíl)íTçı
vector<funcdata> LM::hillclimbing(vector<hist> histvec){
	
	vector<funcdata> fd; //èâä˙íläiî[ópîzóÒ
	
	for(int i = 1 ; i < histvec.size()-1 ; i++){
		int eval = histvec.at(i).freq;
		bool flag = false; //ã…ílåÛï‚îªíËópÉtÉâÉO
		
		if(histvec.at(i-1).freq < eval && histvec.at(i+1).freq < eval){
			flag = true;
			for(int j = -10 ; j <= 10 ; j++)
				if(0 < j+i && j+i < histvec.size()-1)
					if(eval < histvec.at(i+j).freq || eval < 5000) //1000ÇÕåàÇﬂë≈Çø
						flag = false;
		}

		if(flag){
			funcdata fdata;
			fdata.position = i; //ñ{óàCTílÇë„ì¸
			fdata.height = histvec.at(i).freq;
			fdata.width = 10; //ñ{óàFWHM(full width half max)Ç©ÇÁåàíË
			fd.push_back(fdata);
			printf("exValue Position : %d, height : %d, width : %d\n", fdata.position, fdata.height, fdata.width);
		}
	}
	printf("Num. of exValue : %d\n", fd.size());
	return fd;
}

//ã»ê¸ÉtÉBÉbÉeÉBÉìÉO
void LM::curveFitting(string patientName){
	vector<hist> histvec = dataRead(patientName);//ÉqÉXÉgÉOÉâÉÄì¸óÕ
	vector<funcdata> funcVec = hillclimbing(histvec);//ä÷êîÇÃèâä˙ílÇéÊìæ
	int fnum = funcVec.size(); //ìñÇƒÇÕÇﬂÇΩÇ¢ÉKÉEÉXä÷êîÇÃå¬êî

	vector<hist>::iterator histIter = histvec.begin();
	int histnum = histvec.size();

	//Numerical RecipeÇÃå`éÆÇ÷ïœä∑
	float *value = Yvector(1, histnum);
	float *freq = Yvector(1, histnum);
	
	float min = FLT_MAX, max = FLT_MIN;

	for(int i = 1 ; i <= histnum ; i++){
		hist h = histvec.at(i-1);
		//value[i] = h.value; //ÉqÉXÉgÉOÉâÉÄíÜÇÃâÊëfíl
		value[i] = i; //íPèÉÇ»ç¿ïWíl
		freq[i] = h.freq;
		max = fmax(h.freq, max);
		min = fmin(h.freq, min);
	}

	//åvéZópïœêîíËã`
	float *a = Yvector(1, fnum*3);
	int *ia = ivector(1, fnum*3), iter, itst;
	float *sigma = Yvector(1, histnum);
	float **covar = matrix(1, fnum*3, 1, fnum*3), **alpha = matrix(1, fnum*3, 1, fnum*3);
	float chisq, ochisq;

	float alamda = -1;

	//ä÷êîÉpÉâÉÅÅ[É^Ç…èâä˙ílë„ì¸
	for(int i  = 1 ; i <= fnum * 3 ; i += 3){
		a[i] = funcVec.at((int)(i/3+1)-1).height; //çÇÇ≥
		a[i+1] = funcVec.at((int)(i/3+1)-1).position; //à íu
		a[i+2] = funcVec.at((int)(i/3+1)-1).width; //ïù
	}

	for(int i = 1 ; i <= histnum ; i++)
		sigma[i] = 1; //éQçlï∂å£Ç∆ìØÇ∂ÅDèdÇ›ïtÇØÇ»Çµ
			
	for(int i = 1 ; i <= fnum*3 ; i += 3){
		ia[i] = 0; //ëSÇƒÇÃÉpÉâÉÅÅ[É^a[1..ma]ÇïœçXÇ∑ÇÈéwé¶(=0ÇÃèÍçáì¸óÕÇå≈íË)
		ia[i+1] = 0;
		ia[i+2] = 1;
	}
	mrqmin(value, freq, sigma, histnum, a, ia, fnum * 3, covar, alpha,  &chisq, fgauss, &alamda);

	itst = 0;

	for(int iter = 1 ; iter <= 1000 ; iter++){
		ochisq = chisq;
		mrqmin(value, freq, sigma, histnum, a, ia, fnum * 3, covar, alpha,  &chisq, fgauss, &alamda);
		//printf("chisq : %1.0f\n", chisq);
				
		//if(chisq > ochisq)
		if(fabs(ochisq - chisq) > 0.1)
			itst = 0;
		else if(fabs(ochisq - chisq) < 0.1)
			itst++;
		if(itst < 100)
			continue;

		//ÉtÉBÉbÉgäÆóπ?
		alamda = 0.0;
		mrqmin(value, freq, sigma, histnum, a, ia, fnum * 3, covar, alpha,  &chisq, fgauss, &alamda);
		break;
	}
	
	//CSVÇ÷ÇÃèoóÕ
	string filename, str;
    /*
	filename = "C:\\Users\\YASUDA\\Desktop\\";
	filename += "\\";
	filename += "hist_";
	filename += patientName;
	filename += "_LM.csv";
*/
    filename = "¥¥Users¥¥t_yasuda¥¥Desktop¥¥hist_NCCHE_3DABD_1_3_result.csv";
	ofstream ofs(filename);

	for(int i = 1 ; i <= fnum*3 ; i++)
		//printf("covar : %1.0f ", covar[i][i]);
	//printf("\n");

	//ç¿ïWílÇ…ëŒÇ∑ÇÈÉKÉEÉXä÷êîÇÃílÇåvéZ
	for(int i = 0 ; i < histnum ; i++){
		for(int j = 1 ; j <= fnum*3 ; j+=3){
			ofs << calcGauss(i, a, j) << ",";
		}
		ofs << endl;
		//ofs << calcGauss(i, a, fnum) << endl;
	}
	ofs.close();
	
	//ÉÇÉfÉãÉtÉBÉbÉeÉBÉìÉOëOå„ÇèoóÕ
	for(int i  = 1 ; i <= fnum * 3 ; i += 3){
		printf("func no. : %d ", i);
		printf("org.height : %d ", funcVec.at((int)(i/3+1)-1).height);
		printf("a.height : %1.0f \n",  a[i]);

		printf("org.position : %d ", funcVec.at((int)(i/3+1)-1).position);
		printf("a.position : %1.0f \n",  a[i+1]);

		printf("org.width : %d ", funcVec.at((int)(i/3+1)-1).width);
		printf("a.width : %1.0f \n",  a[i+2]);
		printf("\n");
	}
	
	//ä÷êîÉpÉâÉÅÅ[É^ÅCÉKÉEÉXä÷êîÇÃå¬êîÇÉNÉâÉXïœêîÇ÷äiî[
	funcParam = a;
	funcNum = fnum;
}

//ç≈ëÂCTílÇÃï™ïz(nå¬ñ⁄)Ç∆1Ç¬è¨Ç≥Ç¢ï™ïzÇ…ëŒâûÇ∑ÇÈÉKÉEÉXä÷êîÇÃåì_Çì±èo
//KangÇÁÇÃéËñ@Ç…Ç®ÇØÇÈLTÇï‘ñﬂ
short LM::getThreshold(){
	/*
	//ä÷êîÉpÉâÉÅÅ[É^Ç…èâä˙ílë„ì¸
	for(int i  = 1 ; i <= fnum * 3 ; i += 3){
		a[i] = funcVec.at((int)(i/3+1)-1).height; //çÇÇ≥
		a[i+1] = funcVec.at((int)(i/3+1)-1).position; //à íu
		a[i+2] = funcVec.at((int)(i/3+1)-1).width; //ïù
	}
	*/
	
	float g1, g2, e1, e2, b1, b2;
	b1 = funcParam[funcNum*3-5];
	e1 = funcParam[funcNum*3-4];
	g1 = funcParam[funcNum*3-3];
	
	b2 = funcParam[funcNum*3-2];
	e2 = funcParam[funcNum*3-1];
	g2 = funcParam[funcNum*3];

	float g1p2 = pow(g1, 2), g2p2 = pow(g2, 2);
	float e1p2 = pow(e1, 2), e2p2 = pow(e2, 2);
	float a = g1p2 - g2p2;
	float b = 2*(g2p2*e1 - g1p2*e2);
	float c = -g2p2 * e1p2 + g1p2 * e2p2 - g1p2 * g2p2 * log(b2/b1);
	float x1 = (-b + sqrt(pow(b, 2) - 4*a*c)) / (2*a);
	float x2 = (-b - sqrt(pow(b, 2) - 4*a*c)) / (2*a);

	if(e1 < x1 && x1 < e2)
		return x1;
	else
		return x2;

}


float LM::calcGauss(float x, float a[], int fnum){
	float fac, ex, arg, y = 0.0;
	/*
	for(int i = 1 ; i < fnum * 3 ; i += 3){
		arg = ( x-a[i+1]) / a[i+2];
		ex = exp(-arg * arg);
		fac = a[i] * ex * 2.0 * arg;
		y += a[i] * ex;
	}*/
	
	//fnumÇ…ëŒâûÇ∑ÇÈílÇåvéZ
	y = a[fnum] * exp(-pow((x-a[fnum+1])/a[fnum+2],2));
/*
	for(int i = fnum ; i <= fnum ; i++){
		arg = ( x-a[i+1]) / a[i+2];
		ex = exp(-arg * arg);
		fac = a[i] * ex * 2.0 * arg;
		y += a[i] * ex;
	}
	*/
	return y;
}

/*
#define NPT 100
#define MA 6
#define SPREAD 0.001
//êßå‰óp

void curveFitting(){
	long idum=(-911);
	int i,*ia,iter,itst,j,k,mfit=MA;
	float alamda,chisq,ochisq,*x,*y,*sig,**covar,**alpha;
	static float a[MA+1]=
		{0.0,5.0,2.0,3.0,2.0,5.0,3.0};
	static float gues[MA+1]=
		{0.0,4.5,2.2,2.8,2.5,4.9,2.8};

	ia=ivector(1,MA);
	x=Yvector(1,NPT);
	y=Yvector(1,NPT);
	sig=Yvector(1,NPT);
	covar=matrix(1,MA,1,MA);
	alpha=matrix(1,MA,1,MA);
	// First try a sum of two Gaussians
	for (i=1;i<=NPT;i++) {
		x[i]=0.1*i;
		y[i]=0.0;
		for (j=1;j<=MA;j+=3) {
			y[i] += a[j]*exp(-SQR((x[i]-a[j+1])/a[j+2]));
		}
		y[i] *= (1.0+SPREAD*gasdev(&idum));
		sig[i]=SPREAD*y[i];
	}
	for (i=1;i<=mfit;i++) ia[i]=1;
	for (i=1;i<=MA;i++) a[i]=gues[i];
	for (iter=1;iter<=2;iter++) {
		alamda = -1;
		mrqmin(x, y, sig, 100, a, ia, MA, covar, alpha, &chisq,fgauss,&alamda);
		k=1;
		itst=0;
		for (;;) {
			printf("\n%s %2d %17s %10.4f %10s %9.2e\n","Iteration #",k,
				"chi-squared:",chisq,"alamda:",alamda);
			printf("%8s %8s %8s %8s %8s %8s\n",
				"a[1]","a[2]","a[3]","a[4]","a[5]","a[6]");
			for (i=1;i<=6;i++) printf("%9.4f",a[i]);
			printf("\n");
			k++;
			ochisq=chisq;
			mrqmin(x,y,sig,NPT,a,ia,MA,covar,alpha,&chisq,fgauss,&alamda);
			if (chisq > ochisq)
				itst=0;
			else if (fabs(ochisq-chisq) < 0.1)
				itst++;
			if (itst < 4) continue;
			alamda=0.0;
			mrqmin(x,y,sig,NPT,a,ia,MA,covar,alpha,&chisq,fgauss,&alamda);
			printf("\nUncertainties:\n");
			for (i=1;i<=6;i++) printf("%9.4f",sqrt(covar[i][i]));
			printf("\n");
			printf("\nExpected results:\n");
			printf(" %7.2f %8.2f %8.2f %8.2f %8.2f %8.2f\n",
				5.0,2.0,3.0,2.0,5.0,3.0);
			break;
		}
		if (iter == 1) {
			printf("press return to continue with constraint\n");
			(void) getchar();
			printf("holding a[2] and a[5] constant\n");
			for (j=1;j<=MA;j++) a[j] += 0.1;
			a[2]=2.0;
			ia[2]=0;
			a[5]=5.0;
			ia[5]=0;
		}
	}
	free_matrix(alpha,1,MA,1,MA);
	free_matrix(covar,1,MA,1,MA);
	free_vector(sig,1,NPT);
	free_vector(y,1,NPT);
	free_vector(x,1,NPT);
	free_ivector(ia,1,MA);
	return 0;
}
*/
LM::LM(void)
{
}


LM::~LM(void)
{
}

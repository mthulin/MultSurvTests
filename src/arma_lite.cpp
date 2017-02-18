#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

/* Note: this program attempts to lessen
the number of functions exposed to the 
user. */
/*-------------Y_ik-------------*/
double Yik_cpp_arma(arma::mat x, int t, int k)
{
	double length = 0;
	arma::vec x_k = x.col(k-1);
	arma::vec::const_iterator it_end = x_k.end(); 

	/*Arma iterator type*/
	for(arma::vec::const_iterator it = x_k.begin(); it != it_end; ++it)
	{
		/*it is a pointer of const. type. Need to dereference(*) to read
		element it points to*/
		if(*it >= t) 
			length += 1;
	}

	return length;
}

/*-------------summand.G-------------*/
double summandG_cpp_arma(arma::mat x, arma::mat y, int n, 
				 int t, double t_delta, int k)
{
	double G_sum = 0;
	
	G_sum = (t_delta * Yik_cpp_arma(y, t, k) ) / n;

	return G_sum;
}

/*-------------summand.L-------------*/
double summandL_cpp_arma(arma::mat x, arma::mat y, int n, int t, 
	double t_delta, int k)
{
	double L_sum;
	double Yik_1 = Yik_cpp_arma(x,t,k);
	double Yik_2 = Yik_cpp_arma(y,t,k);

	L_sum = t_delta*Yik_2 / (Yik_1 + Yik_2);

	return L_sum;
}

/*-------------gehan.test-------------*/
double GehanTest_cpp_arma(arma::mat x, arma::mat y, int n1, int n2, 
	arma::mat delta_x, arma::mat delta_y, int k)
{
	arma::vec x_k = x.col(k-1), y_k = y.col(k-1),
		delta_x_k = delta_x.col(k-1), delta_y_k = delta_y.col(k-1);
	int n = n1+n2;

	double T = 0;

	for(int i = 0; i < n1; ++i)
	{
		T += summandG_cpp_arma(x, y, n, x_k(i), delta_x_k(i), k);
	}
	for(int j = 0; j < n2; ++j)
	{
		T -= summandG_cpp_arma(y, x, n, y_k(j), delta_y_k(j), k);
	}

	return T / pow(n, 0.5);
}

/*--------------mvlogrank.test----------------*/
double mvlogrankTest_cpp_arma(arma::mat x, arma::mat y, int n1, int n2,
	arma::mat delta_x, arma::mat delta_y, int k)
{
	arma::vec x_k = x.col(k-1), y_k = y.col(k-1),
		delta_x_k = delta_x.col(k-1), delta_y_k = delta_y.col(k-1);
	
	int n = n1+n2;
	double T = 0;

	/*iterators of const type*/
	arma::vec::const_iterator it_end = x_k.end(), it_end2 = y_k.end(),
		it_end3 = delta_x_k.end(), it_end4 = delta_y_k.end();

	for(arma::vec::const_iterator it = x_k.begin(), it3 = delta_x_k.begin(); 
		it != it_end && it3 != it_end3; ++it, ++it3)
	{
		T += summandL_cpp_arma(x, y, n, *it, *it3, k);
	}

	for(arma::vec::const_iterator it2 = y_k.begin(), it4 = delta_y_k.begin();
		it2 != it_end2 && it4 != it_end4; ++it2, ++it4)
	{
		T -= summandL_cpp_arma(y, x , n, *it2, *it4, k);
	}

	return T / pow(n ,0.5);
}

/*------------------------------------------------*/
/*mu.G*/
double muG_cpp_arma(arma::mat x, arma::mat y, int t, int k, int n = 0)
{
	arma::vec x_k = x.col(k-1), y_k = y.col(k-1);
	int n1 = x_k.size() + y_k.size();
	double T;
	
	if(n == 0) {
		T = Yik_cpp_arma(y, t, k)/n1;
	} else {
		T = Yik_cpp_arma(y, t, k)/n;
	}

	return T;
}

/*------------------------------------------------*/
/*mu.L*/
double muL_cpp_arma(arma::mat x, arma::mat y, int t, int k)
{
	double T;

	T = 1-(Yik_cpp_arma(x,t,k)/(Yik_cpp_arma(x,t,k) + Yik_cpp_arma(y, t, k)));
	return T;
}

/*------------------------------------------------*/
/*psi.G*/
double psiG_cpp_arma(arma::mat x, arma::mat y, int t,
	arma::mat delta_x, int k, int n)
{
	arma::vec x_k = x.col(k-1), delta_x_k = delta_x.col(k-1);
	double S = 0;
	
	arma::vec::const_iterator it_end = x_k.end(), it_end2 = delta_x_k.end();

	for(arma::vec::const_iterator it = x_k.begin(), it2 = delta_x_k.begin(); 
		(it != it_end) && (it2 != it_end2); ++it, ++it2)
	{
		if(*it <= t){
			S += *it2*muG_cpp_arma(x,y,*it,k,n)/Yik_cpp_arma(x, *it, k);
		} else{ 
			S += 0;
		}
	}

	return S;
}

/*------------------------------------------------*/
/*psi.L*/
double psiL_cpp_arma(arma::mat x, arma::mat y, int t, 
	arma::mat delta_x, int k)
{
	arma::vec x_k = x.col(k-1), delta_x_k = delta_x.col(k-1);
	double S = 0;
	arma::vec::const_iterator it_end = x_k.end(), 
		it_end2 = delta_x_k.end();
	
	for(arma::vec::const_iterator it = x_k.begin(), it2 = delta_x_k.begin();
		it != it_end && it2 != it_end2; ++it, ++it2)
	{
		if(*it <= t) {
			S += *it2 * muL_cpp_arma(x,y,*it,k) / Yik_cpp_arma(x,*it,k);
		} else {
			S += 0;
		}
	}
	return S; 
}

/*------------------------------------------------*/
/*sigma.ikl.G*/
/*removed int t, delta.y not used*/
double sigma_iklG_cpp_arma(arma::mat x, arma::mat y,
	arma::mat delta_x, arma::mat delta_y, int k, int l, int n)
	{
		arma::vec x_k = x.col(k-1), delta_x_k = delta_x.col(k-1),
			x_l = x.col(l-1), delta_x_l = delta_x.col(l-1);
		double U = 0;
		int nx = x_k.size();

		arma::vec::const_iterator it_end = x_k.end(), it_end2 = x_l.end(),
			it_end3 = delta_x_k.end(), it_end4 = delta_x_l.end();

		for(arma::vec::const_iterator it = x_k.begin(), it2 = x_l.begin(),
			it3 = delta_x_k.begin(), it4 = delta_x_l.begin();
			it != it_end && it2 != it_end2 && it3 != it_end3 && it4 != it_end4;
			++it, ++it2, ++it3, ++it4)
		{
			U += (muG_cpp_arma(x, y, *it, k, n) * *it3 - psiG_cpp_arma(x, y, *it, delta_x, k, n)) *
			(muG_cpp_arma(x, y, *it2, l, n) * *it4 - psiG_cpp_arma(x, y, *it2, delta_x, l, n));
		}

		return U / nx;
	}

/*-----------------------------------------------*/
/*sigma.ikl.L*/
double sigma_iklL_cpp_arma(arma::mat x, arma::mat y,
	arma::mat delta_x, arma::mat delta_y, int k, int l)
	{
		arma::vec x_k = x.col(k-1), delta_x_k = delta_x.col(k-1),
			x_l = x.col(l-1), delta_x_l = delta_x.col(l-1);
		int nx = x_k.size();
		double U = 0;

		arma::vec::const_iterator it_end = x_k.end(), it_end2 = x_l.end(),
			it_end3 = delta_x_k.end(), it_end4 = delta_x_l.end();

		for(arma::vec::const_iterator it = x_k.begin(), it2 = x_l.begin(),
			it3 = delta_x_k.begin(), it4 = delta_x_l.begin();
			it != it_end && it2 != it_end2 && it3 != it_end3 && it4 != it_end4;
			++it, ++it2, ++it3, ++it4)
		{
			U += (muL_cpp_arma(x, y, *it, k) * *it3 - psiL_cpp_arma(x, y, *it, delta_x, k)) *
			(muL_cpp_arma(x, y, *it2, l) * *it4 - psiL_cpp_arma(x, y, *it2, delta_x, l));
		}

		return U / nx;
	}

/*-----------------------------------------------*/
/*sigma.kl.G*/
double sigma_klG_cpp_arma(arma::mat x, arma::mat y,
	arma::mat delta_x, arma::mat delta_y, int k, int l,
	int n1, int n2)
{
	int n = n1+n2;
	double sigma_klg = 0;
	
	sigma_klg = sigma_iklG_cpp_arma(x,y,delta_x, delta_y,k,l,n) * n1 / n + 
	sigma_iklG_cpp_arma(y,x,delta_y,delta_x,k,l,n) * n2 / n;
	
	return sigma_klg;
}

/*-----------------------------------------------*/
/*sigma.kl.L*/
double sigma_klL_cpp_arma(arma::mat x, arma::mat y,
	arma::mat delta_x, arma::mat delta_y, int k, int l,
	int n1, int n2)
{
	int n = n1+n2;
	double sigma_klL = 0;

	sigma_klL = sigma_iklL_cpp_arma(x,y,delta_x,delta_y,k,l) * n1 / n + 
	sigma_iklL_cpp_arma(y,x,delta_y,delta_x,k,l) * n2 / n;

	return sigma_klL;
}

/*-----------------------------------------------*/
/*sigma.G*/
arma::mat sigma_G_cpp_arma(arma::mat x, arma::mat y,
	arma::mat delta_x, arma::mat delta_y, int k, int l,
	int n1, int n2, int p)
{
	arma::mat sigma(p,p);
	sigma.zeros();
	
	for(int i = 0; i < p ; i++){
		for(int j = i ; j < p ; j++){
			sigma (j, i) = sigma_klG_cpp_arma(x, y, delta_x, delta_y, i+1, j+1, n1, n2);
			sigma(i, j) = sigma (j, i);
		}
	}

	return sigma;
}

/*-----------------------------------------------*/
/*Sigma L*/
arma::mat sigma_L_cpp_arma(arma::mat x, arma::mat y, 
	arma::mat delta_x, arma::mat delta_y, int k, int l,
	int n1, int n2, int p)
{
	arma::mat sigma(p,p);
	sigma.zeros();

	for(int i = 0; i < p ; i++){
		for(int j = i; j < p; j++){
			sigma(j, i) = sigma_klL_cpp_arma(x, y, delta_x, delta_y, i+1, j+1, n1, n2);
			sigma(i, j) = sigma(j, i);
		}
	}
	return sigma;
}


//' Gehan test
//'
//' Performs the gehan test
//'
//' @param x Matrix
//' @param y Matrix
//' @param delta.x Matrix
//' @param delta.y Matrix
//' @param n1 Integer. Set as rows in x
//' @param n2 Integer. Set as rows in y
//' @param p Integer. No clue
//'
//' @return 1x1 matrix containing a numeric 
//'
//' @examples
//' data <- as.matrix(wltestdata)
//' x <- data[1:23, c(2, 4, 6, 8)]
//' y <- data[24:47, c(2, 4, 6, 8)]
//' delta.x <- data[1:23, c(3, 5, 7, 9)]
//' delta.y <- data[24:47, c(3, 5, 7, 9)]
//'
//' n1 <- dim(x)[1]
//' n2 <- dim(y)[1]
//' p <- dim(x)[2]
//'
//' gehan(x, y, delta.x, delta.y, n1, n2, p)
//' @export
// [[Rcpp::export]]
arma::mat gehan(arma::mat x, arma::mat y, arma::mat delta_x,
	arma::mat delta_y, int n1, int n2, int p, int k = 1, int l = 1)
{
	arma::Col<double> Ts(p);
	Ts.zeros();

	for(int i = 0; i < p; i++){
		Ts(i) = GehanTest_cpp_arma(x, y, n1, n2, delta_x, delta_y, i+1);
	}

	return Ts.t() * inv(sigma_G_cpp_arma(x, y, delta_x, delta_y, k, l, n1, n2, p)) * Ts;

}

//' Mvlogrank test
//'
//' Performs the mvlogrank test
//'
//' @param x Matrix
//' @param y Matrix
//' @param delta.x Matrix
//' @param delta.y Matrix
//' @param n1 Integer. Set as rows in x
//' @param n2 Integer. Set as rows in y
//' @param p Integer. No clue
//'
//' @return 1x1 matrix containing a numeric 
//'
//' @examples
//' data <- as.matrix(wltestdata)
//' x <- data[1:23, c(2, 4, 6, 8)]
//' y <- data[24:47, c(2, 4, 6, 8)]
//' delta.x <- data[1:23, c(3, 5, 7, 9)]
//' delta.y <- data[24:47, c(3, 5, 7, 9)]
//'
//' n1 <- dim(x)[1]
//' n2 <- dim(y)[1]
//' p <- dim(x)[2]
//'
//' mvlogrank(x, y, delta.x, delta.y, n1, n2, p)
//' @export
// [[Rcpp::export]]
arma::mat mvlogrank(arma::mat x, arma::mat y,
	arma::mat delta_x, arma::mat delta_y, int n1, int n2, 
	int p, int k = 1, int l = 1)
{
	arma::Col<double> Ts(p);
	Ts.zeros();

	for(int i = 0; i < p; i++){
		Ts(i) = mvlogrankTest_cpp_arma(x, y, n1, n2, delta_x, delta_y, i+1);
	}

	return Ts.t() * inv(sigma_L_cpp_arma(x, y, delta_x, delta_y, k, l, n1, n2, p)) * Ts;
}


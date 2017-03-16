#include "algebra.h"
#include "f2c.h"
#include "clapack.h"

vector<double> eigenValues(vector<double>matrix, long n) {
    double *a = &matrix[0];
    double *w=new double[n];
	long lda = n;
    double *work, wkopt;
    long lwork = -1;
    long info;
    dsyev_("N", "Upper", &n, a, &lda, w, &wkopt, &lwork, &info);
    lwork = (int) wkopt;
    work = new double[lwork];
    dsyev_("N", "Upper", &n, a, &lda, w, work, &lwork, &info);
    delete work;
    //cout<<"info"<<info<<endl;
    vector<double> ret;
    ret.assign(w, w + n);
	delete []w;
    return ret;
}
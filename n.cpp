#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <iomanip>

using namespace std;

double func(double x)
{
	return pow(x, 2) - 2*x + exp(-x);
}

void sort_through(double A, double B, double eps)
{
	int n = (B - A)/ eps;	//Число отрезков разбиения
	vector<double> Xvec(n);
	vector<double> Yvec(n);
	double Ymin = Yvec[0];
	int itmin;
	Xvec[0] = -1;
	for(int i = 0; i<n; i++)
	{
		Xvec[i+1] = Xvec[i] + eps;
	}
	cout<<endl;
	for(int i = 0; i<n; i++)
	{
		Yvec[i] = func(Xvec[i]);
		if(Yvec[i]<Ymin)
		{
			Ymin = Yvec[i];
			itmin = i;
		}

	}
	cout << "Метод перебора" <<endl;
	cout << "Координата точки минимума" <<endl;
	cout << fixed << setprecision(3)<<"("<<Xvec[itmin]<<"," << Yvec[itmin]<<")"<<endl<<endl;
}

void dichotomy(double A, double B, double eps)
{
	double a = A;
	double b = B;
	double delta = eps/2;
	double epsn = (b-a)/2;
	int iter = 1;

	cout <<"Метод дихотомии"<<endl;
	cout<<"Результаты вычислений"<<endl<<endl;;
	cout<<"№ Итерации"
		<<setw(5)<<"a"
		<<setw(8)<<"b"
		<<setw(9)<<"(b-a)/2"
		<<setw(5)<<"x1"
		<<setw(7)<<"x2"
		<<setw(8)<<"f(x1)"
		<<setw(8)<<"f(x2)"<<endl;
	while(epsn >= eps)
	{
		double x1 = (b + a - delta)/2;
		double x2 = (b + a + delta)/2;
		double fx1 = func(x1);
		double fx2 = func(x2);
		if (fx1 <= fx2)
		{
			b = x2;
		}
		else if (fx1 > fx2)
		{
			a = x1;
		}
		cout<<fixed <<setprecision(3)
			<<setw(6)<<iter
			<<setw(11)<<a
			<<setw(7)<<b
			<<setw(7)<<epsn
			<<setw(7)<<x1
			<<setw(7)<<x2
			<<setw(7)<<fx1
			<<setw(8)<<fx2<<endl	;
		epsn = (b-a)/2;
		iter++;
	}
	double xmin = (a+b)/2;
	double ymin = func(xmin);
	cout<<endl;
	cout << "Координата точки минимума" <<endl;
	cout<<"("<<xmin<<","<<ymin<<")"<<endl;
}

void golden_ratio(double A, double B, double eps)
{
	double a = A;
	double b = B;
	double t = (sqrt(5)-1)/2;
	double epsn = (b-a)/2;
	int iter = 1;
	cout <<"Метод Золотого сечения"<<endl;
		cout<<"Результаты вычислений"<<endl<<endl;;
		cout<<"№ Итерации"
			<<setw(5)<<"a"
			<<setw(8)<<"b"
			<<setw(9)<<"(b-a)/2"
			<<setw(5)<<"x1"
			<<setw(7)<<"x2"
			<<setw(8)<<"f(x1)"
			<<setw(8)<<"f(x2)"<<endl;
	while(epsn > eps)
	{
		double x1 = a + (3-sqrt(5))/2*(b - a);
		double x2 = a + (sqrt(5)-1)/2*(b - a);
		double fx1 = pow(x1, 2) - 2*x1 + exp(-x1);
		double fx2 = pow(x2, 2) - 2*x2 + exp(-x2);
		if(fx1 <= fx2)
		{
			b = x2;
			x2 = x1;
			fx2 = fx1;
			x1 = b - t*(b-a);
			fx1 = pow(x1, 2) - 2*x1 + exp(-x1);
		}
		else if(fx1 > fx2)
		{
			a = x1;
			x1 = x2;
			fx1 = fx2;
			x2 = b - t*(b-a);
			fx2 = pow(x2, 2) - 2*x2 + exp(-x2);

		}
		epsn *= t;
		cout<<fixed <<setprecision(3)
					<<setw(6)<<iter
					<<setw(11)<<a
					<<setw(7)<<b
					<<setw(7)<<epsn
					<<setw(7)<<x1
					<<setw(7)<<x2
					<<setw(7)<<fx1
					<<setw(8)<<fx2<<endl	;
				epsn = (b-a)/2;

		iter++;
	}
	double xmin = (a+b)/2;
	double ymin = pow(xmin, 2) - 2*xmin + exp(-xmin);
	cout << "Координата точки минимума" <<endl;
	cout<<"("<<xmin<<","<<ymin<<")"<<endl;

}

void bitwise(double A, double B, double eps)
{
	double a = A;
	double b = B;
	double del 	= (B-A)/4;
	double x0 = A;
	double x1 = A+del;
	double f0 = func(x0);
	double xmin;
	double ymin;

		label:
		x1 = x0 + del;
		double f1 = func(x1);
		if(fabs(del)<=eps)
				{
					xmin = x0;
					ymin = func(x0);
				}
		if(f0 >= f1)
		{
			if(a <= x0 && x0 <= b)

			goto label;
		}

		else
		{
			x0=x1;
			f0=f1;
			del=-del/4;
			cout<<del<<endl;
			goto label;
		}
	cout<<xmin<<endl;
	cout<<ymin<<endl;
}



int main()
{
	double A 	= -1;
	double B 	= 1.5;
	double eps 	= 0.001;

	//sort_through(A, B, eps);
	//dichotomy(A, B, eps);
	//golden_ratio(A, B, eps);
	bitwise(A, B, eps);

;
	return 0;
}

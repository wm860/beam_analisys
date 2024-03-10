#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include "winbgi2.h"
#include <time.h>
#include <fstream>
using namespace std;

const int mx = 10;
const int my = 10;

int fix[2 * (mx + 1)*(my + 1)];

#include "MesLib.h"
#include "solvers.h"
//#include "MatrixLib.h" //usunac!!!!!!!!!!!!!

double* alloc_vector(unsigned n)
{
	double *v = (double*)malloc(n * sizeof(double));
	for (unsigned i = 0; i < n; i++)
	{
		v[i] = 0.0;
	}
	return v;
}
double** alloc_matrix(unsigned n)
{
	double **mat = (double**)malloc(n * sizeof(double*));
	for (unsigned i = 0; i < n; i++)
	{
		mat[i] = alloc_vector(n);
	}
	return mat;
}
void free_vector(double* v)
{
	free(v);
}
void free_matrix(double**mat, unsigned n)
{
	for (unsigned i = 0; i < n; i++)
	{
		free(mat[i]);
	}
	free(mat);
}


int main()
{												//liczba wezlow*2
	const unsigned n = 2 * (mx + 1)*(my + 1);  //globalna liczba stopni swobody bo wezlow jest (mx+1)*(my+1) i kazdy wezel moze sie poruszac w 2 kierunkach
	double *d = alloc_vector(n);				//alokacja wektora przemieszczen
	double *F = alloc_vector(n);				//alokacja wektora sil
	double **Kg = alloc_matrix(n);			//alokacja macierzy sztywnosci

	//1 - assembling globalnej macierzy sztywnosci

	for (int b = 0; b < mx; b++)		//przechodzimy po wszystkich kolumnach
	{
		for (int a = 0; a < my; a++)	//przechodzimy po wszystkich wierszach
		{
			int global_index[8];
			for (int i = 0; i < 8; i++)
			{
				global_index[i] = DOF(b, a, i); //mapuje numeracje lokalna stopni swobody na globalna - znam globalny indeks stopnia swobody
			}

			for (int c = 0; c < 8; c++)		//przechodzimy po wszystkich wierszach lokalnej macierzy sztywnosci
			{
				for (int d = 0; d < 8; d++)		//przechodzimy po wszystkich kolumach lokalnej macierzy sztywnosci
				{
					Kg[global_index[c]][global_index[d]] = Kg[global_index[c]][global_index[d]] +Md*K[c][d];
				}
			}
		}
	}

	/*
	//2 - okreslenie wektora prawych stron (wektora sil)

	int g_index;		//zmienna pomocnicza
	g_index = P(mx/2, 0, 1); //globalny indeks przemieszczenia (wiersza w rownaniu) sily przylozonej do punktu mniej wiecej w polowie belki sila przylozona w kierunku pionowym
	F[g_index] = -0.1;		//sila

	for (int i = 0; i < mx+1; i++)		//obciazenie pod katem
	{
		g_index = P(i, my, 1); 
		F[g_index] = 0.1;
		g_index = P(i, my, 0);
		F[g_index] = 0.04;

	}
	
	//3 - warunki brzegowe
	for (int i = 0; i < n; i++)
	{
		fix[i] = 0;			//zerujemy wektor fix, ktory wskazuje nam gdzie znajduja sie stopien swobody,ktory bedzie odebrany 
	}
	//okreslam ktore stopnie swobody sa zamurowane
	
	for (unsigned i = 0; i < my+1; i++)	//funkcja P uzywam do tego
	{
		fix[P(mx, i, 0)] = 1;		//belka zamurowana
		fix[P(mx, i, 1)] = 1;
	}
	fix[P(0, 0, 1)] = 1;		//dodatkowy warunek brzegowy na poczatku belki
	
	//warunki brzegowe
	for (int i = 0; i < n; i++)
	{
		if (fix[i] == 1)
		{
			F[i] = 0;//zerujemy prawa strone itego rownania
			for (int a = 0; a < n; a++)	//zerowanie wszystkiego w danym wierszu
			{
				Kg[i][a] = 0;
			}
			for (int a = 0; a < n; a++)	//zerowanie wszystkiego w danej kolumnie
			{
				Kg[a][i] = 0;
			}
			Kg[i][i] = 1; //wpisanie na diagonali 1 zgodnie z instrukcja
		}
	}
	*/
	//NARZUCENIE OBCI¥¯ENIA
	F[7] = -0.05;

	//WYZEROWANIE WEKTORA STOPNII SWOBODY
	for (unsigned i = 0; i < n; i++)
		fix[i] = 0;

	//ODBIERANIE STOPNII SWOBODY, tj. tam gdzie stoi 1 mamy utwierdzenie - w tym przypadku ca³a lewa strona jest zamurowana 
	for (unsigned i = 0; i < my; i++)
	{
		fix[DOF(0, i, 0)] = 1;
		fix[DOF(0, i, 1)] = 1;
		fix[DOF(0, i, 6)] = 1;
		fix[DOF(0, i, 7)] = 1;
	}

	//WYZEROWANIE ODPOWIEDNICH WIERSZY I KOLUMN W RÓWNANIU MACIERZOWYM, TJ. NARZUCENIE WARUNKÓW BRZEGOWYCH 
	for (unsigned i = 0; i < n; i++)
	{
		if (fix[i] == 1)
			for (unsigned j = 0; j < n; j++)
			{
				if (i == j)
					Kg[i][j] = 1;
				else
					Kg[i][j] = 0;
				F[i] = 0;
			}
	}

	//4 - rozwiazanie
	clock_t start, stop;
	double time;
	double *r = alloc_vector(n);
	double par;
	
	
	ofstream zapis("wyncpp.txt");

	FILE* p; 
	p = fopen("wyniki.txt", "w");
	/*
	residuum(Kg, d, F, n, r);

	cout.precision(10);
	cout<<endl<<"przed gaussem  "<<nor_vector(r, n)<<endl<<endl;
	start = clock();
	Gauss(n, Kg, F, d);	
	stop = clock();
	time = (clock()-start)*1000;
	cout << "czas gaussa " << (clock() - start) << endl;
	residuum(Kg, d, F, n, r);
	cout << endl<< "norma wektora po Gaussie  "<<nor_vector(r, n)<<endl;*/

	
	for (int i = 0; i < n; i++)
	{
		d[i]=0.;
	}
	residuum(Kg, d, F, n, r);
	cout << endl << "norma residuum przed Jacobim\t" << nor_vector(r, n) << endl;
	start = clock();
	Jacob(Kg, F, n, d);
	residuum(Kg, d, F, n, r);
	cout << endl << "norma residuum po Jacobim\t" << nor_vector(r, n)<<endl;
	stop = clock();
	//time = (clock() - start) * 1000;
	cout << "czas po jakobim " << (clock() - start) << endl;
	fprintf(p, "%lf\t%lf\n", nor_vector(r, n), (clock() - start));
	zapis << nor_vector(r, n) << "     " << (clock() - start) << " czas\n";

	
	/*
	

	start = clock();
	Jacobpod(Kg, F, n, d);
	residuum(Kg, d, F, n, r);
	cout << endl << "norma residuum po Jacobimpod\t" << nor_vector(r, n) << endl;
	stop = clock();
	time = (clock() - start) * 1000;
	cout << "czas jakobimpod  " << (clock() - start) << endl;
	fprintf(p, "%lf\t%lf\n", nor_vector(r, n), (clock() - start));
	zapis << nor_vector(r, n) << "     " << (clock() - start) << " czas\n";

	start = clock();
	Gauss_Seidl(Kg, F, n, d);
	residuum(Kg, d, F, n, r);
	cout << endl << "norma residuum po GAUSS SEIDL\t" << nor_vector(r, n) << endl;
	stop = clock();
	time = (clock() - start) * 1000;
	cout << "czas gaussa Seidla   " << (clock() - start) << endl;
	fprintf(p, "%lf\t%lf\n", nor_vector(r, n), (clock() - start));
	zapis << nor_vector(r, n) << "     " << (clock() - start) << " czas\n";

	start = clock();
	Minres(Kg, F, n, d);
	residuum(Kg, d, F, n, r);
	cout << endl << "norma residuum po Minres\t" << nor_vector(r, n) << endl;
	stop = clock();
	time = (clock() - start) * 1000;
	cout << "czas MINRES  " << (clock() - start) << endl;
	fprintf(p, "%lf\t%lf\n", nor_vector(r, n), (clock() - start));
	zapis << nor_vector(r, n) << "     " << (clock() - start) << " czas\n";

	graphics(1000, 1000);
	scale(0.0 - 2, 0.5 * (my - mx - 3), mx + 3, 0.5 * (my + mx + 3));
	title("x", "y", "ugiecie belki");
	draw(d, F);
	wait();

	start = clock();
	Minres2(Kg, F, n, d);
	residuum(Kg, d, F, n, r);
	cout << endl << "norma residuum po Minres2\t" << nor_vector(r, n) << endl;
	stop = clock();
	time = (clock() - start) * 1000;
	cout << "czas MINRES2  " << (clock() - start) << endl;
	fprintf(p, "%lf\t%lf\n", nor_vector(r, n), (clock() - start));
	zapis << nor_vector(r, n) << "     " << (clock() - start) << " czas\n";
	*/

	start = clock();
	Grad(Kg, F, n, d);
	residuum(Kg, d, F, n, r);
	cout << endl << "norma residuum po Grad\t" << nor_vector(r, n) << endl;
	stop = clock();
	time = (clock() - start) * 1000;
	cout << "czas Grad  " << (clock() - start) << endl;
	fprintf(p, "%lf\t%lf\n", nor_vector(r, n), (clock() - start));
	zapis << nor_vector(r, n) << "     " << (clock() - start) << " czas\n";


	start = clock();
	Gradpod(Kg, F, n, d);
	residuum(Kg, d, F, n, r);
	cout << endl << "norma residuum po Gradpod\t" << nor_vector(r, n) << endl;
	stop = clock();
	time = (clock() - start) * 1000;
	cout << "czas Gradpod  " << (clock() - start) << endl;
	fprintf(p, "%lf\t%lf\n", nor_vector(r, n), (clock() - start));
	zapis << nor_vector(r, n) << "     " << (clock() - start) << " czas\n";


	start = clock();
	Gradpod2(Kg, F, n, d);
	residuum(Kg, d, F, n, r);
	cout << endl << "norma residuum po Gradpod2\t" << nor_vector(r, n) << endl;
	stop = clock();
	time = (clock() - start) * 1000;
	cout << "czas Gradpod2  " << (clock() - start) << endl;
	fprintf(p, "%lf\t%lf\n", nor_vector(r, n), (clock() - start));
	zapis << nor_vector(r, n) << "     " << (clock() - start) << " iteracji\n\n\n";

	/*
	for (int i = 0; i < n; i++)
	{
		d[i] = 0.;
	}
	start = clock();
	MetodaJacobiego(Kg, d, F, n);
	residuum(Kg, d, F, n, r);
	cout << endl << "norma residuum po Jacobim od Antka\t" << nor_vector(r, n) << endl;
	stop = clock();
	time = (clock() - start) * 1000;
	cout << "czas Antka  " << (clock() - start) << endl;
	*/

	/* tests
	double **A = alloc_matrix(3);
	double *h=alloc_vector(3);
	double *c = alloc_vector(3);
	h[0] = 3; h[1] = 6; h[2] = 9;
	c[0] = 1; c[1] = 2; c[2] = 3;
	double q = 1.0;
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			A[i][j] = q;
			printf("%lf\t", A[i][j]);
			q++;
		}
		printf("\n");
	}
	double *wyn = alloc_vector(3);
	wyn=sca_vector(h, 5, 3);
	for (int i = 0; i < 3; i++)
	{
		printf("%lf\n", wyn[i]);
	}
	*/

	//5 - rysowanie wynikow
	
	graphics(1000, 1000);
	scale(0.0-2, 0.5 * (my - mx - 3), mx + 3, 0.5 * (my + mx + 3));
	title("x", "y", "ugiecie belki");
	draw(d, F);
	
	
	wait();
	wait();
	wait();
	wait();
	wait();

	//zwolnienie miejsca w pamieci
	free_vector(d);
	free_vector(F);
	free_vector(r);
	free_matrix(Kg, n);
	fclose(p);
	zapis.close(); //obowi¹zkowo nale¿y zamkn¹æ plik
	return 0;
}
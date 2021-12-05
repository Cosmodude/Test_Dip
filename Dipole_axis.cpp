#include<iostream>
#include<fstream>
#include<iomanip>
#include <string>
#include<cstring>
#include <random>
#include <ctime>
#include <vector>

#define PI 3.14159265

using namespace std;
#pragma comment(linker, "/STACK:10000000")

class Galaxy
{
public:
	double ra;
	double dec;
	long long cw, ccw;
	double cosine(double alpha, double dzeta)
	{
		double cosine = cos(dzeta * PI / 180) * cos(dec * PI / 180) * cos((alpha - ra) * PI / 180) + sin(dzeta * PI / 180) * sin(dec * PI / 180);
		return cosine;
	}
};


inline void read_rand(string file, vector <vector<Galaxy>>(&mass))
{
	srand(static_cast<unsigned int>(time(0)));

	mass.clear();
	string str;
	vector <Galaxy> vec;
	Galaxy g;
	for (int j = 0; j < 1000; j++)
	{
		ifstream fp(file.c_str());
		if (!fp.is_open())
		{
			cout << "not opened" << endl;
		}
		vec.clear();
		while (fp >> g.cw)
		{
			fp.get();
			fp >> g.ccw;
			fp.get();
			fp >> g.ra;
			fp.get();
			fp >> g.dec;
			getline(fp, str);

			int val = rand();
			if (val % 2 == 1)
			{
				g.cw = 1;
				g.ccw = 0;
			}
			else
			{
				g.cw = 0;
				g.ccw = 1;
			}
			vec.push_back(g);
		}
		mass.push_back(vec);
	}
}

inline void read_t(vector <Galaxy>(&mass), string file)
{
	mass.clear();
	string str;
	Galaxy g;
	ifstream fp(file.c_str());
	if (!fp.is_open())
	{
		cout << "not opened" << endl;
	}
	while (fp >> g.cw)
	{
		fp.get();
		fp >> g.ccw;
		fp.get();
		fp >> g.ra;
		fp.get();
		fp >> g.dec;
		getline(fp, str);
		mass.push_back(g);
	}
}

double func_real(double alpha, double dzeta, vector <Galaxy>(&mass))
{
	double sum = 0;
	vector<double> cosine;

	for (int i = 0; i < mass.size(); i++)
	{
		cosine.push_back(mass[i].cosine(alpha, dzeta));

		if (abs(cosine[i]) * (mass[i].cw + mass[i].ccw * (-1)) == cosine[i])
		{
			sum++;
		}
	}
	//cout << sum << endl;
	return pow(mass.size() - sum, 2) / mass.size();
}

double func_rand(double alpha, double dzeta, vector <vector<Galaxy>>(&mass), double(&sumr))
{
	double sumra[1000];
	double sigma = 0;
	double sum = 0;
	int t = 0;
	for (int i = 0; i < 1000; i++)
	{
		sumra[i]= func_real(alpha, dzeta, mass[i]);
		sum += sumra[i];
	}
	//cout << sumra[2];
	sumr = sum / 1000;

	double x = 0;
	for (int i = 0; i < 1000; i++)
	{
		x += pow(sumr - sumra[i], 2);
	}

	sigma = sqrt(x / 1000);

	return sigma;
}
int main()
{
	vector<Galaxy> tdata;
	vector <vector <Galaxy>> fdata;
	read_t(tdata, "catalog.csv");
	read_rand("catalog.csv", fdata);
	double sumr;
	ofstream fout("test_dip.txt");
	if (!fout.is_open())
	{
		cout << "not opened" << endl;
	}
	//cout << fdata[5].size()<<endl;
	//cout << func_real(165, 40, tdata) << endl;
	//cout << func_rand(165, 40, fdata,sumr) << endl;
	//cout <<  fdata[2][1].dec << endl;
	//cout << sumr << endl;
	/*cout << func_rand(165, 40, fdata, sumr) << endl;
	cout << data[0].ccw << ", " << data[0].ra << ", " << endl;
	cout << data[77839].ccw << ", " << data[77839].ra << ", " << endl;
	cout << data[0].cosine(165, 40) << endl;*/
	for (double ra =0; ra < 361; ra += 5)
	{	
		for (double dec = -90; dec < 91; dec += 5)
		{
			double temp;
			sumr = 0;
			temp = func_rand(ra, dec, fdata, sumr);
			fout << ra << "\t" << dec << "\t" << abs(sumr - func_real(ra, dec, tdata) )/ temp << endl;
		}
		fout << endl;
	}



	cout << "done" << endl;

	cout << "galaxy" << endl;
	system("pause");
	return 0;
}
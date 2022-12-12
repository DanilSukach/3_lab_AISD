#include <iostream>
#include <locale>
#include <cmath>
#include <iomanip>
#include <complex>
#include <forward_list> 
#include <iterator>
using namespace std;

template<typename T>
class Polynomial {
	forward_list<T> values;
	double eps = 0.0001;
public:
	Polynomial() = default;
	~Polynomial() = default;
	Polynomial(const Polynomial& v) = default;
	Polynomial& operator =(const Polynomial& v) = default;
	Polynomial(const size_t n) {
		for (size_t i = 0; i <= n; i++) {
			values.push_front(1);
		}
	}
	bool operator ==(const Polynomial& v) {
		auto i1 = values.begin();
		auto i2 = v.values.begin();
		auto b1 = values.end();
		auto b2 = v.values.end();
		while (i1 != b1 || i2 != b2) {
			if (i1 != b1 && i2 != b2) {
				if (*i1 == *i2) {
					i1++;
					i2++;
				}
				else return false;
			}
			else {
				if (i1 != b1) {
					if (*i1 == 0)i1++;
					else return false;
				}
				if (i2 != b2) {
					if (*i2 == 0)i2++;
					else return false;
				}
			}
		}
		return true;
	}
	void Set(const T& value, const size_t ind) {
		if (ind < 0) throw "invalid index";
		auto i1 = values.begin();
		auto i2 = values.end();
		for (size_t i = 0; i < ind; i++) {
			i1++;
			if (i1 == i2)throw"invalid index";
		}
		*i1 = value;
	}
	T& operator[](const size_t ind){
		if (ind < 0) throw "invalid index";
		auto i1 = values.begin();
		auto i2 = values.end();
		for (size_t i = 0; i < ind; i++) {
			i1++;
			if (i1 == i2)throw"invalid index";
		}
		return *i1;
	}
	friend ostream& operator<<(ostream& os, const Polynomial& v) {
		size_t ind = 0;
		for (const auto& n : v.values) {
			if (n == 0) {
				ind++;
			}
			else {
				if (n > 0) {
					if (ind == 0) {
						os << n;
					}
					else {
						if (n == 1) {
							os << "+x^" << ind;
						}
						else {
							os << "+" << n << "x^" << ind;
						}
					}
				}
				else {
					if (ind == 0) {
						os << n;
					}
					else {
						if (n == -1) {
							os << "-x^" << ind;
						}
						else {
							os << n << "x^" << ind;
						}
					}
				}
				ind++;
			}
		}
		return os;
	}

	Polynomial operator +(const Polynomial& v) const {
		Polynomial temp;
		T temp_value;
		auto i1 = values.begin();
		auto i2 = v.values.begin();
		auto b1 = values.end();
		auto b2 = v.values.end();
		while (i1 != b1 || i2 != b2) {
			if (i1 != b1 && i2 != b2) {
				temp_value = *i1 + *i2;
				temp.values.push_front(temp_value);
				i1++;
				i2++;
			}
			else {
				if (i1 != b1) {
					temp.values.push_front(*i1);
					i1++;
				}
				if (i2 != b2) {
					temp.values.push_front(*i2);
					i2++;
				}
			}
		}
		temp.values.reverse();
		return temp;
	}
	Polynomial operator -(const Polynomial& v) const {
		Polynomial temp;
		T temp_value;
		auto i1 = values.begin();
		auto i2 = v.values.begin();
		auto b1 = values.end();
		auto b2 = v.values.end();
		while (i1 != b1 || i2 != b2) {
			if (i1 != b1 && i2 != b2) {
				temp_value = *i1 - *i2;
				temp.values.push_front(temp_value);
				i1++;
				i2++;
			}
			else {
				if (i1 != b1) {
					temp.values.push_front(*i1);
					i1++;
				}
				if (i2 != b2) {
					temp.values.push_front(-*i2);
					i2++;
				}
			}
		}
		temp.values.reverse();
		return temp;
	}
	Polynomial& operator *=(const T& value) {
		for (auto& n : values) {
			n *= value;
		}
		return *this;
	}
	Polynomial operator *(const T& value) const {
		Polynomial temp;
		T value_temp;
		for (const auto& n : values) {
			value_temp = n * value;
			temp.values.push_front(value_temp);
		}
		temp.values.reverse();
		return temp;
	}
	friend Polynomial operator *(const T& value, const Polynomial& v) {
		Polynomial temp;
		T value_temp;
		for (const auto& n : v.values) {
			value_temp = n * value;
			temp.values.push_front(value_temp);
		}
		temp.values.reverse();
		return temp;
	}
	double calculation(const T& x) const {
		double result = 0;
		size_t ind = 0;
		for (const auto& n : values){
			result += n * pow(x, ind);
			ind++;
		}
		return result;
	}
	void formula_Kardano(double& x1, double& x2, double& x3, double& Q) const {
		auto i1 = values.begin();
		auto i2 = values.end();
		T d = *i1;
		i1++;
		if (i1 == i2)throw"invalid Polinomial";
		T c = *i1;
		i1++;
		if (i1 == i2)throw"invalid Polinomial";
		T b = *i1;
		i1++;
		if (i1 == i2)throw"invalid Polinomial";
		T a = *i1;
		i1++;
		if (i1 != i2)throw"invalid Polinomial";
		double p = (3 * static_cast<double>(a) * static_cast<double>(c) - pow(b, 2)) / (3 * pow(a, 2));
		double q = (2 * pow(b, 3) - 9 * static_cast<double>(a) * static_cast<double>(b) * static_cast<double>(c) + 27 * pow(a, 2) * d) / (27 * pow(a, 3));
		Q = pow((p / 3), 3) + pow((q / 2), 2);
		cout << "Q = " << Q << endl;
		double u = cbrt((-q / 2) + sqrt((pow(q, 2) / 4) + (pow(p, 3) / 27)));
		double v = cbrt((-q / 2) - sqrt((pow(q, 2) / 4) + (pow(p, 3) / 27)));
		if (Q > 0) {
			double y = u + v;
			x1 = y - (b / (3 * a));
		}
		if (abs(Q) == 0) {
			double y1 = 2 * cbrt(-q / 2);
			double y2 = -cbrt(-q / 2);
			x1 = y1 - (b / (3 * a));
			x2 = x3 = y2 - (b / (3 * a));
		}
		if (Q < 0) {
			double fi = 0;
			if (abs(q) < 0.1) {
				fi = 3.1415926535 / 2;
			}
			if (q > 0) {
				fi = atan(sqrt(-Q) / (-q / 2)) + 3.1415926535;
			}
			if (q < 0) {
				fi = atan(sqrt(-Q) / (-q / 2));
			}
			x1 = (2 * sqrt(-p / 3) * cos(fi / 3)) - (b / (3 * a));
			x2 = (2 * sqrt(-p / 3) * cos((fi + 2 * 3.1415926535) / 3)) - (b / (3 * a));
			x3 = (2 * sqrt(-p / 3) * cos((fi + 4 * 3.1415926535) / 3)) - (b / (3 * a));
		}
	}
};

template <typename T>
class Polynomial<complex<T>> {
	forward_list<complex<T>> values;
	double eps = 0.0001;
public:
	Polynomial() = default;
	~Polynomial() = default;
	Polynomial(const Polynomial & v) = default;
	Polynomial& operator =(const Polynomial & v) = default;
	Polynomial(const size_t n) {
		complex<T> x(1, 1);
		for (size_t i = 0; i <= n; i++) {
			values.push_front(x);
		}
	}
	complex<T>& operator[](const size_t ind){
		if (ind < 0) throw "invalid index";
		auto i1 = values.begin();
		auto i2 = values.end();
		for (size_t i = 0; i < ind; i++) {
			i1++;
			if (i1 == i2)throw"invalid index";
		}
		return *i1;
	}
	void Set(const complex<T>& value, const size_t ind) {
		if (ind < 0) throw "invalid index";
		auto i1 = values.begin();
		auto i2 = values.end();
		for (size_t i = 0; i < ind; i++) {
			i1++;
			if (i1 == i2)throw"invalid index";
		}
		*i1 = value;
	}
	bool operator ==(const Polynomial<complex<T>>& v) {
		Polynomial<complex<T>> temp;
		auto i1 = values.begin();
		auto i2 = v.values.begin();
		auto b1 = values.end();
		auto b2 = v.values.end();
		while (i1 != b1 || i2 != b2) {
			if (i1 != b1 && i2 != b2) {
				if (*i1 == *i2) {
					i1++;
					i2++;
				}
				else return false;
			}
			else {
				if (i1 != b1) {
					complex<T> x(0, 0);
					if (*i1 == x)i1++;
					else return false;
				}
				if (i2 != b2) {
					complex<T> x(0, 0);
					if (*i2 == x)i2++;
					else return false;
				}
			}
		}
		return true;
	}
	friend ostream& operator<<(ostream& os, const Polynomial<complex<T>>& v) {
		size_t ind = 0;
		for(const auto& n : v.values) {
			T tempr = real(n);
			T tempi = imag(n);
			if ((tempr != 0) & (tempi != 0)) {

				if (ind == 0) {
					os << n << "+";
				}
				else  {
					os << n << "x^" << ind << "+";
				}
			}
			ind++;
		}
		return os;
	}
	Polynomial operator +(const Polynomial<complex<T>>& v) const {
		Polynomial<complex<T>> temp;
		complex<T> temp_value;
		auto i1 = values.begin();
		auto i2 = v.values.begin();
		auto b1 = values.end();
		auto b2 = v.values.end();
		while (i1 != b1 || i2 != b2) {
			if (i1 != b1 && i2 != b2) {
				temp_value = *i1 + *i2;
				temp.values.push_front(temp_value);
				i1++;
				i2++;
			}
			else {
				if (i1 != b1) {
					temp.values.push_front(*i1);
					i1++;
				}
				if (i2 != b2) {
					temp.values.push_front(*i2);
					i2++;
				}
			}
		}
		temp.values.reverse();
		return temp;
	}
	Polynomial operator -(const Polynomial<complex<T>>& v) const {
		Polynomial<complex<T>> temp;
		complex<T> temp_value;
		auto i1 = values.begin();
		auto i2 = v.values.begin();
		auto b1 = values.end();
		auto b2 = v.values.end();
		while (i1 != b1 || i2 != b2) {
			if (i1 != b1 && i2 != b2) {
				temp_value = *i1 - *i2;
				temp.values.push_front(temp_value);
				i1++;
				i2++;
			}
			else {
				if (i1 != b1) {
					temp.values.push_front(*i1);
					i1++;
				}
				if (i2 != b2) {
					temp.values.push_front(-*i2);
					i2++;
				}
			}
		}
		temp.values.reverse();
		return temp;
	}
	Polynomial<complex<T>>& operator *=(const complex<T>& value) {
		for (auto& n : values) {
			n *= value;
		}
		return *this;
	}
	Polynomial<complex<T>> operator *(const complex<T>& value) const {
		Polynomial<complex<T>> temp;
		complex<T> value_temp;
		for (const auto& n : values) {
			value_temp = n * value;
			temp.values.push_front(value_temp);
		}
		temp.values.reverse();
		return temp;
	}
	friend Polynomial<complex<T>> operator *(const complex<T>& value, const Polynomial<complex<T>>& v) {
		Polynomial<complex<T>> temp;
		complex<T> value_temp;
		for (const auto& n : v.values) {
			value_temp = n * value;
			temp.values.push_front(value_temp);
		}
		temp.values.reverse();
		return temp;
	}
	complex<T> calculation(const complex<T>& x) const {
		complex<T> result(0, 0);
		size_t ind = 0;
		for (const auto& n : values) {
			complex<T> temp(pow(x, ind));
			result += n * temp;
			ind++;
		}
		return result;
	}
	void formula_Kardano(T& x1, T& x2, T& x3, double& Q, bool& flag) const {
		auto i1 = values.begin();
		auto i2 = values.end();
		complex<T> d = *i1;
		T dr = real(d);
		T di = imag(d);
		i1++;
		if (i1 == i2)throw"invalid Polinomial";
		complex<T> c = *i1;
		T cr = real(c);
		T ci = imag(c);
		i1++;
		if (i1 == i2)throw"invalid Polinomial";
		complex<T> b = *i1;
		T br = real(b);
		T bi = imag(b);
		i1++;
		if (i1 == i2)throw"invalid Polinomial";
		complex<T> a = *i1;
		T ar = real(a);
		T ai = imag(a);
		i1++;
		if (i1 != i2)throw"invalid Polinomial";
		if ((ai == 0) & (bi == 0) & (ci == 0) & (di == 0)) {
			flag = true;
			double p = (3 * static_cast<double>(ar) * static_cast<double>(cr) - pow(br, 2)) / (3 * pow(ar, 2));
			double q = (2 * pow(br, 3) - 9 * static_cast<double>(ar) * static_cast<double>(br) * static_cast<double>(cr) + 27 * pow(ar, 2) * dr) / (27 * pow(ar, 3));
			Q = pow((p / 3), 3) + pow((q / 2), 2);
			cout << "Q = " << Q << endl;
			double u = cbrt((-q / 2) + sqrt((pow(q, 2) / 4) + (pow(p, 3) / 27)));
			double v = cbrt((-q / 2) - sqrt((pow(q, 2) / 4) + (pow(p, 3) / 27)));
			if (Q > 0) {
				double y = u + v;

				x1 = static_cast<T>(y - (br / (3 * ar)));
			}
			if (abs(Q) == 0) {
				double y1 = 2 * cbrt(-q / 2);
				double y2 = -cbrt(-q / 2);
				x1 = static_cast<T>(y1 - (br / (3 * ar)));
				x2 = x3 = static_cast<T>(y2 - (br / (3 * ar)));
			}
			if (Q < 0) {
				double fi = 0;
				if (abs(q) < 0.1) {
					fi = 3.1415926535 / 2;
				}
				if (q > 0) {
					fi = atan(sqrt(-Q) / (-q / 2)) + 3.1415926535;
				}
				if (q < 0) {
					fi = atan(sqrt(-Q) / (-q / 2));
				}
				x1 = static_cast<T>((2 * sqrt(-p / 3) * cos(fi / 3)) - (br / (3 * ar)));
				x2 = static_cast<T>((2 * sqrt(-p / 3) * cos((fi + 2 * 3.1415926535) / 3)) - (br / (3 * ar)));
				x3 = static_cast<T>((2 * sqrt(-p / 3) * cos((fi + 4 * 3.1415926535) / 3)) - (br / (3 * ar)));
			}
		}
		else {
			flag = false;
		}

	}
};



void menu_1() {
	cout << "1.int" << endl;
	cout << "2.float" << endl;
	cout << "3.double" << endl;
	cout << "4.complex<float>" << endl;
	cout << "5.complex<double>" << endl;
	cout << "0.Выход" << endl;
}

void menu_2() {
	cout << "1.Установить коэффициент" << endl;
	cout << "2.Вывести многочлен" << endl;
	cout << "3.Вывести коэффициент при заданной степени" << endl;
	cout << "4.Сложение" << endl;
	cout << "5.Вычитание" << endl;
	cout << "6.Умножение на скаляр" << endl;
	cout << "7.Вычисление значения в Х" << endl;
	cout << "8.Нахождение действительных корней" << endl;
	cout << "9.Проверка на равенство" << endl;
	cout << "0.Выход" << endl;
}
int main() {
	setlocale(LC_ALL, "Russian");
	while (true) {
		system("cls");
		menu_1();
		size_t flag_1;
		cin >> flag_1;
		switch (flag_1)
		{
		case 1:
		{
			bool fg = true;
			system("cls");
			size_t number_1;
			size_t number_2;
			cout << "Введите степень первого многочлена" << endl;
			cin >> number_1;
			cout << "Введите степень второго многочлена" << endl;
			cin >> number_2;
			Polynomial<int> p1(number_1);
			Polynomial <int>p2(number_2);
			while (fg) {
				system("cls");
				menu_2();
				size_t flag_2_1;
				cin >> flag_2_1;
				switch (flag_2_1)
				{
				case 1:
				{
					system("cls");
					size_t flag_2;
					cout << "У какого многочлена?" << endl;
					cin >> flag_2;
					switch (flag_2)
					{
					case 1:
					{
						size_t ind;
						int number;
						size_t flag_3 = 1;
						while (flag_3) {
							system("cls");
							cout << "Введите степень" << endl;
							cin >> ind;
							cout << "Введите значение" << endl;
							cin >> number;
							try {
								p1.Set(number, ind);
							}
							catch (const char* msg) {
								cout << msg << endl;
							}
							cout << "Продолжить 1" << endl;
							cout << "Выход 0" << endl;
							cin >> flag_3;
						}
						break;
					}
					case 2:
					{
						size_t ind;
						int number;
						size_t flag_3 = 1;
						while (flag_3) {
							system("cls");
							cout << "Введите степень" << endl;
							cin >> ind;
							cout << "Введите значение" << endl;
							cin >> number;
							try {
								p2.Set(number, ind);
							}
							catch (const char* msg) {
								cout << msg << endl;
							}
							cout << "Продолжить 1" << endl;
							cout << "Выход 0" << endl;
							cin >> flag_3;
						}
						break;
					}
					}
					break;
				}
				case 2:
				{
					system("cls");
					cout << p1 << endl;
					cout << p2 << endl;
					system("pause");
					break;
				}
				case 3:
				{
					system("cls");
					size_t flag_2;
					cout << "У какого многочлена?" << endl;
					cin >> flag_2;
					switch (flag_2)
					{
					case 1:
					{
						size_t ind;
						cout << "Введите степень" << endl;
						cin >> ind;
						try {
							cout << p1[ind] << endl;
						}
						catch (const char* msg) {
							cout << msg << endl;
						}
						system("pause");
						break;
					}
					case 2:
					{
						size_t ind;
						cout << "Введите степень" << endl;
						cin >> ind;
						try {
							cout << p2[ind] << endl;
						}
						catch (const char* msg) {
							cout << msg << endl;
						}
						system("pause");
						break;
					}
					}
					break;
				}
				case 4:
				{
					system("cls");
					Polynomial <int> p3 = p1 + p2;
					cout << p1 << " + " << p2 << " = " << p3 << endl;
					system("pause");
					break;
				}
				case 5:
				{
					system("cls");
					size_t flag_2;
					cout << "Из какого многочлена вычесть?" << endl;
					cin >> flag_2;
					switch (flag_2)
					{
					case 1:
					{
						Polynomial <int> p3 = p1 - p2;
						cout << p1 << " - " << "(" << p2 << ")" << " = " << p3 << endl;
						system("pause");
						break;
					}
					case 2:
					{
						Polynomial <int> p3 = p2 - p1;
						cout << p2 << " - " << "(" << p1 << ")" << " = " << p3 << endl;
						system("pause");
						break;
					}
					}
					break;
				}
				case 6:
				{
					system("cls");
					int value;
					int flag_3;
					cout << "На что умножить, тип int?" << endl;
					cin >> value;
					cin.ignore();
					cout << "1. Изменить многочлен" << endl;
					cout << "2. Создать новый многочлен" << endl;
					cin >> flag_3;
					switch (flag_3)
					{
					case 1:
					{
						cout << p1 << " * " << value << " = ";
						p1 *= value;
						cout << p1 << endl;
						cout << p2 << " * " << value << " = ";
						p2 *= value;
						cout << p2 << endl;
						break;
					}
					case 2:
					{
						Polynomial <int> p3 = p1 * value;
						Polynomial <int> p4 = value * p2;
						cout << p1 << " * " << value << " = " << p3 << endl;
						cout << value << " * " << p2 << " = " << p4 << endl;
						break;
					}
					}
					system("pause");
					break;
				}
				case 7:
				{
					system("cls");
					int value;
					cout << "Введите Х" << endl;
					cin >> value;
					cout << p1 << endl;
					cout << p2 << endl;
					cout << p1.calculation(value) << endl;
					cout << p2.calculation(value) << endl;
					system("pause");
					break;
				}
				case 8:
				{
					system("cls");
					cout << p1 << " = 0" << endl;
					cout << p2 << " = 0" << endl;
					double x1, x2, x3, Q;
					try {
						p1.formula_Kardano(x1, x2, x3, Q);
						if (Q > 0) {
							cout << "x1 = " << x1 << endl;
						}
						if (abs(Q) < 0.01) {
							cout << "x1 = " << x1 << endl;
							cout << "x2 = " << x2 << endl;
							cout << "x3 = " << x3 << endl;
						}
						if (Q < 0) {
							cout << "x1 = " << x1 << endl;
							cout << "x2 = " << x2 << endl;
							cout << "x3 = " << x3 << endl;
						}
					}
					catch (const char* msg) {
						cout << msg << endl;
					}
					try {
						p2.formula_Kardano(x1, x2, x3, Q);
						if (Q > 0) {
							cout << "x1 = " << x1 << endl;
						}
						if (abs(Q) < 0.01) {
							cout << "x1 = " << x1 << endl;
							cout << "x2 = " << x2 << endl;
							cout << "x3 = " << x3 << endl;
						}
						if (Q < 0) {
							cout << "x1 = " << x1 << endl;
							cout << "x2 = " << x2 << endl;
							cout << "x3 = " << x3 << endl;
						}
					}
					catch (const char* msg) {
						cout << msg << endl;
					}
					system("pause");
					break;
				}
				case 9: {
					system("cls");
					bool f;
					f = p1 == p2;
					cout << f << endl;
					system("pause");
					break;
				}
				case 0:
				{
					fg = false;
				}
				}
			}
			break;
		}
		case 2:
		{
			bool fg = true;
			system("cls");
			size_t number_1;
			size_t number_2;
			cout << "Введите степень первого многочлена" << endl;
			cin >> number_1;
			cout << "Введите степень второго многочлена" << endl;
			cin >> number_2;
			Polynomial<float> p1(number_1);
			Polynomial <float>p2(number_2);
			while (fg) {
				system("cls");
				menu_2();
				size_t flag_2_2;
				cin >> flag_2_2;
				switch (flag_2_2)
				{
				case 1:
				{
					system("cls");
					size_t flag_2;
					cout << "У какого многочлена?" << endl;
					cin >> flag_2;
					switch (flag_2)
					{
					case 1:
					{
						size_t ind;
						float number;
						size_t flag_3 = 1;
						while (flag_3) {
							system("cls");
							cout << "Введите степень" << endl;
							cin >> ind;
							cout << "Введите значение" << endl;
							cin >> number;
							try {
								p1.Set(number, ind);
							}
							catch (const char* msg) {
								cout << msg << endl;
							}
							cout << "Продолжить 1" << endl;
							cout << "Выход 0" << endl;
							cin >> flag_3;
						}
						break;
					}
					case 2:
					{
						size_t ind;
						float number;
						size_t flag_3 = 1;
						while (flag_3) {
							system("cls");
							cout << "Введите степень" << endl;
							cin >> ind;
							cout << "Введите значение" << endl;
							cin >> number;
							try {
								p2.Set(number, ind);
							}
							catch (const char* msg) {
								cout << msg << endl;
							}
							cout << "Продолжить 1" << endl;
							cout << "Выход 0" << endl;
							cin >> flag_3;
						}
						break;
					}
					}
					break;
				}
				case 2:
				{
					system("cls");
					cout << p1 << endl;
					cout << p2 << endl;
					system("pause");
					break;
				}
				case 3:
				{
					system("cls");
					size_t flag_2;
					cout << "У какого многочлена?" << endl;
					cin >> flag_2;
					switch (flag_2)
					{
					case 1:
					{
						size_t ind;
						cout << "Введите степень" << endl;
						cin >> ind;
						try {
							cout << p1[ind] << endl;
						}
						catch (const char* msg) {
							cout << msg << endl;
						}
						system("pause");
						break;
					}
					case 2:
					{
						size_t ind;
						cout << "Введите степень" << endl;
						cin >> ind;
						try {
							cout << p2[ind] << endl;
						}
						catch (const char* msg) {
							cout << msg << endl;
						}
						system("pause");
						break;
					}
					}
					break;
				}
				case 4:
				{
					system("cls");
					Polynomial <float> p3 = p1 + p2;
					cout << p1 << " + " << p2 << " = " << p3 << endl;
					system("pause");
					break;
				}
				case 5:
				{
					system("cls");
					size_t flag_2;
					cout << "Из какого многочлена вычесть?" << endl;
					cin >> flag_2;
					switch (flag_2)
					{
					case 1:
					{
						Polynomial <float> p3 = p1 - p2;
						cout << p1 << " - " << "(" << p2 << ")" << " = " << p3 << endl;
						system("pause");
						break;
					}
					case 2:
					{
						Polynomial <float> p3 = p2 - p1;
						cout << p2 << " - " << "(" << p1 << ")" << " = " << p3 << endl;
						system("pause");
						break;
					}
					}
					break;
				}
				case 6:
				{
					system("cls");
					size_t flag_2;
					float value;
					cout << "На что умножить?" << endl;
					cin >> value;
					cout << "1. Изменить многочлен" << endl;
					cout << "2. Создать новый многочлен" << endl;
					cin >> flag_2;
					switch (flag_2)
					{
					case 1:
					{
						cout << p1 << " * " << value << " = ";
						p1 *= value;
						cout << p1 << endl;
						cout << p2 << " * " << value << " = ";
						p2 *= value;
						cout << p2 << endl;
						break;
					}
					case 2:
					{
						Polynomial <float> p3 = p1 * value;
						Polynomial <float> p4 = value * p2;
						cout << p1 << " * " << value << " = " << p3 << endl;
						cout << value << " * " << p2 << " = " << p4 << endl;
						break;
					}
					}
					system("pause");
					break;
				}
				case 7:
				{
					system("cls");
					float value;
					cout << "Введите Х" << endl;
					cin >> value;
					cout << p1 << endl;
					cout << p2 << endl;
					cout << p1.calculation(value) << endl;
					cout << p2.calculation(value) << endl;
					system("pause");
					break;
				}
				case 8:
				{
					system("cls");
					cout << p1 << " = 0" << endl;
					cout << p2 << " = 0" << endl;
					double x1, x2, x3, Q;
					try {
						p1.formula_Kardano(x1, x2, x3, Q);
						if (Q > 0) {
							cout << "x1 = " << x1 << endl;
						}
						if (abs(Q) < 0.01) {
							cout << "x1 = " << x1 << endl;
							cout << "x2 = " << x2 << endl;
							cout << "x3 = " << x3 << endl;
						}
						if (Q < 0) {
							cout << "x1 = " << x1 << endl;
							cout << "x2 = " << x2 << endl;
							cout << "x3 = " << x3 << endl;
						}
					}
					catch (const char* msg) {
						cout << msg << endl;
					}
					try {
						p2.formula_Kardano(x1, x2, x3, Q);
						if (Q > 0) {
							cout << "x1 = " << x1 << endl;
						}
						if (abs(Q) < 0.01) {
							cout << "x1 = " << x1 << endl;
							cout << "x2 = " << x2 << endl;
							cout << "x3 = " << x3 << endl;
						}
						if (Q < 0) {
							cout << "x1 = " << x1 << endl;
							cout << "x2 = " << x2 << endl;
							cout << "x3 = " << x3 << endl;
						}
					}
					catch (const char* msg) {
						cout << msg << endl;
					}
					system("pause");
					break;
				}
				case 9: {
					system("cls");
					bool f;
					f = p1 == p2;
					cout << f << endl;
					system("pause");
					break;
				}
				case 0:
				{
					fg = false;
				}
				}
			}
			break;
		}
		case 3:
		{
			bool fg = true;
			system("cls");
			size_t number_1;
			size_t number_2;
			cout << "Введите степень первого многочлена" << endl;
			cin >> number_1;
			cout << "Введите степень второго многочлена" << endl;
			cin >> number_2;
			Polynomial<double> p1(number_1);
			Polynomial <double>p2(number_2);
			while (fg) {
				system("cls");
				menu_2();
				size_t flag_2_3;
				cin >> flag_2_3;
				switch (flag_2_3)
				{
				case 1:
				{
					system("cls");
					size_t flag_2;
					cout << "У какого многочлена?" << endl;
					cin >> flag_2;
					switch (flag_2)
					{
					case 1:
					{
						size_t ind;
						double number;
						size_t flag_3 = 1;
						while (flag_3) {
							system("cls");
							cout << "Введите степень" << endl;
							cin >> ind;
							cout << "Введите значение" << endl;
							cin >> number;
							try {
								p1.Set(number, ind);
							}
							catch (const char* msg) {
								cout << msg << endl;
							}
							cout << "Продолжить 1" << endl;
							cout << "Выход 0" << endl;
							cin >> flag_3;
						}
						break;
					}
					case 2:
					{
						size_t ind;
						double number;
						size_t flag_3 = 1;
						while (flag_3) {
							system("cls");
							cout << "Введите степень" << endl;
							cin >> ind;
							cout << "Введите значение" << endl;
							cin >> number;
							try {
								p2.Set(number, ind);
							}
							catch (const char* msg) {
								cout << msg << endl;
							}
							cout << "Продолжить 1" << endl;
							cout << "Выход 0" << endl;
							cin >> flag_3;
						}
						break;
					}
					}
					break;
				}
				case 2:
				{
					system("cls");
					cout << p1 << endl;
					cout << p2 << endl;
					system("pause");
					break;
				}
				case 3:
				{
					system("cls");
					size_t flag_2;
					cout << "У какого многочлена?" << endl;
					cin >> flag_2;
					switch (flag_2)
					{
					case 1:
					{
						size_t ind;
						cout << "Введите степень" << endl;
						cin >> ind;
						try {
							cout << p1[ind] << endl;
						}
						catch (const char* msg) {
							cout << msg << endl;
						}
						system("pause");
						break;
					}
					case 2:
					{
						size_t ind;
						cout << "Введите степень" << endl;
						cin >> ind;
						try {
							cout << p2[ind] << endl;
						}
						catch (const char* msg) {
							cout << msg << endl;
						}
						system("pause");
						break;
					}
					}
					break;
				}
				case 4:
				{
					system("cls");
					Polynomial <double> p3 = p1 + p2;
					cout << p1 << " + " << p2 << " = " << p3 << endl;
					system("pause");
					break;
				}
				case 5:
				{
					system("cls");
					size_t flag_2;
					cout << "Из какого многочлена вычесть?" << endl;
					cin >> flag_2;
					switch (flag_2)
					{
					case 1:
					{
						Polynomial <double> p3 = p1 - p2;
						cout << p1 << " - " << "(" << p2 << ")" << " = " << p3 << endl;
						system("pause");
						break;
					}
					case 2:
					{
						Polynomial <double> p3 = p2 - p1;
						cout << p2 << " - " << "(" << p1 << ")" << " = " << p3 << endl;
						system("pause");
						break;
					}
					}
					break;
				}
				case 6:
				{
					system("cls");
					size_t flag_2;
					double value;
					cout << "На что умножить?" << endl;
					cin >> value;
					cout << "1. Изменить многочлен" << endl;
					cout << "2. Создать новый многочлен" << endl;
					cin >> flag_2;
					switch (flag_2)
					{
					case 1:
					{
						cout << p1 << " * " << value << " = ";
						p1 *= value;
						cout << p1 << endl;
						cout << p2 << " * " << value << " = ";
						p2 *= value;
						cout << p2 << endl;
						break;
					}
					case 2:
					{
						Polynomial <double> p3 = p1 * value;
						Polynomial <double> p4 = value * p2;
						cout << p1 << " * " << value << " = " << p3 << endl;
						cout << value << " * " << p2 << " = " << p4 << endl;
						break;
					}
					}
					system("pause");
					break;
				}
				case 7:
				{
					system("cls");
					double value;
					cout << "Введите Х" << endl;
					cin >> value;
					cout << p1 << endl;
					cout << p2 << endl;
					cout << p1.calculation(value) << endl;
					cout << p2.calculation(value) << endl;
					system("pause");
					break;
				}
				case 8:
				{
					system("cls");
					cout << p1 << " = 0" << endl;
					cout << p2 << " = 0" << endl;
					double x1, x2, x3, Q;
					try {
						p1.formula_Kardano(x1, x2, x3, Q);
						if (Q > 0) {
							cout << "x1 = " << x1 << endl;
						}
						if (abs(Q) < 0.01) {
							cout << "x1 = " << x1 << endl;
							cout << "x2 = " << x2 << endl;
							cout << "x3 = " << x3 << endl;
						}
						if (Q < 0) {
							cout << "x1 = " << x1 << endl;
							cout << "x2 = " << x2 << endl;
							cout << "x3 = " << x3 << endl;
						}
					}
					catch (const char* msg) {
						cout << msg << endl;
					}
					try {
						p2.formula_Kardano(x1, x2, x3, Q);
						if (Q > 0) {
							cout << "x1 = " << x1 << endl;
						}
						if (abs(Q) < 0.01) {
							cout << "x1 = " << x1 << endl;
							cout << "x2 = " << x2 << endl;
							cout << "x3 = " << x3 << endl;
						}
						if (Q < 0) {
							cout << "x1 = " << x1 << endl;
							cout << "x2 = " << x2 << endl;
							cout << "x3 = " << x3 << endl;
						}
					}
					catch (const char* msg) {
						cout << msg << endl;
					}
					system("pause");
					break;
				}
				case 9: {
					system("cls");
					bool f;
					f = p1 == p2;
					cout << f << endl;
					system("pause");
					break;
				}
				case 0:
				{
					fg = false;
				}
				}
			}
			break;
		}
		case 4:
		{
			bool fg = true;
			system("cls");
			size_t number_1;
			size_t number_2;
			cout << "Введите степень первого многочлена" << endl;
			cin >> number_1;
			cout << "Введите степень второго многочлена" << endl;
			cin >> number_2;
			Polynomial <complex<float>> p1(number_1);
			Polynomial <complex<float>> p2(number_2);
			while (fg) {
				system("cls");
				menu_2();
				size_t flag_2;
				cin >> flag_2;
				switch (flag_2)
				{
				case 1:
				{
					system("cls");
					size_t flag_5;
					cout << "У какого многочлена?" << endl;
					cin >> flag_5;
					switch (flag_5)
					{
					case 1:
					{
						size_t ind;
						size_t flag = 1;
						while (flag) {
							system("cls");
							cout << "Введите степень" << endl;
							cin >> ind;
							float real;
							cout << "Введите действительную часть" << endl;
							cin >> real;
							float im;
							cout << "Введите мнимую часть" << endl;
							cin >> im;
							complex<float>x(real, im);
							try {
								p1.Set(x, ind);
							}
							catch (const char* msg) {
								cout << msg << endl;
							}
							cout << "Продолжить 1" << endl;
							cout << "Выход 0" << endl;
							cin >> flag;
						}
						break;
					}
					case 2:
					{
						size_t ind;
						size_t flag = 1;
						while (flag) {
							system("cls");
							cout << "Введите степень" << endl;
							cin >> ind;
							float real;
							cout << "Введите действительную часть" << endl;
							cin >> real;
							float im;
							cout << "Введите мнимую часть" << endl;
							cin >> im;
							complex<float>x(real, im);
							try {
								p2.Set(x, ind);
							}
							catch (const char* msg) {
								cout << msg << endl;
							}
							cout << "Продолжить 1" << endl;
							cout << "Выход 0" << endl;
							cin >> flag;
						}
						break;
					}
					}
					break;
				}
				case 2:
				{
					system("cls");
					cout << p1 << endl;
					cout << p2 << endl;
					system("pause");
					break;
				}
				case 3:
				{
					system("cls");
					size_t flag_0;
					cout << "У какого многочлена?" << endl;
					cin >> flag_0;
					switch (flag_0)
					{
					case 1:
					{
						size_t ind;
						cout << "Введите степень" << endl;
						cin >> ind;
						try {
							cout << p1[ind] << endl;
						}
						catch (const char* msg) {
							cout << msg << endl;
						}
						system("pause");
						break;
					}
					case 2:
					{
						size_t ind;
						cout << "Введите степень" << endl;
						cin >> ind;
						try {
							cout << p2[ind] << endl;
						}
						catch (const char* msg) {
							cout << msg << endl;
						}
						system("pause");
						break;
					}
					}
					break;
				}
				case 4:
				{
					system("cls");
					Polynomial<complex<float>> p3 = p1 + p2;
					cout << p1 << " + " << p2 << " = " << p3 << endl;
					system("pause");
					break;
				}
				case 5:
				{
					system("cls");
					size_t flag_0;
					cout << "Из какого многочлена вычесть?" << endl;
					cin >> flag_0;
					switch (flag_0)
					{
					case 1:
					{
						Polynomial<complex <float>> p3 = p1 - p2;
						cout << p1 << " - " << "(" << p2 << ")" << " = " << p3 << endl;
						system("pause");
						break;
					}
					case 2:
					{
						Polynomial<complex <float>> p3 = p2 - p1;
						cout << p2 << " - " << "(" << p1 << ")" << " = " << p3 << endl;
						system("pause");
						break;
					}
					}
					break;
				}
				case 6:
				{
					system("cls");
					size_t flag_0;
					float value;
					cout << "На что умножить?" << endl;
					cin >> value;
					cout << "1. Изменить многочлен" << endl;
					cout << "2. Создать новый многочлен" << endl;
					cin >> flag_0;
					switch (flag_0)
					{
					case 1:
					{
						cout << p1 << " * " << value << " = ";
						p1 *= value;
						cout << p1 << endl;
						cout << p2 << " * " << value << " = ";
						p2 *= value;
						cout << p2 << endl;
						break;
					}
					case 2:
					{
						Polynomial<complex <float>> p3 = p1 * value;
						Polynomial<complex <float>> p4 = value * p2;
						cout << p1 << " * " << value << " = " << p3 << endl;
						cout << value << " * " << p2 << " = " << p4 << endl;
						break;
					}
					}
					system("pause");
					break;
				}
				case 7:
				{
					system("cls");
					float real;
					cout << "Введите действительную часть" << endl;
					cin >> real;
					float im;
					cout << "Введите мнимую часть" << endl;
					cin >> im;
					complex<float> value(real, im);
					cout << p1 << endl;
					cout << p2 << endl;
					cout << p1.calculation(value) << endl;
					cout << p2.calculation(value) << endl;
					system("pause");
					break;
				}
				case 8:
				{
					system("cls");
					cout << p1 << " = 0" << endl;
					cout << p2 << " = 0" << endl;
					try {
						float x1, x2, x3;
						double Q;
						bool flag;
						p1.formula_Kardano(x1, x2, x3, Q, flag);
						cout << "Первое уравнение:" << endl;
						if (flag) {
							cout << "x1 = " << x1 << endl;
							cout << "x2 = " << x2 << endl;
							cout << "x3 = " << x3 << endl;
						}
						else {
							cout << "Действительных корней нет" << endl;
						}
					}
					catch (const char* msg) {
						cout << msg << endl;
					}
					try {
						float x1, x2, x3;
						double Q;
						bool flag;
						p2.formula_Kardano(x1, x2, x3, Q, flag);
						cout << "Второе уравнение:" << endl;
						if (flag) {
							cout << "x1 = " << x1 << endl;
							cout << "x2 = " << x2 << endl;
							cout << "x3 = " << x3 << endl;
						}
						else {
							cout << "Действительных корней нет" << endl;
						}
					}
					catch (const char* msg) {
						cout << msg << endl;
					}
					system("pause");
					break;

				}
				case 9: {
					system("cls");
					bool f;
					f = p1 == p2;
					cout << f << endl;
					system("pause");
					break;
				}

				case 0:
				{
					fg = false;
				}
				}
			}
			break;
		}
		case 5:
		{
			bool fg = true;
			system("cls");
			size_t number_1;
			size_t number_2;
			cout << "Введите степень первого многочлена" << endl;
			cin >> number_1;
			cout << "Введите степень второго многочлена" << endl;
			cin >> number_2;
			Polynomial <complex<double>> p1(number_1);
			Polynomial <complex<double>> p2(number_2);
			while (fg) {
				system("cls");
				menu_2();
				size_t flag_2;
				cin >> flag_2;
				switch (flag_2)
				{
				case 1:
				{
					system("cls");
					size_t flag_0;
					cout << "У какого многочлена?" << endl;
					cin >> flag_0;
					switch (flag_0)
					{
					case 1:
					{
						size_t ind;
						size_t flag = 1;
						while (flag) {
							system("cls");
							cout << "Введите степень" << endl;
							cin >> ind;
							double real;
							cout << "Введите действительную часть" << endl;
							cin >> real;
							double im;
							cout << "Введите мнимую часть" << endl;
							cin >> im;
							complex<double>x(real, im);
							try {
								p1.Set(x, ind);
							}
							catch (const char* msg) {
								cout << msg << endl;
							}
							cout << "Продолжить 1" << endl;
							cout << "Выход 0" << endl;
							cin >> flag;
						}
						break;
					}
					case 2:
					{
						size_t ind;
						size_t flag = 1;
						while (flag) {
							system("cls");
							cout << "Введите степень" << endl;
							cin >> ind;
							float real;
							cout << "Введите действительную часть" << endl;
							cin >> real;
							float im;
							cout << "Введите мнимую часть" << endl;
							cin >> im;
							complex<double>x(real, im);
							try {
								p2.Set(x, ind);
							}
							catch (const char* msg) {
								cout << msg << endl;
							}
							cout << "Продолжить 1" << endl;
							cout << "Выход 0" << endl;
							cin >> flag;
						}
						break;
					}
					}
					break;
				}
				case 2:
				{
					system("cls");
					cout << p1 << endl;
					cout << p2 << endl;
					system("pause");
					break;
				}
				case 3:
				{
					system("cls");
					size_t flag_0;
					cout << "У какого многочлена?" << endl;
					cin >> flag_0;
					switch (flag_0)
					{
					case 1:
					{
						size_t ind;
						cout << "Введите степень" << endl;
						cin >> ind;
						try {
							cout << p1[ind] << endl;
						}
						catch (const char* msg) {
							cout << msg << endl;
						}
						system("pause");
						break;
					}
					case 2:
					{
						size_t ind;
						cout << "Введите степень" << endl;
						cin >> ind;
						try {
							cout << p2[ind] << endl;
						}
						catch (const char* msg) {
							cout << msg << endl;
						}
						system("pause");
						break;
					}
					}
					break;
				}
				case 4:
				{
					system("cls");
					Polynomial<complex<double>> p3 = p1 + p2;
					cout << p1 << " + " << p2 << " = " << p3 << endl;
					system("pause");
					break;
				}
				case 5:
				{
					system("cls");
					size_t flag_0;
					cout << "Из какого многочлена вычесть?" << endl;
					cin >> flag_0;
					switch (flag_0)
					{
					case 1:
					{
						Polynomial<complex <double>> p3 = p1 - p2;
						cout << p1 << " - " << "(" << p2 << ")" << " = " << p3 << endl;
						system("pause");
						break;
					}
					case 2:
					{
						Polynomial<complex <double>> p3 = p2 - p1;
						cout << p2 << " - " << "(" << p1 << ")" << " = " << p3 << endl;
						system("pause");
						break;
					}
					}
					break;
				}
				case 6:
				{
					system("cls");
					size_t flag_0;
					cout << "1. Изменить многочлен" << endl;
					cout << "2. Создать новый многочлен" << endl;
					cin >> flag_0;
					switch (flag_0)
					{
					case 1:
					{
						double real;
						cout << "Введите действительную часть" << endl;
						cin >> real;
						double im;
						cout << "Введите мнимую часть" << endl;
						cin >> im;
						complex<double> value(real, im);
						cout << p1 << " * " << value << " = ";
						p1 *= value;
						cout << p1 << endl;
						cout << p2 << " * " << value << " = ";
						p2 *= value;
						cout << p2 << endl;
						break;
					}
					case 2:
					{
						double real;
						cout << "Введите действительную часть" << endl;
						cin >> real;
						double im;
						cout << "Введите мнимую часть" << endl;
						cin >> im;
						complex<double> value(real, im);
						Polynomial<complex <double>> p3 = p1 * value;
						Polynomial<complex <double>> p4 = value * p2;
						cout << p1 << " * " << value << " = " << p3 << endl;
						cout << value << " * " << p2 << " = " << p4 << endl;
						break;
					}
					}
					system("pause");
					break;
				}
				case 7:
				{
					system("cls");
					double real;
					cout << "Введите действительную часть" << endl;
					cin >> real;
					double im;
					cout << "Введите мнимую часть" << endl;
					cin >> im;
					complex<double> value(real, im);
					cout << p1 << endl;
					cout << p2 << endl;
					cout << p1.calculation(value) << endl;
					cout << p2.calculation(value) << endl;
					system("pause");
					break;
				}
				case 8:
				{
					system("cls");
					cout << p1 << " = 0" << endl;
					cout << p2 << " = 0" << endl;
					try {
						double x1, x2, x3, Q;
						bool flag;
						p1.formula_Kardano(x1, x2, x3, Q, flag);
						cout << "Первое уравнение:" << endl;
						if (flag) {
							cout << "x1 = " << x1 << endl;
							cout << "x2 = " << x2 << endl;
							cout << "x3 = " << x3 << endl;
						}
						else {
							cout << "Действительных корней нет" << endl;
						}
					}
					catch (const char* msg) {
						cout << msg << endl;
					}
					try {
						double x1, x2, x3, Q;
						bool flag;
						cout << "Второе уравнение:" << endl;
						p2.formula_Kardano(x1, x2, x3, Q, flag);
						if (flag) {
							cout << "x1 = " << x1 << endl;
							cout << "x2 = " << x2 << endl;
							cout << "x3 = " << x3 << endl;
						}
						else {
							cout << "Действительных корней нет" << endl;
						}
					}
					catch (const char* msg) {
						cout << msg << endl;
					}
					system("pause");
					break;

				}
				case 9: {
					system("cls");
					bool f;
					f = p1 == p2;
					cout << f << endl;
					system("pause");
					break;
				}

				case 0:
				{
					fg = false;
				}
				}
			}
			break;
		}
		case 0:
		{
			return 0;
		}
		}
	}
}

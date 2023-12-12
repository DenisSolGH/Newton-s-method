// Differentiation.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//


#include <iostream>
#include <cmath>

using namespace std;

double f1(double x, double y, double z)
{

    //double f1 = x + 3 * log10(x) - y * y;
    //double f2 = 2 * x * x - x * y - 5 * x + 1;
    //  (2*x+7)^(5*x+2)/((2*x-4)^(5*x+2))
    //  (1 / x) + (1 / y) - (3./ 2.)
    // pow(x,2) + sqrt(3)*pow(y,3)+z+56.7

    return  pow(2.718721872, x) - 2.1/pow(y,3) -5.54; //функция нескольких переменных

}

double f2(double x, double y, double z)
{
    // pow(log10(x), 2) +3*pow(y,2) + sqrt(2)*z - 27

    return x * y + pow(y, 2) - z - 7.62; //функция нескольких переменных

}

double f3(double x, double y, double z)
{
   // x + 8.1 * y + pow(z, 2) + 13.05

    return pow((x-z), 3) +y -89.1; //функция нескольких переменных

}

double diff1x(double x, double y, double z, double eps = 0.0001)
{
    return (f1(x + eps, y, z) - f1(x, y, z)) / eps;
}
double diff2x(double x, double y, double z, double eps = 0.0001)
{
    return (f2(x + eps, y, z) - f2(x, y, z)) / eps;
}
double diff3x(double x, double y, double z, double eps = 0.0001)
{
    return (f3(x + eps, y, z) - f3(x, y, z)) / eps;
}

double diff1y(double x, double y, double z, double eps = 0.0001)
{
    return (f1(x, y + eps, z) - f1(x, y, z)) / eps;
}
double diff2y(double x, double y, double z, double eps = 0.0001)
{
    return (f2(x, y + eps, z) - f2(x, y, z)) / eps;
}
double diff3y(double x, double y, double z, double eps = 0.0001)
{
    return (f3(x, y + eps, z) - f3(x, y, z)) / eps;
}

double diff1z(double x, double y, double z, double eps = 0.0001)
{
    return (f1(x, y, z + eps) - f1(x, y, z)) / eps;
}
double diff2z(double x, double y, double z, double eps = 0.0001)
{
    return (f2(x, y, z + eps) - f2(x, y, z)) / eps;
}
double diff3z(double x, double y, double z, double eps = 0.0001)
{
    return (f3(x, y, z + eps) - f3(x, y, z)) / eps;
}

/*void mulMatr(double matrixJacobiReversed[3][3], double funcArr[3], double multMatr[3])
{
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            multMatr[i] += matrixJacobiReversed[i][j] * funcArr[j];
        }
    }
}
*/

int main()
{
    setlocale(LC_ALL, "Rus");

    int n = 50; // количество итераций

    double x = 4;  // Нулевые приближения
    double y = 5;
    double z = 6;

    const int t = 3;

    if ( t == 2)
    { 
        // Случай двух уравнений с двумя неизвестными
        for (int i = 0; i < n; i++)
        {
            double xvect[50][2]; // значения x и y на разных приближениях 

            xvect[0][i] = x; // Запись нулевых приближений в вектор-столдец x
            xvect[1][i] = y;

            //Вывод вектора x на i-м приближении
            cout << "x" << i << " = " << xvect[0][i] << "\t\t\t";
            cout << "y" << i << " = " << xvect[1][i] << endl;

            // Вывод вектор-функции от вектора x на i-м приближении
            cout << "f1(x" << i << ") = " << f1(x, y, z) << "\t\t\t";
            cout << "f2(y" << i << ") = " << f2(x, y, z) << endl;

            // Составление матрицы Якоби второго порядка
            double matrixJacobi[2][2];
            matrixJacobi[0][0] = diff1x(x, y, z);
            matrixJacobi[0][1] = diff1y(x, y, z);
            matrixJacobi[1][0] = diff2x(x, y, z);
            matrixJacobi[1][1] = diff2y(x, y, z);

            cout << "Матрица Якоби: " << endl << diff1x(x, y, z) << "\t" << diff1y(x, y, z) << endl
                << diff2x(x, y, z) << "\t" << diff2y(x, y, z) << endl;

            // Вычисление определителя матрицы Якоби второго порядка
            double delta = diff1x(x, y, z) * diff2y(x, y, z) - diff1y(x, y, z) * diff2x(x, y, z);

            cout << "Определитель delta = " << delta << endl;

            // Вычисление транспонированной матрицы
            double matrixJacobiTransp[2][2];
            matrixJacobiTransp[0][0] = matrixJacobi[0][0];
            matrixJacobiTransp[0][1] = matrixJacobi[1][0];
            matrixJacobiTransp[1][0] = matrixJacobi[0][1];
            matrixJacobiTransp[1][1] = matrixJacobi[1][1];

            cout << "Транспонированная матрица: " << endl
                << matrixJacobiTransp[0][0] << "\t"
                << matrixJacobiTransp[0][1] << endl
                << matrixJacobiTransp[1][0] << "\t"
                << matrixJacobiTransp[1][1] << endl;

            // Вычисление матрицы алгебраических дополнений
            double matrixJacobiTranspAlg[2][2];
            matrixJacobiTranspAlg[0][0] = matrixJacobiTransp[1][1];
            matrixJacobiTranspAlg[0][1] = -matrixJacobiTransp[1][0];
            matrixJacobiTranspAlg[1][0] = -matrixJacobiTransp[0][1];
            matrixJacobiTranspAlg[1][1] = matrixJacobiTransp[0][0];

            cout << "Матрица алгебраических дополнений: " << endl
                << matrixJacobiTranspAlg[0][0] << "\t"
                << matrixJacobiTranspAlg[0][1] << endl
                << matrixJacobiTranspAlg[1][0] << "\t"
                << matrixJacobiTranspAlg[1][1] << endl;

            // Вычисление обратной матрицы для матрицы Якоби
            double matrixJacobiReversed[2][2];
            matrixJacobiReversed[0][0] = matrixJacobiTranspAlg[0][0] / delta;
            matrixJacobiReversed[0][1] = matrixJacobiTranspAlg[0][1] / delta;
            matrixJacobiReversed[1][0] = matrixJacobiTranspAlg[1][0] / delta;
            matrixJacobiReversed[1][1] = matrixJacobiTranspAlg[1][1] / delta;

            cout << "Обратная матрица Якоби: " << endl
                << matrixJacobiReversed[0][0] << "\t"
                << matrixJacobiReversed[0][1] << endl
                << matrixJacobiReversed[1][0] << "\t"
                << matrixJacobiReversed[1][1] << endl;

            // Перемножение матрицы обратной матрице Якоби и матрицы, содержащей функции
            double multMatr[2];
            multMatr[0] = matrixJacobiReversed[0][0] * f1(x, y, z) + matrixJacobiReversed[0][1] * f2(x, y, z);
            multMatr[1] = matrixJacobiReversed[1][0] * f1(x, y, z) + matrixJacobiReversed[1][1] * f2(x, y, z);

            //cout << multMatr[0] << "\t" << multMatr[1];
            // Вычисление значений вектора x на следующем приближении
            x = x - multMatr[0];
            y = y - multMatr[1];



            cout << endl << endl << endl << endl;
        }
    }
    else if (t == 3)
    {
        // Случай двух уравнений с тремя неизвестными
        for (int i = 0; i < n; i++)
        {
            double xvect[50][3]; // значения x и y и z на разных приближениях 

            xvect[0][i] = x; // Запись нулевых приближений в вектор-столдец x
            xvect[1][i] = y;
            xvect[2][i] = z;

            //Вывод вектора x на i-м приближении
            cout << "x" << i << " = " << xvect[0][i] << "\t\t\t";
            cout << "y" << i << " = " << xvect[1][i] << "\t\t\t";
            cout << "z" << i << " = " << xvect[2][i] << endl;

            // Вывод вектор-функции от вектора x на i-м приближении
            cout << "f1(x" << i << ") = " << f1(x, y, z) << "\t\t\t";
            cout << "f2(y" << i << ") = " << f2(x, y, z) << "\t\t\t";
            cout << "f3(y" << i << ") = " << f3(x, y, z) << endl;


            // Составление матрицы Якоби второго порядка
            double matrixJacobi[3][3];
            matrixJacobi[0][0] = diff1x(x, y, z);
            matrixJacobi[0][1] = diff1y(x, y, z);
            matrixJacobi[0][2] = diff1z(x, y, z);
            matrixJacobi[1][0] = diff2x(x, y, z);
            matrixJacobi[1][1] = diff2y(x, y, z);
            matrixJacobi[1][2] = diff2z(x, y, z);
            matrixJacobi[2][0] = diff3x(x, y, z);
            matrixJacobi[2][1] = diff3y(x, y, z);
            matrixJacobi[2][2] = diff3z(x, y, z);

            cout << "Матрица Якоби: " << endl
                << diff1x(x, y, z) << "\t\t"
                << diff1y(x, y, z) << "\t\t" 
                << diff1z(x, y, z) << endl
                << diff2x(x, y, z) << "\t\t"
                << diff2y(x, y, z) << "\t\t" 
                << diff2z(x, y, z) << endl
                << diff3x(x, y, z) << "\t\t"
                << diff3y(x, y, z) << "\t\t" 
                << diff3z(x, y, z) << endl;
        
            // Вычисление определителя матрицы Якоби третьего порядка
            double delta = 0.;

            for (int j = 0; j < 3; j++)
            {
                for (int k = 0; k<3; k++)
                {
                    delta = matrixJacobi[0][0] * matrixJacobi[1][1] * matrixJacobi[2][2] +
                        matrixJacobi[0][1] * matrixJacobi[1][2] * matrixJacobi[2][0] +
                        matrixJacobi[1][0] * matrixJacobi[2][1] * matrixJacobi[0][2] -
                        matrixJacobi[0][2] * matrixJacobi[1][1] * matrixJacobi[2][0] -
                        matrixJacobi[0][1] * matrixJacobi[1][0] * matrixJacobi[2][2] -
                        matrixJacobi[0][0] * matrixJacobi[1][2] * matrixJacobi[2][1];
                }
            }

            cout << "Определитель delta = " << delta << endl;

            // Вычисление транспонированной матрицы
            double matrixJacobiTransp[3][3];

            for (int j = 0; j < 3; j++)
            {
                for (int k = 0; k < 3; k++)
                {
                    matrixJacobiTransp[j][k] = matrixJacobi[k][j];
                }
            }

            cout << "Транспонированная матрица: " << endl;
            for (int j = 0; j < 3; j++)
            {
                for (int k = 0; k < 3; k++)
                {
                    cout << matrixJacobiTransp[j][k] << "\t";
                }
                cout << endl;
            }

            // Вычисление матрицы алгебраических дополнений
            double matrixJacobiTranspAlg[3][3];


            double deltatemp[3][3];
            deltatemp[0][0] = matrixJacobiTransp[1][1] * matrixJacobiTransp[2][2] -
                matrixJacobiTransp[2][1] * matrixJacobiTransp[1][2];
            deltatemp[0][1] = matrixJacobiTransp[1][0] * matrixJacobiTransp[2][2] -
                matrixJacobiTransp[2][0] * matrixJacobiTransp[1][2];
            deltatemp[0][2] = matrixJacobiTransp[1][0] * matrixJacobiTransp[2][1] -
                matrixJacobiTransp[2][0] * matrixJacobiTransp[1][1];
            deltatemp[1][0] = matrixJacobiTransp[0][1] * matrixJacobiTransp[2][2] -
                matrixJacobiTransp[2][1] * matrixJacobiTransp[0][2];
            deltatemp[1][1] = matrixJacobiTransp[0][0] * matrixJacobiTransp[2][2] -
                matrixJacobiTransp[2][0] * matrixJacobiTransp[0][2];
            deltatemp[1][2] = matrixJacobiTransp[0][0] * matrixJacobiTransp[2][1] -
                matrixJacobiTransp[2][0] * matrixJacobiTransp[0][1];
            deltatemp[2][0] = matrixJacobiTransp[0][1] * matrixJacobiTransp[1][2] -
                matrixJacobiTransp[1][1] * matrixJacobiTransp[0][2];
            deltatemp[2][1] = matrixJacobiTransp[0][0] * matrixJacobiTransp[1][2] -
                matrixJacobiTransp[1][0] * matrixJacobiTransp[0][2];
            deltatemp[2][2] = matrixJacobiTransp[0][0] * matrixJacobiTransp[1][1] -
                matrixJacobiTransp[1][0] * matrixJacobiTransp[0][1];

            for (int j = 0; j < 3; j++)
            {
                for (int k = 0; k < 3; k++)
                {
                    matrixJacobiTranspAlg[j][k] = pow(-1,(j+k))* deltatemp[j][k];
                }
            }

            cout << "Матрица алгебраических дополнений: " << endl;
            for (int j = 0; j < t; j++)
            {
                for (int k = 0; k < t; k++)
                {
                    cout << matrixJacobiTranspAlg[j][k] << "\t";
                }
                cout << endl;
            }

            // Вычисление обратной матрицы для матрицы Якоби
            double matrixJacobiReversed[3][3];

            for (int j = 0; j < 3; j++)
            {
                for (int k = 0; k < t; k++)
                {
                    matrixJacobiReversed[j][k] = matrixJacobiTranspAlg[j][k] / delta;
                }
            }

            cout << "Обратная матрица Якоби: " << endl;
            for (int j = 0; j < 3; j++)
            {
                for (int k = 0; k < 3; k++)
                {
                    cout << matrixJacobiReversed[j][k] << "\t";
                }
                cout << endl;
            }
            // ------------------------------------------------------------
            
            // Перемножение матрицы обратной матрице Якоби и матрицы, содержащей функции
            double multMatr[3];
            multMatr[0] = matrixJacobiReversed[0][0] * f1(x, y, z) +
                matrixJacobiReversed[0][1] * f2(x, y, z) +
                matrixJacobiReversed[0][2] * f3(x, y, z);
            multMatr[1] = matrixJacobiReversed[1][0] * f1(x, y, z) +
                matrixJacobiReversed[1][1] * f2(x, y, z) +
                matrixJacobiReversed[1][2] * f3(x, y, z);
            multMatr[2] = matrixJacobiReversed[2][0] * f1(x, y, z) +
                matrixJacobiReversed[2][1] * f2(x, y, z) +
                matrixJacobiReversed[2][2] * f3(x, y, z);


            cout << "После перемножения матрицы функций и обратной матрицы: " << endl;
            for (int i = 0; i < 3; i++)
            {
                cout << multMatr[i] << endl;
            }

            // Вычисление значений вектора x на следующем приближении
            x = x - multMatr[0];
            y = y - multMatr[1];
            z = z - multMatr[2];

          
            cout << endl << endl << endl << endl;
        }
    }
        
    

           

   

}

// Запуск программы: CTRL+F5 или меню "Отладка" > "Запуск без отладки"
// Отладка программы: F5 или меню "Отладка" > "Запустить отладку"

// Советы по началу работы 
//   1. В окне обозревателя решений можно добавлять файлы и управлять ими.
//   2. В окне Team Explorer можно подключиться к системе управления версиями.
//   3. В окне "Выходные данные" можно просматривать выходные данные сборки и другие сообщения.
//   4. В окне "Список ошибок" можно просматривать ошибки.
//   5. Последовательно выберите пункты меню "Проект" > "Добавить новый элемент", чтобы создать файлы кода, или "Проект" > "Добавить существующий элемент", чтобы добавить в проект существующие файлы кода.
//   6. Чтобы снова открыть этот проект позже, выберите пункты меню "Файл" > "Открыть" > "Проект" и выберите SLN-файл.

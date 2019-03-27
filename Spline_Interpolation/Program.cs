using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Spline_Interpolation
{
    static class MatrixOperations
    {
        public static double[] GausMethod(double[,] matrix, double[] b)
        {
            int width = matrix.GetLength(1);
            int height = matrix.GetLength(0);
            for (int j = 0; j < height; j++)
            {
                double cur_diag = matrix[j, j];
                if (cur_diag != 0)
                {
                    for (int i = 0; i < width; i++)
                        matrix[j, i] /= cur_diag;
                    b[j] /= cur_diag;
                }
                for (int k = 0; k < height; k++)
                {
                    if (k != j)
                    {
                        double cur = matrix[k, j];
                        for (int i = 0; i < width; i++)
                            matrix[k, i] -= matrix[j, i] * cur;
                        b[k] -= b[j] * cur;
                    }
                }
            }
            for (int j = height - 1; j > -1; j--)
            {
                double cur_diag = matrix[j, j];
                if (cur_diag != 0)
                {
                    for (int i = 0; i < width; i++)
                        matrix[j, i] /= cur_diag;
                    b[j] /= cur_diag;
                }
                for (int k = height - 1; k > -1; k--)
                {
                    if (k != j)
                    {
                        double cur = matrix[k, j];
                        for (int i = 0; i < width; i++)
                            matrix[k, i] -= matrix[j, i] * cur;
                        b[k] -= b[j] * cur;
                    }
                }
            }
            return b;
        }
        public static void MatrixPrint(double[] matrix)
        {
            Console.WriteLine();
            for (int i = 0; i < matrix.GetLength(0); i++)
                Console.Write(matrix[i] + " ");
            Console.WriteLine();
        }
        public static void MatrixPrint(double[,] matrix)
        {
            Console.WriteLine();
            for (int i = 0; i < matrix.GetLength(0); i++)
            {
                for (int j = 0; j < matrix.GetLength(1); j++)
                    Console.Write($"{matrix[i, j]:0.###}\t");
                Console.WriteLine();
            }
            Console.WriteLine();
        }
        public static double[,] CreateCoefficientMatrix(Condition[,] EqualityOfValues, Condition[,] EqualityOfDerivatives, Condition[,] ZeroingBoundaries, int prerequisites_num)
        {
            double[,] CoefMatrix = new double[prerequisites_num, prerequisites_num];
            int EOVHeight = EqualityOfValues.GetLength(0);
            int EOVWidth = EqualityOfValues.GetLength(1);
            int EODHeight = EqualityOfDerivatives.GetLength(0);
            int EODWidth = EqualityOfDerivatives.GetLength(1);

            for (int i = 0; i < CoefMatrix.GetLength(0); i++)
            {
                if (i < EOVHeight)
                {
                    for (int k = 0; k < CoefMatrix.GetLength(1); k++)
                        if (k >= i * 4 && k < i * 4 + 4)
                            CoefMatrix[i, k] = EqualityOfValues[i, 0].Coeffcients[k - i * 4];
                        else
                            CoefMatrix[i, k] = 0;
                }
                if (EOVHeight <= i && i < EOVHeight * EOVWidth)
                {
                    int curi = i - EOVHeight;
                    for (int k = 0; k < CoefMatrix.GetLength(1); k++)
                        if (k >= curi * 4 && k < curi * 4 + 4)
                            CoefMatrix[i, k] = EqualityOfValues[curi, 1].Coeffcients[k - curi * 4];
                        else
                            CoefMatrix[i, k] = 0;
                }
                if (EOVHeight * EOVWidth <= i && i < EOVHeight * EOVWidth + EODHeight)
                {
                    int curi = i - EOVHeight * EOVWidth;
                    for (int k = 0; k < CoefMatrix.GetLength(1); k++)
                        if (k >= curi * 4 && k < curi * 4 + 4)
                            CoefMatrix[i, k] = EqualityOfDerivatives[curi, 0].Coeffcients[k - curi * 4];
                        else if (k >= curi * 4 + 4 && k < curi * 4 + 8)
                            CoefMatrix[i, k] = -EqualityOfDerivatives[curi, 1].Coeffcients[k - curi * 4 - 4];
                        else
                            CoefMatrix[i, k] = 0;
                }
                if (i >= EOVHeight * EOVWidth + EODHeight && i < EOVHeight * EOVWidth + EODHeight * EODWidth / 2)
                {
                    int curi = i - EOVHeight * EOVWidth - EODHeight;
                    for (int k = 0; k < CoefMatrix.GetLength(1); k++)
                        if (k >= curi * 4 && k < curi * 4 + 4)
                            CoefMatrix[i, k] = EqualityOfDerivatives[curi, 2].Coeffcients[k - curi * 4];
                        else if (k >= curi * 4 + 4 && k < curi * 4 + 8)
                            CoefMatrix[i, k] = -EqualityOfDerivatives[curi, 3].Coeffcients[k - curi * 4 - 4];
                        else
                            CoefMatrix[i, k] = 0;
                }
                if (i == CoefMatrix.GetLength(0) - 2)
                {
                    for (int k = 0; k < 4; k++)
                        CoefMatrix[i, k] = ZeroingBoundaries[0, 0].Coeffcients[k];
                    for (int k = 4; k < CoefMatrix.GetLength(1); k++)
                        CoefMatrix[i, k] = 0;
                }
                if (i == CoefMatrix.GetLength(0) - 1)
                {
                    for (int k = CoefMatrix.GetLength(1) - 4; k < CoefMatrix.GetLength(1); k++)
                        CoefMatrix[i, k] = ZeroingBoundaries[0, 0].Coeffcients[k - CoefMatrix.GetLength(1) + 4];
                    for (int k = 0; k < CoefMatrix.GetLength(1) - 4; k++)
                        CoefMatrix[i, k] = 0;
                }
            }
            return CoefMatrix;
        }
        public static double[] CreateBMatrix(Condition[,] EqualityOfValues, Condition[,] EqualityOfDerivatives, Condition[,] ZeroingBoundaries, int prerequisites_num)
        {
            double[] bmatrix = new double[prerequisites_num];
            int EOVHeight = EqualityOfValues.GetLength(0);
            int EOVWidth = EqualityOfValues.GetLength(1);
            int EODHeight = EqualityOfDerivatives.GetLength(0);
            int EODWidth = EqualityOfDerivatives.GetLength(1);

            for (int i = 0; i < bmatrix.GetLength(0); i++)
            {
                if (i < EOVHeight)
                {
                    bmatrix[i] = EqualityOfValues[i, 0].b;
                } else
                if (EOVHeight <= i && i < EOVHeight * EOVWidth)
                {
                    int curi = i - EOVHeight;
                    bmatrix[i] = EqualityOfValues[curi, 1].b;
                } else 
                    bmatrix[i] = 0;
            }
            return bmatrix;
        }
    }
    class Condition
    {
        public double[] Coeffcients; //aka C
        public double b = 0;
        double point;
        public Condition(int degree, double point, double value)
        {
            this.point = point;
            Coeffcients = new double[degree + 1];
            for (int i = 0; i < Coeffcients.GetLength(0); i++)
                Coeffcients[i] = Math.Pow(point, i);
            //MatrixOperations.MatrixPrint(Coeffcients);
            b = value;
        }
        private Condition(double[] Coeffcients, double point, double b)
        {
            this.Coeffcients = Coeffcients;
            this.point = point;
            this.b = b;
        }
        public Condition Derivative(int derivate_num, double b = 0)
        {
            double[] newcofbody = new double[Coeffcients.GetLength(0)];
            for (int i = 0; i < newcofbody.GetLength(0); i++)
            {
                if (i >= derivate_num)
                {
                    int mult_cof = 1;
                    for (int j = 0; j < derivate_num; j++)
                        mult_cof *= (i) - j;
                    newcofbody[i] = Coeffcients[i] * mult_cof / Math.Pow(point, derivate_num);
                }
                else
                    newcofbody[i] = 0;
            }
            
            //MatrixOperations.MatrixPrint(newcofbody);
            return new Condition(newcofbody, point, b);

        }
    }
    class Spline
    {
        public double[] Coeffcients; // aka C
        public double startx;
        public double endx;
        public Spline(double[] Coeffcients, double a, double b)
        {
            this.Coeffcients = Coeffcients;
            this.startx = a;
            this.endx = b;
        }
        public double ValueInPoint(double point)
        {
            double result = 0;
            for (int i = 0; i < Coeffcients.GetLength(0); i++)
                result += Coeffcients[i] * Math.Pow(point, i); // += C[i]*x^i
            return result;
        }
    }
    class SplinePolynomial
    {
        Spline[] MainSplines;
        double a;
        double b;
        public SplinePolynomial(double[,] table, int splinesdegree) // table of points and function values
        {
            int splines_num = table.GetLength(0) - 1;
            int prerequisites_num = splines_num * (splinesdegree + 1); // number of conditions = number of variables
            Condition[,] EqualityOfValues = new Condition[splines_num, 2];
            Condition[,] EqualityOfDerivatives = new Condition[splines_num - 1, 4];
            Condition[,] ZeroingBoundaries = new Condition[2, 1];
            for (int i = 0; i < EqualityOfValues.GetLength(0); i++)
            {
                EqualityOfValues[i, 0] = new Condition(splinesdegree, table[i, 0], table[i, 1]);
                EqualityOfValues[i, 1] = new Condition(splinesdegree, table[i + 1, 0], table[i + 1, 1]);
            }
            for (int i = 0; i < EqualityOfDerivatives.GetLength(0); i++)
            {
                EqualityOfDerivatives[i, 0] = EqualityOfValues[i, 1].Derivative(1);
                EqualityOfDerivatives[i, 1] = EqualityOfValues[i + 1, 0].Derivative(1);
                EqualityOfDerivatives[i, 2] = EqualityOfValues[i, 1].Derivative(2);
                EqualityOfDerivatives[i, 3] = EqualityOfValues[i + 1, 0].Derivative(2);
            }
            ZeroingBoundaries[0, 0] = EqualityOfValues[0, 0].Derivative(2);
            ZeroingBoundaries[1, 0] = EqualityOfValues[EqualityOfValues.GetLength(0) - 1, 1].Derivative(2);

            double[,] CoefficientMatrix = MatrixOperations.CreateCoefficientMatrix(EqualityOfValues, EqualityOfDerivatives, ZeroingBoundaries, prerequisites_num);
            double[] B = MatrixOperations.CreateBMatrix(EqualityOfValues, EqualityOfDerivatives, ZeroingBoundaries, prerequisites_num);
            //MatrixOperations.MatrixPrint(CoefficientMatrix);
            //MatrixOperations.MatrixPrint(B);
            double[] PolynomialCoefficients = MatrixOperations.GausMethod(CoefficientMatrix, B);
            //MatrixOperations.MatrixPrint(PolynomialCoefficients);
            MainSplines = new Spline[table.GetLength(0) - 1];
            for (int i = 0; i < MainSplines.GetLength(0); i++)
            {
                double[] coef_body = new double[splinesdegree + 1];
                for (int k = 0; k < coef_body.GetLength(0); k++)
                    coef_body[k] = PolynomialCoefficients[i * coef_body.GetLength(0) + k];
                MainSplines[i] = new Spline(coef_body, table[i, 0], table[i + 1, 0]);
            }
            b = table[table.GetLength(0) - 1, 0];
            a = table[0, 0];
        }
        public double ValueInPoint(double x)
        {
            double result = 0;
            if (x < a) return MainSplines[0].ValueInPoint(x);
            if (x >= b) return MainSplines.Last().ValueInPoint(x);
            for (int i = 0; i < MainSplines.GetLength(0); i++)
                if (x >= MainSplines[i].startx && x < MainSplines[i].endx)
                    return MainSplines[i].ValueInPoint(x);
            return result;
        }

    }
    class Program
    {
        public static double function(double x)
        {
            return 3 * x - Math.Cos(x) - 1;
        }
        public static double[,] CreateTable(double a, double b, int numnode)
        {
            double[,] table = new double[numnode + 1, 2];
            double x = a;
            for (int i = 0; i < numnode + 1; i++)
            {
                table[i, 0] = x;
                table[i, 1] = function(x);
                x += (b - a) / numnode;
            }
            return table;
        }
        static int splinesdegree = 3; 
        static void Main(string[] args)
        {
            double a = -3;
            double b = 3;
            int numnode = 3;
            output(a, b, numnode, 40);
            Console.ReadLine();
        }
        public static void output(double a, double b, int numnode, int accuracy)
        {
            double[,] table = CreateTable(a, b, numnode);
            Console.WriteLine("SPLINE INTERPOLATION METHOD");
            Console.WriteLine("###########################");
            Console.WriteLine("x\t\tf(x)\t\tSpl(x)\t\t|f(x)-Spl(x)|");
            SplinePolynomial pol = new SplinePolynomial(table, splinesdegree);
            double x = a;
            int i = 0;
            do
            {
                i++;
                double polvalue = pol.ValueInPoint(x);
                var func = function(x);
                Console.WriteLine($"{x:0.###}\t\t{func:0.###}\t\t{polvalue:0.###}\t\t{(Math.Abs(func - polvalue)):0.###}");
                x += (b - a) / accuracy;
            } while (i <= accuracy);
        }
    }
}

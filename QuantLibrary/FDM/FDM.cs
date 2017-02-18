using System;
using MathNet.Numerics.LinearAlgebra;

namespace QuantLibrary
{
    public class FDM
    {
        IBSPde pde;

        public Vector<double> a, bb, c;
        public Vector<double> rhs;
        public Matrix<double> L;

        double[][] alpha = new double[3][];
        double[][] beta = new double[3][];
        double[][] gamma = new double[3][];
        double[][] delta = new double[3][];

        public Vector<double> vecOld, result;

        //Constructor
        public FDM(IBSPde myPDE)
        {
            pde = myPDE;
        }

        // Initialize the grid
        public void initIC(Vector<double> xarr)
        {
            //Initalize solution at time zero
            vecOld = Vector<double>.Build.Dense(xarr.Count, xarr.MinimumIndex());

            //Initialize boundaries
            vecOld[0] = pde.bcl(0.0);
            vecOld[vecOld.Count - 1] = pde.bcr(0.0);

            //Initialize values in the interior  
            for (int i = 1; i < xarr.Count - 1; i++)
            {
                vecOld[i] = pde.ic(xarr[i]);
            }
            result = vecOld; // Initialize the result
        }

        // Get current solution
        public Vector<double> current()
        {
            return result;
        }

        public void calculateDifferentialCoefficients(Vector<double> xarr)
        {
            beta[0] = new double[xarr.Count];
            beta[1] = new double[xarr.Count];
            beta[2] = new double[xarr.Count];

            delta[0] = new double[xarr.Count];
            delta[1] = new double[xarr.Count];
            delta[2] = new double[xarr.Count];

            double partialXi;
            double partialXiPlus;
            double partialXiPlusPlus;
            double partialXiMinus;

            // Lower Boundary
            // Using forward differeces - gamma
            int j = 0;
            partialXiPlus = xarr[j + 1] - xarr[j];
            partialXiPlusPlus = xarr[j + 2] - xarr[j + 1];

            beta[0][j] = (-2 * partialXiPlus + partialXiPlusPlus) / (partialXiPlus * (partialXiPlus + partialXiPlusPlus));
            beta[1][j] = (partialXiPlus + partialXiPlusPlus) / (partialXiPlus * partialXiPlusPlus);
            beta[2][j] = -partialXiPlus / (partialXiPlusPlus * (partialXiPlus + partialXiPlusPlus));

            // Interior
            // Using central differences - beta
            for (int i = 1; i < xarr.Count - 1; i++)
            {
                partialXi = xarr[i] - xarr[i - 1];
                partialXiPlus = xarr[i + 1] - xarr[i];

                // Beta
                beta[0][i] = -partialXiPlus / (partialXi * (partialXi + partialXiPlus));
                beta[1][i] = (partialXiPlus - partialXi) / (partialXi * partialXiPlus);
                beta[2][i] = partialXi / (partialXiPlus * (partialXi + partialXiPlus));


                // Delta
                delta[0][i] = 2 / (partialXi * (partialXi + partialXiPlus));
                delta[1][i] = -2 / (partialXiPlus * partialXi);
                delta[2][i] = 2 / (partialXiPlus * (partialXi + partialXiPlus));
            }

            // Upper boundary
            // Using backwarddifferences
            j = xarr.Count - 1;
            partialXi = xarr[j] - xarr[j - 1];
            partialXiMinus = xarr[j - 1] - xarr[j - 2];

            beta[0][j] = partialXi / (partialXiMinus * (partialXiMinus + partialXi));
            beta[1][j] = (-partialXiMinus - partialXi) / (partialXiMinus * partialXi);
            beta[2][j] = (partialXiMinus + 2 * partialXi) / (partialXi * (partialXiMinus + partialXi));
        }

        // Calculate the coefficient at given time point
        public void calculateCoefficients(Vector<double> xarr, double tprev, double tnow)
        {
            // Explicit method
            a = Vector<double>.Build.Dense(xarr.Count - 1, 0.0);
            bb = Vector<double>.Build.Dense(xarr.Count, 0.0);
            c = Vector<double>.Build.Dense(xarr.Count - 1, 0.0);
            rhs = Vector<double>.Build.Dense(xarr.Count, 0.0);
            L = Matrix<double>.Build.Dense(xarr.Count, xarr.Count);

            double k = tnow - tprev; //Time step size

            // Interior
            for (int i = 0; i < xarr.Count; i++)
            {

                bb[i] = 1 - k * (pde.mu(xarr[i], tprev) * beta[1][i] + pde.sigma(xarr[i], tprev) * delta[1][i]);

                rhs[i] = k * pde.f(xarr[i], tprev);

                if (i < xarr.Count - 1)
                    L[i, i + 1] = -k * (pde.mu(xarr[i], tprev) * beta[2][i] + pde.sigma(xarr[i], tprev) * delta[2][i]);
                L[i, i] = bb[i];
                if (i > 0)
                    L[i, i - 1] = -k * (pde.mu(xarr[i], tprev) * beta[0][i] + pde.sigma(xarr[i], tprev) * delta[0][i]);
            }
            // Boundaries
            //L[0, 0] = 1 - k * (pde.mu(xarr[i], tprev) * beta[1][i] + pde.sigma(xarr[0], tprev) * delta[1][0]); ;
            Console.WriteLine(L.ToString());
        }
        //// Calculate the coefficient at given time point
        //public void calculateCoefficients(Vector<double> xarr, double tprev, double tnow)
        //{
        //    // Explicit method
        //    a = Vector<double>.Build.Dense(xarr.Count - 1, 0.0);
        //    bb = Vector<double>.Build.Dense(xarr.Count - 1, 0.0);
        //    c = Vector<double>.Build.Dense(xarr.Count - 1, 0.0);
        //    rhs = Vector<double>.Build.Dense(xarr.Count, 0.0);
        //    L = Matrix<double>.Build.Dense(xarr.Count, xarr.Count);

        //    double tmp1, tmp2;
        //    double k = tnow - tprev; //Time step size
        //    double h; // Spatial step size

        //    for (int i = 1; i < xarr.Count - 1; i++)
        //    {
        //        h = xarr[i + 1] - xarr[i];

        //        tmp1 = k * (pde.sigma(xarr[i], tprev) / (h * h));
        //        tmp2 = k * (pde.mu(xarr[i], tprev) * 0.5 / h);

        //        a[i] = tmp1 - tmp2;
        //        bb[i] = 1.0 - (2.0 * tmp1) + (k * pde.b(xarr[i], tprev));
        //        c[i] = tmp1 + tmp2;
        //        rhs[i] = k * pde.f(xarr[i], tprev);



        //        L[i, i + 1] = c[i];
        //        L[i, i] = bb[i];
        //        L[i, i - 1] = a[i];
        //    }
        //    //Console.WriteLine(L.ToString());
        //}

        // Solve the time step given the current coefficients
        public void solve(double tnow)
        {
            result[0] = pde.bcl(tnow);
            result[result.Count - 1] = pde.bcr(tnow);

            //for (int i = 1; i < result.Count - 1; i++) //Loop in spatial (substitute with matrix algebra?)
            //{
            //    result[i] = (a[i] * vecOld[i - 1])
            //                + (bb[i] * vecOld[i])
            //                + (c[i] * vecOld[i + 1])
            //                - rhs[i];
            //}

            result = L.Multiply(vecOld).Subtract(rhs);
            Console.WriteLine(L.ToString());
            //result = L.Solve(vecOld);
            //Console.WriteLine(result.ToString());
            Console.ReadLine();
            vecOld = result;
        }
    }
}


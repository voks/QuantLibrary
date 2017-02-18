using MathNet.Numerics.LinearAlgebra;
using QuantLibrary;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace QuantLib.FDM.OneDimensional
{
    public class ThetaMethod
    {
        IBSPde pde;

        private Vector<double> c, u, l;
        private Vector<double> rhs;
        private Matrix<double> A;
        private Matrix<double> I;

        private Vector<double> vecOld, result;
        private Vector<double> xarr;

        //Constructor
        public ThetaMethod(IBSPde myPDE, Vector<double> xarr)
        {
            this.xarr = xarr;
            I = Matrix<double>.Build.DiagonalIdentity(xarr.Count());
            A = Matrix<double>.Build.Dense(xarr.Count(),xarr.Count());
            c = Vector<double>.Build.Dense(xarr.Count());
            u = Vector<double>.Build.Dense(xarr.Count());
            l = Vector<double>.Build.Dense(xarr.Count());
            pde = myPDE;
        }

        // Initialize the grid
        private void initIC()
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


        // Calculate the coefficient at given time point
        private void calculateCoefficients(double tprev, double tnow)
        {
            double k = tnow - tprev; //Time step size
            double dX = xarr[1] - xarr[0];

            // Interior
            for (int i = 0; i < xarr.Count-1; i++)
            {
                c[i] = -Math.Pow(pde.sigma(xarr[i], tprev), 2) / Math.Pow(dX, 2) - pde.f(xarr[i], tprev);
                A[i, i] = c[i];
                if (i != xarr.Count()-1)
                {
                    u[i] = 0.5 * pde.mu(xarr[i], tprev) / dX + 0.5 * pde.sigma(xarr[i], tprev) / Math.Pow(dX, 2);
                    A[i, i + 1] = u[i];
                }
                if (i != 0)
                {
                    l[i] = -0.5 * pde.mu(xarr[i], tprev) / dX + 0.5 * pde.sigma(xarr[i], tprev) / Math.Pow(dX, 2);
                    A[i, i - 1] = l[i];
                }
            }
        }

        private Vector<double> getBoundaries(double dt)
        {
            Vector<double> boundaries = Vector<double>.Build.Dense(xarr.Count());
            boundaries[0] = pde.bcl(dt);
            boundaries[boundaries.Count()-1] = pde.bcr(dt);
            return boundaries;
        }

        private void solveTridiagonal(double theta, double dt)
        {
            //theta = 1;

            var getB = getBoundaries(dt);
            var tmp = I-dt * I.Multiply(A);
            var tmp2 = vecOld;

            // theta = 0
            var tmp = vecOld;
            var tmp2 = (I + dt * A);


            //var getB = getBoundaries(dt);
            //var tmp = I.Add((1 - theta) * dt).Multiply(A);
            //var tmp2 = (I-theta * dt *A*(vecOld;
            result = tmp.Solve(tmp2);
        }

        public void solve(double[] t, double theta)
        {
            initIC();
            for (int i = 1; i < t.Count(); i++)
            {
                calculateCoefficients(t[i - 1], t[i]);
                solveTridiagonal(theta, t[i] - t[i - 1]);
                vecOld = result;
            }
        }
    }
}

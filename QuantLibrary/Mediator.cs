using QuantLib.FDM.OneDimensional;
using MathNet.Numerics.LinearAlgebra;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace QuantLibrary
{
    class Mediator
    {
        static void Main()
        {
            double spot = 100;
            double expiry = 2;
            double strike = 100;
            double volatility = 0.1;
            double interest = 0.02;
            double dividend = 0.01;
            string type = "CALL";

            QLNet.PlainVanillaPayoff payoff = new QLNet.PlainVanillaPayoff(QLNet.Option.Type.Call, strike);
            QLNet.BlackScholesCalculator bscal = new QLNet.BlackScholesCalculator(payoff, spot, dividend, volatility*expiry, interest);

            Console.WriteLine(bscal.value());
            //BlackScholesClosedForm bs = ConsoleBSOptionFactory.create();
            BlackScholesClosedForm bs = new BlackScholesClosedForm(type, expiry, strike, interest, dividend, volatility);
            Console.WriteLine(bs.Price(spot));
            double[] price = bs.Price(0, 500, 200);
            //for (int i = 0; i < price.Length; i++)
            //{
            //    Console.WriteLine(price[i]);
            //}
            ////Console.ReadLine();

            // Finite Difference
            double truncation = 500;

            PdeBS pdeBS = new PdeBS(expiry, strike, volatility, interest, dividend, truncation, type);
            Mesh xMesh = new Mesh(100, 0, 500);
            Mesh tMesh = new Mesh(1000, 0, 2);

            ThetaMethod test = new ThetaMethod(pdeBS, xMesh.getMesh());
            test.solve(tMesh.getMesh().ToArray(),1);
            test.current();
            //FDM fdm = new FDM(pdeBS);
            //FDMDirector fdmDirector = new FDMDirector(fdm, xMesh.getMesh(), tMesh.getMesh());

            //fdmDirector.Start();
            //Vector<double> result = fdmDirector.current();
            //double valueResult = MathNet.Numerics.Interpolate.Linear(xMesh.getMesh(), result).Interpolate(strike);

            ////double[] test = new double[price.Length];
            ////for (int i = 0; i < price.Length; i++)
            ////{
            ////    test[i] = price[i] - result[i];
            ////    Console.WriteLine(test[i]);
            ////}

            //Console.WriteLine(valueResult);
            //Console.ReadLine();
        }
    }
}

using System;

namespace QuantLibrary
{
    public class ConsoleBSOptionFactory
    {
        public static BlackScholesClosedForm create()
        {
            double r;
            double sig;
            double K;
            double T;
            double d;
            string type;

            Console.Write("Strike ");
            K = Convert.ToDouble(Console.ReadLine());

            Console.Write("Sig ");
            sig = Convert.ToDouble(Console.ReadLine());

            Console.Write("Rate ");
            r = Convert.ToDouble(Console.ReadLine());

            Console.Write("Dividend ");
            d = Convert.ToDouble(Console.ReadLine());

            Console.Write("Expiry ");
            T = Convert.ToDouble(Console.ReadLine());

            Console.Write("Type ");
            type = Convert.ToString(Console.ReadLine());

            BlackScholesClosedForm option = new BlackScholesClosedForm(type, T, K, d, r, sig);
            return option;
        }
    }

    public class BlackScholesClosedForm : IOption
    {
        private double r;
        private double sig;
        private double K;
        private double T;
        private double d;
        private string type;

        //Constructor
        public BlackScholesClosedForm(string optionType, double expiry, double strike,
             double interest, double dividend, double volatility)
        {
            type = optionType;
            T = expiry;
            K = strike;
            d = dividend;
            r = interest;
            sig = volatility;
        }

        private double CallPrice(double S)
        {
            double tmp = sig * Math.Sqrt(T);
            double d1 = (Math.Log(S / K) + (r - d + (sig * sig) * 0.5) * T) / tmp;
            double d2 = d1 - tmp;
                
            return S * Math.Exp(-d * T) * MathNet.Numerics.Distributions.Normal.CDF(0, 1, d1)
                - K * Math.Exp(-r * T) * MathNet.Numerics.Distributions.Normal.CDF(0, 1, d2);
        }

        private double PutPrice(double S)
        {
            double tmp = sig * Math.Sqrt(T);
            double d1 = (Math.Log(S / K) + (r - d + (sig * sig) * 0.5) * T) / tmp;
            double d2 = d1 - tmp;

            return K * Math.Exp(-r * T) * MathNet.Numerics.Distributions.Normal.CDF(0, 1, -d2)
                - S * Math.Exp(-d * T) * MathNet.Numerics.Distributions.Normal.CDF(0, 1, -d1);
        }

        public double CallDelta(double U)
        {
            double tmp = sig * Math.Sqrt(T);
            double d1 = (Math.Log(U / K) + (d + (sig * sig) * 0.5) * T) / tmp;
            return Math.Exp((d - r) * T) * MathNet.Numerics.Distributions.Normal.CDF(d1, 0, 1);
        }

        public double CallGamma(double U)
        {
            double tmp = sig * Math.Sqrt(T);
            double d1 = (Math.Log(U / K) + (d + (sig * sig) * 0.5) * T) / tmp;

            return MathNet.Numerics.Distributions.Normal.CDF(d1, 0, 1) * Math.Exp((d - r) * T) / (U * tmp);
        }

        public double Price(double spot)
        {
            if (type.ToUpper() == "CALL")
            {
                return CallPrice(spot);
            }
            else
            {
                return PutPrice(spot);
            }

        }

        public double Delta(double spot)
        {
            return CallDelta(spot);

        }

        public double Gamma(double spot)
        {
            return CallGamma(spot);

        }
    }
}

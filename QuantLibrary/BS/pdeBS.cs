using System;
namespace QuantLibrary
{
    public class PdeBS : IBSPde
    {
        private double T, K, vol, r, d, Smax;
        private string type;

        //Constructor
        public PdeBS(double expiry, double strike, double volatility, double interest,
            double dividend, double truncation, string type)
        {
            T = expiry;
            K = strike;
            vol = volatility;
            r = interest;
            d = dividend;
            Smax = truncation;
            this.type = type;
        }

        public double sigma(double x, double t)
        {
            double sigmaS = vol * vol;
            return 0.5 * sigmaS * x * x;
        }
        public double mu(double x, double t)
        {
            return (r - d) * x;
        }
        public double b(double x, double t)
        {
            return -r;
        }
        public double f(double x, double t)
        {
            return 0.0;
        }

        // Left boundary
        public double bcl(double tau)
        {
            if (type == "CALL")
            {
                return 0.0;
            }
            else
            {
                return K * Math.Exp(-r * tau);
            }
        }

        // Right boundary
        public double bcr(double tau)
        {
            if (type == "CALL")
            {
                return Smax;
            }
            else
            {
                return 0.0;
            }
        }

        // Initial condition
        public double ic(double x)
        {
            //Put: max(0,K-x)
            if (x < K & type == "PUT")
            {
                return K - x;
            }
            //Call: max(0,x-K)
            if (x > K & type == "CALL")
            {
                return x - K;
            }
            return 0.0;
        }
    }
}


using System;
namespace QuantLibrary
{
	public static class OptionMixins
	{
		public static void Display(this BlackScholesClosedForm option, double S)
		{
			Console.WriteLine(
				"Display Ooption price+greeks in an Extension method");
			Console.WriteLine("Price: {0}", option.Price(S));
			Console.WriteLine("Delta: {0}", option.Delta(S));
			Console.WriteLine("Gamma: {0}", option.Gamma(S));
			//Console.WriteLine("Vega: {0}", option.Vega(S));
			//Console.WriteLine("Theta: {0}", option.Theta(S));
			//Console.WriteLine("Rho: {0}", option.Rho(S));

		}

		public static double[] Price(this BlackScholesClosedForm option, double low, double upper, int NSteps)
		{
			double h = (upper - low) / NSteps;
			double S = low;
			double[] price = new double[NSteps];

			for (int i = 0; i < NSteps; i++)
			{
				price[i] = option.Price(S);
				S += h;
			}
			return price;
		}
	}
}

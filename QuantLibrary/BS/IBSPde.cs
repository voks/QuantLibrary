using System;
namespace QuantLibrary
{
	public interface IBSPde
	{
		//Coefficient of PDE
		double sigma(double x, double t);
		double mu(double x, double t);
		double b(double x, double t);
		double f(double x, double t);

		//Boundary conditions (Dirichlet)
		double bcl(double t);
		double bcr(double t);

		//Initial conditioon
		double ic(double x);

	}
}


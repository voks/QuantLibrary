using QuantLibrary;
using System;
namespace QuantLibrary
{
	public interface IOptionFactory
	{
		IOption create();
	}

	public interface IOption
	{
        double Price(double Spot);
        double Delta(double Spot);
	}
}


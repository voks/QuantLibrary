using MathNet.Numerics.LinearAlgebra;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace QuantLibrary
{
    public class Mesh
    {
        private int numberOfElements;
        private double lowerValue;
        private double upperValue;

        // Constructor
        public Mesh(int numberOfElements, double lowerValue, double upperValue)
        {
            this.numberOfElements = numberOfElements;
            this.lowerValue = lowerValue;
            this.upperValue = upperValue;
        }

        //Public methods
        public Vector<double> getMesh()
        {
            double increments = (upperValue - lowerValue) / (double)numberOfElements;

            Vector<double> xmesh = Vector<double>.Build.Dense(numberOfElements);
            for (int i = 0; i < numberOfElements; i++)
            {
                xmesh[i] = lowerValue + increments *  i;
            }
            return xmesh;
        }
    }
}

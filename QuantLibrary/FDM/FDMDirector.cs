using System;
using MathNet.Numerics.LinearAlgebra;

namespace QuantLibrary
{
	public class FDMDirector
	{
		private double tprev, tnow;
		private Vector<double> xarr;
		private Vector<double> tarr;
		private FDM fdm;

		public FDMDirector(FDM fdScheme, Vector<double> xmesh, Vector<double> tmesh)
		{
			fdm = fdScheme;
			xarr = xmesh;
			tarr = tmesh;
		}

		public Vector<double> current()
		{
			return fdm.current();
		}

		// Run
		public void Start()
		{
			fdm.initIC(xarr);
			doit();
		}

		public void doit()
		{
			tnow = tprev = tarr[tarr.MinimumIndex()];

			for (int n = tarr.MinimumIndex() + 1; n <= tarr.MaximumIndex(); n++)
			{
				tnow = tarr[n];// get time point from mesh
                fdm.calculateDifferentialCoefficients(xarr);
				fdm.calculateCoefficients(xarr, tprev, tnow); //get coefficients
				fdm.solve(tnow); // solve the time step
				tprev = tnow; // update tprev
			}
		}

	}
}


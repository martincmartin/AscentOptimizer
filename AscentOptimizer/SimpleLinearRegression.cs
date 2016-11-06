using System;
namespace AscentOptimizer
{
	public class SimpleLinearRegression
	{
		// We want to choose a, b to minimize:
		//
		// sum( 1/2 * wi * (yi - a - b * xi)^2)
		//
		// diff w.r.t b:
		// 0 = sum( wi * (yi - a - b * xi) * (- xi) )
		// 0 = - sum(wi * yi * xi) + sum( wi * a * xi) + sum(wi * b * xi^2)
		// sum(wi * yi * xi) = a * sum(wi * xi) + b * sum( wi * xi * xi)
		//
		// diff w.r.t a:
		// 0 = sum(wi * (yi - a - b * xi) * -1)
		// 0 = - sum(wi * yi) + sum(wi * a) + sum(wi * b * xi)
		// sum(wi * yi) = a * sum(wi) + b * sum(wi * xi)
		//
		// Solve for a:
		// a = (sum(wi * yi) - b * sum(wi * xi)) / sum(wi)
		// Substitue in the other formula above:
		//
		// sum(wi * yi * xi) = (sum(wi * yi) - b * sum(wi * xi)) / sum(wi) * sum(wi * xi) + b * sum(wi * xi * xi)
		// Multiply through by sum(wi):
		// sum(wi) * sum(wi * yi * xi) = (sum(wi * yi) - b * sum(wi * xi)) * sum(wi * xi) + b * sum(wi * xi * xi) * sum(wi)
		//
		// Solve for b:
		//
		// sum(wi) * sum(wi * xi * yi) - sum(wi * yi) * sum(wi * xi) = b * (sum(wi * xi)^2 + sum(wi * xi * xi) * sum(wi))
		//
		//     sum(wi) * sum(wi * xi * yi) - sum(wi * yi) * sum(wi * xi)
		// b = ---------------------------------------------------------
		//           sum(wi * xi)^2 + sum(wi * xi * xi) * sum(wi)
		//
		// Var(x) = sum(wi * (xi - x_)^2) / sum(wi)
		//        = sum(wi * (xi^2 - 2 * xi * x_ + x_^2)) / sum(wi)
		//        = sum(wi * xi^2) / sum(wi) - 2 * x_ * sum(wi * xi)/sum(wi) + x_^2
		//        = sum(wi * xi^2) / sum(wi) - x_^2
		//
		// What do we do when the denominator is close to zero?  The denominator is essentially the variance of the
		// observed xis, so this happens when we keep the independent variable fixed (e.g. keep the controls on the
		// spaceship fixed).  We can just average in a prior using an Inverse Variance Weighting.  So what's the
		// variance?  The standard error in the estimate of b is average error divided by the variance.
		public double sum_w_x = 0;
		public double sum_w_y = 0;
		public double sum_w_x_x = 0;
		public double sum_w_x_y = 0;
		public double sum_w = 0;

		public SimpleLinearRegression()
		{
		}

		public void observe(double x, double y, double w)
		{
			double w_x = w * x;
			sum_w_x += w_x;
			sum_w_y += w * y;
			sum_w_x_x += w_x * x;
			sum_w_x_y += w_x * y;
			sum_w += w;
		}

		public void decay(double gamma)
		{
			sum_w_x *= gamma;
			sum_w_y *= gamma;
			sum_w_x_x *= gamma;
			sum_w_x_y *= gamma;
			sum_w *= gamma;
		}

		public void solve(out double intercept, out double x_coefficient, out double x_var)
		{
			// What to do when this denominator is near zero?
			//
			// Happens when sum(wi * (xi - x_)^2) == 0, i.e. when all xis are the same (assuming weights are non-zero).
			//
			// Regularize?  Give confience interval?
			double x_mean = sum_w_x / sum_w;
			x_coefficient = (sum_w_x_y * sum_w - sum_w_x * sum_w_y) / (sum_w_x_x * sum_w - sum_w_x * sum_w_x);
			intercept = sum_w_y / sum_w - x_coefficient * x_mean;
			x_var = sum_w_x_x / sum_w - x_mean * x_mean;
		}
	}
}

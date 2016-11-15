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

		public double prior_x_y_covariance;
		public double prior_var_x;

		public SimpleLinearRegression(double prior_x_coefficient, double prior_var_x)
		{
			this.prior_var_x = prior_var_x;
			prior_x_y_covariance = prior_x_coefficient * prior_var_x;
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

		public void solve(out double intercept, out double x_coefficient, out double var_x)
		{
			// What to do when this denominator is near zero?
			//
			// Happens when sum(wi * (xi - x_)^2) == 0, i.e. when all xis are the same (assuming weights are non-zero).
			//
			// We do inverse variance weighting!  The variance of x_coefficient is proportional to 1/var_x.  So,
			// inverse variance weighting is:
			//
			// x_coeff / var(x_coeff) + prior_x_coeff / var(prior_x_coeff)
			// -----------------------------------------------------------
			//     1 / var(x_coeff) + 1 / var(prior_x_coeff)
			//
			//   x_coeff * var_x + prior_x_coeff * prior_var_x
			// = ----------------------------------------------
			//              var_x + prior_var_x
			//
			// Where I kind of glossed over the constant of proportionality (which is the variance of the y error),
			// assuming it was the same for both the observed & the prior.
			//
			// Note that x_coeff * var_x is just x_y_covariance.

			double mean_x = sum_w_x / sum_w;
			double mean_y = sum_w_y / sum_w;
			var_x = sum_w_x_x / sum_w - mean_x * mean_x;
			double x_y_covariance = (sum_w_x_y / sum_w - mean_x * mean_y);
			x_coefficient = (x_y_covariance + prior_x_y_covariance) / (var_x + prior_var_x);
			intercept = mean_y - x_coefficient * mean_x;
		}
	}
}

I also added a file project_template_deriv_onto_gh.h5 that contains the results of: 
Computing a template spectrum for that line region,
Applying a Doppler shift (rv in m/s), 
Subtracting the template,
Projecting the residuals onto the Gauss-Hermite (or truncated Gauss-Hermite) basis functions, 
Fitting a quadratic polynomial to the relationship of the resulting scores versus rv.

The first dimension is the order of the GH coefficients (starting with order 0 through 12).
The second dimension is order of the polynomial fit (so the first value is the zero point (should be zero, but actual values aren't quite zero), the second is the linear term i.e. the important one, and the final is the quadratic term so you can get a sense for how good the linear approximation is).
The third dimension is the index of the (line,order) pair.

I inspected three lines and the results looked reasonable to me.  But it's certainly possible that you'll find something suspicious for some lines.
The main reason that the values in gh_vs_mps_coeffs[2,2,:] aren't the same is because the lines have different depths. 
If you plot median_fit_gauss_depth vs tgh_vs_mps_coeffs[2,2,:], then it looks roughly linear.
There's non-negligible scatter.  There are obviously some outliers that are probably best removed (e.g., lines that are very narrow, have a<<1, have b = -0.7 or 0.7 due to limit of non-linear fitting, etc.).  But once if you exclude those, I don't know what causing the scatter... is it really the slight differences in the line shapes? 
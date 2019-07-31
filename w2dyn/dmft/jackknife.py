import numpy as np
import h5py
import time
import sys


# TODO: Exceptions
# TODO: Documentation
class Jackknife(object):
    """A class implementing a jackknife estimation.

    The jackknife is a resampling method that allows to estimate the
    variance and bias of a parameter. This can further be used to obtain
    the bias-corrected jackknife estimator.

    Given :math:`n` input samples :math:`x_i` and a transformation
    function :math:`f()`, the goal is to estimate :math:`y=f(x)`, with
    :math:`x` and :math:`y` being the true values. A simple estimate
    would just be the transformed input mean :math:`f(\\bar{x})`, with
    :math:`\\bar{x}` being the input sample mean. In general this
    carries a bias of order :math:`1/n`, where :math:`n` is the sample
    size. Using jackknife resampling this leading :math:`1/n` term can
    be estimated and removed, yielding the bias-corrected jackknife
    estimator.

    First, the resampling of the input is done by:

    .. math::
        x_i \\rightarrow x'_i = (n\\bar{x} - x_i) / (n-1),

    where :math:`x'_i` is the new sample. This is then transformed in
    the following way:

    .. math::
        y_i = n f(\\bar{x}) - (n-1) f(x'_i),

    where :math:`y_i` is the transformed, bias-corrected output sample.
    The jackknife estimator :math:`\\hat{y}_{\\textrm{jackknife}}` can
    now be obtained simply by calculating the sample mean of the
    :math:`y_i`. Other statistical quantities like the variance and
    covariance for example can also be calculated from the output
    samples. For more details on that, see :meth:`do_estimation`.
    """

    def __init__(self, input_sample_generator_function, n_samples,
                 transformation_function, correlation_axis=-1,
                 output_file_prefix="", output_samples_should_be_stored=False):
        #: A function that returns a generator yielding the input
        #: samples :math:`x_i`
        self.input_sample_generator_function = input_sample_generator_function
        #: :math:`n`
        self.n_samples = n_samples
        if n_samples < 2:
            raise ValueError("""At least 2 samples are needed for the jackknife
                             estimation.""")
        #: :math:`f()`
        self.transformation_function = transformation_function
        #: Axis along which the covariance and correlation matrices
        #: are calculated
        self.correlation_axis = correlation_axis
        #:
        self.output_file_prefix = output_file_prefix
        self._date_and_time = ""
        self._output_file_name = ""
        self._input_mean = None
        self._transformed_mean = None
        self._output_mean = None
        self._output_variance = None
        self._output_standard_deviation = None
        self._output_standard_error_of_mean = None
        self._output_covariance = None
        self._output_correlation = None
        self._output_samples = []
        #:
        self.should_store_output_samples = output_samples_should_be_stored

    @property
    def output_file_name(self):
        return self._output_file_name

    def do_estimation(self):
        """Estimate :math:`y=f(x)` and some statistical quantities.

        The jackknife estimator :math:`\\hat{y}_{\\textrm{jackknife}}`
        is simply given by the sample mean :math:`\\bar{y}`,

            .. math::
                \\hat{y}_{\\textrm{jackknife}} = \\bar{y} = \\frac{1}{n}
                    \\sum_i y_i,

        with output samples :math:`y_i` and sample size :math:`n`. Other
        important statistical quantities are estimated with the sample
        variance :math:`s^2`, the sample standard deviation :math:`s`,
        the standard error of the mean :math:`\\textrm{SEM}`, the sample
        covariance matrix :math:`Q` and the sample correlation matrix
        :math:`R`. They are calculated using the following formulas:

            .. math::
                s^2 &= \\frac{1}{n-1}\\left(\\sum_i y_i^2 - \\frac{1}{n}
                    \\left(\\sum_i y_i\\right)^2\\right),\\\\
                s& = \\sqrt{s^2},\\\\
                \\textrm{SEM} &= \\frac{s}{\\sqrt{n}},\\\\
                Q &= \\frac{1}{n-1} \\left(\\sum_i y_i\\otimes y_i -
                    \\frac{1}{n}\\left(\\sum_i y_i\\right) \\otimes
                    \\left(\\sum_i y_i\\right)\\right),\\\\
                R_{ij} &= \\frac{Q_{ij}}{s_i s_j}.

        Here, :math:`\\otimes` denotes the outer product along the
        :attr:`correlation_axis`.

        .. note::
            While the sample variance :math:`s^2` is an unbiased
            estimator, the sample standard deviation :math:`s` is not.
            It generally underestimates the standard deviation. Because
            of that, also the standard error of the mean is
            underestimated.
        """

        start = time.time()
        print("-- Begin jackknife estimation")
        self._date_and_time = time.strftime("%Y-%m-%d_%H-%M-%S")

        get_samples = self.input_sample_generator_function
        f = self.transformation_function
        n = self.n_samples

        print("-- Calculating transformed input mean ...")
        self._input_mean = sum(i for i in get_samples()) / n
        self._transformed_mean = f(self._input_mean)

        # To reduce round-off errors when calculating the sample
        # variance, the sample mean should be close to zero. Therefore
        # all samples should be shifted by it. Since the mean is not
        # known thus far, it is estimated with a random (in this case
        # the first) y value.
        print("-- Calculating shift ...")
        samples = get_samples()
        shift = self._resample_and_transform(next(samples))
        sum_ = 0
        sum_of_squares = 0
        sum_of_outer_products = 0
        i = 1
        for sample in samples:
            print("-- Calculating sample {} ...".format(i))
            y = self._resample_and_transform(sample)
            y -= shift
            if self.should_store_output_samples:
                self._output_samples.append(y)
            sum_ += y
            sum_of_squares += np.abs(y)**2
            sum_of_outer_products += self._calculate_outer_product(y, y)
            i += 1

        self._calculate_statistical_quantities(shift, sum_, sum_of_squares,
                                               sum_of_outer_products)

        end = time.time()
        print("\n-- The jackknife estimation took " + str(end - start) + "s.")
        return

    def write_results_to_file(self):
        """Write the results of the last jackknife estimation to a file.

        Write the following quantities of the last jackknife run to an
        HDF5 file:

        - sample mean
        - sample variance
        - sample standard deviation
        - standard error of the mean
        - sample covariance matrix
        - sample correlation matrix
        - transformed input mean
        - number of samples

        The name of the file consists of the :attr:`output_file_prefix`
        as well as the date and time of the start of the last jackknife
        estimation.
        """

        self._output_file_name = self.output_file_prefix + "_"
        self._output_file_name += self._date_and_time + ".hdf5"
        output_file = h5py.File(self._output_file_name, "w")
        output_file.create_group(".config")
        output_file[".config"].attrs.create("n_samples", self.n_samples)
        output_file.create_dataset("mean", data=self._output_mean)
        output_file.create_dataset("variance", data=self._output_variance)
        output_file.create_dataset("standard_deviation",
                                   data=self._output_standard_deviation)
        output_file.create_dataset("standard_error_of_mean",
                                   data=self._output_standard_error_of_mean)
        output_file.create_dataset("transformed_input_mean",
                                   data=self._transformed_mean)
        output_file.create_dataset("covariance",
                                   data=self._output_covariance)
        output_file.create_dataset("correlation",
                                   data=self._output_correlation)
        if self.should_store_output_samples:
            output_file.create_group("output_samples")
            i = 0
            for x in self._output_samples:
                output_file["output_samples"].create_dataset(str(i), data=x)
                i += 1
        output_file.close()
        return

    def _resample_and_transform(self, sample):
        # Do a resampling and a bias-corrected transformation of the
        # given sample and return it.
        n = self.n_samples
        f = self.transformation_function
        resample = (n * self._input_mean - sample) / (n - 1)
        return n * self._transformed_mean - (n - 1) * f(resample)

    def _calculate_outer_product(self, array_1, array_2):
        # Calulate the outer product of the given arrays along the
        # `correlation_axis` according to the following formular:
        #   c_{...,i,j,...}=a_{...,i,...}b_{...,j,...}^*,
        axis = self.correlation_axis
        a = np.moveaxis(array_1, axis, -1)
        b = np.moveaxis(array_2, axis, -1)
        result = a[..., :, np.newaxis] * np.conj(b[..., np.newaxis, :])
        result = np.moveaxis(result, -1, axis)
        result = np.moveaxis(result, -1, axis)
        return result

    def _calculate_statistical_quantities(self, shift, sum_, sum_of_squares,
                                          sum_of_outer_products):
        # Calculate the output mean, variance, standard deviation,
        # standard error of the mean, covariance and correlation.
        n = self.n_samples

        # Add the shift back again to get the correct output mean
        self._output_mean = sum_ / n + shift
        # The sample variance is invariant under shifts -> no extra
        # stuff to do
        self._output_variance = sum_of_squares - np.abs(sum_)**2 / n
        self._output_variance /= (n - 1)
        self._output_standard_deviation = np.sqrt(self._output_variance)
        self._output_standard_error_of_mean = \
            self._output_standard_deviation / np.sqrt(n)
        self._output_covariance = sum_of_outer_products \
            - self._calculate_outer_product(sum_, sum_) / n
        self._output_covariance /= (n - 1)
        self._output_correlation = self._output_covariance \
            / self._calculate_outer_product(self._output_standard_deviation,
                                            self._output_standard_deviation)
        return

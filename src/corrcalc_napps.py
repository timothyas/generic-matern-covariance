
from corrcalc import CorrelationCalculator
from timer import Timer


if __name__ == "__main__":

    stdout = "stdout.correlation-timing.1000samples.02apps.log"

    kw = {  'log10tol'      : -3,
            'n_samples'     : 1000,
            'isoxy'         : False,
            'load_samples'  : True,
            'persist'       : False
            }
    localtime = Timer(filename=stdout)
    walltime = Timer(filename=stdout)

    walltime.start("Starting Job")

    for n_range in [5, 10, 15, 20]:
        for n_applications in [2]:

            localtime.start(f"(n_range, n_applications) = ({n_range}, {n_applications})")

            cc = CorrelationCalculator(n_range=n_range, n_applications=n_applications, **kw)
            cc()
            localtime.stop()

    walltime.stop("Total Walltime")

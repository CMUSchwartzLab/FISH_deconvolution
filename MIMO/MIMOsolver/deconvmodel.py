"""Create model for deconv."""
import numpy as np

from pyscipopt import quicksum
import pyscipopt

# Default upper bound for copy numbers
DEFAULT_MAX_COPY_NUMBER = 10


def _create_freq_matrix(scip_model, msamples, kclones):
    """Initialize the frequency variables."""
    freq_matrix = np.empty((msamples, kclones), dtype=object)
    for i in range(msamples):
        for k in range(kclones):
            freq_matrix[i, k] = scip_model.addVar(
                name="f[%d, %d]" % (i + 1, k + 1), ub=1.0
            )
    return freq_matrix


def _create_copy_matrix(scip_model, kclones, nfeatures, max_copy_number):
    """Initialize a matrix of copy number variables."""
    copy_matrix = np.empty((kclones, nfeatures), dtype=object)
    for k in range(kclones):
        for j in range(nfeatures):
            copy_matrix[k, j] = scip_model.addVar(
                name="c[%d, %d]" % (k + 1, j + 1),
                vtype="I",
                ub=max_copy_number,
            )
    return copy_matrix


def _create_clone_half_ploidy(scip_model, kclones):
    """Initialize a vector representing clone ploidy"""
    clone_half_ploidy = np.empty((kclones,), dtype=object)
    for k in range(kclones):
        clone_half_ploidy[k] = scip_model.addVar(name="ploidy[%d]" % (k + 1,))
    return clone_half_ploidy


class DeconvModel:
    """Model for deconvolution of bulk copy-number data.

    Bulk data consists of the observation of 'nfeatures' copy number
    features for 'msamples' samples, and thus is represented by a
    '(msamples x nclones)' matrix 'bulk_samples'.

    Data is to be deconvolved into 'kclones' clones.  The variables in
    this deconvolution are a '(kclones x nfeatures)' copy number
    matrix 'copy_matrix' and a '(msamples x kclones)' frequency matrix
    'freq_matrix'.  Ideally

        freq_matrix * copy_matrix == bulk_samples

    but in most cases there will be some residual.

    """

    def __init__(
        self, name, bulk_data, kclones, max_copy_number=DEFAULT_MAX_COPY_NUMBER
    ):
        """Constructor"""
        msamples, nfeatures = bulk_data.shape
        self.model = pyscipopt.Model(name)
        self.copy_matrix = _create_copy_matrix(
            self.model, kclones, nfeatures, max_copy_number
        )
        self.clone_half_ploidy = _create_clone_half_ploidy(self.model, kclones)
        self.freq_matrix = _create_freq_matrix(self.model, msamples, kclones)

        self._make_freqs_sum_to_one()
        self._minimize_residual_norm(bulk_data)

        self.edge_matrix = None

    @property
    def msamples(self):
        """The number of samples that have been measured."""
        return self.freq_matrix.shape[0]

    @property
    def nfeatures(self):
        """Return the number of features, e.g. genes, measured."""
        return self.copy_matrix.shape[1]  # cols of copy_matrix

    @property
    def kclones(self):
        """Return the number of clones modeled."""
        return self.copy_matrix.shape[0]  # rows of copy_matrix

    def get_edge_matrix(self):
        """Return the phylogenetic tree for the clone structure.

        If a phylogenetic tree for the clone structure has been
        provided, return its adjacency matrix.  Otherwise return
        None.
        """
        if self.edge_matrix is None:
            return None
        return self.edge_matrix.copy()

    def _minimize_residual_norm(self, bulk_data):
        """Add the norm of B - FC to the objective."""
        copy_matrix = self.copy_matrix
        freq_matrix = self.freq_matrix
        for i in range(self.msamples):
            for j in range(self.nfeatures):
                neg = self.model.addVar(
                    name="n[%d, %d]" % (i + 1, j + 1), obj=1.0
                )
                pos = self.model.addVar(
                    name="p[%d, %d]" % (i + 1, j + 1), obj=1.0
                )
                product = quicksum(
                    freq_matrix[i, k] * copy_matrix[k, j]
                    for k in range(self.kclones)
                )
                self.model.addCons(
                    pos - neg + product == bulk_data[i, j],
                    name="resid[%d, %d]" % (i + 1, j + 1),
                )

    def _make_freqs_sum_to_one(self):
        """Require that the columns of freq_matrix sum to one."""
        freq_matrix = self.freq_matrix
        for i in range(self.msamples):
            row_sum = quicksum(freq_matrix[i, k] for k in range(self.kclones))
            self.model.addCons(row_sum == 1, "sum_freq[%d]" % (i + 1,))

    def get_solution_as_vector(self):
        """Record the best solution as a vector of numbers.

        SCIP stores the best solution to the optimization problem
        internally.  Convert this internal format to a vector of
        numbers, and return the result.
        """
        variables = self.model.getVars()
        solution = np.empty((len(variables),))
        for i, value in enumerate(variables):
            solution[i] = self.model.getVal(value)
        return solution

    def set_initial_point(self, initial_point):
        """Set an initial point based on a previously saved solution."""
        initial_solution = self.model.createSol()
        variables = self.model.getVars()
        for k, val in zip(variables, initial_point):
            self.model.setSolVal(initial_solution, k, val)
        self.model.addSol(initial_solution)

    def set_copy_numbers(self, copy_numbers):
        """Set the values of the copy number matrix."""
        initial_point = self.model.createSol()
        for k in range(self.kclones):
            for j in range(self.nfeatures):
                self.model.setSolVal(
                    initial_point, self.copy_matrix[k, j], copy_numbers[k, j]
                )
        self.model.addSol(initial_point)

    def make_freqs_constant(self, clone_freq):
        """Set the frequencies within the SCIP model to constant values."""
        for i in range(self.msamples):
            for k in range(self.kclones):
                self.model.chgVarLb(self.freq_matrix[i, k], clone_freq[i, k])
                self.model.chgVarUb(self.freq_matrix[i, k], clone_freq[i, k])

    def make_ploidy_constant(self, clone_half_ploidy):
        """Set the frequencies within the SCIP model to constant values."""
        for k in range(self.kclones):
            self.model.chgVarLb(
                self.clone_half_ploidy[k], clone_half_ploidy[k]
            )
            self.model.chgVarUb(
                self.clone_half_ploidy[k], clone_half_ploidy[k]
            )

    def make_cpn_constant(self, copy_numbers):
        """Set the copy numbers within the SCIP model to constants values."""
        for k in range(self.kclones):
            for j in range(self.nfeatures):
                self.model.chgVarLb(self.copy_matrix[k, j], copy_numbers[k, j])
                self.model.chgVarUb(self.copy_matrix[k, j], copy_numbers[k, j])

    def read_copy_numbers(self):
        """Read the copy number data from the best SCIP solution."""
        copy_numbers = np.empty((self.kclones, self.nfeatures))
        for k in range(self.kclones):
            for j in range(self.nfeatures):
                copy_numbers[k, j] = self.model.getVal(self.copy_matrix[k, j])
        return copy_numbers

    def read_clone_freqs(self):
        """Return the clone frequencies in the best solution."""
        clone_freqs = np.zeros((self.msamples, self.kclones))
        for i in range(self.msamples):
            for k in range(self.kclones):
                clone_freqs[i, k] = self.model.getVal(self.freq_matrix[i, k])
        return clone_freqs

    def read_half_ploidy(self):
        """Return the clone frequencies in the best solution."""
        clone_half_ploidy = np.zeros((self.kclones,))
        for k in range(self.kclones):
            clone_half_ploidy[k] = self.model.getVal(self.clone_half_ploidy[k])
        return clone_half_ploidy

    def optimize(self):
        """Call SCIP to optimize the model."""
        self.model.optimize()
        return self.model.getObjVal()

    def free_scip_data(self):
        """Free SCIP understanding of the problem and solution."""
        self.model.freeProb()

    def set_time_limit(self, max_time):
        """Set the solvers time limit"""
        self.model.setRealParam("limits/time", max_time)

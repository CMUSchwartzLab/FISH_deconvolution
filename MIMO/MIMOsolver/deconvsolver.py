"""A solver for the copy number deconvolution problem."""
from math import inf
import pickle

import numpy as np

from deconvmodel import DeconvModel
import treemodel
from treemodel import TreeModel


def initialize_clone_freqs(msamples, kclones, clone_decay):
    """Create initial values for the clonal frequencies.

    Use this function when there is no informed initial guess for the
    clone frequencies.  The function ensures that the frequency matrix
    is nonsigular.
    """
    clone_freqs = np.ones((kclones,))
    for k in range(1, kclones):
        clone_freqs[k] = clone_decay * clone_freqs[k - 1]
    freq_matrix = np.zeros((msamples, kclones))
    for i in range(msamples):
        freq_matrix[i] = clone_freqs.copy()
        if i < kclones:
            freq_matrix[i, i] += 1
        freq_matrix[i] /= freq_matrix[i].sum()

    return freq_matrix


def _initial_copy_numbers(kclones, bulk_data):
    """Create an initial copy-number matrix based on bulk data."""
    msamples, nfeatures = bulk_data.shape
    copy_numbers = np.zeros((kclones, nfeatures))
    for k in range(kclones):
        for j in range(nfeatures):
            copy_numbers[k, j] = int(round(bulk_data[(k % msamples), j]))
    return copy_numbers


class DeconvSolver:
    """A solver for the copy-number deconvolution proble."""

    def __init__(
        self,
        freq_regularization=None,
        fish_regularization=None,
        copy_regularization=None,
        tree_regularization=None,
    ):
        self.copy_regularization = None
        self.freq_regularization = None
        self.tree_regularization = None
        self.fish_regularization = None

        if copy_regularization and copy_regularization.weight > 0:
            self.copy_regularization = copy_regularization
        if freq_regularization and freq_regularization.weight > 0:
            self.freq_regularization = freq_regularization
        if tree_regularization and tree_regularization.weight > 0:
            self.tree_regularization = tree_regularization
        if fish_regularization and fish_regularization.weight > 0:
            self.fish_regularization = fish_regularization

    def _create_model(self, bulk_data, edge_matrix, kclones):
        """Create the SCIP model for these data."""

        bulk_data = np.asarray(bulk_data)

        model = DeconvModel("Example", bulk_data, kclones)

        if self.copy_regularization is not None:
            self.copy_regularization.regularize(model)

        if self.fish_regularization is not None:
            self.fish_regularization.regularize(model)

        if self.freq_regularization is not None:
            self.freq_regularization.regularize(model)

        if self.tree_regularization is not None:
            model.edge_matrix = edge_matrix.copy()
            self.tree_regularization.regularize(model)

        return model

    def _create_phylogenetic_model(self, copy_numbers):
        """Create a model of the phylogenetic tree constraint."""

        fixed_nodes = self.tree_regularization.fixed_nodes

        phylogenetic_model = TreeModel(
            name="phylogenetic_model",
            edge_weights=treemodel.calc_cpn_distance(
                copy_numbers, fixed_nodes
            ),
        )

        return phylogenetic_model

    def _validate(self, bulk_data, initial_freqs):
        """Run basic validation tests on the input data

        Currently makes sure the data has the right shape.
        """
        msamples, nfeatures = bulk_data.shape
        assert msamples > 0 and nfeatures > 0
        assert msamples == initial_freqs.shape[0]

        kclones = initial_freqs.shape[1]
        assert kclones > 0
        if self.copy_regularization:
            k, n = self.copy_regularization.reference_values.shape
            assert k == kclones and n == nfeatures
        if self.freq_regularization:
            m, k = self.freq_regularization.reference_values.shape
            assert m == msamples and k == kclones
        if self.tree_regularization:
            assert nfeatures == self.tree_regularization.fixed_nodes.shape[1]

    def optimize(
        self, bulk_data, initial_freqs, max_iteration=10, solution_tol=1e-4
    ):
        """Create and solve the copy-number deconvolution problem."""
        bulk_data = np.asarray(bulk_data)
        clone_freqs = np.asarray(initial_freqs).copy()
        self._validate(bulk_data, initial_freqs)

        saved_result = None
        objective = inf
        kclones = initial_freqs.shape[1]
        edge_matrix = np.zeros((2 * kclones + 1, 2 * kclones + 1))
        clone_half_ploidy = np.ones((kclones,))
        for iteration_count in range(1, max_iteration + 1):
            model = self._create_model(bulk_data, edge_matrix, kclones)
            model.make_freqs_constant(clone_freqs)
            model.make_ploidy_constant(clone_half_ploidy)

            if saved_result is not None:
                print("====", "Setting initial solution", "====", flush=True)
                model.set_initial_point(saved_result)
            else:
                pass
                # Do not try to initialize the copy number, because this
                # results in an infeasible starting point, which is not used.
                # copy_numbers = _initial_copy_numbers(model.kclones, bulk_data)
                # model.set_copy_numbers(copy_numbers)

            print("====", "Optimizing copy numbers", "====", flush=True)
            model.optimize()

            copy_numbers = np.rint(model.read_copy_numbers())

            model.free_scip_data()
            model = None

            if (
                self.tree_regularization
                and self.tree_regularization.weight > 0
            ):
                model = self._create_phylogenetic_model(copy_numbers)
                print("====", "Optimizing pylogeny", "====", flush=True)
                model.optimize()
                edge_matrix = np.rint(model.read_edges())
                model.free_scip_data()
                model = None

            model = self._create_model(bulk_data, edge_matrix, kclones)
            model.make_cpn_constant(copy_numbers)

            print("====", "Optimizing frequencies", "====", flush=True)
            new_objective = model.optimize()

            clone_freqs = model.read_clone_freqs()
            clone_half_ploidy = model.read_half_ploidy()

            if abs(objective - new_objective) <= solution_tol:
                return model, iteration_count

            if iteration_count == max_iteration:
                return model, max_iteration + 1

            objective = new_objective
            saved_result = model.get_solution_as_vector()
            model.free_scip_data()

        return None, max_iteration


class Simulation:
    """Instances of Simulation represent simulated tests of DeconvSolver."""

    def __init__(self, data):
        """Constructor"""
        self.data = data

    @classmethod
    def read(cls, filename):
        """Read a simulation from the named pickle file."""
        with open(filename, "rb") as fileh:
            return cls(pickle.load(fileh))

    @property
    def bulk_data(self):
        """Bulk tumor data"""
        return self.data["TumorSample"].T

    @property
    def reference_frequencies(self):
        """Frequencies for frequency normalization"""
        return self.data["FRefer"].T

    @property
    def reference_cells(self):
        """Counts for copy-number regularization"""
        return self.data["CRefer"].T

    @property
    def reference_fish(self):
        """FISH counts for FISH regularization"""
        return np.delete(self.data["FISHRefer"], -1, axis=0).T

    @property
    def reference_fish_cell_counts(self):
        return self.data["FISHRefer"][-1].tolist()

    @property
    def true_cells(self):
        """Major cells simulated bulk, used to calculate correctness"""
        return self.data["CTrue"]

    @property
    def true_frequencies(self):
        """True frequency of major cells in simulated bulk

        Used to calculate correctness
        """
        return self.data["FTrue"]

    @property
    def reference_cells_ploidy(self):
        """Ploidy for reference cell used as rescaling factor"""
        return self.data.get("ReferPloidy")

    @property
    def true_cells_ploidy(self):
        """Ploidy for major cells in the simulated bulk

        Used as rescaling factor.
        """
        return self.data.get("TruePloidy")

    @property
    def reference_fish_ploidy(self):
        """ploidy for FISH"""
        return self.data.get("referFISHploidy")

    @property
    def expand_FISH(self):
        """ploidy for FISH"""
        return self.data["EnrichReferFISH"].T

    @property
    def correlate_pos(self):
        """get the correlated position for each probes"""
        return self.data.get("correlatePositions")
    
    @property
    def expand_uFISH(self):
            """ploidy for FISH"""
            return self.data["uEnrichReferFISH"].T


class MarkerPositions:
    """Positions of the FISH and copy number markers"""

    def __init__(self, data):
        self.data = data

    @classmethod
    def read(cls, filename):
        """Read positions from a pickle file."""
        with open(filename, "rb") as fileh:
            data = pickle.load(fileh)
            data["FISHInterval"].set_index("probe")
            del data["OriginalSCSInterval"]
            return cls(data)

    @property
    def fish_probes(self):
        """Genetic interval of the gene containing the FISH probe"""
        return self.data["FISHInterval"]

    @property
    def copy_numbers(self):
        """The positions of the copy number intervals"""
        return self.data["CompressedSCSInterval"]





class Real:
    """Instances of Real Data represent simulated tests of DeconvSolver."""

    def __init__(self, data):
        """Constructor"""
        self.data = data

    @classmethod
    def read(cls, filename):
        """Read a simulation from the named pickle file."""
        with open(filename, "rb") as fileh:
            return cls(pickle.load(fileh))

    @property
    def bulk_data(self):
        """Bulk tumor data"""
        return self.data["TumorSample"].T

    @property
    def reference_frequencies(self):
        """Frequencies for frequency normalization"""
        return self.data["FRefer"].T

    @property
    def reference_cells(self):
        """Counts for copy-number regularization"""
        return self.data["CRefer"].T

    @property
    def reference_fish_ploidy(self):
        """reference ploidy from FISH"""
        return self.data.get("referFISHploidy")

    @property
    def expand_uFISH(self):
            """ploidy for FISH"""
            return self.data["uEnrichReferFISH"].T

    @property
    def correlate_pos(self):
        """get the correlated position for each probes"""
        return self.data.get("correlatePositions")
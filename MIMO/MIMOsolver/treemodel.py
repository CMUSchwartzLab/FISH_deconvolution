"""A phylogentic tree model for copy-number deconvolution"""

import numpy as np

import pyscipopt
from pyscipopt import quicksum


def _create_edge_matrix(model, edge_weights):
    """Create a matrix of edges with the specified weights."""
    nnodes = edge_weights.shape[0]
    edge_matrix = np.empty((nnodes, nnodes), dtype=object)
    for u in range(nnodes):
        for v in range(nnodes):
            # The objective will be to minimize the sum of edge_weights
            # over the edges in the tree.
            edge_matrix[u, v] = model.addVar(
                name="edge[%d, %d]" % (u + 1, v + 1),
                vtype="B",
                obj=edge_weights[u, v],
            )
    return edge_matrix


def _create_flow_variables(model, nnodes):
    """Create the variables representing flow.

    The element flow[u, v, t] represents the flow directed towards t
    flowing through edge[u,v].
    """
    flow = np.empty(shape=(nnodes, nnodes, nnodes), dtype=object)
    for u in range(nnodes):
        for v in range(nnodes):
            for t in range(nnodes):
                flow[u, v, t] = model.addVar(
                    name="flow[%d, %d, %d]" % (u + 1, v + 1, t + 1), vtype="B"
                )
    return flow


class TreeModel:
    """A weighted phylogentic tree model for copy-number deconvolution"""

    def __init__(self, name, edge_weights, tree_root=0):
        self.tree_root = tree_root
        self.model = pyscipopt.Model(name)
        self.edge_weights = np.copy(edge_weights)
        self.edge_matrix = _create_edge_matrix(self.model, edge_weights)
        self.flow = _create_flow_variables(self.model, self.nnodes)

        self._no_self_edges()
        self._add_flow_source()
        self._add_flow_sinks()
        self._make_flow_conserved_at_nodes()
        self._no_flow_through_unselected_edges()
        self._ignore_unweighted_edges()

    @property
    def nnodes(self):
        """The number of nodes in the tree"""
        return self.edge_matrix.shape[0]

    def _no_self_edges(self):
        """There is no edge from a node to itself -- constraint 16."""
        for u in range(self.nnodes):
            self.model.chgVarUb(self.edge_matrix[u, u], 0)
            # As a consequence, there is no flow from a node to itself.
            for t in range(self.nnodes):
                self.model.chgVarUb(self.flow[u, u, t], 0)

    def _add_flow_source(self):
        """For each node in the tree, one unit of flow originates at root.

        Constraint 15.
        """
        nnodes = self.nnodes
        tree_root = self.tree_root
        flow = self.flow
        for t in range(nnodes):
            # Neither edges, nor flow, enter root.
            self.model.chgVarUb(self.edge_matrix[t, tree_root], 0)
            for u in range(nnodes):
                self.model.chgVarUb(self.flow[u, tree_root, t], 0)
            if t != tree_root:
                # For each terminal t (except root), one unit of flow
                # directed towards t originates at root.
                source = quicksum(flow[tree_root, u, t] for u in range(nnodes))
                self.model.addCons(source == 1, name="source[%s]" % (t + 1,))

    def _add_flow_sinks(self):
        """One unit of flow terminates in each node in the tree.

        Constraint 17.
        """
        nnodes = self.nnodes
        flow = self.flow
        for t in range(nnodes):
            if t != self.tree_root:
                # One unit of flow directed towards t enters t, except
                # when t is root (in which case no flow enters).
                sink = quicksum(flow[u, t, t] for u in range(nnodes))
                self.model.addCons(sink == 1, name="sink[%s]" % (t + 1,))
            for u in range(nnodes):
                # No flow towards t exits t (even if t is root).
                self.model.chgVarUb(flow[t, u, t], 0)

    def _make_flow_conserved_at_nodes(self):
        """Flow at a node is conserved, unless the node is a source or sink.

        The amount of flow directed towards terminal t that enters a
        node u must equal the amount of flow directed towards t that
        exits node u, unless node u is the source of flow (the root
        node) or a sink (the terminal t itself).

        Constraint 14.
        """
        nnodes = self.nnodes
        flow = self.flow
        for t in range(nnodes):
            for u in range(nnodes):
                if u not in [self.tree_root, t]:
                    # The vertex is not a source or sink
                    expr = quicksum(
                        flow[u, v, t] for v in range(nnodes)
                    ) - quicksum(flow[v, u, t] for v in range(nnodes))
                    self.model.addCons(
                        expr == 0, name="conserv[%d, %d]" % (u + 1, t + 1)
                    )

    def _no_flow_through_unselected_edges(self):
        """There can be flow only through edges that are part of the tree.

        Constraint 18."""
        nnodes = self.nnodes
        flow = self.flow
        edge = self.edge_matrix
        for u in range(nnodes):
            for v in range(nnodes):
                for t in range(nnodes):
                    self.model.addCons(
                        edge[u, v] >= flow[u, v, t],
                        name="present[%s, %s, %s]" % (u + 1, v + 1, t + 1),
                    )

    def _ignore_unweighted_edges(self):
        """Set edges having a weight of zero to a fixed, harmless, value.

        For most edges that have weight zero, setting the edge to 1
        has no affect on the solution of the optimization problem.
        The edge variables that cannot be set to 1 are those that must
        be set to 0 -- self edges and edges entering the tree root.

        The method self.read_edges() reads the flow variables to find
        appropriate values for the ignored edges.
        """
        for u in range(self.nnodes):
            for v in range(self.nnodes):
                if (
                    u != v
                    and v != self.tree_root
                    and self.edge_weights[u, v] == 0
                ):
                    self.model.chgVarLb(self.edge_matrix[u, v], 1)

    def read_edges(self):
        """Read the edge data from the best solution."""
        edges = np.empty((self.nnodes, self.nnodes))
        for u in range(self.nnodes):
            for v in range(self.nnodes):
                if self.edge_weights[u, v] != 0:
                    edges[u, v] = self.model.getVal(self.edge_matrix[u, v])
                else:
                    # There was no weight on the edge, so the value of
                    # the edge_matrix variable stored in the solution
                    # is meaningless.  Find a good edge assignment by
                    # looking at the flow.
                    edges[u, v] = int(
                        any(
                            round(self.model.getVal(self.flow[u, v, t]))
                            for t in range(self.nnodes)
                        )
                    )
        return edges

    def optimize(self):
        """Call SCIP to optimize the model."""
        self.model.optimize()

    def free_scip_data(self):
        """Free SCIP's understanding of the problem and solution."""
        self.model.freeProb()


def calc_cpn_distance(cpn, fixed_nodes):
    """Calculate the L1 distance between each pair of rows in [cpn;
    fixed_nodes].
    """
    kclones, _ = cpn.shape
    nnodes = 2 * kclones + 1
    edge_weights = np.empty((nnodes, nnodes))
    for u in range(nnodes):
        for v in range(nnodes):
            if u < kclones:
                c_u = cpn[u, :]
            else:
                c_u = fixed_nodes[u - kclones, :]
            if v < kclones:
                c_v = cpn[v, :]
            else:
                c_v = fixed_nodes[v - kclones, :]

            edge_weights[u, v] = np.linalg.norm(c_u - c_v, ord=1)
    return edge_weights

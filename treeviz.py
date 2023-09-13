import math

import msprime
import tskit
import yaml
from graphviz import Digraph


def vars():
    with open("_variables.yml") as fh:
        try:
            data = yaml.safe_load(fh)
        except yaml.YAMLError as e:
            print(e)
            raise
    return data


##############################
# 3d-effect tree example from
# https://tskit.dev/tutorials/viz.html#d-effects
##############################
def make_7_tree_4_tip_ts():
    ts = msprime.sim_ancestry(
        4, ploidy=1, random_seed=889, sequence_length=1000, recombination_rate=0.001
    )
    ts = msprime.sim_mutations(ts, rate=2e-3, random_seed=123)

    # Check we have picked a random seed that gives a nice plot of 7 trees
    tip_orders = {
        tuple(u for u in t.nodes(order="minlex_postorder") if t.is_sample(u))
        for t in ts.trees()
    }
    topologies = {tree.rank() for tree in ts.trees()}
    assert tip_orders == {(0, 1, 2, 3)} and len(topologies) > 1 and ts.num_trees == 7

    return ts


def make_3d_tree(*, style=None, tree_width=100, y_step=40, lmargin=20, rmargin=20):
    ts = make_7_tree_4_tip_ts()

    # Set some parameters: these can be adjusted to your liking
    tree_width = tree_width
    height = 200  # Normal height for tree + x-axis
    y_step = y_step  # Stagger between trees (i.e. 0 for all trees in a horizontal line)
    skew = 0.6  # How skewed the trees are, in radians

    width = (
        tree_width * ts.num_trees + lmargin + rmargin
    )  # L & R margins in draw_svg = 20px
    angle = math.atan(y_step / tree_width)
    ax_mv = y_step, (ts.num_trees - 1) * y_step + math.tan(skew) * (tree_width * 0.9)

    if style is None:
        # CSS transforms used to skew the axis and stagger + skew the
        # trees
        style = (
            f".x-axis {{ transform: translate({ax_mv[0]}px,"
            f"{ax_mv[1]}px) skewY(-{angle}rad)}}"
        )
        for i in range(ts.num_trees):
            # Stagger each tree vertically by y_step, transforming the
            # "plotbox" tree container
            style += (
                f".tree.t{i} > .plotbox "
                + "{ transform:"
                + f"translateY({(ts.num_trees - i - 1) * y_step}px) skewY({skew}rad)"
                + "}"
            )

        style = "#threedtree { " + style + "}"
    # Define a bigger canvas size so we don't crop the moved trees from the drawing
    size = (width, height)
    canvas_size = (
        width + y_step,
        height + ts.num_trees * y_step + math.tan(skew) * tree_width,
    )

    return ts.draw_svg(
        size=size,
        x_scale="treewise",
        style=style,
        canvas_size=canvas_size,
        root_svg_attributes={"id": "threedtree"},
    )


##############################
# cf https://github.com/tskit-dev/tutorials/issues/43
##############################
def ts_to_arg_dot(ts, *, size="4,6"):
    """Convert ts to arg image. Code from  https://github.com/tskit-dev/tutorials/issues/43"""
    dot = Digraph(strict=False)
    dot.attr(size=size)

    # Add sample nodes
    with dot.subgraph() as s:
        s.attr(rank="same")
        for n in ts.samples():
            s.node(str(n), **{"shape": "doublecircle"})

    # Internal nodes
    itr = iter(range(ts.num_samples, ts.num_nodes))
    for i in itr:
        if ts.node(i).flags == msprime.NODE_IS_RE_EVENT:
            # Only add one of the recombination nodes and make the other one invisible
            with dot.subgraph() as s:
                s.attr(rank="same")
                s.node(str(i), **{"shape": "rect"}, label=str(i) + "/" + str(i + 1))
                s.node(str(i + 1), **{"style": "invis"})
                next(itr, None)
        else:
            dot.node(str(i), **{"shape": "circle"})
        # Add invisible edges to fix the ranks
        dot.edge(str(i), str(i - 1), **{"style": "invis"})

    # Add edges
    check = [0, 0]
    for e in ts.edges():
        ch = e.child
        pa = e.parent

        # Preventing duplicate edges
        if check == [pa, ch]:
            continue
        check = [pa, ch]

        # Check which edge to draw if there is a recombination
        if (
            e.parent >= 1
            and ts.node(pa).time == ts.node(pa - 1).time
            and ts.node(pa).flags == ts.node(pa - 1).flags == msprime.NODE_IS_RE_EVENT
        ):
            ch = pa = -1
        elif (
            e.child >= 1
            and ts.node(ch).time == ts.node(ch - 1).time
            and ts.node(ch).flags == ts.node(ch - 1).flags == msprime.NODE_IS_RE_EVENT
        ):
            ch -= 1

        lab = ""
        if ch >= 0 and ts.node(ch).flags == msprime.NODE_IS_RE_EVENT:
            lab = "[" + str(int(e.left)) + "," + str(int(e.right)) + ")"

        if ch >= 0 and pa >= 0:
            dot.edge(str(pa), str(ch), label=lab)

    return dot.render("arg", format="svg")


##############################
# Custom functions
##############################
def recombination_trees():
    """Plot multiple trees with colored egdes to highlight recombination"""
    # 4-color viridis: "#440154FF" "#31688EFF" "#35B779FF" "#FDE725FF"
    css_style = "#recombination { .edge {stroke-width: 5px}"
    css_style += " .n0 .edge {stroke: #31688EFF}"
    css_style += " .n1 .edge {stroke: #31688EFF}"
    css_style += " .n6 .edge {stroke: #31688EFF}"
    css_style += " .n10 .edge {stroke: #31688EFF}"
    css_style += " .n11 .edge {stroke: #31688EFF}"
    css_style += " .n2 .edge {stroke: #35B779FF}"
    css_style += " .n7 .edge {stroke: #35B779FF}"
    css_style += " .n3 .edge {stroke: #35B779FF}"
    css_style += " .n5 .edge {stroke: #35B779FF}"
    css_style += " .n8 .edge {stroke: #35B779FF}"
    css_style += " .n9 .edge {stroke: #35B779FF}"
    css_style += " .n4 .edge {stroke: #FDE725FF}"
    css_style += "}"
    kwargs = dict(
        size=(1000, 300),
        x_axis=False,
        node_labels={},
        symbol_size=0,
        style=css_style,
        root_svg_attributes={"id": "recombination"},
    )
    ts = msprime.sim_ancestry(
        sequence_length=1e6, recombination_rate=1e-6, samples=4, random_seed=19
    )
    return ts.draw_svg(**kwargs)


def tree_topology(
    model,
    *,
    svgid,
    size=(300, 500),
    x_axis=False,
    node_labels={},
    symbol_size=0,
    style=".edge {stroke-width: 2px}",
):
    """Plot tree topology under different evolutionary scenarios"""
    kwargs = dict(
        size=size,
        x_axis=x_axis,
        node_labels=node_labels,
        symbol_size=symbol_size,
        style=f"#{svgid} {{ {style} }}",
        root_svg_attributes={"id": svgid},
    )
    if model == "neutral":
        ts = msprime.sim_ancestry(10, random_seed=12)
    elif model == "expansion":
        demography = msprime.Demography()
        demography.add_population(name="A", initial_size=10_000, growth_rate=0.1)
        ts = msprime.sim_ancestry(
            samples={"A": 10}, demography=demography, random_seed=12
        )
    elif model == "bottleneck":
        demography = msprime.Demography()
        demography.add_population(name="A", initial_size=1_000)
        demography.add_instantaneous_bottleneck(time=100, strength=1000, population=0)
        ts = msprime.sim_ancestry(
            samples={"A": 10}, demography=demography, random_seed=12
        )
    elif model == "selection":
        Ne = 1_000
        L = 1e6
        sweep_model = msprime.SweepGenicSelection(
            position=L / 2,
            start_frequency=1.0 / (2 * Ne),
            end_frequency=1.0 - (1.0 / (2 * Ne)),
            s=0.25,
            dt=1e-6,
        )
        ts = msprime.sim_ancestry(
            10,
            model=[sweep_model, msprime.StandardCoalescent()],
            population_size=Ne,
            sequence_length=L,
            random_seed=119,
        )
    return ts.draw_svg(**kwargs)


def basics_tree():
    ts = tskit.load("data/basics.trees")
    ts = ts.delete_sites(0)
    return ts


def treemut():
    ts = basics_tree()
    tsmut = msprime.sim_mutations(ts, rate=1e-5, random_seed=228)
    return tsmut

#!/usr/bin/env python

from sys import argv, stderr, exit

try:
    from igraph import *
except ImportError:
    raise "Please install the Python IGraph library: `pip install python-igraph`"

def main():
    if len(argv) != 2:
        print >> stderr, "Usage: ./visualize_segmented_trace_graph.py graph.gml"
        exit(1)

    graph_fname = argv[1]
    assert graph_fname.endswith(".gml"), "This script only plots GML files (for now)"

    st_graph = Graph.Read_GML(graph_fname)

    visual_style = {
        "vertex_label": [ "(%i, %i) (%i, %i)" %
            (v["track"], v["chr"], v["start"], v["stop"])
            for v in st_graph.vs ],
        # python-igraph does not support labelling edges!!!
        # they recommend changing the width of the edge to represent weight
        # http://lists.gnu.org/archive/html/igraph-help/2011-01/msg00020.html
        "vertex_label_dist": 2,
        "layout": st_graph.layout("fr"),
        "bbox": (800, 800),
        "margin": 50
    }

    plot(st_graph, **visual_style)

if __name__ == '__main__':
    main()

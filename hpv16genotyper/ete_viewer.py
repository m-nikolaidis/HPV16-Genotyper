import re
import sys

from ete3 import Tree, TreeStyle, TextFace


def main():
    tree_path = sys.argv[1]
    selected_seq = sys.argv[2]
    gene = sys.argv[3]

    t = Tree(tree_path)
    t.ladderize(direction=1)

    ts = TreeStyle()
    ts.title.add_face(
        TextFace(gene + " Gene", fsize=13),
        column=1,
    )
    ts.show_branch_support = True
    ts.show_leaf_name = False

    for leaf in t.iter_leaves():
        color = "Black"

        if leaf.name == selected_seq:
            color = "Gray"
        elif re.match(r"^A\d+_\S\d", leaf.name):
            color = "Green"
        elif re.match(r"^B\d+_\S\d", leaf.name):
            color = "SteelBlue"
        elif re.match(r"^C\d+_\S\d", leaf.name):
            color = "Orange"
        elif re.match(r"^D\d+_\S\d", leaf.name):
            color = "FireBrick"

        leaf.add_face(
            TextFace(leaf.name, fgcolor=color),
            column=0,
            position="branch-right",
        )

    # This is now correct: ETE owns this process/event loop.
    t.show(tree_style=ts)


if __name__ == "__main__":
    main()

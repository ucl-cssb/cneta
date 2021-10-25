#!/usr/bin/env python

'''
This script is used to convert NEWICK trees into the format that can be read by program svtree*.
The naming coventions for trees accepted by svtree*:
Suppose there are Ns disease samples, named from 1 to Ns, and one normal sample named by Ns + 1.
The root is named Ns + 2.
The internal nodes are named from Ns + 3 to Ns + Ns + 1.
Two additional edges (Ns+2, Ns+1) and (Ns+2, Ns+3) are added to the original tree due to the addition of a node representing unaltered genome.
The length of edge (Ns+2, Ns+1) is 0.
The length of edge (Ns+2, 2*Ns+1) has to be specified by the user.

Sample command to run the program:
To convert a NEXUS tree file (which may contains multiple trees) to a list of TXT trees (in the format of edge lists):
    python newick2elist.py -f 0 -t AllTreesNr5.txt -b 5.25
To convert a NEWICK tree file (which only contains one tree) to a TXT tree:
    python newick2elist.py -f 1 -t tree1.nwk

'''

import networkx as nx
import os
import sys
import re
import argparse


class NewickError(Exception):
    pass


def write_edge_list(tree_file, tree_renamed, nleaf, scaling, i = -1):
    dir = os.path.dirname(tree_file)
    bname = os.path.basename(tree_file)
    if i > 0:
        fname = os.path.splitext(bname)[0] + "_" + str(i) + ".txt"
    else:
        fname = os.path.splitext(bname)[0] + ".txt"
    fname_full = os.path.join(dir, fname)
    # print("write edge list to file {}".format(fname_full))

    header="start\tend\tlength\teid\n"
    with open(fname_full, "w") as fout:
        fout.write(header)
        j = 1
        for s, e in tree_renamed.edges:
            l = tree_renamed.nodes[e]['length']
            if scaling > 0:
                l = l / scaling
            if e != nleaf and l <= 0:  # avoid length 0 in branches other than that leads to normal node
                l = 1
            line = "\t".join([str(s), str(e), str(l), str(j)])
            line += "\n"
            fout.write(line)
            j += 1


'''
build a tree from newick file without using other libraries
'''
def newick_to_tree(newick):
    tree = nx.DiGraph()
    tokens = [
    (r"\(",                                  'open parens'),
    (r"\)",                                  'close parens'),
    (r"[^\s\(\)\[\]\'\:\;\,]+",              'unquoted node label'),
    (r"\:[0-9]*\.?[0-9]+([eE][+-]?[0-9]+)?", 'edge length'),
    (r"\,",                                  'comma'),
    (r"\[(\\.|[^\]])*\]",                    'comment'),
    (r"\'(\\.|[^\'])*\'",                    'quoted node label'),
    (r"\;",                                  'semicolon'),
    (r"\n",                                  'newline'),
    ]
    tokenizer = re.compile('(%s)' % '|'.join(token[0] for token in tokens))
    tokens = re.finditer(tokenizer, newick.strip())
    parents = {} # to record the parent of current node which is far away
    prev_node = ''
    current_node = ''
    lp_count = 0
    rp_count = 0

    for match in tokens:
        token = match.group()
        # print 'token', token
#         print 'prev_node', prev_node
#         if current_node!='':
#             print 'current node %s' % (current_node)
#         if current_node in parents:
#             print 'parent node %s' % (parents[current_node])
        if token.startswith("'"):
                # quoted label; add characters to clade name
            prev_node = current_node
            current_node = token[1:-1]
            tree.add_node(current_node)
            parents[current_node] = prev_node

        elif token.startswith('['):
                # comment
                pass
#                 comment = token[1:-1]
#                 tree[current_node]['comment'] = comment
#                 if comments_are_confidence:
#                     # Try to use this comment as a numeric support value
#                     tree[current_node]['confidence'] = _parse_confidence(comment)

        elif token == '(':
            # start a new clade, which is a child of the current clade
            # current_clade = new_clade(current_clade)
            lp_count += 1
            prev_node = current_node
            # It will ignore current name of internal nodes
            current_node = 't' + str(lp_count)
            tree.add_node(current_node)
            if prev_node != '':
                parents[current_node] = prev_node

        elif token == ',':
#                 # if the current clade is the root, then the external parentheses
#                 # are missing and a new root should be created
#                 if current_clade is root_clade:
#                     root_clade = new_clade()
#                     current_clade.parent = root_clade
                # start a new child clade at the same level as the current clade
#                 if current_node not in tree:
#                     tree.add_node(current_node)
                tree.add_edge(parents[current_node], current_node)
                current_node = parents[current_node]

        elif token == ')':
                # done adding children for this parent clade
                if current_node not in parents:
                    raise NewickError('Parenthesis mismatch.')
#                     print 'Parenthesis mismatch.'
#                     break
                rp_count += 1
                next_token = next(tokens)
                nval = next_token.group()
                # print("next token after ) is {}".format(nval))
                # Skip internal node name
                if not nval.startswith(":"):
                    pass
#                 if current_node not in tree:
#                     tree.add_node(current_node)
                tree.add_edge(parents[current_node], current_node)
                current_node = parents[current_node]
                if nval.startswith(":"):
                    value = float(nval[1:])
                    # print("Add length {} to {}".format(value, current_node))
                    tree.nodes[current_node]['length'] = value

        elif token == ';':
                break
        elif token.startswith(':'):
                # branch length or confidence
                value = float(token[1:])
                # print("Add length {} to {}".format(value, current_node))
                tree.nodes[current_node]['length'] = value
                # tree.edges[parents[current_node], current_node]['length'] = value
        elif token == '\n':
                pass
        else:
                # unquoted node label
                prev_node = current_node
                if token.startswith('t'):
                    current_node = token[1:]
                else:
                    current_node = token
                tree.add_node(current_node)
                parents[current_node] = prev_node

    if not lp_count == rp_count:
            raise NewickError('Number of open/close parentheses do not match.')
    try:
            next_token = next(tokens)
            raise NewickError('Text after semicolon in Newick tree: %s'
                              % next_token.group())
    except StopIteration:
            pass

    return tree


def get_num_leaf(tree):
    num_leaf = 0
    for n in tree.nodes:
        if tree.out_degree(n) == 0:
            num_leaf += 1
    return num_leaf


def get_root(tree):
    root = -1
    for n in tree.nodes:
        if tree.in_degree(n) == 0 and tree.out_degree(n) > 0:
            root = n
            break
    return root



# Rename all node names to integer, order as that is defined in sveta
# num_leaf is the number of leaves in the final tree (plus one when the normal sample is not included)
def rename_tree(tree, num_leaf, incl_normal, node_prefix=""):
    mapping = {}
    # print(tree.edges())
    k = 0
    for n in nx.dfs_postorder_nodes(tree):
        n = str(n)
        # print("visit node {}".format(n))
        # print(tree.nodes[e]['length']
        if n.isdigit() and int(n) <= num_leaf:
            mapping[n] = int(n)

        if n.startswith("t"):   # Internal nodes got names bottom up with prefix t when traversing the tree
            k = k + 1
            if k == num_leaf - 1 and incl_normal:   # root when the normal sample is originally included in the tree
                # print("root node is {}".format(n))
                mapping[n] = num_leaf + 1
            else:
                mapping[n] = k + num_leaf + 1
            # print("{}th node {} mapped to {}".format(k, n, mapping[n]))

        if node_prefix != '' and n.startswith(node_prefix):
            new_name = n.replace(node_prefix, "")
            new_name = new_name.lstrip("0")
            mapping[n] = new_name

    # print(mapping)

    tree_renamed = nx.relabel_nodes(tree, mapping)

    return tree_renamed


# Convert the nexus trees to edge list as it is
def convert_nexus_file(tree_file, node_prefix, incl_normal, scaling, root_blen = 0):
    with open(tree_file, "r") as fin:
        dir = os.path.dirname(tree_file)
        i = 0 # TREE ID
        for line in fin:
            line = line.strip()
            if not (line.startswith("tree") or line.startswith("TREE")):
                continue
            # print(line)
            line = line.split("=")[-1]
            newick = line.strip()
            # print("{}th tree {}".format(i, newick))

            tree_renamed, nleaf = newick_to_edgelist(newick, node_prefix, incl_normal, root_blen)

            # print(tree)
            i += 1
            write_edge_list(tree_file, tree_renamed, nleaf, scaling, i)



def get_newick_from_file(tree_file):
    with open(tree_file, "r") as fin:
        lines = fin.readlines()
        assert(len(lines) == 1)
        newick = lines[0].strip()

    return newick


def newick_to_edgelist(newick, node_prefix, incl_normal, root_blen = 0):
    tree = newick_to_tree(newick)
    # print("edges in original tree: ")
    # print(tree.edges())

    # Find the root of the tree
    num_leaf = get_num_leaf(tree)
    root = get_root(tree)
    if not incl_normal:
        print("adding normal edge to root {}".format(root))
        s = num_leaf + 2  # new root
        tree.add_edge(s, root)
        tree.nodes[root]['length'] = root_blen
        # print("blen above {} is {}".format(root, tree.nodes[root]['length']))
        e2 = num_leaf + 1 # normal sample
        tree.add_edge(s, e2)
        tree.nodes[e2]['length'] = 0
        num_leaf = num_leaf + 1

    # print("edges in original tree after adding normal edges when required: ")
    # print(tree.edges())

    tree_renamed = rename_tree(tree, num_leaf, incl_normal, node_prefix)

    # print("edges in renamed tree: ")
    # print(tree_renamed.edges())

    if incl_normal:
        root_id = num_leaf + 1
        norm_id = num_leaf
        # print("root {}, normal {}".format(root_id, norm_id))
        tree_renamed.nodes[norm_id]['length'] = 0

    return tree_renamed, num_leaf


# Convert the newick tree to edge list as it is
def convert_newick_file(tree_file, node_prefix, incl_normal, scaling, root_blen = 0):
    newick = get_newick_from_file(tree_file)
    # print("original tree is {}".format(newick))
    tree_renamed, nleaf = newick_to_edgelist(newick, node_prefix, incl_normal, root_blen)
    write_edge_list(tree_file, tree_renamed, nleaf, scaling)



if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Converting NEWICK format to EGDE list')
    parser.add_argument('-f','--format', dest='format', metavar='N', type=int, default=1,
                        help='The format of input file (0: NEXUS (which may contain multiple NEWICK strings), 1: NEWICK)')
    parser.add_argument('-t','--tree_file', dest='tree_file', type=str, required=True,
                        help='The input tree file')
    parser.add_argument('-p','--node_prefix', dest='node_prefix', type=str, default="",
                        help='The prefix of the node to be replaced by space so that node name is an integer')
    parser.add_argument('-b','--root_blen', dest='root_blen', type=float, default=0,
                        help='The additional branch length connecting LUCA to MRCA if normal sample is not in the original tree')
    parser.add_argument('-s','--scaling', dest='scaling', type=float, default=1,
                        help='Scaling the branch length')
    parser.add_argument('-n','--incl_normal', action='store_true', help='The normal sample is included in the tree')

    args = parser.parse_args()

    if args.format == 0:    # Convert all trees in a NEXUS file
        print("convert nexus tree files")
        convert_nexus_file(args.tree_file, args.node_prefix, args.incl_normal, args.scaling, args.root_blen)
    else:
        print("convert newick tree file")
        convert_newick_file(args.tree_file, args.node_prefix, args.incl_normal, args.scaling, args.root_blen)

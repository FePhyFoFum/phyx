#!/usr/bin/env python3

# generates all (2 * (n + 1) - 5)!! rooted trees for n taxa
# hacked version of https://stackoverflow.com/a/14901512/3389183

import sys
import string

# A very simple representation for Nodes.
# Leaves are anything which is not a Node.
class Node(object):
  def __init__(self, left, right):
    self.left = left
    self.right = right

  def __repr__(self):
    return '(%s %s)' % (self.left, self.right)

# Given a tree and a label, yields every possible augmentation of the tree by
# adding a new node with the label as a child "above" some existing Node or Leaf.
def add_leaf(tree, label):
  yield Node(label, tree)
  if isinstance(tree, Node):
    for left in add_leaf(tree.left, label):
      yield Node(left, tree.right)
    for right in add_leaf(tree.right, label):
      yield Node(tree.left, right)

# Given a list of labels, yield each rooted, unordered full binary tree with
# the specified labels.
def enum_unordered(labels):
  if len(labels) == 1:
    yield labels[0]
  else:
    for tree in enum_unordered(labels[1:]):
      for new_tree in add_leaf(tree, labels[0]):
        yield new_tree

def dfactorial(n):
     if n == 0 or n == 1:
         return 1
     else:
         return n * dfactorial(n-2)

def get_num_trees(n):
  # (2 * (n + 1) - 5)!!
  N = 2 * (n + 1) - 5
  return dfactorial(N)

def usage():
  print("Usage: python3 python_trees.py num_tips")

def main(labels):
  i = 1
  for tree in enum_unordered(labels):
    phy = str(tree).replace(" ", ",")
    print(str(i) + ". " + phy + ";")
    i += 1

if __name__ == "__main__":
  n = None
  labels = None
  #labels = tuple(string.ascii_uppercase[:n])
  
  if len(sys.argv) == 2:
    nt = sys.argv[1]
    if not nt.isdigit():
      usage()
      exit()
    else:
      n = int(nt)
      if n < 3:
        print("Minimum of 3 taxa.")
        exit()
  else:
    usage()
    exit()
  #labels = tuple(string.ascii_uppercase[:n]) # simple letters (limited to 26)
  labels = tuple("taxon_" + str(i) for i in range(1, (n + 1))) # match pxbdsim labels
  print("Generating " + str(get_num_trees(n)) + " trees for " + str(n) + " taxa.")
  main(labels)

#!/usr/bin/env python
# encoding: utf-8

from treelib import Node, Tree
tree = Tree()
tree.create_node("n1", 1)
tree.create_node("n2", 2, parent=1)
tree.create_node("n3", 3, parent=1)
tree.create_node("n4", 4, parent=2)
tree.create_node("n5", 5, parent=4)
tree.create_node("n6", 6, parent=3)
tree.create_node("n7", 7, parent=3)
tree.create_node("n9", 9, parent=3)
tree.create_node("n8", 8, parent=6)
tree.create_node("n10", 10, parent=8)
tree.show()


dd = []
for x in tree.all_nodes():
    if x.is_leaf():
        d = []
        print x.identifier
        d.append(x)
        nid = x.bpointer
        print 'dd', nid
        xx = tree.get_node(nid)
        while not xx.is_root():
            d.append(xx)
            nid = xx.bpointer
            xx = tree.get_node(nid)
        dd.append(d)

print len(dd)
for x in dd:
    print len(x)
    print 'P',
    for xx in x:
        print xx.identifier,
    print

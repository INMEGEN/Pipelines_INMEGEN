#! /usr/bin/env python2


import sys
import csv
from lxml import etree


__author__ = 'Adam Caldwell'


root = etree.Element('krona')

# print attributes list, currently only count
attributes = etree.SubElement(root, "attributes", magnitude="count")
attribute = etree.SubElement(attributes, "attribute", display="Count")
attribute.text = "count"

datasets = etree.SubElement(root, "datasets")

# read header, print mothur groups as krona datasets
with open(sys.argv[1], 'r') as taxonomyFile:
    line = taxonomyFile.readline()

# mothur prints a spurious extra tab at the end of the line, so only read to -1, as of mothur 1.37 the spurious tab is gone
    for group in line.split('\t')[5:]:
        dataset = etree.SubElement(datasets, "dataset")
        dataset.text = group

# read through taxonomy summary and print node counts
taxonomyReader = csv.reader(open(sys.argv[1]), delimiter='\t')


def writeCount(row):
    for group in row[5:]:
            val = etree.SubElement(count, "val")
            val.text = group

for row in taxonomyReader:
    if row[0] == "0":
        L0node = etree.SubElement(root, "node", name=row[2])
        count = etree.SubElement(L0node, "count")
        writeCount(row)
    if row[0] == "1":
        L1node = etree.SubElement(L0node, "node", name=row[2])
        count = etree.SubElement(L1node, "count")
        writeCount(row)
    if row[0] == "2":
        L2node = etree.SubElement(L1node, "node", name=row[2])
        count = etree.SubElement(L2node, "count")
        writeCount(row)
    if row[0] == "3":
        L3node = etree.SubElement(L2node, "node", name=row[2])
        count = etree.SubElement(L3node, "count")
        writeCount(row)
    if row[0] == "4":
        L4node = etree.SubElement(L3node, "node", name=row[2])
        count = etree.SubElement(L4node, "count")
        writeCount(row)
    if row[0] == "5":
        L5node = etree.SubElement(L4node, "node", name=row[2])
        count = etree.SubElement(L5node, "count")
        writeCount(row)
    if row[0] == "6":
        L6node = etree.SubElement(L5node, "node", name=row[2])
        count = etree.SubElement(L6node, "count")
        writeCount(row)
    if row[0] == "7":
        L7node = etree.SubElement(L6node, "node", name=row[2])
        count = etree.SubElement(L7node, "count")
        writeCount(row)
    if row[0] == "8":
        L8node = etree.SubElement(L7node, "node", name=row[2])
        count = etree.SubElement(L8node, "count")
        writeCount(row)
    if row[0] == "9":
        L9node = etree.SubElement(L8node, "node", name=row[2])
        count = etree.SubElement(L9node, "count")
        writeCount(row)

print etree.tostring(root, pretty_print=True, xml_declaration=True)

#!/usr/bin/python3.3

import argparse
#from xml.etree.ElementTree import ElementTree, fromstring, tostring, parse, Element, register_namespace
from xml.etree.ElementTree import ElementTree, fromstring, tostring, parse, Element
# apparently register_namespace() is not actually contained within the ElementTree class
#from xml.etree.ElementTree import ElementTree, Element
#import xml.etree.ElementTree 
import sys
import csv
from collections import defaultdict
import urllib.request
import shutil
import gzip
import pickle
import os
import time

## for me
import inspect
def lineno():
    """Returns the current line number in our program."""
    return "line " + str(inspect.currentframe().f_back.f_lineno) + ": "
###

def parseArgs():

    desc = ("CHAT identifies nodes in a network that are significantly more "
            "connected to contextually relevant nodes (e.g. differentially "
            "expressed genes) than is expected by chance.\n\nTo run CHAT:\n\n"
            "Provide an XGMML formatted network file (see chat.primesdb.eu "
            "for format requirements) and specify the attribute containing "
            "context:"
            "\nchat.py -x filename -a attributeName\n\nOR\n\n"
            "A tab delimited text file of human/mouse Ensembl Gene IDs (first column, " 
            "labeled with \"ens\") plus contextual information (second column, "
            "labeled \"pVal\"). CHAT will build a network and then run the analysis:"
            "\nchat.py -n filename\n")
    usg = "%(prog)s [-n listOfGenes | -x xgmmlFile] [optional flags]" 
    cite = "Please cite: [paper details here]"

    parser = argparse.ArgumentParser(usage = usg,
                                     add_help=False, 
                                     description=desc, 
                                     epilog=cite,
                                     formatter_class=argparse.RawTextHelpFormatter)

    inputDataGroup = parser.add_mutually_exclusive_group()
    inputDataGroup.add_argument('-x', help='Input xgmml file', metavar='xgmml')
    inputDataGroup.add_argument('-n', help='List of genes and context data for network', metavar='genes')

    parser.add_argument('-h', action='help', help='Show this help message and exit')
    parser.add_argument('-a', help='Attribute that contains context (case sensitive)', metavar='attName')
    parser.add_argument('-e', help='Attribute that contains ENS ID (case sensitive)', metavar='ensAtt')
    parser.add_argument('-aFC', help=argparse.SUPPRESS, metavar='attFC') # specify attribute that has fold-change data
    parser.add_argument('-aPV', help=argparse.SUPPRESS, metavar='attPV') # specify attribute that has pValue data

    parser.add_argument('-pV', help='P-value threshold, default 0.1', metavar='pThresh')
    parser.add_argument('-fc', help='Fold-change boundaries, e.g. 2', metavar='foldChangeBounds')
    #sekrit range option
    parser.add_argument('-r', help=argparse.SUPPRESS)

    parser.add_argument('-g', help='Number of connections to classify a hub (default n=5)', metavar='hubDeGree', default = 5, type=int)
    parser.add_argument('-m', help='MITAB file of interactions', metavar='mitab')
    parser.add_argument('-p', help='Specify a prefix to be appended to output files.', metavar='prefix', default='')
    parser.add_argument('-v', help='Print all hub parameters', action='store_true')

    # secret arguments
    #parser.add_argument('-U', help='URL for mitab file', metavar='URL', default="http://innatedb.com/download/interactions/all.mitab.gz")
    parser.add_argument('-U', help=argparse.SUPPRESS)
    parser.add_argument('-t', help=argparse.SUPPRESS) # print timing information to err.txt
    #parser.add_argument('-d', help='Delimiter for context information file (default: tab)', metavar='delimiter', default = '\t')
    parser.add_argument('-d', help=argparse.SUPPRESS, metavar='delimiter', default = '\t')
    parser.add_argument('-ve', help=argparse.SUPPRESS, action='store_true') # print verbose errors to err.txt
    parser.add_argument('-spoke', help=argparse.SUPPRESS, action='store_true') # don't decompose complexes
    parser.add_argument('-interactorList', help=argparse.SUPPRESS) # just list the interactors for this node (idTag)
    parser.add_argument('-networkOnly', help=argparse.SUPPRESS, action='store_true') # get the network and STOP
    parser.add_argument('-nx', help=argparse.SUPPRESS, action='store_true') # use the network generated from PSICQUIC #HAX!
    # Include the upper limit for context range, ie: x <= upper
    parser.add_argument('-u', help=argparse.SUPPRESS, default=False, action='store_true') 
    # include the lower limit for specified context range, ie lower <= x
    parser.add_argument('-l', help=argparse.SUPPRESS, default=False, action='store_true')
    parser.add_argument('-printURL', help=argparse.SUPPRESS, action='store_true') # print the network URL generated to stdout
    parser.add_argument('-debug', help=argparse.SUPPRESS, action='store_true') # print random shite
    parser.add_argument('-testRange', help=argparse.SUPPRESS, type=int) # is this value contextual?

    # in case the user has a network and a CSV of data. Consider taking this out...
    contextDataGroup = parser.add_mutually_exclusive_group()
    contextDataGroup.add_argument('-c', help=argparse.SUPPRESS) # file of contextual data, tab delimited by default (use -d to specify delimiter)
    contextDataGroup.add_argument('-i', help=argparse.SUPPRESS, action='store_true') #File contains contextual information written by Cerebral or InnateDB

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    # input checking
    try: 
        args = parser.parse_args()
    except NameError:
        print("One of your options was specified incorrectly, check your input and try again.")
        exit()
    #if args.nx:
    #    args.x = args.p + "-network.xgmml" # hax!
    #    args.a = "context1"
    if args.x:  # for pre-generated xgmml networks
        if args.i:
            args.a = 'cond_1_pvals'
        if args.a: # if user supplies an attribute, obviously the context is embedded
            args.i = True
            if args.fc and args.pV: # with -a can only specify one test type
                choosePlx("-fc", "-pV", parser)
        if not (args.c or args.i): # where is the context? 
            #choosePlx("-c", "-i", parser)
            if not (args.a or (args.aFC and args.aPV)):
                print("Please specify -a")
                exit()
    elif args.n:
        if args.fc and args.pV:
            args.aFC = 'context2'
            args.aPV = 'context1'
        elif args.fc: 
            args.a = 'context2'
        else: 
            args.a = 'context1'
    else: # where is the network?
        choosePlx("-n", "-x", parser)
    if args.p:
        args.p = args.p + "-"
    if not args.pV:
        if not args.fc:
            args.pV = 0.1

    return vars(args)

def choosePlx(a, b, parser):
    parser.print_help()
    print("\nError: Please specify either " + a + " or " + b)
    sys.exit(1)

class Node():
    """ Contains node data """
    def __init__(self, name):
        self.name = name
        self.idTag = ""
        self.context1 = 0  #pValue by default
        self.context2 = 0  #foldchange by default, any secondary attribute
        self.interactors = [] # list of unique interacting nodes
        self.xref = {}  # dictionary of cross reference IDs
        self.selfInteractions = 0
        self.isComplex = False
        self.sameAs = "" # this should be an idTag
        self.isContextual = False

    def toString(self):
        """ What has it got in its pocketses?"""
        return self.name + " " +  str(self.idTag) + " " +  str(self.context1) + " " +  str(self.context2) + " " +  str(len(self.interactors)) + " " +  str(self.isContextual) + " " +  str(self.xref)

    def sigInteractors(self):
        """ How many contextual interactors do I have?"""
        sigs = []
        for node in self.interactors:
            if node.isContextual:
                sigs.append(node)
        pairs = 0
        for sigNode in sigs:
            if any(x for x in sigs if x.idTag == sigNode.sameAs):
                pairs += 1
        if pairs % 2 != 0:
            print(lineno(), pairs)
            print("Something went wrong with counting significant interactors for " + self.name + ". Check your output for discrepancies.")
        return len(sigs) - int(pairs/2)

    def getXref(self, xrefString):
        """ returns xrefs: this may be a string or a tuple of strings """
        if xrefString in self.xref:
            return self.xref[xrefString]
        else:
            return "Xref not found"

############### DIY network stuff

### csv loading for DIY network
def loadNetworkCSV(filename, delim):
    """ takes a filename and returns a list of node data
    [ (ens, context1, context2), (ens, context1, context2)...] """
    dataList = []
    with open(filename, 'rt') as f:
        reader = csv.DictReader(f, delimiter=delim)
        for row in reader:
            dataList.append(row)
    return dataList

### node maker for DIY network
def makeLittleNode(data):
    """ data contains: (ens, context1 (should be a pValue), context2 (should be a fold-change)) """
    # this may need to take a unique identifier. For now try using n.name = ens until it breaks.
    if getName(data[0]):
        node = Node(getName(data[0]))
        #print(lineno(), getName(data[0]))
    else:
        node = Node(data[0])
        #print(lineno(), data[0])
    node.xref['ens'] = data[0]
    node.context1 = data[1]
    node.context2 = data[2]
    return node
def makeAllNodes(dataList):
    """ takes a dictionary of node data {ens:ens, att1:val1, att2:val2} and returns a 
        dictionary of label -> node 
        In this case, label == ens """
    nodes = {}
    for item in dataList:
        if 'ens' in item:
            if 'pVal' and 'FC' in item: # this should be the most common case
                nodes[item['ens']] = makeLittleNode((item['ens'], item['pVal'], item['FC']))
            elif 'pVal' in item:
                nodes[item['ens']] = makeLittleNode((item['ens'], item['pVal'], 0))
            elif 'FC' in item:
                nodes[item['ens']] = makeLittleNode((item['ens'], 0, item['FC']))
            else: # if there's only ens and no pVal or FC, all genes are significant
                print("Found no pValue or FC data in CSV; all genes are being assigned pValue = 0.1")
                nodes[item['ens']] = makeLittleNode((item['ens'], 0.1, 0))
        else:
            print("Input CSV file must have Ensembl IDs in a tab delimited (or other specified using -d) column labeled 'ens'.")
            exit()
        #nodes[item[0]] = makeNode(item)
    return nodes

### network builder
def buildNetwork(nodesIn, interactionDict, err=False, timing=False):
    """ take a dictionary of little nodes (ens: node(name == ENS, 
        context1, context2)) and a processed mitab file and make a 
        CHAT-style network dictionary out of it."""

    print("Building network...")
    lookupStart = time.clock()
    nodeDict = {} # contains uniqueID -> node
    newNodes = {}
    i = 1
    for node in nodesIn:
        # get interactors
        if nodesIn[node].xref['ens'] in interactionDict:
            nodesIn[node].interactors = interactionDict[nodesIn[node].xref['ens']]
            nodesIn[node].idTag = i
            nodeDict[i] = nodesIn[node]
            i += 1
        else: 
            if err:
                with open ('err.txt', 'a') as e:
                    if 'ens' in nodesIn[node].xref:
                        e.write(nodesIn[node].xref['ens'] + " not found in mitab.\n")
                    else: 
                        e.write(nodesIn[node].name + " not found in mitab.\n")
            continue

        # if interactors not in network already, make it and set it aside. 
        for ens in nodesIn[node].interactors:
            if ens not in nodesIn and ens not in newNodes: 
                newNode = makeLittleNode((ens, 0, 0))
                newNode.interactors = interactionDict[ens]
                newNodes[ens] = newNode
    allEns = []
    for item in newNodes:
        allEns.append(item)
    for thingy in nodesIn:
        allEns.append(thingy)

    if timing: 
        with open('err.txt', 'a') as e:
            e.write("Network lookup took " + str((time.clock() - lookupStart)) + "s\n")


    start = time.clock()
    for node in newNodes:
        for tag in newNodes[node].interactors:
            if tag in allEns: 
                pass
            else:
                newNodes[node].interactors = [x for x in newNodes[node].interactors if x != tag]
        newNodes[node].idTag = i
        nodeDict[i] = newNodes[node]
        i += 1

    # right now node.interactors is an array of ENS IDs, not an array of nodes.
    # We need to make this be an array of nodes so node.isContextual works correctly. 
    # nodesIn and newNodes should, together, contain every node in our network mapped
    # to their ENS id (ens -> node). 
    # so I combine those dictionaries, pull the nodes by ENS ID, replace those 
    # arrays, and bob's your uncle

    allNodes = dict(list(nodesIn.items()) + list(newNodes.items()))
    for j in nodeDict:
        nodeList = []
        for tag in nodeDict[j].interactors:
            nodeList.append(allNodes[tag])
        nodeDict[j].interactors = nodeList
        #print(lineno(), str(len(nodeDict[j].interactors)))

    if timing: 
        with open('err.txt', 'a') as e:
            e.write("Network building took " + str((time.clock() - start)) + "s\n")

    return nodeDict

def writeNetwork(nodeDict, name):
    if not name:
        name = 'CHAT'
    print("Writing network to " + name + ".xgmml")
    xgmml = "http://www.cs.rpi.edu/XGMML"

    ensLookup = {} #make this so we can search it for edge creation
    for label in nodeDict:
        ensLookup[nodeDict[label].name] = label

    nodes = []
    edges = {} # { (source, target) : elSourceTarget} in order to not make duplicates

    #for idTag in nodeDict:
    #    print(idTag, nodeDict[idTag].name)

    for idTag in nodeDict:
        #print(lineno(), nodeDict[idTag].toString())
        # build node element
        el = createNodeEl(str(idTag), nodeDict[idTag].name)
        el.append(createAttEl('context1', str(nodeDict[idTag].context1), 'real'))
        el.append(createAttEl('context2', str(nodeDict[idTag].context2), 'real'))
        if 'ens' in nodeDict[idTag].xref:
            ensTag = nodeDict[idTag].xref['ens']
        else: ensTag = "None"
        el.append(createAttEl('ens', ensTag, 'string'))
        nodes.append(el)

        dupes = 0
        # build edge elements 
        for n in nodeDict[idTag].interactors:
            source = idTag
            #print(lineno(), n.name)
            try: target = ensLookup[n.name]
            except KeyError as e:
                try: os.remove('net.pkl')
                except FileNotFoundError: pass
                print(lineno(), ensLookup)
                exit()
            if (target, source) in edges: dupes += 1
            else: 
                #print(lineno(), nodeDict[idTag].name)
                edges[(source, target)] = createEdgeEl((nodeDict[idTag].name + " interacts with " + n.name), str(source), str(target))

    # build a tree structure
    root = Element('graph')
    root.attrib['directed'] = '1'
    root.attrib['label'] = name
    root.attrib['xmlns'] = xgmml
    
    # populate the network
    for node in nodes:
        root.append(node)
    for pair in edges:
        root.append(edges[pair])

    # lay out the network

    # wrap in an ElementTree instance and write out
    tree = ElementTree(root)
    indent(tree.getroot())
    outfile = name + ".xgmml"
    tree.write(outfile, encoding="utf-8", xml_declaration=True)

def createAttEl(name, value, attType):
    return Element('att', {'name': name, 'value': value, 'type':attType})
def createNodeEl(idTag, label):
    return Element('node', {'id': idTag, 'label': label})
def createEdgeEl(label, source, target):
    return Element('edge', {'label': label, 'source': source, 'target': target})
########################### end DIY network stuff

def makeNodeDict(xml, args):
    # fuck me this is too complicated. 
    if args['aFC']:
        if not args['aPV']:
            print("Please supply -aPV as well as -aFC, or just use -a to test one attribute only.")
        else: 
            return extractNodes(xml, args, args['aPV'], False, args['aFC']) 

    if args['i'] and not args['a']:
        return extractNodes(xml, args, findContextValue(root), True) # deduce context information from xgmml
    elif args['a']: # use user-provided context attribute
        if args['fc']: # store context value in node.context2 and node.context1
            return extractNodes(xml, args, args['a'], False, args['a']) 
        else: 
            return extractNodes(xml, args, args['a']) # store value in node.context1 only
    elif args['n'] or args['nx']:
        return extractPsicquicNodes(xml)
    else:
        noContextDict = extractNodes(xml, args) # this expects an InnateDB-generated file, contains node.idTag -> Node
        context = loadCSV(args['c'], args['d'])
        return addContext(context, noContextDict) 

def extractNodes(xmlNodes, args, contextValue = "", cerebral=False, secondaryAttName=""):
    nodeDict = {} # dict of nodes: node.idTag -> Node
    if args['ve']: err = True
    else: err = False
    if args['e']: ensName = args['e']
    else: ensName = ""

    for el in xmlNodes:
        n = Node(el.attrib["label"])
        n.idTag = el.attrib["id"]
        xrefs = 0

        for att in el.findall("{http://www.cs.rpi.edu/XGMML}att"):
            # xref extraction
            if ensName: 
                if att.get("name") == ensName:
                    n.xref["ens"] = att.get("value")
                    xrefs += 1
            else: 
                try:
                    if att.get("name") == "Cross-references":
                        n.xref["ens"] = att.get("value").split("|")[2]
                        n.xref["idb"] = att.get("value").split("|")[1]
                        xrefs += 1
                except IndexError:
                    if err:
                        print("err exists!")
                        with open('err.txt', 'a') as e:
                            e.write("Cannot parse Cross-reference data for " + el.attrib['label'] + ": " + att.get('value'))

            # context value extraction
            if cerebral:
                contextValName = contextValue + "_pvals" #old style: *_sig
                if att.get("name") == contextValue:
                    n.context2 = Decimal(att.get("value"))
            else: 
                contextValName = contextValue

            if att.get("name") == contextValName:
                n.context1 = Decimal(att.get("value"))
            if att.get("name") == secondaryAttName:
                n.context2 = Decimal(att.get("value"))

        if "complex" in n.name: # this is for InnateDB files, may not work as desired for other databases
            n.isComplex = True
        nodeDict[n.idTag] = n

    if xrefs < 1: # There needs to be some kind of error printed if nothing is cross-referenced.
        if args['ve']:
            with open('err.txt', 'a') as f:
                f.write("No ENS IDs collected from nodes; specify -e for cross-referencing.\n")
        print("No ENS IDs collected from nodes; specify -e for cross-referencing.")

    return nodeDict

def findContextValue(xml): 
    atts = xml.findall("{http://www.cs.rpi.edu/XGMML}att")
    for att in atts:
        if att.get("name") == "cerebralPCAttributes":
            return att.get("value").split('=')[1]
    return "none"

def addEdges(xmlEdges, nodeDict, spoke=False):
    complexesToDecompose = defaultdict(list)
    
    for edge in xmlEdges:
        sourceId = edge.attrib["source"]
        targetId = edge.attrib["target"]

        if sourceId == targetId:          # self interaction
            #print(nodeDict[sourceId].name, sourceId)
            nodeDict[sourceId].selfInteractions += 1
        else: 
            if spoke == False: # decompose complexes
                if nodeDict[sourceId].isComplex:    # complex (source)
                    complexesToDecompose[sourceId].append(targetId)
                elif nodeDict[targetId].isComplex:   # complex (target)
                    complexesToDecompose[targetId].append(sourceId)
                else: 
                    addNonRedundantEdge((sourceId, targetId), nodeDict)
                
            else:                               # everything else
                addNonRedundantEdge((sourceId, targetId), nodeDict)

    # turned off by default
    # nb: decomposing complexes increases the # of edges by A LOT
    if spoke == False: 
        for key in complexesToDecompose:
            for i, node in enumerate(complexesToDecompose[key]): 
                for j, node in enumerate(complexesToDecompose[key]):
                    if (complexesToDecompose[key][i] == complexesToDecompose[key][j]):
                        pass
                    else: 
                        addNonRedundantEdge((complexesToDecompose[key][i], complexesToDecompose[key][j]), nodeDict)
            del nodeDict[key]

    return nodeDict

def addNonRedundantEdge(edge, nodeDict): #edge is a (source, target) tuple
    """ For a pair of nodes, add each node object to the other's node.interactors[] list """
    sourceId = edge[0]
    targetId = edge[1]
    # add target to source's interactors list if it's not already there
    if nodeDict[targetId] not in nodeDict[sourceId].interactors:
        nodeDict[sourceId].interactors.append(nodeDict[targetId])
    # add source to target's interactors list if it's not already there
    if nodeDict[sourceId] not in nodeDict[targetId].interactors:
        nodeDict[targetId].interactors.append(nodeDict[sourceId])

def combineGenes(nodeDict, err=False):
    # find gene pairs, store in a list of tuples
    geneDict = {}
    nameDict = {}
    toCombine = [] # list of tuples (node1, node2)
    for key, n in nodeDict.items(): 
        nameDict[n.name] = n
        if 'gene' in n.name:
            #print(n.name)
            geneDict[n.name] = n
    for key, n in geneDict.items():
        smallerName = key[:-5]
        #print(smallerName)
        #n.name = smallerName
        if smallerName in nameDict.keys():
            toCombine.append((n, nameDict[smallerName]))

    for tup in toCombine:

        # update interactors (UNIQUE interactors), 
        if err: 
            with open('debug.txt', 'a') as e:
                e.write(str(lineno()) + "Combining " + tup[0].name + " (" + tup[0].idTag + ") and " + tup[1].name +" (" + tup[1].idTag + ")\n" + tup[0].name + " (" + tup[0].idTag + "): ")
                for interactor in tup[0].interactors:
                    e.write(interactor.name + ", ")
                e.write("\n")
                e.write(tup[1].name + " (" + tup[1].idTag + "): ")
                for interactor in tup[1].interactors:
                    e.write(interactor.name + ", ")
                e.write("\n")
        
        tup[0].interactors = list(set(tup[0].interactors) | set(tup[1].interactors))
        tup[1].interactors = tup[0].interactors
        if err:
            with open('debug.txt', 'a') as e:
                e.write(str(lineno()) + "After merging:\n" + tup[0].name + " (" + tup[0].idTag + "): ")
                for interactor in tup[0].interactors:
                    e.write(interactor.name + ", ")
                e.write("\n" + tup[1].name + " (" + tup[1].idTag + "): ")
                for interactor in tup[1].interactors:
                    e.write(interactor.name + ", ")
                e.write("\n")
        
        # context pValue  (take the most significant one, probably only one value between the two)
        if tup[0].context1 == 0 or tup[1].context1 == 0:
            context = max(tup[0].context1, tup[1].context1)
        else: 
            context = min(tup[0].context1, tup[1].context1)
        """
        if err: 
            with open('debug.txt', 'a') as e:
                e.write(str(lineno()) + str(tup[0].context1) + "\t" + str(tup[1].context1) + " -> ")
        """
        tup[0].context1 = tup[1].context1 = context
        """
        if err: 
            with open('debug.txt', 'a') as e:
                e.write(str(lineno()) + str(tup[0].context1) + "\t" + str(tup[1].context1) + "\n")
        """

        # expression value 
        """
        if err: 
            with open('debug.txt', 'a') as e:
                e.write (str(lineno()) + str(tup[0].context2) + "\t" + str(tup[1].context2) + " -> ")
        """
        tup[0].context2 = tup[1].context2 = max(tup[0].context2, tup[1].context2)
        """
        if err: 
            with open('debug.txt', 'a') as e:
                e.write(str(lineno()) + str(tup[0].context2) + "\t" + str(tup[1].context2) + "\n")
        """

        # and n.sameAs for each
        tup[0].sameAs = tup[1].idTag
        tup[1].sameAs = tup[0].idTag
        if err:
            with open('debug.txt', 'a') as e:
                e.write(str(lineno()) + tup[0].name + " (" + tup[0].idTag + ") is the same as " + nodeDict[tup[0].sameAs].name + " (" + nodeDict[tup[0].sameAs].idTag + ")\n")
                e.write(tup[1].name + " (" + tup[1].idTag + ") is the same as " + nodeDict[tup[1].sameAs].name + " (" + nodeDict[tup[1].sameAs].idTag + ")\n")

    # now go through and take out all the "NAME gene" nodes
    # #TODO FIRST
    newDict = {}
    for tag in nodeDict:
        if 'gene' in nodeDict[tag].name and nodeDict[tag].sameAs:
            if nodeDict[tag].sameAs:
            #print("remove " + nodeDict[tag].xref['ens'] + " sameAs " + nodeDict[tag].sameAs)
                pass #for now
        else: 
            newDict[tag] = nodeDict[tag]
    return newDict

## these are for getting and storing name data from the mitab file, for 
## the DIY network results (want to be able to print the name as well as
## the ENS IDs.)
def parseName(xrefs):
    """Take the 'alias' field from the mitab file and return a name for the ENS id""" 
    if '(display_short)' in xrefs:
        for term in xrefs.split('|'):
            if 'display_short' in term:
                return term.split(':')[1].split('(')[0]

def getName(ens):
    pklFile = 'names.pkl'
    if (checkForFile(pklFile)): # check for pickled mitab
        names = unpickleFile(pklFile)
        if ens in names:
            return names[ens]
        else:
            return 'None'
    else:
        return 'None'
###


# new all.mitab has an issue on lines 362783 and 362784, the word "destruction" is 
# surrounded by offensive characters. Should change this method to deal with 
# offensive characters. TODO
# but for right now, just take those lines out and use manually modified mitab.
def processMitab(filename):
    """ Takes a filename and returns a dict of interactions """
    """ It will ALSO create a .pkl file with ens -> name in the working directory. """
    if '.pkl' in filename: #hax!
        return unpickleFile(filename)
    pklFile = filename + '.pkl'
    if (checkForFile(pklFile)): # check for pickled mitab
        print("Using previously processed background mitab file. Delete " + pklFile + " to force reprocessing.")
        return unpickleFile(pklFile)

    if not checkForFile(filename): # can we even find it?
        print("Cannot open " + filename)
        exit()
    print("Processing mitab file...")

    # contains ens -> [ens1, ens2, ens3...] (these xrefs are stored in Node.xref["ens"]
    allInteractions = defaultdict(list)
    allNames = {}
    with open(filename, 'rt') as f:
        reader = csv.DictReader(f, delimiter="\t")
        try: 
            for row in reader:
                try: 
                    ensA = row['alt_identifier_A'].split(":")[1] # new 30/7/13
                    ensB = row['alt_identifier_B'].split(":")[1]
                    aliasA = row['alias_A']
                    aliasB = row['alias_B']
                    if ensA != ensB: # don't keep self interactions
                        if ensB not in allInteractions[ensA]: #don't keep duplicate interactions
                            allInteractions[ensA].append(ensB)
                        if ensA not in allInteractions[ensB]: #same here
                            allInteractions[ensB].append(ensA)
                    ## pull name data
                    allNames[ensA] = parseName(aliasA)
                    allNames[ensB] = parseName(aliasB)
                except IndexError:
                    print("Something is up with your mitab file in line " + str(reader.line_num))
                    exit()
        except UnicodeDecodeError as e:
            print(e + "\nUnable to process mitab file. Is it a normal text file? (Not zipped?)")
            exit()

    if len(allInteractions) < 1:
        print("Unable to process mitab file. Does it contain data?")
        exit()

    # pickle this dictionary so we don't have to parse that mitab file every time
    pickleIt(filename, allInteractions)
    # save the names for later
    pickleIt('names', allNames)

    return allInteractions

def getMitab(url):
    print("Retrieving mitab file...")
    outFile = url.split('/')[-1]
    try: 
        with urllib.request.urlopen(url) as response, open(outFile, 'wb') as out:
            shutil.copyfileobj(response, out)
        return outFile
    except URLError as e:
        print(e.args)
        exit()

def unzip(inFileName, outFileName):
    zippedFile = gzip.open(inFileName, 'rb')
    outFile = open(outFileName, 'wb')
    try: 
        outFile.write(zippedFile.read())
    finally:
        outFile.close()
        zippedFile.close()

def getIDBMitab():
    url='http://innatedb.com/download/interactions/all.mitab.gz'
    unzipped = "all.mitab"
    zipped = "all.mitab.gz"

    # check for pickled mitab
    if (checkForFile(unzipped + '.pkl')):
        print("Using previously processed background mitab file. Delete " + unzipped + ".pkl to force reprocessing.")
        #print(lineno())
        return unzipped + ".pkl"

    # is the file already here?
    if checkForFile(unzipped):
        return unzipped
    # is the .gz here?
    elif not checkForFile(zipped):
        getMitab(url)
    unzip(zipped, unzipped)
    return unzipped

def loadCSV(csvFile, delimiter='\t'): 
    contextInfo = {} # contains tag -> value
    print("Processing " + csvFile + "...")
    try: 
        with open(csvFile, 'r') as f:
            reader = csv.reader(f, delimiter=delimiter, skipinitialspace=True) 
            for row in reader:
                try:  
                    contextInfo[row[0]] = row[1]
                except IndexError: 
                    contextInfo[row[0]] = "999"

    except FileNotFoundError:
        print(csvFile + " not found.")
        exit()

    # if that file only had one column, assume all genes are significant
    if all(val == "999" for val in contextInfo.values()) :
        #set all values to 0.01
        print
        c = dict.fromkeys(contextInfo.keys(), "0.01")
        return c
    return contextInfo

def addContext(csvDict, nodeDict):
    """ csvDict contains ens -> value
        nodeDict contains idTag -> node """
    for tag in nodeDict:
        if 'ens' in nodeDict[tag].xref:
            if nodeDict[tag].xref['ens'] in csvDict:
                nodeDict[tag].context1 = csvDict[nodeDict[tag].xref['ens']]

    return nodeDict

def findHubs(nodeDict, hubThresh):
    hubList = []

    for tag in nodeDict:
        if len(nodeDict[tag].interactors) >= hubThresh:
            hubList.append(nodeDict[tag])

    if len(hubList) == 0:
        print("No nodes with degree >= " + str(hubThresh) + " found.")
        exit()
    return hubList

def printHubs(hubList, prefix):
    hubList.sort(key=lambda n: len(n.interactors), reverse=True) 
    outFile = prefix + "topHubs.txt"
    out = open(outFile, 'w')
    out.write("Top hubs by degree\n")
    for node in hubList:
        #print(node.name)
        if 'gene' in node.name: continue # don't print "NAME gene", it is identical to "NAME"
        if 'ens' in node.xref:
            out.write(node.xref['ens'])
        else: # just ignore complexes
            continue
            out.write(node.name)
        out.write("\t" + str(len(node.interactors)) + "\n")
    out.write("\n")

    out.write("Node\tName\tTotal Interactors\tContext Relevant Interactors\n")
    for node in hubList:
        if 'ens' in node.xref: 
            out.write(node.xref['ens'] + "\t" + node.name + "\t")
        else: 
            out.write(node.name)
        out.write(str(len(node.interactors)) + "\t" + str(node.sigInteractors()) + "\n")
    out.close()

def getDEnodes(nodeDict): #dict of node objects
    """ Given a dictionary of node objects, find how many context-relevant nodes are in the dict """
    n = 0
    for key in nodeDict:
        #print(lineno(), key, nodeDict[key].isContextual)
        if nodeDict[key].isContextual:
            n += 1
    if n == 0:
        print("No contextually relevant nodes found.")
        exit()
    return n

#returns a list of tuples: (idTag, M, n, N, x, prb, overrep) overrep is a boolean, false if underrepresented
def getProbs(hubList, M, n, allNodeDict):
    print("Calculating overrepresentation...")
    allProbs = []
    ratio = n/M
    for hub in hubList:
        x = hub.sigInteractors() # x = DE interactors for the node
        #print(x)
        if 'ens' in hub.xref:
            N = len(allNodeDict[hub.xref['ens']])
            if N == 0:  # no interactors? move on. 
                if args['ve']:
                    with open('err.txt', 'a') as e:
                        e.write(hub.xref['ens'] + " (" + hub.name + ") not found in universe mitab\n")
                continue
                # how can this happen? If an interaction is predicted by orthology by InnateDB or
                # created through decomposition of a complex, it may not have any "actual" 
                # interactors in the universe mitab file. 
                # These cases should not have very many significant interactors, and so should
                # not be that big a deal. 
        else: # no xref means we can't look it up in the mitab file, so ditch it
            if args['debug']:
                with open('debug.txt', 'a') as e:
                    e.write(str(lineno()) + hub.name + " has no ens tag\n")
            if args['ve']:
                with open('err.txt', 'a') as e:
                    e.write(hub.name + " has no ens tag for cross-referencing.\n")
            continue

        prb = AliceHyper(x, M, n, N)
        overrep = True
        if x/N < ratio: # if we're actually seeing FEWER sig nodes than expected by chance
            overrep = False
        if x > N:
            if args['ve']:
                with open('err.txt', 'a') as e:
                    e.write(hub.xref['ens'] + " (" + hub.name + ") has more interactors in the network than in the background mitab file. Check results carefully.\n")

        allProbs.append((hub.idTag, M, n, N, x, prb, overrep))

    return allProbs

from decimal import Decimal

def AliceHyper(x, M, n, N):
    if x > N:
        return 999
    prb = Decimal(Binom(n, x)) * Decimal(Binom((M - n), (N - x))) / Decimal((Binom(M, N)))
    return prb

def Binom(a, b):
    top = 1
    bottom = 1
    for i in range(0, b): #(0, 1, 2... b-1)
        top = Decimal(top * (a - i))
    for j in range(b, 0, -1): # (b, b-1, b-2...1)
        bottom = Decimal(bottom * j)
    return top/bottom

def corrProbs(probsList, alpha=0.05): # probsList contains (idTag, M, n, N, x, prb, overrep)
    """ Correct for multiple testing """
    allProbs = [] #list of tuples: (idTag, M, n, N, x, prb, overrep, BH, BF)
    numTests = len(probsList)

    corrected = {} # contains initial rank -> (BH, Bon)
    probsList.sort(key=lambda tup: tup[5])
    for i, tup in reversed(list(enumerate(probsList, 1))):
        bon = min(tup[5] * numTests, 1) # Bonferroni correction

        #q = ((i+1)*alpha)/numTests 
        adjBH = min(Decimal((numTests/i)) * tup[5], 1)
        if i < len(probsList):
            corrected[i] = (min(adjBH, corrected[(i+1)][0]), bon)
        else:
            corrected[i] = (adjBH, bon)
        #print(i, tup[5], adjBH)
    for key in corrected:
        allProbs.append(probsList[key-1] + corrected[key])
    return allProbs 

import operator

def printAllProbs(hubProbs, nodeDict, prefix, verbose, sort=""): #hubProbs contains (idTag, M, n, N, x, prb, overrep, BH, bon)
    outfile = prefix + "results.txt"
    # M: universe
    # n: DE nodes in universe
    # N: all idb interactors for the node
    # x: DE nodes for the hub
    sort_on = 7 #default: sort on benjamini-hochberg corrected values 
    if sort == "bon":
        sort_on = 8
    elif sort == "prb":
        sort_on = 5
    hubProbs.sort(key=lambda tup: Decimal(tup[sort_on]))

    # print out a "bottom line" with ENS prb
    bottomLine = open(outfile, 'w')
    bottomLine.write("ENS\tname\tp\tp-corr\n")
    for tup in hubProbs:
        #print(tup)
        if tup[6]: # if the value is overrepresented
            #print(lineno(), tup)
            bottomLine.write(nodeDict[tup[0]].getXref('ens') + "\t")
            bottomLine.write(nodeDict[tup[0]].name + "\t")
            bottomLine.write('{:.4e}'.format(tup[5], 'g') + "\t")
            bottomLine.write('{:.4e}'.format(tup[sort_on])) 
            bottomLine.write("\n")
    bottomLine.close()
    
    # print "verbose" output with all parameters to file
    if verbose:
        extras = [] #print the underrepresented stuff afterwards
        vFile = prefix + "allNumbers.txt"
        verbose = open(vFile, 'w')
        # how it's notated in the tuple
        #verbose.write("ENS\tname\tprb\tBH\tbon\tallGenes (M)\tsigGenes (n)\tallInteractors (N)\tsigInteractors (x)\n")
        # how it's notated around the web 
        verbose.write("ENS\tname\tprb\tBH\tbon\tallGenes (N)\tsigGenes (M)\tallInteractors (n)\tsigInteractors (x)\n")
        for tup in hubProbs:
            if tup[6]:
                writeStuff(tup, verbose, nodeDict[tup[0]])
            else:
                extras.append(tup)
        verbose.write("\n\n\n")
        for tup in extras:
            writeStuff(tup, verbose, nodeDict[tup[0]])
    
    #print to stdout for my benefit
    """
    if args['debug']:
        print("id\tM\tn\tN\tx\tprb\tBH\tbon")
        for hub in hubProbs:
            if hub[6]: # this tells us whether the value is overrepresented or not
                for x in range(0, 5):
                    print(str(hub[x]), end="\t")
                for y in range(5,9):
                    if y != 6:
                        print('{:.2e}'.format(hub[y]), end="\t")
                print()
    """

### write out probabilities and hub info
def writeStuff(tup, writer, node):
    writer.write(node.getXref('ens') + "\t")
    writer.write(node.name + "\t")
    for y in range(5, 9):
        if y != 6:
            writer.write('{:.4e}'.format(tup[y], 'g') + "\t")
    for x in range(1, 4):
        writer.write(str(tup[x]) + "\t")
    writer.write(str(tup[4])) # number of significant interactors
    writer.write("\n")

# given an identifier for a node, tell me which things it interacts with.
#
# nb: this only works with the identifiers. You HAVE to look at the xgmml 
# file to figure out what the idTag is (in <node label='someName' id='youWantThisOne'>)
def printInteractors(idTag, nodeDict):
    nodeDict[idTag].interactors.sort(key=lambda n: n.name)
    
    print(nodeDict[idTag].name + " interacts with:")
    #print(nodeDict[idTag].name + " (" + idTag + ") interacts with:")
    sig = 0
    for node in nodeDict[idTag].interactors:
        print(node.name, end="")
        #print(node.name + " (" + idTag + ")", end="")
        #print(node.name, node.xref['ens'])
        #print(node.xref['ens'])
        #if 0 < float(node.context1) <= 0.1:
        if node.isContextual:
            print("*", end="")
            sig += 1
        print(", ", end="")
    
    print("which is " + str(len(nodeDict[idTag].interactors)) + " interactors. (" + str(sig) + " contextual.)")

def openNetwork(filename):
    print("Loading network...")
    try:
        with open(filename, 'rt') as f:
            return parse(f)
    except IOError:
        print("Couldn't open " + filename)
        exit()

#from http://effbot.org/zone/element-lib.htm#prettyprint
def indent(elem, level=0):
    i = "\n" + level*"  "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "  "
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for elem in elem:
            indent(elem, level+1)
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = i

# I don't think I use this anymore...
"""
def printNodeElements(nodeList):
    for node in nodeList:
        print(node.attrib['label'], str(len(node.findall("{http://www.cs.rpi.edu/XGMML}att"))))
"""

## random utilities
def checkForFile(filename):
    try:
        with open(filename, 'rb'):
            return True
    except IOError:
        return False

def unpickleFile(filename):
    with open(filename, 'rb') as f:
        return pickle.load(f)

def pickleIt(filename, thing):
    try:
        out = open(filename + '.pkl', 'wb')
        pickle.dump(thing, out, -1)
    except IOError: # oh well
        pass
###

# return "psicquic" if it has that source attribute
# return  "cytoscape" if it has that source attribute
# return "other" otherwise
def checkNetworkFormat(xml):
    """ Maybe someday. """
    return "other"

### range evaluation
def evaluateAllNodes(nodeDict, rangeInfo):
    """ rangeInfo contains:
        attName : ( (lower, upper, linclude, uinclude), type) """
    # this would be a lot smarter if the context attributes were stored
    # in a dictionary rather than pulled out and named. 
    # Instead: node.context1 is a pValue (pV or r)
    # node.context2 is a fold change value (fc)
    for tag in nodeDict:
        tests = []
        for att in rangeInfo:
            if rangeInfo[att][1] == 'fc':
                tests.append(evaluateContext(nodeDict[tag].context2, rangeInfo[att]))
            else:  # pValues and anything else should be stored in context1. 
                tests.append(evaluateContext(nodeDict[tag].context1, rangeInfo[att]))
        # if they're all true, it's true. If not, False.
        nodeDict[tag].isContextual = all(tests)
        #nodeDict[tag].printAtts()

        """
        nodeDict[tag].isContextual = evaluateContext(nodeDict[tag].context1, rangeInfo[1])
        if rangeInfo2 and nodeDict[tag].isContextual:
            nodeDict[tag].isContextual = evaluateContext(nodeDict[tag].context2, rangeInfo2)
        #print(nodeDict[tag].context1, nodeDict[tag].isContextual)
        """
    return nodeDict

def evaluateContext(value, rangeInfo):
    """ rangeInfo contains ( (lower, upper, lInclude, uInclude), type) """
    #print(rangeInfo)
    comparisonType = rangeInfo[1]
    limits = (rangeInfo[0][0], rangeInfo[0][1])
    lInclude = rangeInfo[0][2]
    uInclude = rangeInfo[0][3]
    value = float(value)

    if comparisonType == 'fc':
        return not limits[0] < value < limits[1]
    elif comparisonType == 'pV':
        return 0 < value <= limits[1]
    elif comparisonType == 'r':
        if lInclude == uInclude == False:
            return limits[0] < value < limits[1]
        if lInclude and uInclude == False:
            return limits[0] <= value < limits[1]
        if lInclude == False and uInclude:
            return limits[0] <= value < limits[1]
        if lInclude and uInclude:
            return limits[0] <= value <= limits[1]
    return False

def getRange(rangeString):
    """ Expects a range in the format 'lower...upper'
        Should be perfectly general. """
    lower = float(rangeString.split('...')[0])
    upper = float(rangeString.split('...')[1])

    if upper < lower:
        print("Range specified should be in the format: lower...upper (i.e. 0...1)")
        exit()
    return (lower, upper)

def processRange(args):
    """ From the arguments given on the command line, return a dictionary with range info:
        attName -> ( (lower, upper, lInclude, uInclude), type) """
    if args['a']:
        if args['fc']:
            return {args['a']: ( ( (-1 * float(args['fc'])), float(args['fc']), True, True), 'fc')}
        elif args['pV']:
            return {args['a']: ( (0, float(args['pV']), False, True), 'pV')}
        elif args['r']:
            userRange = getRange(args['r'])
            return {args['a']: ( (userRange[0], userRange[1], args['l'], args['u']), 'r')}
        else: # default behaviour: look for a pValue
            return {args['a']: ( (0, 0.1, False, True), 'pV')}
    elif args['aFC'] and args['aPV']:
        ret = {}
        ret[args['aFC']] = ( ( (-1 * float(args['fc'])), float(args['fc']), True, True), 'fc')
        ret[args['aPV']] = ( (0, float(args['pV']), False, True), 'pV')
        return ret

    else:
        if args['d']:
            with open('debug.txt', 'a') as f:
                f.write("Range information not specified correctly:")
                f.write(str(lineno()) +  " " + args + "\n")
                f.write(str(lineno()) + " " + str(args['pV']) + " " + str(args['fc']) + " " + str(args['r']) + " " + args['a'])
        exit()

def setUpLogs(args):
    if args['ve']:
        try: os.remove('err.txt')
        except (FileNotFoundError, OSError): pass
        with open('err.txt', 'w') as f:
            f.write('Verbose errors:\n')
    if args['debug']:
        try: os.remove('debug.txt')
        except FileNotFoundError: pass
        with open('debug.txt', 'w') as f:
            f.write('Debugging output:\n')


if __name__ == '__main__':
    args = parseArgs()
    setUpLogs(args)

    # get range test information
    rangeInfoDict = processRange(args)
    
    # handle mitab file
    if (args['m']):
        if checkForFile(args['m']):
            allInteractions = processMitab(args['m'])
        else: 
            print(args['m'] + " not found.")
            exit()
    elif (args['U']): # fetch mitab from given URL
        allInteractions = processMitab(getMitab(args['U']))
    else: # fetch from default URL
        allInteractions = processMitab(getIDBMitab())
    M = len(allInteractions) # nodes in the interactome

    # sort out whether we're building our own network. 
    if args['n']:
        nodeDict = buildNetwork(makeAllNodes(loadNetworkCSV(args['n'], args['d'])), allInteractions, args['ve'], args['t'])
        writeNetwork(nodeDict, args['p'][0:-1]) # strip that dash I added to the prefix
        if args['networkOnly']:
            exit()
    else: # for -x 
        root = openNetwork(args['x'])
        xmlNodes = root.findall("{http://www.cs.rpi.edu/XGMML}node")
        nodeDict = makeNodeDict(xmlNodes, args)
        # add interactions to node list
        xmlEdges = root.findall("{http://www.cs.rpi.edu/XGMML}edge")
        # add edges, decompose complexes, and combine nodes that are "name" and "name gene"
        nodeDict = combineGenes(addEdges(xmlEdges, nodeDict, args['spoke']), args['debug'])

    # tag contextual nodes 
    nodeDict = evaluateAllNodes(nodeDict, rangeInfoDict)

    hubList = findHubs(nodeDict, args['g']) # specifies threshold for what defines a hub
    n = getDEnodes(nodeDict) # DE nodes in the network

    #sekrit option
    if args['interactorList']:
        printInteractors(args['interactorList'], nodeDict)
        exit()

    # do the stuff
    printHubs(hubList, args['p']) # could maybe make this optional
    hubProbs = getProbs(hubList, M, n, allInteractions) #contains (idTag, M, n, N, x, prb)
    corrected = corrProbs(hubProbs) # correct for multiple testing
    printAllProbs(corrected, nodeDict, args['p'], args['v']) #specify prefix, verbose flag


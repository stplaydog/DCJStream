JFLAGS = -O -Xlint
JC = javac
.SUFFIXES: .java .class
.java.class:
	$(JC) $(JFLAGS) $*.java

CLASSES = \
	  detector/Detector.java \
	  detector/DetectorCir.java \
	  detector/DetectorLin.java \
	  graphs/CList.java \
	  graphs/Graph.java \
	  graphs/GraphCir.java \
	  graphs/GraphLin.java \
	  order/GeneOrder.java \
	  order/GeneOrderCir.java \
	  order/GeneOrderLin.java \
	  scheduler/LoadBalancer.java \
	  solvers/ASMSolver.java \
	  solvers/ExactSolver.java \
	  solvers/ExactThread.java \
	  solvers/HeuSolver.java \
	  structs/Elem.java \
	  structs/ElemNoBuffer.java \
	  structs/ElemWithBuffer.java \
	  structs/SearchList.java \
	  tools/Const.java \
	  tools/Func.java \
	  tools/Info.java \
	  tools/Params.java \
	  gasts/Solution.java \
	  gasts/GraphXu.java \
	  gasts/GASTS.java \
	  gasts/GAS_Phylogeny.java \
	  gasts/Heuristic_min.java \
	  gasts/Linearalization.java \
	  gasts/Newick.java \
	  gasts/Node.java \
	  gasts/NodePair.java \
	  gasts/Poisson.java \
	  gasts/ReversalDistance.java \
	  gasts/SimulatedGenome.java \
	  gasts/Simulation.java \
	  gasts/SimuNode.java \
	  gasts/Solution.java \
	  gasts/WeightedASM.java \
	  gasts/Genome.java \
	  gasts/EF.java \
	  gasts/Constant.java \
	  gasts/ASM.java \
	  gasts/Adjacency.java \
	  gasts/AdjNode.java \
	  gasts/Adequate.java \
	  DCJ.java

default: classes

classes: $(CLASSES:.java=.class)

clean:
	        $(RM) *.class
	        $(RM) detector/*.class
	        $(RM) graphs/*.class
	        $(RM) order/*.class
	        $(RM) scheduler/*.class
	        $(RM) solvers/*.class
	        $(RM) structs/*.class
	        $(RM) tools/*.class
	        $(RM) gasts/*.class


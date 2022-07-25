

# Initial Setup
CXX=g++
MKDIR=mkdir -p

CND_DISTDIR=bin
CND_BUILDDIR=Build

OBJECTDIR=${CND_BUILDDIR}
DISTDIR=${CND_DISTDIR}

SRC_FOLDER=.#current directory of the makefile

# Required Included Libraries and includes 
LOCAL_DLIBLIBS=-lpthread -lnsl 
LOCAL_DLIBDIR=-I${SRC_FOLDER}/lib/dlib-19.4

# Setup variables for make to reuse existing make functionality
CXXFLAGS=-std=c++11 -O3 ${LOCAL_DLIBDIR}
LDFLAGS=-L/usr/local/lib  
LDLIBS=-lm -ldl -lgmp ${LOCAL_DLIBLIBS}
#-lboost_atomic -lcolamd 

VPATH=src:${OBJECTDIR}


# OBJECTS
OBJECTS= \
	${OBJECTDIR}/dataset.o \
	${OBJECTDIR}/edge.o \
	${OBJECTDIR}/experiment.o \
	${OBJECTDIR}/functionCoupling_MI.o \
	${OBJECTDIR}/graph.o \
	${OBJECTDIR}/groupwiseExperiment.o \
	${OBJECTDIR}/linearAssignmentSolver.o \
	${OBJECTDIR}/main.o \
	${OBJECTDIR}/matcher.o \
	${OBJECTDIR}/matrixSimilarity.o \
	${OBJECTDIR}/node.o \
	${OBJECTDIR}/normalizer.o \
	${OBJECTDIR}/parameterSet.o \
	${OBJECTDIR}/strFuncCoupling_MI.o \
	${OBJECTDIR}/structureCoupling_MI.o \
	${OBJECTDIR}/subjectwiseExperiment.o \
	${OBJECTDIR}/test.o \
	${OBJECTDIR}/tools.o

# Recipes
.PHONY : clean

${DISTDIR}/brainMatch : ${OBJECTS}
	${MKDIR} ${DISTDIR}
	${LINK.cc} -o $@ $^ ${LDLIBS}

${OBJECTS} : | ${OBJECTDIR}

${OBJECTDIR}/%.o : %.cpp 
	${RM} "$@.d"
	${COMPILE.cc} -MMD -MP -MF "$@.d" -o $@ $<


${OBJECTDIR} :
	${MKDIR} ${OBJECTDIR}

clean :
	${RM} ${OBJECTDIR}/*.d ${OBJECTDIR}/*.o


# DO NOT DELETE
build/Debug_Special/GNU-Linux-x86/dataset.o: \
	dataset.cpp dataset.h \
	graph.h edge.h node.h geometry.h utility.h \
	normalizer.h parameterSet.h initializerParameters.h
build/Debug_Special/GNU-Linux-x86/edge.o: \
	edge.cpp edge.h node.h \
	utility.h geometry.h 
build/Debug_Special/GNU-Linux-x86/experiment.o: \
	experiment.cpp \
	experiment.h dataset.h graph.h edge.h node.h \
	geometry.h utility.h normalizer.h parameterSet.h \
	initializerParameters.h matcher.h \
	matrixInitializer.h tools.h
build/Debug_Special/GNU-Linux-x86/functionCoupling_MI.o: \
	functionCoupling_MI.cpp functionCoupling_MI.h matrixInitializer.h \
	graph.h edge.h node.h \
	geometry.h utility.h normalizer.h \
	initializerParameters.h parameterSet.h
build/Debug_Special/GNU-Linux-x86/graph.o: \
	graph.cpp graph.h \
	edge.h node.h \
	utility.h geometry.h 
build/Debug_Special/GNU-Linux-x86/groupwiseExperiment.o: \
    groupwiseExperiment.cpp groupwiseExperiment.h experiment.h \
    dataset.h graph.h edge.h node.h \
    geometry.h utility.h \
    normalizer.h parameterSet.h initializerParameters.h \
    matcher.h matrixInitializer.h matrixSimilarity.h
build/Debug_Special/GNU-Linux-x86/linearAssignmentSolver.o: \
    linearAssignmentSolver.cpp linearAssignmentSolver.h graph.h \
    edge.h node.h geometry.h utility.h \
    matrixInitializer.h normalizer.h \
    initializerParameters.h parameterSet.h \
	${SRC_FOLDER}/lib/dlib-19.4/dlib/optimization/max_cost_assignment.h \
	${SRC_FOLDER}/lib/dlib-19.4/dlib/optimization/max_cost_assignment_abstract.h 
build/Debug_Special/GNU-Linux-x86/main.o: \
	main.cpp graph.h edge.h \
    node.h geometry.h utility.h dataset.h normalizer.h \
    parameterSet.h initializerParameters.h matcher.h \
    matrixInitializer.h experiment.h \
    subjectwiseExperiment.h groupwiseExperiment.h \
    tools.h test.h
build/Debug_Special/GNU-Linux-x86/matcher.o: \
	matcher.cpp matcher.h \
	dataset.h graph.h edge.h node.h geometry.h \
	utility.h normalizer.h parameterSet.h \
	initializerParameters.h matrixInitializer.h \
	strFuncCoupling_MI.h functionCoupling_MI.h \
	structureCoupling_MI.h
build/Debug_Special/GNU-Linux-x86/matrixSimilarity.o: \
	matrixSimilarity.cpp matrixSimilarity.h \
	dataset.h graph.h edge.h node.h geometry.h \
	utility.h normalizer.h parameterSet.h \
	initializerParameters.h
build/Debug_Special/GNU-Linux-x86/node.o: \
	node.cpp  node.h \
	geometry.h utility.h
build/Debug_Special/GNU-Linux-x86/normalizer.o: \
	normalizer.cpp normalizer.h \
	graph.h edge.h node.h \
	utility.h geometry.h 
build/Debug_Special/GNU-Linux-x86/parameterSet.o: \
	parameterSet.cpp parameterSet.h \
	utility.h geometry.h 
build/Debug_Special/GNU-Linux-x86/strFuncCoupling_MI.o: \
    strFuncCoupling_MI.cpp strFuncCoupling_MI.h matrixInitializer.h \
    graph.h edge.h node.h geometry.h utility.h \
    normalizer.h initializerParameters.h parameterSet.h
build/Debug_Special/GNU-Linux-x86/structureCoupling_MI.o: \
    structureCoupling_MI.cpp structureCoupling_MI.h matrixInitializer.h \
    graph.h edge.h node.h geometry.h utility.h \
    normalizer.h initializerParameters.h parameterSet.h
build/Debug_Special/GNU-Linux-x86/subjectwiseExperiment.o: \
    subjectwiseExperiment.cpp subjectwiseExperiment.h \
    experiment.h dataset.h graph.h edge.h node.h \
    geometry.h utility.h normalizer.h parameterSet.h \
    initializerParameters.h   \
    matcher.h matrixSimilarity.h matrixInitializer.h 
build/Debug_Special/GNU-Linux-x86/test.o: \
	test.cpp test.h \
	graph.h edge.h node.h \
	utility.h geometry.h 
build/Debug_Special/GNU-Linux-x86/tools.o: \
	tools.cpp tools.h \
	utility.h geometry.h dataset.h graph.h edge.h \
	node.h normalizer.h parameterSet.h \
	initializerParameters.h

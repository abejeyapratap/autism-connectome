/* 
 * File:   normalizer.cpp
 * Author: yusuf
 *
 * Generated on October 26, 2016, 2:58 PM
 */

#include "normalizer.h"
#include <vector>
#include <iostream>
#include "node.h"
#include "utility.h"

using std::cout;
using std::endl;

//<editor-fold defaultstate="collapsed" desc=" constructors X 2">
Normalizer::Normalizer()
{
   //initialize variables to default values which will not affect the actual values of the 
   //features if we do not wish to use normalization.
   minSpatial=0;
   maxSpatial=-1;
   spatialRange=1;
   
   minStrNodeStrength=0;
   maxStrNodeStrength=-1;
   strNodeStrengthRange=1;
   
   for(int i=0;i<NUM_FUNC_NODE_FEATURES;i++)
   {
      minFuncNodeStrength[i]=0;
      maxFuncNodeStrength[i]=-1;
      funcNodeStrengthRange[i]=1;
   }

   minStructuralConnectedness=0;
   maxStructuralConnectedness=-1;
   structuralConnectednessRange=1;
   
   minSpatialDistance=0;
   maxSpatialDistance=-1;
   spatialDistanceRange=1;
   
   minFunctionalCorrelation = 0;
   maxFunctionalCorrelation = -1;
   functionalCorrelationRange = 1;
}

Normalizer::Normalizer(const Normalizer& other)
{
   //initialize variables to default values which will not affect the actual values of the 
   //features if we do not wish to use normalization.   
   minSpatial=other.minSpatial;
   maxSpatial=other.maxSpatial;
   spatialRange=other.spatialRange;
   
   minStrNodeStrength=other.minStrNodeStrength;
   maxStrNodeStrength=other.maxStrNodeStrength;
   strNodeStrengthRange=other.strNodeStrengthRange;
   
   for(int i=0;i<NUM_FUNC_NODE_FEATURES;i++)
   {
      minFuncNodeStrength[i]=other.minFuncNodeStrength[i];
      maxFuncNodeStrength[i]=other.maxFuncNodeStrength[i];
      funcNodeStrengthRange[i]=other.funcNodeStrengthRange[i];
   }

   minStructuralConnectedness=other.minStructuralConnectedness;
   maxStructuralConnectedness=other.maxStructuralConnectedness;
   structuralConnectednessRange=other.structuralConnectednessRange;
   
   minSpatialDistance=other.minSpatialDistance;
   maxSpatialDistance=other.maxSpatialDistance;
   spatialDistanceRange=other.spatialDistanceRange;
   
   minFunctionalCorrelation = other.minFunctionalCorrelation;
   maxFunctionalCorrelation = other.maxFunctionalCorrelation;
   functionalCorrelationRange = other.functionalCorrelationRange;
}

//</editor-fold>

//This class keeps record of the data range for features of a set of graphs
//Note that, this class does not modify the edges or he nodes of the graphs.
//It is programmers responsibility to make such adjustments later on based on the 
//data ranges that are calculated by this class.
void Normalizer::normalizeData(std::map<int,Graph>& graphs)
{
   //normalize node features
//   normalizeVolume(graphs);
//   normalizeSpatial(graphs);
//   normalizeStructuralNodeStrength(graphs);
//   normalizeFunctionalNodeStrength(graphs);
   normalizeNodeFeatures(graphs);

   //normalize edge features
   normalizeSpatialDistance(graphs);
   normalizeStructuralConnectedness(graphs);
   normalizeFunctionalCorrelation(graphs);
   
//   normalizeEdgeWeightPairwise(graphs);
//   normalizeEdgeDistancePairwise(graphs);
}

//<editor-fold defaultstate="collapsed" desc=" normalizers for a set of graphs">

//<editor-fold defaultstate="collapsed" desc=" normalizers for node features over a set of graphs: normalizeNodeFeatures(), normalizeVolume(), normalizeSpatial(), normalizeFiberCount() ">

void Normalizer::normalizeSpatial(std::map<int,Graph>& graphs)
{
   float spatial;
   std::map<int,Graph>::iterator graphsBegin = graphs.begin();
   std::map<int,Graph>::iterator graphsEnd = graphs.end();
   std::vector<Node>::iterator nodePtr1 = graphsBegin->second.getNodeIteratorBegin();
   std::vector<Node>::iterator nodePtr2 = (std::next(graphsBegin,1))->second.getNodeIteratorBegin();
   
   minSpatial = nodePtr1->calculateDistance(*nodePtr2);
   maxSpatial = minSpatial;
   for(std::map<int,Graph>::iterator g1=graphsBegin; g1!=graphsEnd; g1++)
   {
      std::vector<Node>::iterator nodePtr1Begin = g1->second.getNodeIteratorBegin();
      std::vector<Node>::iterator nodePtr1End = g1->second.getNodeIteratorEnd();
      
      for(std::map<int,Graph>::iterator g2=std::next(g1,1); g2!=graphsEnd; g2++)
      {
         std::vector<Node>::iterator nodePtr2Begin = g2->second.getNodeIteratorBegin();
         std::vector<Node>::iterator nodePtr2End = g2->second.getNodeIteratorEnd();
         
         for(nodePtr1=nodePtr1Begin;nodePtr1!=nodePtr1End;nodePtr1++)
         {
            for(nodePtr2=nodePtr2Begin;nodePtr2!=nodePtr2End;nodePtr2++)
            {
               spatial = nodePtr1->calculateDistance(*nodePtr2);
               if(spatial>maxSpatial)
                  maxSpatial = spatial;
               else if(spatial<minSpatial)
                  minSpatial = spatial;
            }
         }
      }
   }
   spatialRange = maxSpatial - minSpatial;
}

void Normalizer::normalizeStructuralNodeStrength(std::map<int,Graph>& graphs)
{
   float fiberCount;
   std::vector<Node>::iterator nodePtr = graphs.begin()->second.getNodeIteratorBegin();
   std::map<int,Graph>::iterator graphsBegin = graphs.begin();
   std::map<int,Graph>::iterator graphsEnd = graphs.end();
   
   minStrNodeStrength = nodePtr->getFeature(1);
   maxStrNodeStrength = minStrNodeStrength;
   for(std::map<int,Graph>::iterator g1=graphsBegin; g1!=graphsEnd; g1++)
   {
      std::vector<Node>::iterator g1NodesBegin = g1->second.getNodeIteratorBegin();
      std::vector<Node>::iterator g1NodesEnd = g1->second.getNodeIteratorEnd();
      
      for(nodePtr=g1NodesBegin;nodePtr!=g1NodesEnd;nodePtr++)
      {
         fiberCount = nodePtr->getFeature(1);
         if(fiberCount>maxStrNodeStrength)
            maxStrNodeStrength = fiberCount;
         else if(fiberCount<minStrNodeStrength)
            minStrNodeStrength = fiberCount;
      }
   }
   strNodeStrengthRange = maxStrNodeStrength - minStrNodeStrength;
}

void Normalizer::normalizeFunctionalNodeStrength(std::map<int,Graph>& graphs)
{
   float funcNodeStrength[NUM_FUNC_NODE_FEATURES];
   std::map<int,Graph>::iterator graphsBegin = graphs.begin();
   std::map<int,Graph>::iterator graphsEnd = graphs.end();
   std::vector<Node>::iterator nodePtr1 = graphsBegin->second.getNodeIteratorBegin();
   std::vector<Node>::iterator nodePtr2 = (std::next(graphsBegin,1))->second.getNodeIteratorBegin();
   
   for(int i=0;i<NUM_FUNC_NODE_FEATURES;i++)
   {
      minFuncNodeStrength[i] = Utility::absoluteValue(nodePtr1->getFeature(Node::FUNC_NODE_STRENGTH_NEG+i)
                                    - nodePtr2->getFeature(Node::FUNC_NODE_STRENGTH_NEG+i));
      maxFuncNodeStrength[i] = minFuncNodeStrength[i];
   }
   
   for(std::map<int,Graph>::iterator g1=graphsBegin; g1!=graphsEnd; g1++)
   {
      std::vector<Node>::iterator nodePtr1Begin = g1->second.getNodeIteratorBegin();
      std::vector<Node>::iterator nodePtr1End = g1->second.getNodeIteratorEnd();
      
      for(std::map<int,Graph>::iterator g2=std::next(g1,1); g2!=graphsEnd; g2++)
      {
         std::vector<Node>::iterator nodePtr2Begin = g2->second.getNodeIteratorBegin();
         std::vector<Node>::iterator nodePtr2End = g2->second.getNodeIteratorEnd();
         
         for(nodePtr1=nodePtr1Begin;nodePtr1!=nodePtr1End;nodePtr1++)
         {
            for(nodePtr2=nodePtr2Begin;nodePtr2!=nodePtr2End;nodePtr2++)
            {
               for(int i=0;i<NUM_FUNC_NODE_FEATURES;i++)
               {
                  funcNodeStrength[i] = Utility::absoluteValue(nodePtr1->getFeature(Node::FUNC_NODE_STRENGTH_NEG+i)
                                             - nodePtr2->getFeature(Node::FUNC_NODE_STRENGTH_NEG+i));
                  if(funcNodeStrength[i]>maxFuncNodeStrength[i])
                     maxFuncNodeStrength[i] = funcNodeStrength[i];
                  else if(funcNodeStrength[i]<minFuncNodeStrength[i])
                     minFuncNodeStrength[i] = funcNodeStrength[i];
               }
            }
         }
      }
   }
   
   for(int i=0;i<NUM_FUNC_NODE_FEATURES;i++)
      funcNodeStrengthRange[i] = maxFuncNodeStrength[i] - minFuncNodeStrength[i];
}

//function that normalizes all node features at one place (to avoid repetition of the for loops)
//over a set of graphs
void Normalizer::normalizeNodeFeatures(std::map<int,Graph>& graphs)
{
   float fiberCount,spatial,volume,funcNodeStrength[NUM_FUNC_NODE_FEATURES];
   
   std::map<int,Graph>::iterator graphsBegin = graphs.begin();
   std::map<int,Graph>::iterator graphsEnd = graphs.end();
   std::vector<Node>::iterator nodePtr1 = graphsBegin->second.getNodeIteratorBegin();
   std::vector<Node>::iterator nodePtr2 = (std::next(graphsBegin,1))->second.getNodeIteratorBegin();//first node of the secondGraph

   minSpatial = nodePtr1->calculateDistance(*nodePtr2);
   maxSpatial = minSpatial;
   
   minStrNodeStrength = Utility::absoluteValue(nodePtr1->getFeature(Node::STR_NODE_STRENGTH)-nodePtr2->getFeature(Node::STR_NODE_STRENGTH));
   maxStrNodeStrength = minStrNodeStrength;
   
   for(int i=0;i<NUM_FUNC_NODE_FEATURES;i++)
   {
      minFuncNodeStrength[i] = Utility::absoluteValue(nodePtr1->getFeature(Node::FUNC_NODE_STRENGTH_NEG+i)
                                    - nodePtr2->getFeature(Node::FUNC_NODE_STRENGTH_NEG+i));
      maxFuncNodeStrength[i] = minFuncNodeStrength[i];
   }

   for(std::map<int,Graph>::iterator g1=graphsBegin; g1!=graphsEnd; g1++)
   {
      std::vector<Node>::iterator nodePtr1Begin = g1->second.getNodeIteratorBegin();
      std::vector<Node>::iterator nodePtr1End = g1->second.getNodeIteratorEnd();
      
      for(std::map<int,Graph>::iterator g2=std::next(g1,1); g2!=graphsEnd; g2++)
      {
         std::vector<Node>::iterator nodePtr2Begin = g2->second.getNodeIteratorBegin();
         std::vector<Node>::iterator nodePtr2End = g2->second.getNodeIteratorEnd();
         
         for(nodePtr1=nodePtr1Begin;nodePtr1!=nodePtr1End;nodePtr1++)
         {
            for(nodePtr2=nodePtr2Begin;nodePtr2!=nodePtr2End;nodePtr2++)
            {
               spatial = nodePtr1->calculateDistance(*nodePtr2);
               if(spatial>maxSpatial)
                  maxSpatial = spatial;
               else if(spatial<minSpatial)
                  minSpatial = spatial;
               
               fiberCount = Utility::absoluteValue(nodePtr1->getFeature(Node::STR_NODE_STRENGTH)-nodePtr2->getFeature(Node::STR_NODE_STRENGTH));
               if(fiberCount>maxStrNodeStrength)
                  maxStrNodeStrength = fiberCount;
               else if(fiberCount<minStrNodeStrength)
                  minStrNodeStrength = fiberCount;
               
               for(int i=0;i<NUM_FUNC_NODE_FEATURES;i++)
               {
                  funcNodeStrength[i] = Utility::absoluteValue(nodePtr1->getFeature(Node::FUNC_NODE_STRENGTH_NEG+i)
                                             - nodePtr2->getFeature(Node::FUNC_NODE_STRENGTH_NEG+i));
                  if(funcNodeStrength[i]>maxFuncNodeStrength[i])
                     maxFuncNodeStrength[i] = funcNodeStrength[i];
                  else if(funcNodeStrength[i]<minFuncNodeStrength[i])
                     minFuncNodeStrength[i] = funcNodeStrength[i];
               }
            }
         }
      }
   }

   spatialRange = maxSpatial - minSpatial;
   strNodeStrengthRange = maxStrNodeStrength - minStrNodeStrength;
   for(int i=0;i<NUM_FUNC_NODE_FEATURES;i++)
      funcNodeStrengthRange[i] = maxFuncNodeStrength[i] - minFuncNodeStrength[i];
}

//</editor-fold>

//<editor-fold defaultstate="collapsed" desc=" normalizers for edges over a set of graphs: normalizeStructuralConnectedness(), normalizeSpatialDistance(), normalizeFunctionalCorrelation() ">

//NOTE: we are interested in the minimum nonzero feature in each set
//However, while determining the range, one might also be interested in taking the range to be [0,maxValue]
//in cases where minimum nonzero value will be quite small.

void Normalizer::normalizeStructuralConnectedness(std::map<int,Graph>& graphs)
{
   float edgeWeight;
   std::map<int,Graph>::iterator graphsBegin = graphs.begin();
   std::map<int,Graph>::iterator graphsEnd = graphs.end();
   
   minStructuralConnectedness = Utility::INF;
   maxStructuralConnectedness = Utility::N_INF;
   
   for(std::map<int,Graph>::iterator g1=graphsBegin; g1!=graphsEnd; g1++)
   {
      std::vector<Edge>::iterator g1EdgesBegin = g1->second.getEdgeIteratorBegin();
      std::vector<Edge>::iterator g1EdgesEnd = g1->second.getEdgeIteratorEnd();
      
      for(std::vector<Edge>::iterator edgePtr=g1EdgesBegin;edgePtr!=g1EdgesEnd;edgePtr++)
      {
         edgeWeight = edgePtr->getFeature(Edge::STRUCTURAL_CONNECTIVITY);
         if(edgeWeight>maxStructuralConnectedness)
            maxStructuralConnectedness = edgeWeight;
         else if(edgeWeight<minStructuralConnectedness && edgeWeight>Utility::EPSILON)
            minStructuralConnectedness = edgeWeight;
      }
   }

   minStructuralConnectedness = (minStructuralConnectedness < Utility::EPSILON ? Utility::EPSILON : minStructuralConnectedness);
   maxStructuralConnectedness = (maxStructuralConnectedness < Utility::EPSILON ? Utility::EPSILON : maxStructuralConnectedness);
   
   structuralConnectednessRange = maxStructuralConnectedness - minStructuralConnectedness;
}

void Normalizer::normalizeSpatialDistance(std::map<int,Graph>& graphs)
{
   float edgeDistance;
   std::map<int,Graph>::iterator graphsBegin = graphs.begin();
   std::map<int,Graph>::iterator graphsEnd = graphs.end();
   
   minSpatialDistance = Utility::INF;
   maxSpatialDistance = Utility::N_INF;
   
   for(std::map<int,Graph>::iterator g1=graphsBegin; g1!=graphsEnd; g1++)
   {
      std::vector<Edge>::iterator g1EdgesBegin = g1->second.getEdgeIteratorBegin();
      std::vector<Edge>::iterator g1EdgesEnd = g1->second.getEdgeIteratorEnd();
      
      for(std::vector<Edge>::iterator edgePtr=g1EdgesBegin;edgePtr!=g1EdgesEnd;edgePtr++)
      {
         edgeDistance = edgePtr->getFeature(Edge::SPATIAL_DISTANCE);
         if(edgeDistance>maxSpatialDistance)
            maxSpatialDistance = edgeDistance;
         else if(edgeDistance<minSpatialDistance && edgeDistance>Utility::EPSILON)
            minSpatialDistance = edgeDistance;
      }
   }

   minSpatialDistance = (minSpatialDistance < Utility::EPSILON ? Utility::EPSILON : minSpatialDistance);
   maxSpatialDistance = (maxSpatialDistance < Utility::EPSILON ? Utility::EPSILON : maxSpatialDistance);
   
   spatialDistanceRange = maxSpatialDistance - minSpatialDistance;
}

void Normalizer::normalizeFunctionalCorrelation(std::map<int,Graph>& graphs)
{
   float functionalCorrelation;
   std::map<int,Graph>::iterator graphsBegin = graphs.begin();
   std::map<int,Graph>::iterator graphsEnd = graphs.end();
   
   minFunctionalCorrelation = Utility::INF;
   maxFunctionalCorrelation = Utility::N_INF;
   
   for(std::map<int,Graph>::iterator g1=graphsBegin; g1!=graphsEnd; g1++)
   {
      std::vector<Edge>::iterator g1EdgesBegin = g1->second.getEdgeIteratorBegin();
      std::vector<Edge>::iterator g1EdgesEnd = g1->second.getEdgeIteratorEnd();
      
      for(std::vector<Edge>::iterator edgePtr=g1EdgesBegin;edgePtr!=g1EdgesEnd;edgePtr++)
      {
         functionalCorrelation = edgePtr->getFeature(Edge::FUNCTIONAL_CONNECTIVITY);
         if(functionalCorrelation>maxFunctionalCorrelation)
            maxFunctionalCorrelation = functionalCorrelation;
         else if(functionalCorrelation<minFunctionalCorrelation && functionalCorrelation>Utility::EPSILON)
            minFunctionalCorrelation = functionalCorrelation;
      }
   }

   minFunctionalCorrelation = (minFunctionalCorrelation < Utility::EPSILON ? Utility::EPSILON : minFunctionalCorrelation);
   maxFunctionalCorrelation = (maxFunctionalCorrelation < Utility::EPSILON ? Utility::EPSILON : maxFunctionalCorrelation);
   
   functionalCorrelationRange = maxFunctionalCorrelation - minFunctionalCorrelation;
}

//</editor-fold>
//</editor-fold>

//<editor-fold defaultstate="collapsed" desc=" normalizers for a single graph">

//<editor-fold defaultstate="collapsed" desc=" normalizers for edges over a single graph: normalizeStructuralConnectedness(), normalizeSpatialDistance(), normalizeFunctionalCorrelation() ">

void Normalizer::normalizeStructuralConnectedness(Graph &graph)
{
   float edgeWeight;

   minStructuralConnectedness = Utility::INF;
   maxStructuralConnectedness = Utility::N_INF;
   
   for(std::vector<Edge>::iterator edgePtr=graph.getEdgeIteratorBegin();edgePtr!=graph.getEdgeIteratorEnd();edgePtr++)
   {
      edgeWeight = edgePtr->getFeature(Edge::STRUCTURAL_CONNECTIVITY);
      if(edgeWeight>maxStructuralConnectedness)
         maxStructuralConnectedness = edgeWeight;
      else if(edgeWeight<minStructuralConnectedness && edgeWeight>0)
         minStructuralConnectedness = edgeWeight;
   }

   minStructuralConnectedness = (minStructuralConnectedness < Utility::EPSILON ? Utility::EPSILON : minStructuralConnectedness);
   maxStructuralConnectedness = (maxStructuralConnectedness < Utility::EPSILON ? Utility::EPSILON : maxStructuralConnectedness);
   
   structuralConnectednessRange = maxStructuralConnectedness - minStructuralConnectedness;
}

void Normalizer::normalizeSpatialDistance(Graph &graph)
{
   float edgeDistance;
   
   minSpatialDistance = Utility::INF;
   maxSpatialDistance = Utility::N_INF;
   
   for(std::vector<Edge>::iterator edgePtr=graph.getEdgeIteratorBegin();edgePtr!=graph.getEdgeIteratorEnd();edgePtr++)
   {
      edgeDistance = edgePtr->getFeature(Edge::SPATIAL_DISTANCE);
      if(edgeDistance>maxSpatialDistance)
         maxSpatialDistance = edgeDistance;
      else if(edgeDistance<minSpatialDistance && edgeDistance>0)
         minSpatialDistance = edgeDistance;
   }

   minSpatialDistance = (minSpatialDistance < Utility::EPSILON ? Utility::EPSILON : minSpatialDistance);
   maxSpatialDistance = (maxSpatialDistance < Utility::EPSILON ? Utility::EPSILON : maxSpatialDistance);
   
   spatialDistanceRange = maxSpatialDistance - minSpatialDistance;
}

void Normalizer::normalizeFunctionalCorrelation(Graph &graph)
{
   float functionalCorrelation;
   
   minFunctionalCorrelation = Utility::INF;
   maxFunctionalCorrelation = Utility::N_INF;
   
   for(std::vector<Edge>::iterator edgePtr=graph.getEdgeIteratorBegin();edgePtr!=graph.getEdgeIteratorEnd();edgePtr++)
   {
      functionalCorrelation = edgePtr->getFeature(Edge::FUNCTIONAL_CONNECTIVITY);
      if(functionalCorrelation>maxFunctionalCorrelation)
         maxFunctionalCorrelation = functionalCorrelation;
      else if(functionalCorrelation<minFunctionalCorrelation && functionalCorrelation>0)
         minFunctionalCorrelation = functionalCorrelation;
   }

   minFunctionalCorrelation = (minFunctionalCorrelation < Utility::EPSILON ? Utility::EPSILON : minFunctionalCorrelation);
   maxFunctionalCorrelation = (maxFunctionalCorrelation < Utility::EPSILON ? Utility::EPSILON : maxFunctionalCorrelation);
   
   functionalCorrelationRange = maxFunctionalCorrelation - minFunctionalCorrelation;
}

//</editor-fold>

//</editor-fold>

//<editor-fold defaultstate="collapsed" desc=" pairwise normalizers for edges over a set of graphs: normalizeEdgeWeightPairwise(), normalizeEdgeDistancePairwise() ">

//this function calculates the range for the w_{p,q}*d_{a,b}
//it runs over all the pairs of graphs and pairs of edges in corresponding 
//graphs to find the min and max of this multiplication term.
//it takes long to run it and the result is not too much different than the 
//less robust version of this function, which calculates the range based on the
//min and max edge weight only, and then calculates their multiplication for
//determining the range of w_{p,q}*d_{a,b}
void Normalizer::normalizeStructuralConnectednessPairwise(std::map<int,Graph>& graphs)
{
   float edgeWeight;
   
   std::map<int,Graph>::iterator graphsBegin = graphs.begin();
   std::map<int,Graph>::iterator graphsEnd = graphs.end();
   std::vector<Edge>::iterator edgePtr1 = graphsBegin->second.getEdgeIteratorBegin();
   std::vector<Edge>::iterator edgePtr2 = std::next(graphsBegin,1)->second.getEdgeIteratorBegin();
   
   float weight1=edgePtr1->getFeature(Edge::STRUCTURAL_CONNECTIVITY);
   float weight2=edgePtr2->getFeature(Edge::STRUCTURAL_CONNECTIVITY);
   
   weight1 = (weight1 < 1 ? 1 : weight1);
   weight2 = (weight2 < 1 ? 1 : weight2);
   
   minStructuralConnectedness = weight1/weight2;
   maxStructuralConnectedness = minStructuralConnectedness;
   for(std::map<int,Graph>::iterator g1=graphsBegin; g1!=graphsEnd; g1++)
   {
      std::vector<Edge>::iterator g1EdgesBegin = g1->second.getEdgeIteratorBegin();
      std::vector<Edge>::iterator g1EdgesEnd = g1->second.getEdgeIteratorEnd();
         
      for(std::map<int,Graph>::iterator g2=std::next(g1,1); g2!=graphsEnd; g2++)
      {
         std::vector<Edge>::iterator g2EdgesBegin = g2->second.getEdgeIteratorBegin();
         std::vector<Edge>::iterator g2EdgesEnd = g2->second.getEdgeIteratorEnd();

         for(edgePtr1=g1EdgesBegin;edgePtr1!=g1EdgesEnd;edgePtr1++)
         {
            for(edgePtr2=g2EdgesBegin;edgePtr2!=g2EdgesEnd;edgePtr2++)
            {
               weight1 = (edgePtr1->getFeature(Edge::STRUCTURAL_CONNECTIVITY) < 1 ? 1 : edgePtr1->getFeature(Edge::STRUCTURAL_CONNECTIVITY));
               weight2 = (edgePtr2->getFeature(Edge::STRUCTURAL_CONNECTIVITY) < 1 ? 1 : edgePtr2->getFeature(Edge::STRUCTURAL_CONNECTIVITY));
   
               edgeWeight = weight1/weight2;
               if(edgeWeight>maxStructuralConnectedness)
                  maxStructuralConnectedness = edgeWeight;
               else if(edgeWeight<minStructuralConnectedness)
                  minStructuralConnectedness = edgeWeight;
            }
         }
      }
   }
   structuralConnectednessRange = maxStructuralConnectedness - minStructuralConnectedness;
}

void Normalizer::normalizeSpatialDistancePairwise(std::map<int,Graph>& graphs)
{
   float edgeDistance;
   
   std::map<int,Graph>::iterator graphsBegin = graphs.begin();
   std::map<int,Graph>::iterator graphsEnd = graphs.end();
   std::vector<Edge>::iterator edgePtr1 = graphsBegin->second.getEdgeIteratorBegin();
   std::vector<Edge>::iterator edgePtr2 = std::next(graphsBegin,1)->second.getEdgeIteratorBegin();
   
   float distance1=edgePtr1->getFeature(Edge::SPATIAL_DISTANCE);
   float distance2=edgePtr2->getFeature(Edge::SPATIAL_DISTANCE);
   
   distance1 = (distance1 < 1 ? 1 : distance1);
   distance2 = (distance2 < 1 ? 1 : distance2);
   
   minStructuralConnectedness = distance1/distance2;
   maxSpatialDistance = minSpatialDistance;
   for(std::map<int,Graph>::iterator g1=graphsBegin; g1!=graphsEnd; g1++)
   {
      std::vector<Edge>::iterator g1EdgesBegin = g1->second.getEdgeIteratorBegin();
      std::vector<Edge>::iterator g1EdgesEnd = g1->second.getEdgeIteratorEnd();
         
      for(std::map<int,Graph>::iterator g2=std::next(g1,1); g2!=graphsEnd; g2++)
      {
         std::vector<Edge>::iterator g2EdgesBegin = g2->second.getEdgeIteratorBegin();
         std::vector<Edge>::iterator g2EdgesEnd = g2->second.getEdgeIteratorEnd();

         for(edgePtr1=g1EdgesBegin;edgePtr1!=g1EdgesEnd;edgePtr1++)
         {
            for(edgePtr2=g2EdgesBegin;edgePtr2!=g2EdgesEnd;edgePtr2++)
            {
               distance1 = (edgePtr1->getFeature(Edge::SPATIAL_DISTANCE) < 1 ? 1 : edgePtr1->getFeature(Edge::SPATIAL_DISTANCE));
               distance2 = (edgePtr2->getFeature(Edge::SPATIAL_DISTANCE) < 1 ? 1 : edgePtr2->getFeature(Edge::SPATIAL_DISTANCE));
   
               edgeDistance = distance1/distance2;
               if(edgeDistance>maxSpatialDistance)
                  maxSpatialDistance = edgeDistance;
               else if(edgeDistance<minSpatialDistance)
                  minSpatialDistance = edgeDistance;
            }
         }
      }
   }
   spatialDistanceRange = maxSpatialDistance - minSpatialDistance;
}

//</editor-fold>

//<editor-fold defaultstate="collapsed" desc=" load/save/print functions">
void Normalizer::print()
{
   cout<<"Spatial:("<<minSpatial<<","<<maxSpatial<<","<<spatialRange<<")"<<endl;
   cout<<"Structural node strength:("<<minStrNodeStrength<<","<<maxStrNodeStrength<<","<<strNodeStrengthRange<<")"<<endl;
   cout<<"Functional node strength:("<<minFuncNodeStrength<<","<<maxFuncNodeStrength<<","<<funcNodeStrengthRange<<")"<<endl;
   cout<<"Edge Weight:("<<minStructuralConnectedness<<","<<maxStructuralConnectedness<<","<<structuralConnectednessRange<<")"<<endl;
   cout<<"Edge Distance:("<<minSpatialDistance<<","<<maxSpatialDistance<<","<<spatialDistanceRange<<")"<<endl;
   cout<<"BrainVolume:("<<minBrainVolume<<","<<maxBrainVolume<<","<<brainVolumeRange<<")"<<endl;
}  

void Normalizer::save(std::string filename)
{
   std::ofstream file;

   file.open(filename.c_str());

   file.setf(std::ios::fixed,std::ios::floatfield);
   
   //ranges for node features
	file<<"#minSpatial, maxSpatial, spatialRange"<<endl;
   file<<minSpatial<<"\t"<<maxSpatial<<"\t"<<spatialRange<<endl;
   
   file<<"#minStrNodeStrength, maxStrNodeStrength, strNodeStrengthRange"<<endl;
	file<<minStrNodeStrength<<"\t"<<maxStrNodeStrength<<"\t"<<strNodeStrengthRange<<endl;
   
   file<<"#minFuncNodeStrength, maxFuncNodeStrength, funcNodeStrengthRange -- 1)FUNC_NODE_STRENGTH_NEG_SUM, 2)FUNC_NODE_STRENGTH_POS_SUM, 3)FUNC_NODE_STRENGTH_NEG_CNT, 4)FUNC_NODE_STRENGTH_POS_CNT"<<endl;
   for(int i=0; i<NUM_FUNC_NODE_FEATURES;i++)
   {
      file<<minFuncNodeStrength[i]<<"\t"<<maxFuncNodeStrength[i]<<"\t"<<funcNodeStrengthRange[i]<<endl;
   }
		
   //ranges for edge features
   file<<"#minSpatialDistance, maxSpatialDistance, spatialDistanceRange"<<endl;
   file<<minSpatialDistance<<"\t"<<maxSpatialDistance<<"\t"<<spatialDistanceRange<<endl;
   
   file<<"#minStructuralConnectedness, maxStructuralConnectedness, structuralConnectednessRange"<<endl;
   file<<minStructuralConnectedness<<"\t"<<maxStructuralConnectedness<<"\t"<<structuralConnectednessRange<<endl;
   
   file<<"#minFunctionalCorrelation, maxFunctionalCorrelation, functionalCorrelationRange"<<endl;
   file<<minFunctionalCorrelation<<"\t"<<maxFunctionalCorrelation<<"\t"<<functionalCorrelationRange<<endl;
   
   file<<"#minBrainVolume, maxBrainVolume, brainVolumeRange"<<endl;
   file<<minBrainVolume<<"\t"<<maxBrainVolume<<"\t"<<brainVolumeRange<<endl;
   
   file.close();
}

void Normalizer::load(std::string filename)
{
   std::ifstream file; 
   std::istringstream is,iss; 
   std::string line;  

   file.open(filename.c_str());
   if(file.fail())
   {
     cout << "File not found!!!!   " << filename << "\n";
     return;
   }	

   //read node feature ranges
   getline(file, line);	//ignore comment line
   getline(file, line);	
   is.clear();
   is.str(line);
   is>>minSpatial>>maxSpatial>>spatialRange;
   
   getline(file, line);	//ignore comment line
   getline(file, line);	
   is.clear();
   is.str(line);
   is>>minStrNodeStrength>>maxStrNodeStrength>>strNodeStrengthRange;
   
   getline(file, line);	//ignore comment line
   for(int i=0;i<NUM_FUNC_NODE_FEATURES;i++)
   {
      getline(file, line);	
      is.clear();
      is.str(line);
      is>>minFuncNodeStrength[i]>>maxFuncNodeStrength[i]>>funcNodeStrengthRange[i];
   }
   
   //read edge feature ranges
   getline(file, line);	//ignore comment line
   getline(file, line);	
   is.clear();
   is.str(line);
   is>>minSpatialDistance>>maxSpatialDistance>>spatialDistanceRange;
   
   getline(file, line);	//ignore comment line
   getline(file, line);	
   is.clear();
   is.str(line);
   is>>minStructuralConnectedness>>maxStructuralConnectedness>>structuralConnectednessRange;
   
   getline(file, line);	//ignore comment line
   getline(file, line);	
   is.clear();
   is.str(line);
   is>>minFunctionalCorrelation>>maxFunctionalCorrelation>>functionalCorrelationRange;
   
   //read brain feature ranges
   getline(file, line);	//ignore comment line
   getline(file, line);	
   is.clear();
   is.str(line);
   is>>minBrainVolume>>maxBrainVolume>>brainVolumeRange;
   
   file.close();
}

//</editor-fold>
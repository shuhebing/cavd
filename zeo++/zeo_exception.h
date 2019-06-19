/**
  * This file defines some exception classes for handling
  * exceptions that occur during zeo++ calculations.
  *
  */
#ifndef ZEO_EXCEPTION_H
#define ZEO_EXCEPTION_H

#include <exception>

class InvalidParticlesNumException : public exception{
    public:
        const char * what () const throw (){
            return "Error: Invalid number of particles provided for Voronoi decomposition.\n";
        }
};
class InvalidBoxDimException : public exception{
    public:
        const char * what () const throw (){
            return "Error: Ivalid box dimensions calculated for Voronoi decomposition, please check the input lattices.\n";
        }
};
class HugeGridException : public exception{
   public:
    const char * what () const throw (){
            return "Error: voro++: The number of computational blocks exceeds the maximum value.\n";
        }
};
class AttemptException : public exception{
    public:    
    const char * what () const throw (){
            return "Error: Software did not pass the volume check after Voronoi decomposition.\n";
        }
};
class VoronoiDecompException : public exception{
    public:
        const char * what () const throw (){
            return "Error: Unable to begin Voronoi decomposition.\n";
        }
};
class CoordNumException : public exception{
    public:
        const char * what () const throw (){
            return "Error: Improper number of node coordinates in Voronoi decomposition.\n";
        }
};

class MatrixException : public exception{
    public:
        const char * what () const throw (){
            return "Error: Determinant of provided matrix is 0. Matrix is not invertible.\n";
        }
};

class TRIPLETOperatorException : public exception{
    public:
        const char * what () const throw (){
            return "Error: Invalid index to [] operator for TRIPLET instance.\n";
        }
};

class XYZOperatorException : public exception{
    public:
        const char * what () const throw (){
            return "Error: Invalid index to [] operator for XYZ instance.\n";
        }
};

class AssignDummyException : public exception{
    public:
        const char * what () const throw (){
            return "ERROR: should not call assign_dummy_site() with num_sites!=2.\n";
        }
};

class Zero_NonzeroException : public exception{
    public:
        const char * what () const throw (){
            return "ERROR: Passed trivial equivalence filter but both loops are completely zero (no ratios between elements could be found).\n";
        }
};

class ParseFilenameException : public exception{
    public:
        const char * what () const throw (){
            return "ERROR: Improper input filename No support extension found.\n";
        }
};

class CalculateUnitEdgException : public exception{
    public:
        const char * what () const throw (){
            return "ERROR: found a basic edge length is sufficiently different to the previous length; at the moment, nets with more than one edge length are not handled.\n";
        }
};

class CreateUnitCellException : public exception{
    public:
        const char * what () const throw (){
            return "ERROR: create_unit_cell_from_vectors() called with provided vectors.\n";
        }
};

class zAxisFoundException : public exception{
    public:
        const char * what () const throw (){
            return "ERROR: after setting x and y vectors, more than one vector remains to be assigned to z.\n";
        }
};

class ConclassNetException : public exception{
    public:
        const char * what () const throw (){
            return "ERROR: could not conclass net: no vertex could be found which overlaps periodically with provided vertex and edge.\n";
        }
};

class ConnectionsException : public exception{
    public:
        const char * what () const throw (){
            return "ERROR: there are an odd number of one-way connections between vertices - this should not be the case because connections are expressed redundantly.\n";
        }
};

class VertexEdgeOverlapException : public exception{
    public:
        const char * what () const throw (){
            return "ERROR: determined that vertex %d edge %d overlaps with more than one vertex!\n";
        }
};

class VertexConnectionException : public exception{
    public:
        const char * what () const throw (){
            return "ERROR: could not find corresponding connection for vertex %d edge %d overlapping with vertex!\n";
        }
};

class FitMoleculeToVertexException : public exception{
    public:
        const char * what () const throw (){
            return "ERROR: cannot fit molecule with sites and dummy sites to a vertex with sites and dummy sites!\n";
        }
};

class OutputRotatedMoleculeException : public exception{
    public:
        const char * what () const throw (){
            return "ERROR: could not open output rotated molecule file with provided name!\n";
        }
};

class FindAllUCVectorsException : public exception{
    public:
        const char * what () const throw (){
            return "ERROR: did not find all uc vectors from loop analysis!\n";
        }
};

class OverrideException : public exception{
    public:
        const char * what () const throw (){
            return "ERROR: was expecting to override exactly 1 cell side length value, but were overwritten!\n";
        }
};

class NETException : public exception{
    public:
        const char * what () const throw (){
            return "ERROR: detected that both atom and node flags are used in the input net file - this is currently assumed to indicate an invalid input file\n";
        }
};

#endif
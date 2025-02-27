//
// Created by Shang Zhang on 5/7/17.
//

#ifndef PEBBLEGAMETEST_SITERP_CONT_H
#define PEBBLEGAMETEST_SITERP_CONT_H

#include <vector>
#include <stack>
#include <iostream>
#include <fstream>
#include <math.h>				// Basic math functions
#include <string>
#include "bond.h"
#include "site.h"


class SiteRP_cont {
protected:
    int SIZE; // The number of vertices in the graph
    std::string FileName;
    float BoxLength;
    std::string ReadPATH;
    std::string OutPATH;
public:
    short* pc;                // Creates the pebble count at each vertex.

    int numparts;           // the number of particles (not pebbles) present in the system
    int numbonds;                    // The number of non-redundant bonds (original bonds and crossbraces) in the system
    int rbonds;                        // The number of redundant bonds in the system
    int giantsize_site;                  // The size of the giant rigid cluster
    int giantsize_bond;                  // The size of the giant rigid cluster
    int giantindex;                 // The index for the giant rigid cluster

    std::vector<site> vertices;        // the vertices information, store all the vertices
    //int XSpanLastStatus;             // the last status of whether to have a spanning rigid cluster
    //int YSpanLastStatus;             // the last status of whether to have a spanning rigid cluster


    std::vector<int>* rcluster_site;     // Store all the information about rigid cluster decomposition in sites
    std::vector<int>* thegraph;        // thegraph is the graph of all loaded edges
    std::vector<int>* rgraph;        // rgraph is the graph of redundant bonds that don't take up any edges
    std::vector<int>* undirectedgraph;     // undirected graph for calculating the coordination number
    std::vector<Bond> edges;             //bonds only contains loaded edges
    std::vector<int>* giantrigidcluster;    //giantrigidcluster is the graph for the giant rigid cluster
    std::vector<int>* giantrigidclusterOBC;    //giantrigidcluster is the graph for the giant rigid cluster with open boundary condition
    std::stack<int> placesbeen;            // The list of places been while looking for a pebble
    std::ofstream myfile;                  // The file stream to output the mainly wanted info
    std::ofstream rclusterfile;            // the file stream to output the rigid cluster decomposition info
    std::ofstream coordnumfile;             // the file stream to output the coordination number info
    std::ofstream coordnum_hist_file;       // the file stream to output the coordination number info in histogram
    std::ofstream bondorderparafile;        // the file stream to output the bond order parameter (phi-6)

    
    double* beenthere;    // The array to mark if a vertex is visited in DFS
    
    //std::ifstream InputDataFile;            // read the data from Monte-Carlo data files
    int EMPTY;

    int* CoordDist;        // the coordination number distribution

    double BondCriteria;

    //Constructor and destructor
    //SiteRP_cont(std::vector<site>& VerticesInput);
    //SiteRP_cont(int a) {numparts = a;};
    SiteRP_cont(int sizeIN, std::string FileNameIN, std::string ReadPATHIN, std::string OutPATHIN, double BondLEN){
        SIZE = sizeIN;
        pc = new short[sizeIN];
        rcluster_site = new std::vector<int>[sizeIN];
        thegraph = new std::vector<int>[sizeIN];
        rgraph = new std::vector<int>[sizeIN];
        undirectedgraph = new std::vector<int>[sizeIN];
        giantrigidcluster = new std::vector<int>[sizeIN];
        giantrigidclusterOBC = new std::vector<int>[sizeIN];
        beenthere = new double[sizeIN];
        CoordDist = new int[sizeIN];
        EMPTY = -sizeIN-1;
        FileName = FileNameIN;
        ReadPATH = ReadPATHIN;
        OutPATH = OutPATHIN;
        BondCriteria = BondLEN;
    };
    
    ~SiteRP_cont(){
        delete[] pc;
        delete[] rcluster_site;
        delete[] thegraph;
        delete[] rgraph;
        delete[] undirectedgraph;
        delete[] giantrigidcluster;
        delete[] beenthere;
        delete[] CoordDist;
    };

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // ADDING, REMOVING, AND CHECKING FOR EDGES

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // addedge adds an edge pointing from vertex i to vertex j in thegraph
    void addedge(int i, int j);

    // addredundant adds an edge between i and j in rgraph, the separate graph of redundant edges
    void addredundant(int i, int j);

    // badremoveedge removes an edge pointing from i to j in the graph.
    // I call it "bad" because it looks through all of i's elements instead of
    // using which one might be open in a tree search.
    // This should return an error if there is no edge pointing from i to j.
    int badremoveedge(int i, int j);

    // contains returns 1 if there is a non-redundant edge pointing from i to j (but doesn't check j to i) and 0 otherwise
    bool contains(int i, int j);

    // rcontains returns 1 if there is a redundant edge pointing from i to j (but doesn't check j to i) and 0 otherwise
    bool rcontains(int i, int j);

    // isempty returns 0 if there is any kind of redundant or nonredundant brace pointing in either direction between i and j, and 1 otherwise
    // Modification required: isempty only looks at one type of bond, not the six (or at least, three) kinds of the triangular lattice
    bool isempty(int i, int j);
    bool isempty(int i);

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// ADDING BONDS

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// This takes places been and reverses all the edges on the path.
// So if placesbeen is from 2 to 12 to 5
// It should remove edges from 2 to 12 and from 12 to 5
// and add edges from 5 to 12 and 12 to 2
    void reversepath();

// findpebble executes a depth first search for pebbles in the lattice, starting at i
// and searching all the other vertices
// placesbeen must be empty before findpebble is called
// if findpebble finds a pebble, it should set placesbeen to a path from i to the site with the pebble
// Otherwise, it should leave placesbeen blank.
// It returns 1 if a pebble was found and 0 if it wasn't or if the path wasn't empty to start
    bool findpebble(int i);


// When called with two arguments, findpebble skips over the second site to avoid infinitely swapping pebble
// between the two sites the brace connects, by marking skip as a place that we've already been
    bool findpebble(int i, int skip);


// loadsite looks for pebbles and moves them onto i until i has two pebbles or it stops finding pebbles
    bool loadsite(int i);

// loadsite with two arguments does the same thing, but won't try to take pebbles from skip to move onto i
// Which is to keep the two sites from swapping the three pebbles between themselves endlessly
    bool loadsite(int i, int skip);
// loadsites tries to move pebbles until there are two on both sites i and j
    bool loadsites(int i, int j);

// addbond tries to load the sites. If it succeeds, it adds an edge from i to j and takes a pebble from i. Otherwise, it adds a redundant edge
    void addbond(int i, int j);

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// REPORTING TO SCREEN

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    void setfilestream_rigid_decomposition();
    void setfilestream_coordination_number();
    void setfilestream_bond_order_parameter();
    void log(std::pair<int,int> span);

// listedges lists the edges from site i
    void listedges(int i);

// listalledges lists all the edges from the sites
    void listalledges();

// initgiantrigidcluster() initializes the giantrigidcluster graph
    void initgiantrigidcluster();

// initemptytrigraph() updates numbonds, rbonds, thegraph, rgraph, placesbeen to a triangular graph with no particles or bonds
    int initemptytrigraph();

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// rigid cluster

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    bool isredundant(int i, int j);  //see if the test bond between (i,j) is redundant(dependent)
    bool isbondrigid(Bond &a, Bond &b); //see if the two bonds a and b are rigid to each other
    void rigidcluster(); // mark the rigid clusters


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// THE RIGID CLUSTER INFO (FOR CLUSTER STATISTICS INFO)

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // store the rigid cluster decomposition information in sites (initial rcluster_site before the running of this function)
    void StoreRigidInfoOfSite();

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// SPANNING RIGID CLUSTER

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// We need to pick out the giant rigid cluster from the network and then determine if it is the spanning cluster.

    // a function to calculate the bond length between two particles (without the effect of PBC)
    float BondCoordShift(int prosp, int cl, int CoordDirection);
    
    std::pair<int, int> spanningrcluster();
    
    std::pair<int, int> spanningrclusterNEW();

    // Key difference for SiteRP to SiteRP_cont

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// BUILDING THE NETWORK

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    void ReadVerticeInfo();
    
    void ReadVerticeInfofromdump();
    
    void ReadVerticeInfoFromCSV();

    float distance(site &a, site &b);

    void BuildNetwork(); //build the network with site information
    
    void RigidClusterDecomposition(); //run pebble game and identify the spanning rigid cluster
    
    void StudyArtificialNetwork();

    void RigidAtomWriteBack();
    
    void RigidAtomWriteBack_dump();

    void CoordNumber();
    
    void BondOrderParameter();

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// END OF THE CLASS

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

};


#endif //PEBBLEGAMETEST_SITERP_CONT_H

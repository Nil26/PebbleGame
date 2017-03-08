//
// Created by Zeb & Shang.
//

#include <iostream>
#include <vector>
#include <stack>
#include <fstream>
#include <math.h>				// Basic math functions
#include <time.h>				// Get clock time, to measure total run time

using namespace std;

class Bond{
public:
    pair<int, int> vertices; // the connected vertices
    int RigidIndex;         //  the rigid cluster index
    Bond (int x_input, int y_input)
    {
        vertices = make_pair(x_input,y_input);
        RigidIndex = 0;
    }
    void initBondRigidIndex ()
    {
        RigidIndex = 0;
    }
    pair<int,int> vertex() const { return vertices;};
};

class SiteRP {
    static const int ll = 64;                                                                                    // The number of vertices on a side of the lattice
private:
    static const int size = ll *
                            ll;                                                                                // The number of vertices in the graph
public:
    short pc[size];                // Creates the pebble count at each vertex.
    short occ[size];             // Says whether the site is occupied with a particle
    int numparts;           // the number of particles (not pebbles) present in the system
    int numbonds;                    // The number of non-redundant bonds (original bonds and crossbraces) in the system
    int rbonds;                        // The number of redundant bonds in the system
    float correlation;              // the correlation constant
    int giantsize_site;                  // The size of the giant rigid cluster
    int giantsize_bond;                  // The size of the giant rigid cluster
    int giantindex;                 // The index for the giant rigid cluster
    vector<int> thegraph[size];        // thegraph is the graph of all loaded edges
    vector<int> rgraph[size];        // rgraph is the graph of redundant bonds that don't take up any edges
    vector<Bond> edges;             //bonds only contains loaded edges
    vector<int> giantrigidcluster[size];    //giantrigidcluster is the graph for the giant rigid cluster
    stack<int> placesbeen;            // The list of places been while looking for a pebble
    ofstream myfile;                                                                                    // The file stream
private:
    const int EMPTY = -size - 1;
public:
    //SiteRP (int LL) : ll(LL){}

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // ADDING, REMOVING, AND CHECKING FOR EDGES

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // addedge adds an edge pointing from vertex i to vertex j in thegraph
    void addedge(int i, int j) {
        thegraph[i].push_back(j);
    }

    // addredundant adds an edge between i and j in rgraph, the separate graph of redundant edges
    void addredundant(int i, int j) {
        rgraph[i].push_back(j);
    }

    // badremoveedge removes an edge pointing from i to j in the graph.
// I call it "bad" because it looks through all of i's elements instead of
// using which one might be open in a tree search.
// This should return an error if there is no edge pointing from i to j.
    int badremoveedge(int i, int j) {
        for (int k = 0; k < thegraph[i].size(); k++) {
            if (thegraph[i].at(k) == j) {
                thegraph[i].erase(thegraph[i].begin() + k);
                return 0;
            }
        }
        cout << "I tried to remove an edge pointing from " << i << " to " << j << " but I couldn't find one.\n";
        return 0;
    }

// contains returns 1 if there is a non-redundant edge pointing from i to j (but doesn't check j to i) and 0 otherwise
    bool contains(int i, int j) {
        for (int k = 0; k < thegraph[i].size(); k++) {
            if (thegraph[i].at(k) == j) {
                return 1;
            }
        }
        return 0;
    }


// rcontains returns 1 if there is a redundant edge pointing from i to j (but doesn't check j to i) and 0 otherwise
    bool rcontains(int i, int j) {
        for (int k = 0; k < rgraph[i].size(); k++) {
            if (rgraph[i].at(k) == j) {
                return 1;
            }
        }


        return 0;
    }


// isempty returns 0 if there is any kind of redundant or nonredundant brace pointing in either direction between i and j, and 1 otherwise
// Modification required: isempty only looks at one type of bond, not the six (or at least, three) kinds of the triangular lattice
    bool isempty(int i, int j) {
        if (contains(i, j) || contains(j, i) || rcontains(i, j) || rcontains(j, i)) {
            return 0;
        }
        else {
            return 1;
        }
    }

    bool isempty(int i) {
        return isempty(i, i + ll + 1);
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// FINDING NEIGHBORS ON THE TRIANGULAR LATTICE WITH PBC

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    int moveright(int site) {
        if (site % ll == ll - 1) {
            return site - ll + 1;
        }
        else {
            return site + 1;
        }
    }

    int moveleft(int site) {
        if (site % ll == 0) {
            return site + ll - 1;
        }
        else {
            return site - 1;
        }
    }

    int moveup(int site) {
        if (site >= ll * ll - ll) {
            return site + ll - ll * ll;
        }
        else {
            return site + ll;
        }
    }

    int movedown(int site) {
        if (site <= ll - 1) {
            return site - ll + ll * ll;
        }
        else {
            return site - ll;
        }
    }


// moving to the right
    int dir1(int site) {
        return moveright(site);
    }


// moving up (geometrically, up and to the right)
    int dir2(int site) {
        return moveup(site);
    }


// moving up and to the left (the moves should commute)
    int dir3(int site) {
        return moveleft(moveup(site));
    }

    int dir4(int site) {
        return moveleft(site);
    }

    int dir5(int site) {
        return movedown(site);
    }

    int dir6(int site) {
        return moveright(movedown(site));
    }

    int choosedir(int site, int d) {
        switch (d) {
            case 1:
                return dir1(site);
            case 2:
                return dir2(site);
            case 3:
                return dir3(site);
            case 4:
                return dir4(site);
            case 5:
                return dir5(site);
            case 6:
                return dir6(site);
        }
        return false;
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// RANDOM FUNCTIONS

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// randprob() just returns a random number between zero and one.
    float randprob() {
        float f = rand();
        return f / RAND_MAX;
    }

// randsite0 finds a random plaquette that may or may not be occupied
    int randsite0() {
        return (rand() % size);
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// ADDING BONDS

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// This takes places been and reverses all the edges on the path.
// So if placesbeen is from 2 to 12 to 5
// It should remove edges from 2 to 12 and from 12 to 5
// and add edges from 5 to 12 and 12 to 2
    void reversepath() {
        int starter = placesbeen.top();                // We start at the last place in the path, the site where we found a pebble
        pc[starter]--;                                // We remove a pebble from this site

        placesbeen.pop();                            // We remove this site from our path, but it is still stored in starter

        int ender;

        while (placesbeen.size() > 0)                    // While there are still sites in the path...
        {
            ender = placesbeen.top();                // ender becomes the most recent site
            addedge(starter, ender);                // and we point an edge from starter to ender
            badremoveedge(ender, starter);            // and remove an edge from ender to starter
            starter = ender;                        // then ender becomes the new starter
            placesbeen.pop();                        // and is removed from the path
        }
        pc[ender]++; // Finally, we add a pebble to the first site in the path, we've moved a pebble from the end to the start, reversing edges along the way
    }


// findpebble executes a depth first search for pebbles in the lattice, starting at i
// and searching all the other vertices
// placesbeen must be empty before findpebble is called
// if findpebble finds a pebble, it should set placesbeen to a path from i to the site with the pebble
// Otherwise, it should leave placesbeen blank.
// It returns 1 if a pebble was found and 0 if it wasn't or if the path wasn't empty to start
    bool findpebble(int i) {
        if (placesbeen.size() >
            0)                        // We shouldn't ever call findpebble unless our path has been cleared, or something didn't close like it should
        {
            cout << "I tried to find a pebble, but when I started the path stack wasn't empty, it had size " <<
            placesbeen.size() << endl;
            return 0;
        }
        else {
            bool beenthere[size] = {};            // We create an array to say whether we've been to these sites before (or should we create a global array and just set it to zero here?)
            beenthere[i] = 1;                    // which we haven't, to start, except for our starting site (and any skip sites)
            placesbeen.push(i);                    // start our path at the starting site
            // cout << "Current location: " << placesbeen.top() << endl;

            int cl;                                // cl is the current location
            int prosp;                            // prosp is a vertex connected to cl that we might prospectively move to
            //cout << "Current location is " << placesbeen.top() << endl;

            while (placesbeen.size() >
                   0)                                        // Until we are forced to retreat all the way back to the first vertex...
            {
                cl = placesbeen.top();                                            // Our current location is the last place in the path
                // cout << "From the top, our current location is " << cl << endl;
                for (int index1 = 0;
                     index1 < thegraph[cl].size(); index1++) // for each place we can go from our current location...
                {
                    prosp = thegraph[cl].at(
                            index1);                            // Our prospective location is one of the places we can go to from cl
                    //cout << "The prospective vertex we consider is " << prosp << endl;
                    if (beenthere[prosp] == 0)                                    // if we haven't been there before...
                    {
                        // cl = prosp; // move our current location to there
                        //cout << "We are moving to " << prosp << endl;

                        placesbeen.push(prosp);                                    // Then we add it to the path
                        //cout << "Current location after pushing: " << placesbeen.top() << endl;


                        if (pc[placesbeen.top()] >
                            0) { return 1; }                // If our new site has a pebble, quit looking for pebbles and say we found one

                        beenthere[prosp] = 1;                                    // Otherwise mark it as having been visited, but keep looking for a pebble
                        //cout << "Current location is, in the for loop " << placesbeen.top() << endl;

                        break;                                                    // And break the for loop-- no point in continuing to explore cl's neighbors
                        // Once we've found one to move onto
                    }
                }


                if (cl == placesbeen.top() || placesbeen.size() ==
                                              0) // If, after the for loop, the new top of the path is the same as the old one,
                {                            // Then we didn't move anywhere, so we need to pop off the last vertex and retreat
                    placesbeen.pop();
                    /*
                    if (placesbeen.size() == 0){ cout << "The stack is empty." << endl; }
                    else
                    {
                    cout << "Current location after popping: " << placesbeen.top() << endl;
                    }
                    */
                }
            }
            return 0;
        }
    }


// When called with two arguments, findpebble skips over the second site to avoid infinitely swapping pebble
// between the two sites the brace connects, by marking skip as a place that we've already been
    bool findpebble(int i, int skip) {
        if (placesbeen.size() > 0) {
            cout << "I tried to find a pebble, but when I started the path stack wasn't empty, it had size " <<
            placesbeen.size() << endl;
            return 0;
        }
        else {
            bool beenthere[size] = {};            // We create an array to say whether we've been to these sites before
            beenthere[i] = 1;                    // which we haven't, to start, except for our starting site (and any skip sites)
            beenthere[skip] = 1;
            placesbeen.push(i);                    // start our path at the starting site
            // cout << "Current location: " << placesbeen.top() << endl;


            int cl;                                // cl is the current location
            int prosp;                            // prosp is a vertex connected to cl that we might prospectively move to
            //cout << "Current location is " << placesbeen.top() << endl;

            while (placesbeen.size() > 0) // Until we are forced to retreat all the way back to the first vertex...
            {
                cl = placesbeen.top();
                // cout << "From the top, our current location is " << cl << endl;
                for (int index1 = 0;
                     index1 < thegraph[cl].size(); index1++) // for each place we can go from our current location...
                {
                    prosp = thegraph[cl].at(index1);
                    //cout << "The prospective vertex we consider is " << prosp << endl;
                    if (beenthere[prosp] == 0) // if we haven't been there before...
                    {
                        // cl = prosp; // move our current location to there
                        //cout << "We are moving to " << prosp << endl;

                        placesbeen.push(prosp);      // and add it to the path
                        //cout << "Current location after pushing: " << placesbeen.top() << endl;
                        if (pc[placesbeen.top()] >
                            0) { return 1; }        // If our new site has a pebble, quit looking for pebbles
                        beenthere[prosp] = 1;
                        //cout << "Current location is, in the for loop " << placesbeen.top() << endl;
                        break;
                        /*
                        if (pc[cl]>0)			// Then, if there is a pebble at our new location, stop looking for pebbles
                        {
                        return 0;
                        }
                        else
                        {
                        break;					  // Otherwise, just quit moving through the old location's neighbors
                        }
                        */
                    }
                }
                if (cl == placesbeen.top() || placesbeen.size() ==
                                              0) // If, after the for loop, the new top of the path is the same as the old one,
                {                            // Then we didn't move anywhere, so we need to pop off the last vertex and retreat
                    placesbeen.pop();
                    /*
                    if (placesbeen.size() == 0){ cout << "The stack is empty." << endl; }
                    else
                    {
                    cout << "Current location after popping: " << placesbeen.top() << endl;
                    }
                    */
                }
            }

            return 0;
        }
    }


// loadsite looks for pebbles and moves them onto i until i has two pebbles or it stops finding pebbles
    bool loadsite(int i) {
        if (placesbeen.size() == 0) {
            while (pc[i] < 2 && findpebble(i)) // while the site is not loaded and you are finding pebbles...
                // c++ documentation says && short circuits, so you shouldn't
                // even look for a pebble if the site is loaded.
            {
                reversepath();                    // reverse the path, which also shifts the pebbles
            }
        }
        else {
            cout << "I tried to find a pebble, but when I started the path stack wasn't empty, it had size " <<
            placesbeen.size() << endl;
        }

        if (pc[i] == 2) { return 1; }        // return 1 if the site loaded successfully.
        else { return 0; }
    }

// loadsite with two arguments does the same thing, but won't try to take pebbles from skip to move onto i
// Which is to keep the two sites from swapping the three pebbles between themselves endlessly
    bool loadsite(int i, int skip) {
        if (placesbeen.size() == 0) {
            while (pc[i] < 2 && findpebble(i, skip)) // while the site is not loaded and you are finding pebbles...
                // C++ documentation says && short-circuits, so this shouldn't even
                // look for a pebble if the site is loaded.
            {
                reversepath();                    // reverse the path, which also shifts the pebbles
            }
        }
        else {
            cout << "I tried to load site " << i << " but the placesbeen stack wasn't empty.\n";
        }
        if (pc[i] == 2) { return 1; }        // return 1 if the site loaded successfully.
        else { return 0; }
    }

// loadsites tries to move pebbles until there are two on both sites i and j
    bool loadsites(int i, int j) {

        while (pc[j] < 2 && findpebble(j)) {
            reversepath();
        }

        while (pc[i] < 2 &&
               findpebble(i, j))    // while the first site is not loaded and you are finding pebbles skipping
        {                                        // second site, load the first site, then load the second.
            reversepath();
            while (pc[j] < 2 && findpebble(j)) {
                reversepath();
            }
        }
        if (pc[i] == 2 && pc[j] == 2) { return 1; }
        else { return 0; }
    }

// addbond tries to load the sites. If it succeeds, it adds an edge from i to j and takes a pebble from i. Otherwise, it adds a redundant edge
    void addbond(int i, int j) {
        if (numbonds < 2 * size - 3 &&
            loadsites(i, j))            // If there are at least four pebbles left, we try to load the sites
        {
            addedge(i, j);                // if we succeed, we add the edge from i to j and remove a pebble from i
            numbonds++;
            edges.push_back(Bond(i, j));
            pc[i]--;
        }
        else {
            addredundant(i, j);            // otherwise, we leave the pebbles where we shuffled them and
            rbonds++;
        }                                // place only a redundant bond
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// REPORTING TO SCREEN

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    void setfilestream(float cval, int tval) {


        char mychar[] = "data/cxxxtxxxx.txt";
        int start = sizeof(mychar) - 13;
        cval = cval * 100;
        int intc = cval;
        mychar[start + 4] = '0' + (tval % 10000) / 1000;
        mychar[start + 5] = '0' + (tval % 1000) / 100;
        mychar[start + 6] = '0' + (tval % 100) / 10;
        mychar[start + 7] = '0' + (tval % 10);

        mychar[start] = '0' + (intc % 1000) / 100;
        mychar[start + 1] = '0' + (intc % 100) / 10;
        mychar[start + 2] = '0' + (intc % 10);

        myfile.close();
        myfile.open(mychar);

        //cout << mychar;
    }

    void log(int span = -1) {
        myfile << ll << "\t" << correlation << "\t" << numparts << "\t" << numbonds << "\t" << rbonds << "\t" <<
        giantsize_bond << "\t" << giantsize_site << "\t" << span << "\n";
    }

// listedges lists the edges from site i
    void listedges(int i) {
        cout << "\nvertex number " << i << " points towards the following vertices:\n";

        for (int index = 0; index < thegraph[i].size(); index++) {
            cout << thegraph[i].at(index) << " ";
        }

        cout << endl;
    }

// listalledges lists all the edges from the sites
    void listalledges() {
        for (int index = 0; index < size; index++) {
            listedges(index);
        }
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// FIXED-PROBABILITY TRIALS AND ATTENDANT FUNCTIONS

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// BOND DILUTED NETWORK TO CONSIDER RIGIDITY PERCOLATION

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// functions in original codes by Zeb


    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// SITE DILUTED NETWORK TO CONSIDER RIGIDITY PERCOLATION

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// initgiantrigidcluster() initializes the giantrigidcluster graph
    void initgiantrigidcluster() {
        giantsize_site = 0;
        giantsize_bond = 0;

        // rewrite the RigidIndex for edges
        for (vector<Bond>::iterator it = edges.begin(); it != edges.end(); ++it) {
            it->initBondRigidIndex();
        }

        for (int bondindex = 0; bondindex < size; bondindex++)  // Clear the graphs of giantrigidcluster
        {
            giantrigidcluster[bondindex].clear();
        }
    }

// initemptytrigraph() updates numbonds, rbonds, thegraph, rgraph, placesbeen to a triangular graph with no particles or bonds
    int initemptytrigraph() {
        numbonds = 0; // the number of bonds
        edges.clear();
        numparts = 0; // the number of particles
        rbonds = 0;   // the number of redundant bonds

        for (int pcindex = 0; pcindex < size; pcindex++) // Just setting the pebble count to 2 everywhere.
        {                                               // and setting which sites are occupied
            pc[pcindex] = 2;
            occ[pcindex] = 0;
        }

        while (placesbeen.size() > 0)                    // Clear the places been stack
        {
            placesbeen.pop();
        }

        for (int bondindex = 0; bondindex < size; bondindex++)  // Clear the graphs of redundant and nonredundant bonds
        {
            rgraph[bondindex].clear();
            thegraph[bondindex].clear();
            giantrigidcluster[bondindex].clear();

        }

        initgiantrigidcluster();

        return numparts;
    }

    void addtricluster2(int site,
                        float c) // Has already added the rigidcluster function, as well as the spanning cluster
    {

        if (occ[site] == 0) {
            int newsite;
            int numneighbors = occ[choosedir(site, 1)] + occ[choosedir(site, 2)] + occ[choosedir(site, 3)] +
                               occ[choosedir(site, 4)] + occ[choosedir(site, 5)] + occ[choosedir(site, 6)];

            if (randprob() < pow(1. - c, 6 - numneighbors)) {
                occ[site] = 1;

                numparts++;
                for (int k = 1; k <= 6; k++) {
                    newsite = choosedir(site, k);
                    if (occ[newsite] == 1 && isempty(site, newsite)) {
                        addbond(site, newsite);
                    }
                }
                // choose some densities for the rigid cluster
                if (numparts % (ll / 2) == 0) {
                    rigidcluster();

                    int span = spanningrcluster();
                    log(span);
                }
                else
                    log();

            }
        }
    }

    void onetritrial2(int max_out_time, float c) {
        int numattempts = 0;
        int maxout = max_out_time * ll * ll;
        initemptytrigraph();

        while (numattempts < maxout && numparts < ll * ll) {
            numattempts++;
            addtricluster2(randsite0(), c);
            //plot();
        }

        myfile.close();
    }

    void multictrial(int maxout, float c1, float c2, float dc, int numtrials) {
        for (correlation = c1; correlation <= c2; correlation += dc) {
            for (int mtc = 14; mtc <= numtrials; mtc++) {
                int numattempts = 0;
                initemptytrigraph();

                setfilestream(correlation, mtc);

                while (numattempts < maxout && numparts < ll * ll) {
                    numattempts++;

                    addtricluster2(randsite0(), correlation);
                }
            }
        }
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// rigid cluster

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    bool isredundant(int i, int j)  //see if the test bond between (i,j) is redundant(dependent)
    {
        if (numbonds < 2 * size - 3 &&
            loadsites(i, j))            // if there are at least four pebbles left, we try to load the sites
        {
            return 0;                // if we succeed, then the edge from i to j is independent
        }
        else {
            return 1;            // otherwise, the edge is redundant(dependent)
        }
    }

    bool isbondrigid(Bond &a, Bond &b) //see if the two bonds a and b are rigid to each other
    {
        if (isredundant(a.vertices.first, b.vertices.first) && isredundant(a.vertices.first, b.vertices.second) &&
            isredundant(a.vertices.second, b.vertices.first) && isredundant(a.vertices.second, b.vertices.second)) {
            return 1; // if all the four test bonds between a1,b1; a1,b2; a2,b1; a2,b2 are redundant
        }
        else {
            return 0;
        }
    }

    void rigidcluster() // mark the rigid clusters and put into rcluster, and return the size of the giant rigid cluster.
    {
        initgiantrigidcluster(); //empty the vector array first

        int rcnum = 0; //index of the rigid cluster
        giantsize_bond = 1; //if there exists any bond, the smallest possible giant rigid cluster bond size is 1
        giantindex = 0; //the index for the giant rigid cluster (in this function)

        for (vector<Bond>::iterator refBond = edges.begin(); refBond != edges.end(); ++refBond) {
            if (refBond->RigidIndex == 0) {
                int rclustersize_bond = 1;//the bond size of this rigid cluster
                rcnum++;
                refBond->RigidIndex = rcnum;
                for (vector<Bond>::iterator testBond = edges.begin(); testBond != edges.end(); ++testBond) {
                    if (testBond->RigidIndex == 0) {
                        //the test and ref bonds are not in some rigid clusters
                        if (isbondrigid(*refBond, *testBond) == true) { //if test bond is rigid with respect to refbond
                            testBond->RigidIndex = rcnum;
                            rclustersize_bond++;
                        }
                    }
                }
                if (rclustersize_bond >= giantsize_bond) { //update the size of the giant rigid cluster
                    giantsize_bond = rclustersize_bond;
                    giantindex = rcnum;
                }
            }
        }

        // pick out the giant rigid cluster and store it in the vector "giantrigidcluster" (!!! we need it to become a undirected adjacent list)

        for (vector<Bond>::iterator it = edges.begin(); it != edges.end(); ++it) {
            if (it->RigidIndex == giantindex) {
                int site_I = it->vertices.first;
                int site_J = it->vertices.second; //the two sites of the rigid bond
                giantrigidcluster[site_I].push_back(site_J);
                giantrigidcluster[site_J].push_back(site_I);
            }
        }

        for (int i = 0; i <= size - 1; i++) {
            if (!giantrigidcluster[i].empty()) {
                giantsize_site++;
            }
        }
    }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// SPANNING RIGID CLUSTER

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// We need to pick out the giant rigid cluster from the network and then determine if it is the spanning cluster.

    bool spanningrcluster() {
        // Do the DFS in the giant rigid cluster
        int beenthere[size] = {};
        std::fill_n(beenthere, size, EMPTY);
        stack<int> DFS_rcluster;
        // Find a starting point
        int v_start;

        for (v_start = 0; v_start <= size - 1; v_start++) {
            if (!giantrigidcluster[v_start].empty()) {
                break;
            }
        }

        //v_start is the starting point for the DFS
        beenthere[v_start] = 0;
        DFS_rcluster.push(v_start);
        int cl; //cl is the current location
        int prosp; //prosp is the prospective location we're going to move to
        while (DFS_rcluster.size() > 0) // make sure the DFS is in the giant rigid cluster
        {
            cl = DFS_rcluster.top();
            for (int index1 = 0; index1 < giantrigidcluster[cl].size(); index1++) {
                prosp = giantrigidcluster[cl].at(index1);
                if (beenthere[prosp] == EMPTY) {
                    DFS_rcluster.push(prosp);
                    //Do sth here
                    //beenthere[prosp]=1;


                    //Get the x-displacement for prosp from the x-displacement for cl (x_prosp=x_cl +1 OR +0 OR -1)
                    if (prosp == choosedir(cl, 2) || prosp == choosedir(cl, 5)) {
                        beenthere[prosp] = beenthere[cl];
                    } else if (prosp == choosedir(cl, 1) || prosp == choosedir(cl, 6)) {
                        beenthere[prosp] = beenthere[cl] + 1;
                    } else if (prosp == choosedir(cl, 3) || prosp == choosedir(cl, 4)) {
                        beenthere[prosp] = beenthere[cl] - 1;
                    }
                    else {
                        cout << "Whoooops, sth goes wrong." << endl;
                    }

                    break;
                }
                else //Existing marked displacement can be compared to current place to see if the spanning cluster exists.
                {
                    if (abs(beenthere[cl] - beenthere[prosp]) == ll - 1) {
                        return true;
                    }
                }
            }
            if (cl == DFS_rcluster.top() || DFS_rcluster.size() == 0) {
                DFS_rcluster.pop();
            }
        }
        return false;// After the DFS still no finding about the spanning cluster, then it is not spanning.
    }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// PLOT THE NETWORK

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    void onetritrial2_plot(float p, float c) {
        initemptytrigraph();
        //while(numparts <= p*ll*ll)
        while (numparts <= p * ll * ll) {
            addtricluster2_withoutRIGID(randsite0(), c);
        }
        myfile.close();
    }

    void addtricluster2_withoutRIGID(int site, float c) // Has not added the rigidcluster function, as well as the spanning cluster
    {

        if (occ[site] == 0) {
            int newsite;
            int numneighbors = occ[choosedir(site, 1)] + occ[choosedir(site, 2)] + occ[choosedir(site, 3)] +
                               occ[choosedir(site, 4)] + occ[choosedir(site, 5)] + occ[choosedir(site, 6)];

            if (randprob() < pow(1. - c, 6 - numneighbors)) {
                occ[site] = 1;

                numparts++;
                for (int k = 1; k <= 6; k++) {
                    newsite = choosedir(site, k);
                    if (occ[newsite] == 1 && isempty(site, newsite)) {
                        addbond(site, newsite);
                    }
                }
                //log();
            }
        }
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// UNI-TEST & APPLICATIONS

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    void OneTrialTest(float cfor, int trial) //Generate one-time trial for triangular lattice (site RP)
    {
        cfor = 0.0; // correlation constant
        trial = 1; // trial counting

        setfilestream(cfor,trial);
        onetritrial2(300000,cfor);
    }

    void MultiTrialTest() //Generate multiple-times trial for triangular lattice (site RP)
    {
        multictrial(ll*ll*100000,0.0,0.1,.2,20); //That's it!
    }

    void PlotNetworkTest() //Generate network plot file
    {
        float cfor;
        float p;
        cout<<"Type in the number for constant c: \n";
        cin>>cfor;
        cout<<"Type in network density p: \n";
        cin>>p;


        onetritrial2_plot(p,cfor);
        rigidcluster();

        int span = spanningrcluster();

        cout<<"The size of the giant rigid cluster is "<< giantsize_site <<" with the spanning status in "<<span<<"\n";

        ofstream outfile;
        outfile.open("rclusterout.txt");
        for(int i=0; i < size; i++)
        {
            if(occ[i] != 0)
                outfile << i%ll << "\t" << i/ll << "\t" << !giantrigidcluster[i].empty() << "\n";
        }
        outfile.close();
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// END OF THE CLASS

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// MAIN FUNCTION

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main()
{
    clock_t t1, t2;													// Creates clock variables to track the runtime
    t1 = clock();
    srand(time(NULL));

    SiteRP TriLattice;
    TriLattice.PlotNetworkTest();

    t2 = clock();
    float clocktime((float)t2 - (float)t1);
    cout << "\n The total run time was " << clocktime / CLOCKS_PER_SEC << endl;
    cin.ignore();
    return 0;
}

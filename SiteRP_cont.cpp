//
// Created by Shang Zhang on 5/7/17.
//

//
// Created by Zeb & Shang.
//

#include "SiteRP_cont.h"
#include <boost/tokenizer.hpp>
#include <boost/foreach.hpp>
#include <sstream>

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// ADDING, REMOVING, AND CHECKING FOR EDGES

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// addedge adds an edge pointing from vertex i to vertex j in thegraph
void SiteRP_cont::addedge(int i, int j) {
    thegraph[i].push_back(j);
}

// addredundant adds an edge between i and j in rgraph, the separate graph of redundant edges
void SiteRP_cont::addredundant(int i, int j) {
    rgraph[i].push_back(j);
}

// badremoveedge removes an edge pointing from i to j in the graph.
// I call it "bad" because it looks through all of i's elements instead of
// using which one might be open in a tree search.
// This should return an error if there is no edge pointing from i to j.
int SiteRP_cont::badremoveedge(int i, int j) {
    for (int k = 0; k < thegraph[i].size(); k++) {
        if (thegraph[i].at(k) == j) {
            thegraph[i].erase(thegraph[i].begin() + k);
            return 0;
        }
    }
    std::cout << "I tried to remove an edge pointing from " << i << " to " << j << " but I couldn't find one.\n";
    return 0;
}

// contains returns 1 if there is a non-redundant edge pointing from i to j (but doesn't check j to i) and 0 otherwise
bool SiteRP_cont::contains(int i, int j) {
    for (int k = 0; k < thegraph[i].size(); k++) {
        if (thegraph[i].at(k) == j) {
            return 1;
        }
    }
    return 0;
}


// rcontains returns 1 if there is a redundant edge pointing from i to j (but doesn't check j to i) and 0 otherwise
bool SiteRP_cont::rcontains(int i, int j) {
    for (int k = 0; k < rgraph[i].size(); k++) {
        if (rgraph[i].at(k) == j) {
            return 1;
        }
    }


    return 0;
}


// isempty returns 0 if there is any kind of redundant or nonredundant brace pointing in either direction between i and j, and 1 otherwise
// Modification required: isempty only looks at one type of bond, not the six (or at least, three) kinds of the triangular lattice
bool SiteRP_cont::isempty(int i, int j) {
    if (contains(i, j) || contains(j, i) || rcontains(i, j) || rcontains(j, i)) {
        return 0;
    }
    else {
        return 1;
    }
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
void SiteRP_cont::reversepath() {
    int starter = placesbeen.top();                // We start at the last place in the path, the site where we found a pebble
    pc[starter]--;                                // We remove a pebble from this site

    placesbeen.pop();                            // We remove this site from our path, but it is still stored in starter

    int ender = 0;

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
bool SiteRP_cont::findpebble(int i) {
    if (placesbeen.size() > 0)                        // We shouldn't ever call findpebble unless our path has been cleared, or something didn't close like it should
    {
        std::cout << "I tried to find a pebble, but when I started the path stack wasn't empty, it had size " << placesbeen.size() << std::endl;
        return 0;
    }
    else {
        bool *beenthere = new bool[size]; // We create an array to say whether we've been to these sites before (or should we create a global array and just set it to zero here?)
        for (int ii = 0; ii < size; ii++)
            beenthere[ii] = false;
        
        beenthere[i] = 1;                    // which we haven't, to start, except for our starting site (and any skip sites)
        placesbeen.push(i);                    // start our path at the starting site
        // cout << "Current location: " << placesbeen.top() << endl;

        int cl;                                // cl is the current location
        int prosp;                            // prosp is a vertex connected to cl that we might prospectively move to
        //cout << "Current location is " << placesbeen.top() << endl;

        while (placesbeen.size() > 0)        // Until we are forced to retreat all the way back to the first vertex...
        {
            cl = placesbeen.top();         // Our current location is the last place in the path
            // cout << "From the top, our current location is " << cl << endl;
            for (int index1 = 0; index1 < thegraph[cl].size(); index1++) // for each place we can go from our current location...
            {
                prosp = thegraph[cl].at(index1);         // Our prospective location is one of the places we can go to from cl
                //cout << "The prospective vertex we consider is " << prosp << endl;
                if (beenthere[prosp] == 0)          // if we haven't been there before...
                {
                    // cl = prosp; // move our current location to there
                    //cout << "We are moving to " << prosp << endl;

                    placesbeen.push(prosp);         // Then we add it to the path
                    //cout << "Current location after pushing: " << placesbeen.top() << endl;


                    if (pc[placesbeen.top()] > 0) {
                        delete[] beenthere;
                        return 1; }                // If our new site has a pebble, quit looking for pebbles and say we found one

                    beenthere[prosp] = 1;            // Otherwise mark it as having been visited, but keep looking for a pebble
                    //cout << "Current location is, in the for loop " << placesbeen.top() << endl;

                    break;                                                    // And break the for loop-- no point in continuing to explore cl's neighbors
                    // Once we've found one to move onto
                }
            }


            if (cl == placesbeen.top() || placesbeen.size() == 0) // If, after the for loop, the new top of the path is the same as the old one,
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
        delete[] beenthere;
        return 0;
    }
}


// When called with two arguments, findpebble skips over the second site to avoid infinitely swapping pebble
// between the two sites the brace connects, by marking skip as a place that we've already been
bool SiteRP_cont::findpebble(int i, int skip) {
    if (placesbeen.size() > 0) {
        std::cout << "I tried to find a pebble, but when I started the path stack wasn't empty, it had size " << placesbeen.size() << std::endl;
        return 0;
    }
    else {
        bool *beenthere = new bool[size]; // We create an array to say whether we've been to these sites before (or should we create a global array and just set it to zero here?)
        for (int ii = 0; ii < size; ii++)
            beenthere[ii] = false;
        
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
                        0) { delete[] beenthere; return 1; }        // If our new site has a pebble, quit looking for pebbles
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
            if (cl == placesbeen.top() || placesbeen.size() == 0) // If, after the for loop, the new top of the path is the same as the old one,
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
        
        delete[] beenthere;
        return 0;
    }
}


// loadsite looks for pebbles and moves them onto i until i has two pebbles or it stops finding pebbles
bool SiteRP_cont::loadsite(int i) {
    if (placesbeen.size() == 0) {
        while (pc[i] < 2 && findpebble(i)) // while the site is not loaded and you are finding pebbles...
            // c++ documentation says && short circuits, so you shouldn't
            // even look for a pebble if the site is loaded.
        {
            reversepath();                    // reverse the path, which also shifts the pebbles
        }
    }
    else {
        std::cout << "I tried to find a pebble, but when I started the path stack wasn't empty, it had size " <<
                  placesbeen.size() << std::endl;
    }

    if (pc[i] == 2) { return 1; }        // return 1 if the site loaded successfully.
    else { return 0; }
}

// loadsite with two arguments does the same thing, but won't try to take pebbles from skip to move onto i
// Which is to keep the two sites from swapping the three pebbles between themselves endlessly
bool SiteRP_cont::loadsite(int i, int skip) {
    if (placesbeen.size() == 0) {
        while (pc[i] < 2 && findpebble(i, skip)) // while the site is not loaded and you are finding pebbles...
            // C++ documentation says && short-circuits, so this shouldn't even
            // look for a pebble if the site is loaded.
        {
            reversepath();                    // reverse the path, which also shifts the pebbles
        }
    }
    else {
        std::cout << "I tried to load site " << i << " but the placesbeen stack wasn't empty.\n";
    }
    if (pc[i] == 2) { return 1; }        // return 1 if the site loaded successfully.
    else { return 0; }
}

// loadsites tries to move pebbles until there are two on both sites i and j
bool SiteRP_cont::loadsites(int i, int j) {

    while (pc[j] < 2 && findpebble(j)) {
        reversepath();
    }

    while (pc[i] < 2 && findpebble(i, j))    // while the first site is not loaded and you are finding pebbles skipping
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
void SiteRP_cont::addbond(int i, int j) {
    if (numbonds < 2 * size - 3 && loadsites(i, j))            // If there are at least four pebbles left, we try to load the sites
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

void SiteRP_cont::setfilestream() {

    std::string mychar = "./data/";
    std::string mychar_rcluster = "./data/rcluster_";
    std::string FileEnding = ".txt";

    mychar.append(FileName);
    mychar.append(".txt");

    mychar_rcluster.append(FileName);
    mychar_rcluster.append(".txt");

    myfile.close();
    myfile.open(mychar);

    rclusterfile.close();
    rclusterfile.open(mychar_rcluster);

    //cout << mychar;
}

void SiteRP_cont::log(int span) {
    myfile <<  numparts << "\t" << numbonds << "\t" << rbonds << "\t" << giantsize_bond << "\t" << giantsize_site << "\t" << span << "\n";
}

// listedges lists the edges from site i
void SiteRP_cont::listedges(int i) {
    std::cout << "\nvertex number " << i << " points towards the following vertices:\n";

    for (int index = 0; index < thegraph[i].size(); index++) {
        std::cout << thegraph[i].at(index) << " ";
    }

    std::cout << std::endl;
}

// listalledges lists all the edges from the sites
void SiteRP_cont::listalledges() {
    for (int index = 0; index < size; index++) {
        listedges(index);
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// SITE DILUTED NETWORK TO CONSIDER RIGIDITY PERCOLATION

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// initgiantrigidcluster() initializes the giantrigidcluster graph
void SiteRP_cont::initgiantrigidcluster() {
    giantsize_site = 0;
    giantsize_bond = 0;

    // rewrite the RigidIndex for edges
    for (std::vector<Bond>::iterator it = edges.begin(); it != edges.end(); ++it) {
        it->initBondRigidIndex();
    }

    for (int bondindex = 0; bondindex < size; bondindex++)  // Clear the graphs of giantrigidcluster
    {
        giantrigidcluster[bondindex].clear();
    }
}

// initemptytrigraph() updates numbonds, rbonds, thegraph, rgraph, placesbeen to a triangular graph with no particles or bonds
int SiteRP_cont::initemptytrigraph() {

    numbonds = 0; // the number of bonds
    edges.clear();
    numparts = 0; // the number of particles
    rbonds = 0;   // the number of redundant bonds

    for (int pcindex = 0; pcindex < size; pcindex++) // Just setting the pebble count to 2 everywhere.
    {                                               // and setting which sites are occupied
        pc[pcindex] = 2;
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
        undirectedgraph[bondindex].clear();
    }

    initgiantrigidcluster();

    //XSpanLastStatus = 0;
    //YSpanLastStatus = 0;

    return numparts;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// rigid cluster

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool SiteRP_cont::isredundant(int i, int j)  //see if the test bond between (i,j) is redundant(dependent)
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

bool SiteRP_cont::isbondrigid(Bond &a, Bond &b) //see if the two bonds a and b are rigid to each other
{
    if (isredundant(a.vertices.first, b.vertices.first) && isredundant(a.vertices.first, b.vertices.second) &&
        isredundant(a.vertices.second, b.vertices.first) && isredundant(a.vertices.second, b.vertices.second)) {
        return 1; // if all the four test bonds between a1,b1; a1,b2; a2,b1; a2,b2 are redundant
    }
    else {
        return 0;
    }
}

void SiteRP_cont::rigidcluster() // mark the rigid clusters
{
    initgiantrigidcluster(); //empty the vector array first

    int rcnum = 0; //index of the rigid cluster
    giantsize_bond = 1; //if there exists any bond, the smallest possible giant rigid cluster bond size is 1
    giantindex = 0; //the index for the giant rigid cluster (in this function)

    for (std::vector<Bond>::iterator refBond = edges.begin(); refBond != edges.end(); ++refBond) {
        if (refBond->RigidIndex == 0) {
            int rclustersize_bond = 1;//the bond size of this rigid cluster
            rcnum++;
            refBond->RigidIndex = rcnum;
            for (std::vector<Bond>::iterator testBond = edges.begin(); testBond != edges.end(); ++testBond) {
                if (testBond->RigidIndex == 0) {
                    //the test and ref bonds are not in some rigid clusters
                    if (isbondrigid(*refBond, *testBond) == true) { //if test bond is rigid with respect to refbond
                        testBond->RigidIndex = rcnum;
                        rclustersize_bond++;
                    }
                }
            }
            //std::cout << "Hello: " << rclustersize_bond << std::endl;
            if (rclustersize_bond >= giantsize_bond) { //update the size of the giant rigid cluster
                giantsize_bond = rclustersize_bond;
                giantindex = rcnum;
            }
        }
        //std::cout << rcnum << std::endl;
    }

    // pick out the giant rigid cluster and store it in the vector "giantrigidcluster" (!!! we need it to become a undirected adjacent list)

    for (std::vector<Bond>::iterator it = edges.begin(); it != edges.end(); ++it) {
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
            //std::cout << i << std::endl;
        }
    }

}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// THE RIGID CLUSTER INFO (FOR CLUSTER STATISTICS INFO)

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// store the rigid cluster decomposition information in sites (initial rcluster_site before the running of this function)
void SiteRP_cont::StoreRigidInfoOfSite(){
    // clear the stored rigid cluster info, which is already not useful
    for (int i = 0; i <= size - 1; ++i) {
        rcluster_site[i].clear();
    }

    // store the rigid cluster info to rcluster_site
    for (std::vector<Bond>::iterator it = edges.begin(); it != edges.end() ; ++it) {
        int site_I = it->vertices.first;
        int site_J = it->vertices.second; //the two sites of the rigid bond

        // if rcluster_site has not stored the RigidIndex, store it
        if (find(rcluster_site[site_I].begin(),rcluster_site[site_I].end(),it->RigidIndex) == rcluster_site[site_I].end()){
            rcluster_site[site_I].push_back(it->RigidIndex);
        }
        if (find(rcluster_site[site_J].begin(),rcluster_site[site_J].end(),it->RigidIndex) == rcluster_site[site_J].end()){
            rcluster_site[site_J].push_back(it->RigidIndex);
        }
        //push the rigid indices to the rcluster_site in order to store the info
    }

    // print out the rigid cluster decomposition info
    for (int st = 0; st <= size - 1; ++st) {
        if (rcluster_site[st].empty()){
            rclusterfile << "0" << "\n";
        }
        else{
            for (std::vector<int>::iterator it = rcluster_site[st].begin(); it != rcluster_site[st].end(); ++it) {
                rclusterfile << *it << "\t";
            }
            rclusterfile << "\n";
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// SPANNING RIGID CLUSTER

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// We need to pick out the giant rigid cluster from the network and then determine if it is the spanning cluster.

bool SiteRP_cont::spanningrcluster() {
    // Do the DFS in the giant rigid cluster
    double *beenthere = new double[size]; // We create an array to say whether we've been to these sites before (or should we create a global array and just set it to zero here?)
    for (int ii = 0; ii < size; ii++)
        beenthere[ii] = EMPTY;
    //std::fill_n(beenthere, size, EMPTY);
    std::stack<int> DFS_rcluster;
    // Find a starting point
    int v_start;

    for (v_start = 0; v_start < size - 1; v_start++) {
        if (!giantrigidcluster[v_start].empty()) {
            break;
        }
    }

    //std::cout << v_start << std::endl;

    //v_start is the starting point for the DFS
    beenthere[v_start] = 0.0;
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

                //Get the x-displacement for prosp from the x-displacement for cl
                beenthere[prosp] = beenthere[cl] + vertices[prosp].XCoordinate - vertices[cl].XCoordinate;

                break;
            }
            else //Existing marked displacement can be compared to current place to see if the spanning cluster exists.
            {
                //std::cout << "The distance for these two atoms in the same cluster: " << fabs(beenthere[cl] - beenthere[prosp]) << std::endl;
                if (fabs(fabs(beenthere[cl] - beenthere[prosp]) - BoxLength) < 1.3) {
                    std::cout << "The distance for these two atoms in the same cluster: " << fabs(beenthere[cl] - beenthere[prosp]) << std::endl;
                    delete[] beenthere;
                    return true;
                }
            }
        }
        if (cl == DFS_rcluster.top() || DFS_rcluster.size() == 0) {
            DFS_rcluster.pop();
        }
    }
    delete[] beenthere;
    return false;
}

// Key difference for SiteRP to SiteRP_cont

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// BUILDING THE NETWORK

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// Do the periodic boundary condition for particle b
float SiteRP_cont::distance(site &a, site &b) {
    float dis1 = sqrt( pow(a.XCoordinate - b.XCoordinate,2 ) + pow(a.YCoordinate - b.YCoordinate,2 ) );
    float dist = dis1;
    float dis2 = sqrt( pow(a.XCoordinate- b.XCoordinate - BoxLength,2 ) + pow(a.YCoordinate - b.YCoordinate,2 ) );
    if (dis2 <= dist) {
        dist = dis2;
    }
    float dis3 = sqrt( pow(a.XCoordinate - b.XCoordinate,2 ) + pow(a.YCoordinate - b.YCoordinate - BoxLength,2 ) );
    if (dis3 <= dist) {
        dist = dis3;
    }
    float dis4 = sqrt( pow(a.XCoordinate - b.XCoordinate + BoxLength,2 ) + pow(a.YCoordinate - b.YCoordinate,2 ) );
    if (dis4 <= dist) {
        dist = dis4;
    }
    float dis5 = sqrt( pow(a.XCoordinate - b.XCoordinate,2 ) + pow(a.YCoordinate - b.YCoordinate + BoxLength,2 ) );
    if (dis5 <= dist) {
        dist = dis5;
    }
    float dis6 = sqrt( pow(a.XCoordinate - b.XCoordinate - BoxLength,2 ) + pow(a.YCoordinate - b.YCoordinate - BoxLength,2 ) );
    if (dis6 <= dist) {
        dist = dis6;
    }
    float dis7 = sqrt( pow(a.XCoordinate - b.XCoordinate - BoxLength,2 ) + pow(a.YCoordinate - b.YCoordinate + BoxLength,2 ) );
    if (dis7 <= dist) {
        dist = dis7;
    }
    float dis8 = sqrt( pow(a.XCoordinate - b.XCoordinate + BoxLength,2 ) + pow(a.YCoordinate - b.YCoordinate - BoxLength,2 ) );
    if (dis8 <= dist) {
        dist = dis8;
    }
    float dis9 = sqrt( pow(a.XCoordinate - b.XCoordinate + BoxLength,2 ) + pow(a.YCoordinate - b.YCoordinate + BoxLength,2 ) );
    if (dis9 <= dist) {
        dist = dis9;
    }
    return dist;
};


void SiteRP_cont::BuildNetwork() // Has already added the rigidcluster function, as well as the spanning cluster
{
    initemptytrigraph();
    
    std::cout << "####################################################################################################" << std::endl;
    std::cout << " Building the Network... " << std::endl;
    std::cout << "####################################################################################################" << std::endl;

    int SiteNum = vertices.size();
    numparts = SiteNum;
    for (int i = 0; i < SiteNum; ++i) {
        for (int j = 0; j < SiteNum; ++j) {
            float AtomDistance = distance(vertices[i],vertices[j]);
            if (AtomDistance <= 1.3 && i != j){ // when the distance between two particles is less than 1.3d (d=1), the bond exists.
                undirectedgraph[i].push_back(j);
                if (isempty(i, j))
                {
                    addbond(i, j);
                }
            }
        }
    }
    
    std::cout << "####################################################################################################" << std::endl;
    std::cout << " Done! " << std::endl;
    std::cout << "####################################################################################################" << std::endl;


    setfilestream();
    
    std::cout << "####################################################################################################" << std::endl;
    std::cout << " Identify the Rigid Cluster... " << std::endl;
    std::cout << "####################################################################################################" << std::endl;

    rigidcluster();

    std::cout << "####################################################################################################" << std::endl;
    std::cout << " Done! " << std::endl;
    std::cout << "####################################################################################################" << std::endl;
    
    std::cout << "####################################################################################################" << std::endl;
    std::cout << " Identify the Spanning Cluster... " << std::endl;
    std::cout << "####################################################################################################" << std::endl;

    
    int span = spanningrcluster();
    log(span);
    
    std::cout << "####################################################################################################" << std::endl;
    std::cout << " Done! " << std::endl;
    std::cout << "####################################################################################################" << std::endl;

    std::cout << "####################################################################################################" << std::endl;
    std::cout << " Store the Rigid Infomation... " << std::endl;
    std::cout << "####################################################################################################" << std::endl;

    StoreRigidInfoOfSite();
    
    std::cout << "####################################################################################################" << std::endl;
    std::cout << " Done! " << std::endl;
    std::cout << "####################################################################################################" << std::endl;

    
    std::cout << "####################################################################################################" << std::endl;
    std::cout << " Write Back the Atoms... " << std::endl;
    std::cout << "####################################################################################################" << std::endl;

    
    RigidAtomWriteBack();
    
    std::cout << "####################################################################################################" << std::endl;
    std::cout << " Done! " << std::endl;
    std::cout << "####################################################################################################" << std::endl;

}

void SiteRP_cont::StudyArtificialNetwork()
{
    std::string mychar = "Data_Detect.txt";
    std::string mychar_rcluster = "Data_Detect_rcluster.txt";
    
    myfile.close();
    myfile.open(mychar);
    
    rclusterfile.close();
    rclusterfile.open(mychar_rcluster);
    
    rigidcluster();
    
    int span = spanningrcluster();
    log();
    
    StoreRigidInfoOfSite();
    
}


void SiteRP_cont::ReadVerticeInfoFromCSV() {
    
    std::string ReadFileSite = "./MCdata/SitePos.csv";
    std::string ReadFileBond = "./MCdata/Bonds.csv";
    
    std::ifstream InputDataFileSite(ReadFileSite);
    std::ifstream InputDataFileBond(ReadFileBond);
    
    std::string line;
    std::string cell;
    int AtomIndex = 0;
    while (std::getline(InputDataFileSite, line)) {
        std::stringstream linestream(line);
        double AtomXCoord = 0.0;
        double AtomYCoord = 0.0;
        ++AtomIndex;
        int flag = 0;
        while (std::getline(linestream, cell, ',')) {
            flag++;
            switch (flag) {
                case 1:
                    AtomXCoord = std::stod(cell);
                    break;
                    
                case 2:
                    AtomYCoord = std::stod(cell);
                    break;
            }
        }
        site Newsite(AtomXCoord,AtomYCoord,AtomIndex);
        vertices.push_back(Newsite);
    }
    
    // then build the network with the bond information
    
    initemptytrigraph();
    
    int SiteNum = vertices.size();
    numparts = SiteNum;
    
    std::string lineBond;
    std::string cellBond;
    while (std::getline(InputDataFileBond, lineBond)) {
        std::stringstream linestreamBond(lineBond);
        int IAtom = 0;
        int JAtom = 0;
        int flagBond = 0;
        while (std::getline(linestreamBond, cellBond, ',')) {
            flagBond++;
            switch (flagBond) {
                case 1:
                    IAtom = std::stoi(cellBond)-1;
                    break;
                    
                case 2:
                    JAtom = std::stoi(cellBond)-1;
                    break;
            }
        }
        if (isempty(IAtom, JAtom) && IAtom != JAtom) {
            addbond(IAtom, JAtom);
            
            undirectedgraph[IAtom].push_back(JAtom);
            undirectedgraph[JAtom].push_back(IAtom);
        }
    }
    InputDataFileSite.close();
    InputDataFileBond.close();
    
}

void SiteRP_cont::ReadVerticeInfo() {
    
    std::cout << "####################################################################################################" << std::endl;
    std::cout << " Read Vertice Information... " << std::endl;
    std::cout << "####################################################################################################" << std::endl;

    
    std::string ReadFile = "./MCdata/";
    ReadFile.append(FileName);
    ReadFile.append(".data");

    std::ifstream InputDataFile(ReadFile);

    // read the data to vertices
    std::string line;

    unsigned int count = 0;
    int LineInitialNum = 15;
    int LineEndNum = LineInitialNum + size;

    while (std::getline(InputDataFile, line))
    {
        ++count;
        if (count > LineEndNum) { break; }    // done
        if (count == 6)
        {
            boost::char_separator<char> sep(" ");
            boost::tokenizer<boost::char_separator<char>> tokens(line, sep);
            int flag = 0;
            BOOST_FOREACH (const auto& t, tokens) {
                ++flag;
                if (flag == 2)
                {
                    BoxLength = std::stod(t);
                }
            }
        }
        if (count > LineInitialNum)
        {
            //std::cout << line << std::endl;
            boost::char_separator<char> sep(" ");
            boost::tokenizer<boost::char_separator<char>> tokens(line, sep);
            int flag = 0;
            int AtomIndex = 0;
            double AtomXCoord = 0.0;
            double AtomYCoord = 0.0;
            BOOST_FOREACH (const auto& t, tokens) {
                ++flag;
                switch (flag){
                    case 1: AtomIndex = std::stod(t);
                    case 3: AtomXCoord = std::stod(t);
                    case 4: AtomYCoord = std::stod(t);
                }
            }
            site Newsite(AtomXCoord,AtomYCoord,AtomIndex);
            vertices.push_back(Newsite);
        }
    }

    InputDataFile.close();
    InputDataFile.close();
    
    std::cout << "####################################################################################################" << std::endl;
    std::cout << " Done! " << std::endl;
    std::cout << "####################################################################################################" << std::endl;

}

void SiteRP_cont::RigidAtomWriteBack() {
    
    std::string ReadFile = "./MCdata/";
    ReadFile.append(FileName);
    ReadFile.append(".data");

    std::string WriteFile = "./MCdata/";
    WriteFile.append(FileName);
    WriteFile.append("OUT.data");

    std::ifstream ReadInputFile(ReadFile);
    std::ofstream WriteBackDataFile;
    WriteBackDataFile.open(WriteFile);

    // read the data to vertices
    std::string line;

    unsigned int count = 0;
    int LineInitialNum = 15;
    int LineEndNum = LineInitialNum + size;

    while (std::getline(ReadInputFile, line))
    {
        ++count;
        std::string newline;

        if (count > LineInitialNum && count <= LineEndNum)
        {
            //std::cout << line << std::endl;
            boost::char_separator<char> sep(" ");
            boost::tokenizer<boost::char_separator<char>> tokens(line, sep);
            int flag = 0;

            for (const auto& t : tokens) {
                ++flag;
                if (!giantrigidcluster[count-LineInitialNum-1].empty() && flag == 2){
                    newline.append("2 ");
                } else {
                    newline.append(t);
                    newline.append(" ");
                }
            }
            newline.append("\n");

            WriteBackDataFile << newline;
            //std::cout << newline << std::endl;

        } else if (count == 4){
            line[0] = '2';
            WriteBackDataFile << line << "\n";
        } else {
            WriteBackDataFile << line << "\n";
        }
    }

    ReadInputFile.close();
    WriteBackDataFile.close();
}

void SiteRP_cont::CoordNumber() {
    
    std::cout << "####################################################################################################" << std::endl;
    std::cout << " Write the CoordNumber... " << std::endl;
    std::cout << "####################################################################################################" << std::endl;

    
    std::fill_n(CoordDist, size, 0);
    for (int i = 0; i < size; ++i) {
        int CoordNumIndex = undirectedgraph[i].size();
        CoordDist[CoordNumIndex]++;
    }

    double sum = 0;

    for (int j = 0; j < size; ++j) {
        std::cout << "Coord " << j << " : " << CoordDist[j] << std::endl;
        sum += j*CoordDist[j];
    }

    double CoordNum = sum/size * 2;

    std::cout << "Coordination Number: " << CoordNum << std::endl;
    
    std::cout << "####################################################################################################" << std::endl;
    std::cout << " Done! " << std::endl;
    std::cout << "####################################################################################################" << std::endl;
    
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// READING FILE TO GET VERTICES

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// END OF THE CLASS

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
